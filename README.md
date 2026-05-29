# Tn5 Insertion Mapping Pipeline

![Language](https://img.shields.io/badge/language-Bash%20%7C%20Python-blue)

This pipeline takes raw paired-end FASTQ files from a transposon insertion enrichment sequencing experiment and produces a complete map of insertion sites across the *C. necator* H16 genome. Each read pair contains the junction between the transposon end sequence and flanking genomic DNA; the pipeline trims away the known transposon and adaptor sequences, aligns the remaining genomic sequence to the reference, and extracts the precise insertion coordinate from the 5′ end of R1.

Jump to [How to Run](#how-to-run)

## Overview

- **Organism:** *Cupriavidus necator* H16
- **Reference:** `GCF_000009285.1` (chromosome + megaplasmid pHG1)
- **Input:** paired-end FASTQ files from Tn5 insertion-junction sequencing
- **Key coordinate definition:** the **insertion site** is the **first genomic base of trimmed R1**

---

## Experimental Context

The library was prepared as follows:

1. Tn5 transposon (with cargo) randomly inserted into the *C. necator* genome
2. gDNA purified and enzymatically fragmented
3. Illumina adaptors ligated to fragment ends
4. **Enrichment PCR:** forward primer binds to inside the transposon, reverse primer binds to the Illumina adaptor
5. PE150 sequencing

This means every read pair that passes enrichment contains one transposon–genome junction, and the genomic sequence flanking the junction is what reveals the insertion site.

---

## Read Structure

```
R1  5'→  [variable primer: cgttcttctgctgctcggggatctgg]
         [Tn5 mosaic end:   TATAAGAGACAG               ]  ← trimmed
         [genomic DNA →                                 ]  ← retained, used for alignment
         [RC of rev primer, if insert short             ]  ← trimmed from 3' end

R2  5'→  [rev enrichment primer: GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC]  ← trimmed (anchored)
         [genomic DNA ←                                             ]  ← retained
         [RC of Tn5 ME, if insert short: CTGTCTTATATACA             ]  ← trimmed from 3' end
```

The **insertion site** is the first genomic base retained in R1 after trimming — i.e., the base immediately 3′ of `TATAAGAGACAG`.

---

## Pipeline Scripts

| File | Role |
|---|---|
| `setup_reference.sh` | One-time: download *C. necator* H16 reference and build Bowtie2 index |
| `tn5_insertion_mapping.sh` | Main per-sample pipeline (all 7 steps) |
| `annotate_insertions.py` | Annotates insertion sites against GFF3 genes; called by the main script |
| `qc_report.py` | Generates QC plots and summary table; called by the main script |
| `tn5_slurm_array.sh` | SLURM array wrapper to process multiple samples in parallel on a cluster |

---

## Inputs

### Required per sample

| Input | Description |
|---|---|
| `sample_R1.fastq.gz` | R1 reads (gzipped FASTQ) |
| `sample_R2.fastq.gz` | R2 reads (gzipped FASTQ) |
| Bowtie2 index | Built from the *C. necator* H16 FASTA by `setup_reference.sh` |
| GFF3 annotation | *C. necator* H16 gene annotation from NCBI (downloaded by `setup_reference.sh`) |

### Sample sheet (for multi-sample SLURM runs)

A plain tab-separated file, no header, one sample per line:

```
lib1    /data/fastq/lib1_R1.fastq.gz    /data/fastq/lib1_R2.fastq.gz
lib2    /data/fastq/lib2_R1.fastq.gz    /data/fastq/lib2_R2.fastq.gz
lib3    /data/fastq/lib3_R1.fastq.gz    /data/fastq/lib3_R2.fastq.gz
```

---

## Outputs

All outputs are written to a per-sample subdirectory under `OUTDIR_ROOT/SAMPLE_NAME/`.

### Alignment files (`aligned/`)

| File | Description |
|---|---|
| `SAMPLE_sorted.bam` + `.bai` | All mapped read pairs, coordinate-sorted and indexed |
| `SAMPLE_dedup.bam` + `.bai` | Coordinate-pair-deduplicated BAM (R1-5′ + R2-5′ position unique) |

### Insertion site files (`sites/`)

| File | Description |
|---|---|
| `SAMPLE_raw_insertion_sites.bed` | BED6: all unique insertion sites, score = raw read count |
| `SAMPLE_dedup_insertion_sites.bed` | BED6: insertion sites from dedup reads only |
| `SAMPLE_raw_5prime_plus.bedgraph` | Plus-strand 5′ pileup (IGV-ready) |
| `SAMPLE_raw_5prime_minus.bedgraph` | Minus-strand 5′ pileup (IGV-ready) |

### Annotation and count files (`counts/`)

| File | Description |
|---|---|
| `SAMPLE_insertions_annotated.tsv` | Per-site table with gene annotation, raw and dedup counts, relative position within gene |
| `SAMPLE_per_gene_counts.tsv` | Per-gene aggregated table: unique sites, total reads, sites-per-kb, mean relative position |

### QC files (`qc/`)

| File | Description |
|---|---|
| `SAMPLE_qc_report.pdf` | 3-page QC report (see QC section below) |
| `SAMPLE_qc_summary.tsv` | Flat TSV of all QC metrics (for aggregation across samples) |

### Stats and logs (`stats/`, `logs/`)

| File | Description |
|---|---|
| `SAMPLE_cutadapt.json` | Trimming statistics (structured JSON) |
| `SAMPLE_bowtie2_stats.txt` | Bowtie2 alignment summary |
| `SAMPLE_samtools_stats.txt` | Full samtools stats output (includes insert size histogram) |
| `SAMPLE_dedup_stats.txt` | Coordinate-pair deduplication metrics |
| `SAMPLE_summary_stats.tsv` | High-level summary: sites, genes hit, duplication rate |
| `SAMPLE.pipeline.log` | Full stderr/stdout log of the entire run |

---

## Pipeline Steps in Detail

### Step 1 — Adapter and transposon trimming (cutadapt)

Four sequences are trimmed simultaneously in a single cutadapt pass:

| Direction | Sequence | What it removes |
|---|---|---|
| R1 5′ (unanchored) | `TATAAGAGACAG` | Tn5 ME end; everything upstream (variable primer) is also removed |
| R1 3′ | RC of rev primer | Bleedthrough when the genomic insert is shorter than 150 bp |
| R2 5′ (anchored `^`) | `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC` | Rev enrichment primer; always at position 0 of R2 |
| R2 3′ | RC of Tn5 ME | Bleedthrough for short inserts |

Read pairs where `TATAAGAGACAG` is **not found** in R1 are discarded (`--discard-untrimmed`). These are background non-enriched fragments. Pairs shorter than 30 bp after trimming are also discarded.

### Step 2 — Alignment (Bowtie2)

Trimmed pairs are aligned to the *C. necator* H16 reference (chromosome NC_008313.2 + megaplasmid pHG1 NC_005241.1) in paired, concordant-only mode. Read pairs that do not align concordantly as a proper pair are discarded (`--no-mixed --no-discordant`). Reads with mapping quality < 20 (probability of mismapping > 1%) are filtered out.

### Step 3 — Coordinate-pair deduplication

Without molecular barcodes (UMIs), true PCR duplicates cannot be distinguished from two independent insertions at the same genomic position. The pipeline handles this by collapsing read pairs with identical **R1 5′ coordinate + R2 5′ coordinate** to a single count. The 5′ coordinate is strand-aware:
- Forward strand: `POS` (leftmost alignment base)
- Reverse strand: `POS + CIGAR_ref_length − 1` (rightmost alignment base = read 5′ end)

Both the raw BAM (all reads) and deduplicated BAM are retained, and all downstream outputs report **both raw and dedup counts** so you can assess the impact of duplication on your data.

### Step 4 — Insertion site extraction

The insertion site is extracted as the **5′ coordinate of each aligned R1 read**. This is the first genomic base immediately downstream of the trimmed Tn5 ME sequence — the exact position where the transposon integrated. Sites within 1 bp on the same strand are merged to account for minor end-repair positional jitter (a known ±1 bp artefact of enzymatic fragmentation). Output is a BED6 file with read counts per site.

### Step 5 — Strand-specific bedGraphs

Per-base 5′-end pileup tracks are written for both strands of the raw (all-reads) BAM. These can be loaded directly into IGV alongside the BAM to visually inspect insertion site calls.

### Step 6 — Gene annotation (annotate_insertions.py)

Each insertion site is intersected against the GFF3 gene models using `bedtools intersect`. For each site the output records:
- Whether it falls inside a gene (intragenic) or between genes (intergenic)
- Gene ID, gene name, and locus tag from the GFF3
- Relative position within the gene body (0.0 = transcription start, 1.0 = stop codon), in the direction of transcription

A per-gene aggregated table is produced with unique site counts, total read support, and sites-per-kilobase (a library saturation proxy) reported separately for raw and dedup counts.

### Step 7 — QC report (qc_report.py)

A 3-page PDF and a flat TSV of all QC metrics are generated. See the QC section below.

---

## QC Report

The QC report (`SAMPLE_qc_report.pdf`) contains 8 plots across 3 pages.

### Page 1 — Read processing

| Plot | What to look for |
|---|---|
| **Read fate funnel** | Gradual step-down is normal. A large drop at trimming means few reads contained the Tn5 ME (enrichment failed or wrong sequence). A large drop at alignment means mapping problems (wrong reference, poor trim). |
| **Trimming outcome** | The majority should be "Passed". High "No Tn5 ME found" fraction = enrichment or sequence problem. High "Too short" = over-fragmentation or insert size issue. |
| **Insert size distribution** | Expect a peak between 100–500 bp. Very short peak (<50 bp) = over-fragmentation. Very broad or flat distribution = under-fragmentation or ligation artefacts. |
| **Alignment categories** | Unique concordant should dominate (>70% is good). High multi-mapping (>20%) can indicate repetitive regions or index contamination. |

### Page 2 — Library complexity and insertion sites

| Plot | What to look for |
|---|---|
| **Duplication rate** | Some duplication is expected without UMIs. >80% duplication suggests the library was under-complex (few independent insertion events) or over-amplified. |
| **Site read depth histogram** | A right-skewed distribution (most sites covered by few reads, a tail of highly covered sites) is normal. A very narrow spike suggests low complexity. |
| **Cumulative site coverage** | Shows what fraction of sites have ≥N reads. A steep drop-off at low N is expected; you want most sites to have ≥2–3 reads for confident calling. |
| **Genomic feature distribution** | Tn5 inserts semi-randomly; expect 70–90% intragenic given *C. necator* has a compact, gene-dense genome (~87% coding). Very low intragenic fraction is unexpected. Also shows how many of the ~7,400 annotated genes were hit. |

### Page 3 — Summary metrics table

All key metrics in a single table with pass (green) / fail (red) highlighting based on the following thresholds:

| Metric | Acceptable range |
|---|---|
| % reads passing trimming | > 50% |
| % reads mapped | > 50% |
| % properly paired | > 50% |
| Duplication rate | < 80% |
| Median insert size | 50–900 bp |
| % intragenic insertions | > 20% |

---

## Deduplication Note

Because this library has **no UMIs or barcodes**, it is impossible to definitively distinguish a PCR duplicate from two genuine independent insertions at the same genomic coordinate. The pipeline reports both:

- **Raw counts** — all reads supporting a site (may include PCR duplicates)
- **Dedup counts** — coordinate-pair-deduplicated reads (may under-count sites in PCR hot spots)

For most analyses (identifying insertion site positions, asking which genes are hit) use **raw counts**. For quantitative analyses comparing insertion frequency between sites or samples, **dedup counts** are more conservative. If duplication rate is low (<30%), the two counts will agree closely.

---

## How to Run

### Step 0 — Install dependencies (once)

The recommended approach is a conda environment:

```bash
conda create -n tn5_mapping -c bioconda -c conda-forge \
    cutadapt bowtie2 samtools bedtools \
    pandas matplotlib seaborn python=3.10
conda activate tn5_mapping
```

Minimum versions: cutadapt ≥ 3.0, bowtie2 ≥ 2.4, samtools ≥ 1.15, bedtools ≥ 2.30.

---

### Step 1 — Build the reference (once)

Edit the `REFS_DIR` path in `setup_reference.sh`, then run it:

```bash
# Edit REFS_DIR="/path/to/refs/cnecator_H16" inside setup_reference.sh first
bash setup_reference.sh
```

This will:
1. Download the *C. necator* H16 FASTA and GFF3 from NCBI (accession GCF_000009285.1)
2. Build the Bowtie2 index (covers both the chromosome and the pHG1 megaplasmid)
3. Index the FASTA with `samtools faidx`

**Output paths** (update these in `tn5_slurm_array.sh`):
```
BOWTIE2_INDEX = <REFS_DIR>/cnecator_H16
GFF           = <REFS_DIR>/cnecator_H16_genomic.gff
```

---

### Step 2 — Create your sample sheet

Create a tab-separated file (no header) listing all samples:

```
lib1    /abs/path/to/lib1_R1.fastq.gz    /abs/path/to/lib1_R2.fastq.gz
lib2    /abs/path/to/lib2_R1.fastq.gz    /abs/path/to/lib2_R2.fastq.gz
```

Use **absolute paths**. Save as e.g. `/path/to/samples.tsv`.

---

### Step 3a — Run a single sample (test / local)

```bash
bash tn5_insertion_mapping.sh \
    -1 /path/to/sample_R1.fastq.gz \
    -2 /path/to/sample_R2.fastq.gz \
    -s SAMPLE_NAME \
    -r /path/to/refs/cnecator_H16/cnecator_H16 \
    -g /path/to/refs/cnecator_H16/cnecator_H16_genomic.gff \
    -o /path/to/output/SAMPLE_NAME \
    -t 8
```

For a quick test, first subset your FASTQ to 100,000 read pairs:

```bash
zcat sample_R1.fastq.gz | head -400000 | gzip > test_R1.fastq.gz
zcat sample_R2.fastq.gz | head -400000 | gzip > test_R2.fastq.gz
```

---

### Step 3b — Run multiple samples on a SLURM cluster

Edit the five paths at the top of `tn5_slurm_array.sh`:

```bash
SCRIPT_DIR="/path/to/pipeline/scripts"
BOWTIE2_INDEX="/path/to/refs/cnecator_H16/cnecator_H16"
GFF="/path/to/refs/cnecator_H16/cnecator_H16_genomic.gff"
SAMPLE_SHEET="/path/to/samples.tsv"
OUTDIR_ROOT="/path/to/output"
```

Set the array range to match your sample count:
```bash
#SBATCH --array=1-N    # N = number of lines in samples.tsv
```

Update the `module load` lines to match your cluster's module names (or switch to the conda block). Then submit:

```bash
sbatch tn5_slurm_array.sh
```

Monitor progress:
```bash
squeue -u $USER
sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS
```

---

### Step 4 — Check outputs

After each sample completes, review:

1. `logs/SAMPLE.pipeline.log` — confirm each step completed without errors
2. `qc/SAMPLE_qc_report.pdf` — check all QC metrics (see QC section above)
3. `stats/SAMPLE_summary_stats.tsv` — quick numeric summary
4. `counts/SAMPLE_per_gene_counts.tsv` — your primary analysis output

Load into IGV for visual inspection:
- Open `aligned/SAMPLE_dedup.bam`
- Load `sites/SAMPLE_raw_5prime_plus.bedgraph` and `sites/SAMPLE_raw_5prime_minus.bedgraph`
- Set reference genome to *C. necator* H16

---

## Resource Requirements

| Resource | Recommendation | Notes |
|---|---|---|
| CPUs | 16 | Used by bowtie2 and cutadapt; diminishing returns above 16 |
| Memory | 32 GB | Bowtie2 index for *C. necator* is small (~200 MB); headroom is for samtools sort |
| Wall time | 4 hours | Typical for 20–50M read pairs PE150; adjust down for smaller libraries |
| Disk (per sample) | ~20 GB | BAMs dominate; trimmed FASTQs can be deleted after alignment if space is tight |

---

## Dependencies Summary

| Tool | Version | Purpose |
|---|---|---|
| cutadapt | ≥ 3.0 | Adapter and transposon end trimming |
| bowtie2 | ≥ 2.4 | Short read alignment |
| samtools | ≥ 1.15 | BAM sorting, indexing, flagstat, stats |
| bedtools | ≥ 2.30 | Genomic interval operations, coverage |
| python3 | ≥ 3.8 | Deduplication, annotation, QC plotting |
| pandas | ≥ 1.3 | Data tables |
| matplotlib | ≥ 3.5 | QC plots |
| seaborn | ≥ 0.11 | QC plot styling |
| NCBI datasets CLI | optional | Faster reference download in setup_reference.sh |
