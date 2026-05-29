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

## Pipeline

1. **Trim reads with cutadapt**
   - Remove the first 8 bp from R1 with `--cut 8`
   - Trim anchored primer from R1: `^TTCTTCTGCTGCTCGGGGATCTGG`
   - Trim reverse-complement Tn5 ME sequence from the 3' end of R2
2. **Align with Bowtie2**
   - Paired-end alignment
   - Concordant pairs only
   - Filter to `MAPQ >= 20`
3. **Deduplicate read pairs**
   - Coordinate-pair deduplication using **R1 5' + R2 5' positions**
   - **No UMIs are used**
4. **Extract insertion sites** from the R1 5' coordinate
5. **Generate strand-specific bedGraphs**
6. **Annotate insertions** against GFF3 genes with `bedtools intersect`
7. **Build QC PDF report**
8. **Create Excel workbook and genome plots**

---

## Quick Start

```bash
# 1. Create and activate the conda environment
conda env create -f environment.yml
conda activate tn5_mapping

# 2. One-time reference setup (edit REFS_DIR inside the script first)
bash setup_reference.sh

# 3. Create samples.tsv (tab-separated: SAMPLE  R1_path  R2_path)

# 4. Edit paths and set --array=1-N in tn5_slurm_array.sh, then submit
sbatch tn5_slurm_array.sh

# 5. Generate Excel workbook and genome figures
bash run_excel_and_plots.sh
```

---

## Outputs Per Sample

| Output | Description |
|---|---|
| `*.sorted.bam` | coordinate-sorted aligned reads |
| `*.dedup.bam` | deduplicated BAM |
| `*.bed` | insertion site calls (raw and dedup) |
| `*.bedGraph` | strand-specific 5′-end pileup tracks (IGV-ready) |
| `*.annotated.tsv` | insertion sites annotated to genes with raw and dedup counts |
| `*_per_gene_counts.tsv` | insertion counts aggregated by gene |
| `*_qc_report.pdf` | 3-page QC report (read funnel, insert size, duplication, site depth) |
| `*_summary_stats.tsv` | sample-level metrics |

### Summary table (`insertion_sites.xlsx`)

A multi-sheet Excel workbook consolidating all results across samples:

| Sheet | Contents |
|---|---|
| `Summary` | QC metrics for all samples side-by-side (reads, alignment rate, duplication rate, sites called, genes hit) |
| `[sample]_sites` | Every unique insertion site with genomic coordinates, raw and dedup read counts, gene annotation, and fractional position within the gene body |
| `[sample]_per_gene` | Per-gene aggregated counts: unique sites, total reads, sites per kb, mean insertion position |
| `combined_genes` | All samples merged on gene ID — side-by-side dedup read counts per gene, sorted by total support |

### Genome-wide visualizations

**Insertion map** (`genome_insertions_[sample].pdf`)
One panel per replicon (chromosome and megaplasmid pHG1). Each insertion site is drawn as a vertical line at its genomic coordinate — upward for the + strand, downward for the − strand. Line height is proportional to log(dedup read count). A combined figure overlaying all samples is also produced.

**Insertion density** (`insertion_density_[sample].pdf`)
Sliding-window plot (default 10 kb bins) showing the number of unique insertion sites per genomic window. Regions with consistently zero coverage across all samples are candidate essential gene regions. A combined figure overlaying all samples is produced for direct comparison.

---

## Repository Structure

| File | Description |
|---|---|
| `setup_reference.sh` | Download the *C. necator* H16 reference and build the Bowtie2 index |
| `tn5_insertion_mapping.sh` | Main per-sample pipeline script |
| `annotate_insertions.py` | Annotate insertion sites against GFF3 genes |
| `qc_report.py` | Generate a 3-page QC PDF report |
| `tn5_slurm_array.sh` | Run multiple samples as a SLURM array job |
| `make_excel_and_plots.py` | Generate Excel workbook and genome-wide figures |
| `run_excel_and_plots.sh` | Cluster wrapper for Excel/plot generation |
| `environment.yml` | Conda environment specification |
| `README.md` | Full user guide and detailed usage notes |

---

## Dependencies

- `cutadapt`
- `bowtie2`
- `samtools`
- `bedtools`
- `python3`
  - `pandas`
  - `matplotlib`
  - `seaborn`
  - `openpyxl`

---

## Deduplication Note

This pipeline performs **coordinate-pair deduplication without UMIs**. That is appropriate for this library design, but duplicate removal is based on identical mapped read-pair coordinates rather than molecule barcodes.
