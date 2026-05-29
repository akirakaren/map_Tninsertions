#!/usr/bin/env bash
# =============================================================================
# Tn5 Transposon Insertion Site Mapping Pipeline
# Organism: Cupriavidus necator H16 (GCF_000009285.1)
# Library:  PE150, no barcodes, no UMIs
#
# Read structure (confirmed from raw read inspection):
#
#   R1: 5'→ [variable prefix: ~7-8 bp junk]
#           [constant junction: TTCTTCTGCTGCTCGGGGATCTGGTATAAGAGACAG] ← trimmed
#           [genomic DNA →]                                            ← retained
#
#   R2: 5'→ [genomic DNA ←]                                           ← retained
#           [RC of Tn5 junction, if insert is short: CTGTCTTATATACA]  ← trimmed from 3'
#
# Trimming strategy:
#   R1 5': unanchored -g TTCTTCTGCTGCTCGGGGATCTGGTATAAGAGACAG
#   R2 3': -A CTGTCTTATATACA (short insert bleedthrough only)
#   --discard-untrimmed --pair-filter=first
#
# Deduplication:
#   No UMIs. Coordinate-pair dedup: collapse pairs with identical
#   R1-5'-coordinate + R2-5'-coordinate. Both raw and dedup counts reported.
#
# Outputs per sample (in OUTDIR/):
#   trimmed/  — trimmed FASTQs
#   aligned/  — sorted BAM (all reads), dedup BAM
#   sites/    — insertion site BEDs, strand bedGraphs
#   counts/   — annotated site TSV, per-gene count table
#   stats/    — cutadapt report, bowtie2 stats, dedup stats
#   qc/       — QC plots (PDF) and QC summary TSV
#   logs/     — full pipeline log
#
# Dependencies (on PATH or loaded via modules):
#   cutadapt  >= 3.0
#   bowtie2   >= 2.4
#   samtools  >= 1.15
#   bedtools  >= 2.30
#   python3   with pandas, matplotlib, seaborn
#
# Usage:
#   bash tn5_insertion_mapping.sh \
#       -1 sample_R1.fastq.gz \
#       -2 sample_R2.fastq.gz \
#       -s SAMPLE_NAME \
#       -r /path/to/bowtie2_index/cnecator_H16 \
#       -g /path/to/cnecator_H16.gff \
#       -o /path/to/output_dir \
#       [-t THREADS]
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# SEQUENCE CONFIGURATION
# ---------------------------------------------------------------------------

# Tn5 mosaic end — found after the variable primer at 5' of R1
TN5_END="TTCTTCTGCTGCTCGGGGATCTGG"
# RC of Tn5 ME — appears at 3' of R2 when insert is short
TN5_END_RC="CTGTCTTATATACA"

# Reverse enrichment primer — always at the very 5' of R2 (anchored)
REV_PRIMER="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"
# RC of reverse primer — appears at 3' of R1 when insert is short

# Minimum overlap for adapter matching (bp)
TN5_OVERLAP=24
PRIMER_OVERLAP=15

# Minimum read length after trimming
MIN_LEN=30

# Bowtie2 settings
ALIGN_MODE="end-to-end"
BOWTIE2_PRESET="--sensitive"
MAX_INSERT=1000

# Mapping quality filter
MIN_MAPQ=20

# Default threads
THREADS=8

# ---------------------------------------------------------------------------
# ARGUMENT PARSING
# ---------------------------------------------------------------------------
usage() {
    echo ""
    echo "Usage: $0 -1 R1.fastq.gz -2 R2.fastq.gz -s SAMPLE -r BOWTIE2_INDEX -g GFF -o OUTDIR [-t THREADS]"
    echo ""
    echo "  -1   R1 FASTQ (gzipped)"
    echo "  -2   R2 FASTQ (gzipped)"
    echo "  -s   Sample name"
    echo "  -r   Bowtie2 index prefix"
    echo "  -g   GFF3 annotation file"
    echo "  -o   Output directory"
    echo "  -t   Threads [default: 8]"
    echo ""
    exit 1
}

while getopts "1:2:s:r:g:o:t:" opt; do
    case $opt in
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        s) SAMPLE="$OPTARG" ;;
        r) BOWTIE2_INDEX="$OPTARG" ;;
        g) GFF="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

for VAR in R1 R2 SAMPLE BOWTIE2_INDEX GFF OUTDIR; do
    [[ -z "${!VAR:-}" ]] && { echo "ERROR: missing required argument for ${VAR}"; usage; }
done

# ---------------------------------------------------------------------------
# DIRECTORY SETUP AND LOGGING
# ---------------------------------------------------------------------------
mkdir -p "$OUTDIR/logs" "$OUTDIR/trimmed" "$OUTDIR/aligned" \
         "$OUTDIR/sites" "$OUTDIR/counts" "$OUTDIR/stats" "$OUTDIR/qc"

LOG="$OUTDIR/logs/${SAMPLE}.pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "======================================================"
echo " Tn5 Insertion Mapping Pipeline"
echo " Sample     : $SAMPLE"
echo " R1         : $R1"
echo " R2         : $R2"
echo " Reference  : $BOWTIE2_INDEX"
echo " Annotation : $GFF"
echo " Output     : $OUTDIR"
echo " Threads    : $THREADS"
echo " Started    : $(date)"
echo "======================================================"

# ---------------------------------------------------------------------------
# STEP 1: TRIMMING WITH CUTADAPT
# ---------------------------------------------------------------------------
# R1 5' (-g): unanchored — find TATAAGAGACAG after the variable primer region.
#             Removes primer + ME; first retained base = genomic insertion site.
# R1 3' (-a): remove RC of reverse primer if R1 reads through short insert.
# R2 5' (-G): anchored (^) — reverse primer is ALWAYS at position 0 of R2.
#             First retained base = first genomic base from the fragmentation end.
# R2 3' (-A): remove RC of Tn5 ME if R2 reads through short insert.
#
# --discard-untrimmed: discard pair if Tn5 ME not found in R1.
#   This is the enrichment filter — any read that doesn't contain the ME
#   is a background fragment, not a transposon-junction read.
# --pair-filter=first: apply length filter to either read in the pair.
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 1] Trimming adapters and transposon sequences..."

TRIM_R1="$OUTDIR/trimmed/${SAMPLE}_R1_trimmed.fastq.gz"
TRIM_R2="$OUTDIR/trimmed/${SAMPLE}_R2_trimmed.fastq.gz"
TRIM_REPORT="$OUTDIR/stats/${SAMPLE}_cutadapt.json"

cutadapt \
    --cut 8 \
    -g "^${TN5_END}" \
    -A "${TN5_END_RC}" \
    --overlap ${TN5_OVERLAP} \
    --discard-untrimmed \
    --pair-filter=first \
    -m ${MIN_LEN} \
    -e 0.1 \
    -q 20,20 \
    -j ${THREADS} \
    -o "${TRIM_R1}" \
    -p "${TRIM_R2}" \
    --json "${TRIM_REPORT}" \
    "${R1}" "${R2}" \
    2>&1 | tee "$OUTDIR/stats/${SAMPLE}_cutadapt.txt"

echo "  Trimming complete. JSON report: ${TRIM_REPORT}"

# ---------------------------------------------------------------------------
# STEP 2: ALIGN TO REFERENCE WITH BOWTIE2
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 2] Aligning to reference genome..."

SORTED_BAM="$OUTDIR/aligned/${SAMPLE}_sorted.bam"
ALIGN_STATS="$OUTDIR/stats/${SAMPLE}_bowtie2_stats.txt"

bowtie2 \
    -x "${BOWTIE2_INDEX}" \
    -1 "${TRIM_R1}" \
    -2 "${TRIM_R2}" \
    --${ALIGN_MODE} \
    ${BOWTIE2_PRESET} \
    -X ${MAX_INSERT} \
    --no-mixed \
    --no-discordant \
    -p ${THREADS} \
    2>"${ALIGN_STATS}" \
| samtools view -bS -F 4 -q ${MIN_MAPQ} - \
| samtools sort -@ ${THREADS} -o "${SORTED_BAM}"

samtools index "${SORTED_BAM}"

echo "  Alignment stats:"
cat "${ALIGN_STATS}"
echo "  BAM: ${SORTED_BAM}"

# Compute insert size metrics from BAM for QC (samtools stats)
SAMTOOLS_STATS="$OUTDIR/stats/${SAMPLE}_samtools_stats.txt"
samtools stats "${SORTED_BAM}" > "${SAMTOOLS_STATS}"
echo "  samtools stats: ${SAMTOOLS_STATS}"

# ---------------------------------------------------------------------------
# STEP 3: COORDINATE-PAIR DEDUPLICATION
# ---------------------------------------------------------------------------
# Without UMIs, collapse read pairs with identical R1-5' + R2-5' coordinates.
# 5' coordinate:
#   + strand: POS - 1 (0-based leftmost base)
#   - strand: POS + cigar_ref_length - 1 (0-based rightmost base = read 5' in genome)
# Outputs both a dedup BAM and dedup statistics for QC.
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 3] Coordinate-pair deduplication..."

DEDUP_BAM="$OUTDIR/aligned/${SAMPLE}_dedup.bam"
DEDUP_STATS="$OUTDIR/stats/${SAMPLE}_dedup_stats.txt"
DEDUP_READNAMES="$OUTDIR/aligned/${SAMPLE}_dedup_readnames.txt"

python3 - "${SORTED_BAM}" "${DEDUP_READNAMES}" "${DEDUP_STATS}" << 'PYEOF'
import sys, re, subprocess

bam_path  = sys.argv[1]
out_names = sys.argv[2]
out_stats = sys.argv[3]

def cigar_ref_length(cigar):
    return sum(int(n) for n, op in re.findall(r'(\d+)([MIDNX=])', cigar))

def five_prime_coord(pos0, cigar, flag):
    if flag & 0x10:
        return pos0 + cigar_ref_length(cigar) - 1
    return pos0

reads = {}
proc = subprocess.Popen(
    ["samtools", "view", "-f", "0x2", bam_path],
    stdout=subprocess.PIPE, text=True
)
for line in proc.stdout:
    f = line.split("\t")
    qname = f[0]; flag = int(f[1]); chrom = f[2]
    pos0  = int(f[3]) - 1; cigar = f[5]
    mate  = 1 if (flag & 0x40) else 2
    strand = "-" if (flag & 0x10) else "+"
    coord  = five_prime_coord(pos0, cigar, flag)
    key    = f"{chrom}:{coord}:{strand}"
    if qname not in reads:
        reads[qname] = {}
    reads[qname][mate] = key
proc.wait()

seen_pairs = {}
total = dedup = no_mate = 0
for qname, mates in reads.items():
    if 1 not in mates or 2 not in mates:
        no_mate += 1; continue
    total += 1
    pk = (mates[1], mates[2])
    if pk not in seen_pairs:
        seen_pairs[pk] = qname; dedup += 1

with open(out_names, "w") as fh:
    for name in seen_pairs.values():
        fh.write(name + "\n")

dups = total - dedup
rate = dups / total if total > 0 else 0.0
with open(out_stats, "w") as fh:
    fh.write(f"total_properly_paired\t{total}\n")
    fh.write(f"unique_coordinate_pairs\t{dedup}\n")
    fh.write(f"duplicate_pairs_removed\t{dups}\n")
    fh.write(f"duplication_rate\t{rate:.4f}\n")
    fh.write(f"pairs_missing_mate\t{no_mate}\n")

print(f"  Total properly paired : {total}")
print(f"  Unique coord pairs    : {dedup}")
print(f"  Duplicates removed    : {dups}")
print(f"  Duplication rate      : {rate:.2%}")
PYEOF

# Filter BAM to dedup read names
samtools view -h "${SORTED_BAM}" \
| awk -v nf="${DEDUP_READNAMES}" '
    BEGIN { while ((getline line < nf) > 0) keep[line]=1 }
    /^@/ { print; next }
    { if ($1 in keep) print }
' \
| samtools view -bS - \
| samtools sort -@ ${THREADS} -o "${DEDUP_BAM}"

samtools index "${DEDUP_BAM}"
echo "  Dedup BAM: ${DEDUP_BAM}"

# ---------------------------------------------------------------------------
# STEP 4: EXTRACT INSERTION SITES
# ---------------------------------------------------------------------------
# Insertion site = 5' end of aligned R1 = first genomic base after Tn5 ME.
# Two site BEDs produced: from all reads (raw) and from dedup reads.
# Sites within 1bp on the same strand are merged (end-repair jitter).
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 4] Extracting insertion sites from R1 5' coordinates..."

extract_sites() {
    local bam="$1"
    local out_raw="$2"
    local out_merged="$3"

    samtools view -f 0x40 -F 0x4 "${bam}" \
    | awk 'BEGIN{OFS="\t"} {
        flag=int($2); chrom=$3; pos0=$4-1; cigar=$6; mapq=$5;
        ref_len=0; s=cigar;
        while (match(s,/[0-9]+[MIDNX=]/)) {
            tok=substr(s,RSTART,RLENGTH);
            ref_len+=substr(tok,1,length(tok)-1)+0;
            s=substr(s,RSTART+RLENGTH);
        }
        if (and(flag,16)) { site=pos0+ref_len-1; strand="-"; }
        else              { site=pos0;            strand="+"; }
        print chrom, site, site+1, "INS", mapq, strand;
    }' > "${out_raw}"

    bedtools sort -i "${out_raw}" \
    | bedtools merge -i stdin -s -d 1 -c 4,6 -o count,distinct \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"INS_"NR,$4,$5}' \
    > "${out_merged}"
}

RAW_SITES_BED="$OUTDIR/sites/${SAMPLE}_raw_read_sites.bed"
RAW_MERGED="$OUTDIR/sites/${SAMPLE}_raw_insertion_sites.bed"
DEDUP_SITES_BED="$OUTDIR/sites/${SAMPLE}_dedup_read_sites.bed"
DEDUP_MERGED="$OUTDIR/sites/${SAMPLE}_dedup_insertion_sites.bed"

extract_sites "${SORTED_BAM}" "${RAW_SITES_BED}"   "${RAW_MERGED}"
extract_sites "${DEDUP_BAM}"  "${DEDUP_SITES_BED}" "${DEDUP_MERGED}"

echo "  Raw unique sites  : $(wc -l < "${RAW_MERGED}")"
echo "  Dedup unique sites: $(wc -l < "${DEDUP_MERGED}")"

# ---------------------------------------------------------------------------
# STEP 5: STRAND-SPECIFIC 5'-END COVERAGE (bedGraph for IGV)
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 5] Generating strand-specific 5'-end coverage bedGraphs..."

BG_PLUS="$OUTDIR/sites/${SAMPLE}_raw_5prime_plus.bedgraph"
BG_MINUS="$OUTDIR/sites/${SAMPLE}_raw_5prime_minus.bedgraph"

# Plus-strand R1: first in pair (0x40), NOT reverse (exclude 0x10), NOT unmapped (exclude 0x4)
samtools view -b -f 0x40 -F 0x14 "${SORTED_BAM}" \
| bedtools genomecov -ibam stdin -bg -5 > "${BG_PLUS}"

# Minus-strand R1: first in pair AND reverse (0x50 = 0x40|0x10), NOT unmapped
samtools view -b -f 0x50 -F 0x04 "${SORTED_BAM}" \
| bedtools genomecov -ibam stdin -bg -5 > "${BG_MINUS}"

echo "  Plus-strand bedGraph : ${BG_PLUS}"
echo "  Minus-strand bedGraph: ${BG_MINUS}"

# ---------------------------------------------------------------------------
# STEP 6: ANNOTATE INSERTIONS
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 6] Annotating insertion sites against gene features..."

ANNOTATED="$OUTDIR/counts/${SAMPLE}_insertions_annotated.tsv"
GENE_COUNTS="$OUTDIR/counts/${SAMPLE}_per_gene_counts.tsv"
STATS_REPORT="$OUTDIR/stats/${SAMPLE}_summary_stats.tsv"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "${SCRIPT_DIR}/annotate_insertions.py" \
    --raw-sites    "${RAW_MERGED}" \
    --dedup-sites  "${DEDUP_MERGED}" \
    --gff          "${GFF}" \
    --sample       "${SAMPLE}" \
    --outfile      "${ANNOTATED}" \
    --per-gene-out "${GENE_COUNTS}" \
    --stats-out    "${STATS_REPORT}"

echo "  Annotated sites : ${ANNOTATED}"
echo "  Per-gene counts : ${GENE_COUNTS}"
echo "  Summary stats   : ${STATS_REPORT}"

# ---------------------------------------------------------------------------
# STEP 7: QC REPORT
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 7] Generating QC report..."

QC_OUT="$OUTDIR/qc"

python3 "${SCRIPT_DIR}/qc_report.py" \
    --sample         "${SAMPLE}" \
    --cutadapt-json  "${TRIM_REPORT}" \
    --bowtie2-stats  "${ALIGN_STATS}" \
    --samtools-stats "${SAMTOOLS_STATS}" \
    --dedup-stats    "${DEDUP_STATS}" \
    --raw-sites      "${RAW_MERGED}" \
    --dedup-sites    "${DEDUP_MERGED}" \
    --summary-stats  "${STATS_REPORT}" \
    --outdir         "${QC_OUT}"

echo "  QC report: ${QC_OUT}/${SAMPLE}_qc_report.pdf"
echo "  QC table : ${QC_OUT}/${SAMPLE}_qc_summary.tsv"

# ---------------------------------------------------------------------------
# DONE
# ---------------------------------------------------------------------------
echo ""
echo "======================================================"
echo " Pipeline complete: ${SAMPLE}"
echo " Finished: $(date)"
echo "======================================================"
echo ""
echo " Key outputs:"
printf "   %-38s %s\n" "All-reads BAM:"          "${SORTED_BAM}"
printf "   %-38s %s\n" "Dedup BAM:"               "${DEDUP_BAM}"
printf "   %-38s %s\n" "Raw insertion sites BED:" "${RAW_MERGED}"
printf "   %-38s %s\n" "Dedup insertion sites BED:" "${DEDUP_MERGED}"
printf "   %-38s %s\n" "+ strand bedGraph:"       "${BG_PLUS}"
printf "   %-38s %s\n" "- strand bedGraph:"       "${BG_MINUS}"
printf "   %-38s %s\n" "Annotated sites TSV:"     "${ANNOTATED}"
printf "   %-38s %s\n" "Per-gene counts TSV:"     "${GENE_COUNTS}"
printf "   %-38s %s\n" "QC report PDF:"           "${QC_OUT}/${SAMPLE}_qc_report.pdf"
echo "======================================================"
