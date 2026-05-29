#!/usr/bin/env bash
# =============================================================================
# SLURM Array Job — Tn5 Insertion Mapping Pipeline
#
# Submits one job per sample using a tab-delimited sample sheet.
# Each job runs the full pipeline: trimming → alignment → dedup →
# insertion site extraction → annotation → QC report.
#
# SAMPLE SHEET FORMAT (no header, tab-separated, one sample per line):
#   SAMPLE_NAME <TAB> /abs/path/to/R1.fastq.gz <TAB> /abs/path/to/R2.fastq.gz
#
#   Example (samples.tsv):
#   lib1    /data/fastq/lib1_R1.fastq.gz    /data/fastq/lib1_R2.fastq.gz
#   lib2    /data/fastq/lib2_R1.fastq.gz    /data/fastq/lib2_R2.fastq.gz
#   lib3    /data/fastq/lib3_R1.fastq.gz    /data/fastq/lib3_R2.fastq.gz
#
# SETUP (one time):
#   1. Run setup_reference.sh to download the C. necator H16 reference
#      and build the Bowtie2 index.
#   2. Edit the USER CONFIGURATION section below.
#   3. Set --array=1-N where N = number of lines in your sample sheet.
#
# SUBMIT:
#   sbatch tn5_slurm_array.sh
#
# MONITOR:
#   squeue -u $USER
#   sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS,CPUTime
#
# AGGREGATE QC ACROSS SAMPLES (after all jobs finish):
#   python3 aggregate_qc.py --qc-dir /path/to/output --outdir /path/to/output/multiqc
# =============================================================================

#SBATCH --job-name=tn5_map
#SBATCH --output=slurm_logs/tn5_%A_%a.out      # %A = array job ID, %a = task index
#SBATCH --error=slurm_logs/tn5_%A_%a.err
#SBATCH --array=1-2                             # <<<< SET TO: 1-N  (N = number of samples)
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00                         # PE150 ~10–50M reads typically finishes in 1–2h
#SBATCH --partition=day                    # <<<< SET TO your cluster partition name
## Uncomment for email notification
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=akira.nakamura@yale.edu

set -euo pipefail

# =============================================================================
# USER CONFIGURATION — edit all paths marked <<<< before submitting
# =============================================================================

# Directory containing the pipeline scripts
# (tn5_insertion_mapping.sh, annotate_insertions.py, qc_report.py)
SCRIPT_DIR="/home/akn27/project/Fulk2026_revisions/NGS/1_analysis/beakr_pipeline/scripts"                        # <<<< EDIT

# Bowtie2 index prefix for C. necator H16
# Built by setup_reference.sh as: <REFS_DIR>/cnecator_H16
BOWTIE2_INDEX="/home/akn27/project/Fulk2026_revisions/NGS/1_analysis/beakr_pipeline/reference/cnecator_H16"

# GFF3 annotation for C. necator H16
# Downloaded by setup_reference.sh as: <REFS_DIR>/cnecator_H16_genomic.gff
GFF="/home/akn27/project/Fulk2026_revisions/NGS/1_analysis/beakr_pipeline/reference/cnecator_H16_genomic.gff"

# Sample sheet: tab-separated, no header, columns: SAMPLE  R1  R2
SAMPLE_SHEET="/home/akn27/project/Fulk2026_revisions/NGS/1_analysis/beakr_pipeline/samples.tsv"                            # <<<< EDIT

# Root output directory (per-sample subdirectories are created automatically)
OUTDIR_ROOT="/home/akn27/project/Fulk2026_revisions/NGS/1_analysis/beakr_pipeline/20260506_run2"                                  # <<<< EDIT

# =============================================================================
# MODULE LOADING
# =============================================================================
# Edit module names to match your cluster's module system.
# To use conda instead, comment out the module lines and uncomment the conda block.

#module purge
module load miniconda
#module load cutadapt/3.7
#module load bowtie2/2.5.1
#module load samtools/1.17
#module load bedtools/2.31.0
#module load python/3.10

# --- Conda alternative ---
# source /path/to/conda/etc/profile.d/conda.sh
conda activate tn5_mapping
# -------------------------

# =============================================================================
# PRE-FLIGHT CHECKS
# =============================================================================

mkdir -p slurm_logs

# Verify all required tools are available
echo "[preflight] Checking tool availability..."
for TOOL in cutadapt bowtie2 samtools bedtools python3; do
    if ! command -v "$TOOL" &>/dev/null; then
        echo "ERROR: $TOOL not found on PATH. Check module load lines."
        exit 1
    fi
    echo "  OK: $TOOL ($(command -v $TOOL))"
done

# Verify Python dependencies
for PKG in pandas matplotlib seaborn; do
    python3 -c "import ${PKG}" 2>/dev/null \
        || { echo "ERROR: Python package '${PKG}' not available."; exit 1; }
    echo "  OK: python3 ${PKG}"
done

# Verify required files exist
for F in "${SCRIPT_DIR}/tn5_insertion_mapping.sh" \
         "${SCRIPT_DIR}/annotate_insertions.py" \
         "${SCRIPT_DIR}/qc_report.py" \
         "${SAMPLE_SHEET}" \
         "${GFF}"; do
    [[ -f "$F" ]] || { echo "ERROR: Required file not found: $F"; exit 1; }
done

# Verify Bowtie2 index
[[ -f "${BOWTIE2_INDEX}.1.bt2" || -f "${BOWTIE2_INDEX}.1.bt2l" ]] \
    || { echo "ERROR: Bowtie2 index not found at ${BOWTIE2_INDEX}. Run setup_reference.sh first."; exit 1; }

echo "[preflight] All checks passed."

# =============================================================================
# RESOLVE SAMPLE FOR THIS ARRAY TASK
# =============================================================================

SAMPLE=$(awk -v idx="${SLURM_ARRAY_TASK_ID}" 'NR==idx {print $1}' "${SAMPLE_SHEET}")
R1=$(    awk -v idx="${SLURM_ARRAY_TASK_ID}" 'NR==idx {print $2}' "${SAMPLE_SHEET}")
R2=$(    awk -v idx="${SLURM_ARRAY_TASK_ID}" 'NR==idx {print $3}' "${SAMPLE_SHEET}")

if [[ -z "${SAMPLE}" || -z "${R1}" || -z "${R2}" ]]; then
    echo "ERROR: Could not parse sample at array index ${SLURM_ARRAY_TASK_ID} from ${SAMPLE_SHEET}"
    echo "  Check that --array=1-N matches the number of lines in your sample sheet."
    exit 1
fi

[[ -f "${R1}" ]] || { echo "ERROR: R1 file not found: ${R1}"; exit 1; }
[[ -f "${R2}" ]] || { echo "ERROR: R2 file not found: ${R2}"; exit 1; }

echo ""
echo "======================================================"
echo " SLURM Array Task"
echo " Array job ID : ${SLURM_ARRAY_JOB_ID}"
echo " Array index  : ${SLURM_ARRAY_TASK_ID}"
echo " Sample       : ${SAMPLE}"
echo " R1           : ${R1}"
echo " R2           : ${R2}"
echo " Node         : $(hostname)"
echo " CPUs         : ${SLURM_CPUS_PER_TASK}"
echo " Memory       : ${SLURM_MEM_PER_NODE:-not set} MB"
echo " Started      : $(date)"
echo "======================================================"

# =============================================================================
# RUN PIPELINE
# =============================================================================

SAMPLE_OUTDIR="${OUTDIR_ROOT}/${SAMPLE}"
mkdir -p "${SAMPLE_OUTDIR}"

bash "${SCRIPT_DIR}/tn5_insertion_mapping.sh" \
    -1 "${R1}" \
    -2 "${R2}" \
    -s "${SAMPLE}" \
    -r "${BOWTIE2_INDEX}" \
    -g "${GFF}" \
    -o "${SAMPLE_OUTDIR}" \
    -t "${SLURM_CPUS_PER_TASK}"

echo ""
echo "Job finished: $(date)"
echo "Output: ${SAMPLE_OUTDIR}"
