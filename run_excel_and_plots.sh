#!/usr/bin/env bash
# =============================================================================
# run_excel_and_plots.sh
# ----------------------
# Generate the Excel workbook and genome-wide insertion figures from pipeline
# output using make_excel_and_plots.py.
#
# Can be run interactively on a login node (fast, <5 min) or submitted to
# SLURM for large numbers of samples.
#
# Usage (interactive):
#   conda activate tn5_mapping
#   bash run_excel_and_plots.sh
#
# Usage (SLURM):
#   sbatch run_excel_and_plots.sh
# =============================================================================

#SBATCH --job-name=tn5_plots
#SBATCH --output=slurm_logs/tn5_plots_%j.out
#SBATCH --error=slurm_logs/tn5_plots_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --partition=day          # <<<< SET to your cluster partition
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=akira.nakamura@yale.edu

set -euo pipefail

# =============================================================================
# USER CONFIGURATION — edit these before running
# =============================================================================

# Directory containing the pipeline scripts (make_excel_and_plots.py lives here)
SCRIPT_DIR="/home/akn27/project/Fulk2026_revisions/NGS/1_analysis/beakr_pipeline/scripts"                  # <<<< EDIT

# Pipeline run output directory (contains per-sample subdirectories)
RUN_DIR="/home/akn27/project/Fulk2026_revisions/NGS/1_analysis/beakr_pipeline/20260506_run2"               # <<<< EDIT

# Sample names — must match subdirectory names inside RUN_DIR
SAMPLES="lib023 lib123"                                        # <<<< EDIT (space-separated)

# Output directory for Excel file and figures
OUTDIR="/home/akn27/project/Fulk2026_revisions/NGS/1_analysis/beakr_pipeline/20260506_run2/insertion_sites"                      # <<<< EDIT

# Bin size in bp for the density plot (default: 10000 = 10 kb)
BIN_SIZE=5000

# =============================================================================
# CONDA / MODULE SETUP
# =============================================================================
# Option A: conda (recommended — matches the tn5_mapping environment)
module load miniconda
conda activate tn5_mapping

which python3
python3 -c "import sys; print(sys.executable)"

# Option B: modules — uncomment and edit if using cluster modules instead
# module purge
# module load python/3.10

# =============================================================================
# DEPENDENCY CHECK
# =============================================================================
echo "[check] Verifying Python dependencies..."

# Auto-install openpyxl if missing
python3 -c "import openpyxl" 2>/dev/null || {
    echo "[setup] openpyxl not found — installing via pip..."
    pip install openpyxl --quiet
    echo "[setup] openpyxl installed."
}

# Final check — all packages must be present
python3 -c "import pandas, matplotlib, seaborn, openpyxl" 2>/dev/null || {
    echo "ERROR: One or more required Python packages are missing."
    echo "  Try: pip install pandas matplotlib seaborn openpyxl"
    exit 1
}
echo "[check] All dependencies OK."

# =============================================================================
# VALIDATE INPUTS
# =============================================================================
[[ -f "${SCRIPT_DIR}/make_excel_and_plots.py" ]] || {
    echo "ERROR: make_excel_and_plots.py not found in ${SCRIPT_DIR}"
    exit 1
}

[[ -d "${RUN_DIR}" ]] || {
    echo "ERROR: Run directory not found: ${RUN_DIR}"
    exit 1
}

for SAMPLE in ${SAMPLES}; do
    SITES="${RUN_DIR}/${SAMPLE}/counts/${SAMPLE}_insertions_annotated.tsv"
    [[ -f "${SITES}" ]] || {
        echo "ERROR: Annotated sites file not found: ${SITES}"
        echo "  Make sure the pipeline completed successfully for sample ${SAMPLE}"
        exit 1
    }
    echo "[check] Found: ${SITES}"
done

mkdir -p "${OUTDIR}"
mkdir -p slurm_logs

# =============================================================================
# RUN
# =============================================================================
echo ""
echo "======================================================"
echo " Generating Excel workbook and insertion figures"
echo " Run dir : ${RUN_DIR}"
echo " Samples : ${SAMPLES}"
echo " Output  : ${OUTDIR}"
echo " Started : $(date)"
echo "======================================================"

python3 "${SCRIPT_DIR}/make_excel_and_plots.py" \
    --run-dir  "${RUN_DIR}" \
    --samples  ${SAMPLES} \
    --outdir   "${OUTDIR}" \
    --bin-size ${BIN_SIZE}

echo ""
echo "======================================================"
echo " Done: $(date)"
echo " Outputs:"
for F in "${OUTDIR}"/*; do
    SIZE=$(du -sh "$F" 2>/dev/null | cut -f1)
    printf "   %-50s %s\n" "$(basename $F)" "${SIZE}"
done
echo "======================================================"

# =============================================================================
# COPY TO LOCAL (optional — uncomment and edit to auto-scp results)
# =============================================================================
# LOCAL_DEST="you@yourmachine.local:/path/to/local/results/"
# scp -r "${OUTDIR}"/* "${LOCAL_DEST}"
