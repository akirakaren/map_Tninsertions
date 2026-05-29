#!/usr/bin/env bash
# =============================================================================
# One-time reference setup for C. necator H16 (GCF_000009285.1)
# Run this ONCE before submitting the SLURM array job.
#
# Downloads the reference FASTA + GFF3 from NCBI and builds the Bowtie2 index.
# =============================================================================

set -euo pipefail

REFS_DIR="/home/akn27/project/Fulk2026_revisions/NGS/1_analysis/20260506_beakr_analysis_run/reference"   # <<<< EDIT: where to store reference files
THREADS=8

mkdir -p "${REFS_DIR}"
cd "${REFS_DIR}"

echo "[1/3] Downloading C. necator H16 genome (FASTA + GFF3)..."

# NCBI Datasets CLI (recommended — install from: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
# If datasets CLI is not available, use the wget fallback below.

if command -v datasets &>/dev/null; then
    datasets download genome accession GCF_000009285.1 \
        --include genome,gff3 \
        --filename cnecator_H16_ncbi.zip
    unzip -o cnecator_H16_ncbi.zip
    # Move files to expected locations
    find ncbi_dataset/data/ -name "*.fna"  -exec cp {} cnecator_H16_genomic.fna \;
    find ncbi_dataset/data/ -name "*.gff"  -exec cp {} cnecator_H16_genomic.gff \;
else
    echo "  datasets CLI not found — using direct FTP download..."
    BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/285/GCF_000009285.1_ASM928v2"
    wget -nc "${BASE_URL}/GCF_000009285.1_ASM928v2_genomic.fna.gz"
    wget -nc "${BASE_URL}/GCF_000009285.1_ASM928v2_genomic.gff.gz"
    gunzip -k GCF_000009285.1_ASM928v2_genomic.fna.gz
    gunzip -k GCF_000009285.1_ASM928v2_genomic.gff.gz
    mv GCF_000009285.1_ASM928v2_genomic.fna cnecator_H16_genomic.fna
    mv GCF_000009285.1_ASM928v2_genomic.gff cnecator_H16_genomic.gff
fi

echo "[2/3] Building Bowtie2 index (chromosome + pHG1 megaplasmid included)..."
bowtie2-build --threads ${THREADS} cnecator_H16_genomic.fna cnecator_H16

echo "[3/3] Indexing FASTA with samtools fai..."
samtools faidx cnecator_H16_genomic.fna

echo ""
echo "Reference setup complete."
echo "  FASTA   : ${REFS_DIR}/cnecator_H16_genomic.fna"
echo "  GFF3    : ${REFS_DIR}/cnecator_H16_genomic.gff"
echo "  BT2 idx : ${REFS_DIR}/cnecator_H16"
echo ""
echo "Update these paths in tn5_slurm_array.sh:"
echo "  BOWTIE2_INDEX=${REFS_DIR}/cnecator_H16"
echo "  GFF=${REFS_DIR}/cnecator_H16_genomic.gff"
