#!/usr/bin/env python3
"""
make_excel_and_plots.py
-----------------------
Reads pipeline output for all samples and produces:
  1. insertion_sites.xlsx  — multi-sheet Excel workbook
       Sheet "Summary"            : per-sample QC metrics side-by-side
       Sheet "lib023_sites"       : all dedup insertion sites with annotation
       Sheet "lib123_sites"       : all dedup insertion sites with annotation
       Sheet "lib023_per_gene"    : per-gene counts for lib023
       Sheet "lib123_per_gene"    : per-gene counts for lib123
       Sheet "combined_genes"     : union of both samples, side-by-side read counts

  2. genome_insertions_lib023.pdf  — genome-wide insertion map for lib023
  3. genome_insertions_lib123.pdf  — genome-wide insertion map for lib123
  4. genome_insertions_combined.pdf — both samples overlaid

Usage:
    python3 make_excel_and_plots.py \
        --run-dir /path/to/20260506_run2 \
        --samples lib023 lib123 \
        --outdir  /path/to/output

    # Or with explicit per-sample paths:
    python3 make_excel_and_plots.py \
        --run-dir /home/akn27/project/Fulk2026_revisions/NGS/1_analysis/beakr_pipeline/20260506_run2 \
        --samples lib023 lib123 \
        --outdir  ./results

Requires: pandas, matplotlib, seaborn, openpyxl
    pip install openpyxl   (if not already installed)
"""

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# STYLE
# ---------------------------------------------------------------------------
sns.set_theme(style="whitegrid", font_scale=1.1)
COLORS = {
    "lib023":      "#3498db",   # blue
    "lib123":      "#e74c3c",   # red
    "intragenic":  "#2ecc71",
    "intergenic":  "#95a5a6",
    "chromosome":  "#2c3e50",
    "megaplasmid": "#8e44ad",
}

# C. necator H16 replicon sizes (bp) — used for axis scaling
# Chromosome NC_008313.2 and megaplasmid pHG1 NC_005241.1
REPLICON_SIZES = {
    "NC_008313.2": 7416678,   # chromosome
    "NC_005241.1": 452156,    # megaplasmid pHG1
}
# Accept common alternate accession formats
REPLICON_ALIASES = {
    "chromosome": "NC_008313.2",
    "chr":        "NC_008313.2",
    "phg1":       "NC_005241.1",
    "megaplasmid":"NC_005241.1",
}
REPLICON_LABELS = {
    "NC_008313.2": "Chromosome (7.4 Mb)",
    "NC_005241.1": "Megaplasmid pHG1 (452 kb)",
}


# ---------------------------------------------------------------------------
# LOADERS
# ---------------------------------------------------------------------------

def load_dedup_sites(run_dir, sample):
    """Load the dedup insertion sites BED for a sample."""
    path = Path(run_dir) / sample / "sites" / f"{sample}_dedup_insertion_sites.bed"
    if not path.exists():
        sys.exit(f"ERROR: File not found: {path}")
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "name", "dedup_count", "strand"])
    df["sample"] = sample
    return df


def load_annotated_sites(run_dir, sample):
    """Load the annotated insertion sites TSV for a sample."""
    path = Path(run_dir) / sample / "counts" / f"{sample}_insertions_annotated.tsv"
    if not path.exists():
        sys.exit(f"ERROR: File not found: {path}")
    df = pd.read_csv(path, sep="\t")
    df["sample"] = sample
    return df


def load_per_gene(run_dir, sample):
    """Load the per-gene count table for a sample."""
    path = Path(run_dir) / sample / "counts" / f"{sample}_per_gene_counts.tsv"
    if not path.exists():
        sys.exit(f"ERROR: File not found: {path}")
    df = pd.read_csv(path, sep="\t")
    df["sample"] = sample
    return df


def load_summary_stats(run_dir, sample):
    """Load summary stats TSV and return as a dict."""
    path = Path(run_dir) / sample / "stats" / f"{sample}_summary_stats.tsv"
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["metric"], df["value"]))


# ---------------------------------------------------------------------------
# EXCEL BUILDER
# ---------------------------------------------------------------------------

def build_excel(run_dir, samples, outdir):
    """Build the multi-sheet Excel workbook."""
    outpath = Path(outdir) / "insertion_sites.xlsx"

    # Load all data
    annotated  = {s: load_annotated_sites(run_dir, s) for s in samples}
    per_gene   = {s: load_per_gene(run_dir, s)        for s in samples}
    stats      = {s: load_summary_stats(run_dir, s)   for s in samples}

    with pd.ExcelWriter(outpath, engine="openpyxl") as writer:

        # ------------------------------------------------------------------
        # Sheet 1: Summary
        # ------------------------------------------------------------------
        summary_rows = [
            "sample",
            "total_raw_reads_at_sites",
            "total_dedup_reads_at_sites",
            "duplication_pct_at_sites",
            "total_unique_insertion_sites",
            "intragenic_sites",
            "intergenic_sites",
            "fraction_intragenic",
            "genes_with_insertions",
            "total_annotated_genes",
            "fraction_genes_hit",
            "median_raw_reads_per_site",
            "max_raw_reads_per_site",
            "median_dedup_reads_per_site",
            "max_dedup_reads_per_site",
        ]
        summary_data = {"Metric": summary_rows}
        for s in samples:
            summary_data[s] = [stats[s].get(r, "") for r in summary_rows]

        summary_df = pd.DataFrame(summary_data)
        summary_df.to_excel(writer, sheet_name="Summary", index=False)
        _format_sheet(writer, "Summary", summary_df, col_widths=[35] + [18]*len(samples))

        # ------------------------------------------------------------------
        # Sheets 2-3: Per-sample annotated insertion sites
        # ------------------------------------------------------------------
        site_cols = [
            "chrom", "ins_start", "ins_end", "site_strand",
            "raw_read_count", "dedup_read_count",
            "feature_type", "gene_id", "gene_name", "locus_tag",
            "gene_start", "gene_end", "gene_strand",
            "relative_pos_in_gene"
        ]
        for s in samples:
            df = annotated[s][[c for c in site_cols if c in annotated[s].columns]].copy()
            df = df.sort_values(["chrom", "ins_start"])
            df.to_excel(writer, sheet_name=f"{s}_sites", index=False)
            _format_sheet(writer, f"{s}_sites", df,
                          col_widths=[14,10,10,8,12,14,12,25,20,15,10,10,10,18])

        # ------------------------------------------------------------------
        # Sheets 4-5: Per-gene counts
        # ------------------------------------------------------------------
        gene_cols = [
            "gene_id", "gene_name", "locus_tag", "chrom",
            "gene_start", "gene_end", "gene_strand", "gene_length_bp",
            "n_unique_sites_raw", "total_raw_reads",
            "n_unique_sites_dedup", "total_dedup_reads",
            "unique_sites_per_kb_raw", "unique_sites_per_kb_dedup",
            "mean_relative_pos"
        ]
        for s in samples:
            df = per_gene[s][[c for c in gene_cols if c in per_gene[s].columns]].copy()
            df = df.sort_values("total_dedup_reads", ascending=False)
            df.to_excel(writer, sheet_name=f"{s}_per_gene", index=False)
            _format_sheet(writer, f"{s}_per_gene", df,
                          col_widths=[25,20,15,14,10,10,8,12,14,12,14,14,16,18,16])

        # ------------------------------------------------------------------
        # Sheet 6: Combined gene table (both samples side by side)
        # ------------------------------------------------------------------
        if len(samples) >= 2:
            s0, s1 = samples[0], samples[1]
            g0 = per_gene[s0][["gene_id","gene_name","locus_tag","chrom",
                                "gene_start","gene_end","gene_strand","gene_length_bp",
                                "n_unique_sites_dedup","total_dedup_reads"]].copy()
            g1 = per_gene[s1][["gene_id",
                                "n_unique_sites_dedup","total_dedup_reads"]].copy()

            g0 = g0.rename(columns={
                "n_unique_sites_dedup": f"unique_sites_{s0}",
                "total_dedup_reads":    f"dedup_reads_{s0}",
            })
            g1 = g1.rename(columns={
                "n_unique_sites_dedup": f"unique_sites_{s1}",
                "total_dedup_reads":    f"dedup_reads_{s1}",
            })

            combined = g0.merge(g1, on="gene_id", how="outer")
            combined = combined.fillna(0)
            combined[f"dedup_reads_{s0}"] = combined[f"dedup_reads_{s0}"].astype(int)
            combined[f"dedup_reads_{s1}"] = combined[f"dedup_reads_{s1}"].astype(int)
            combined["total_dedup_reads_combined"] = (
                combined[f"dedup_reads_{s0}"] + combined[f"dedup_reads_{s1}"]
            )
            combined = combined.sort_values("total_dedup_reads_combined", ascending=False)

            combined.to_excel(writer, sheet_name="combined_genes", index=False)
            _format_sheet(writer, "combined_genes", combined,
                          col_widths=[25,20,15,14,10,10,8,12,14,12,14,12,18])

    print(f"  Excel workbook: {outpath}")
    return outpath


def _format_sheet(writer, sheet_name, df, col_widths=None):
    """Apply column widths and freeze the header row."""
    ws = writer.sheets[sheet_name]
    # Freeze header row
    ws.freeze_panes = "A2"
    # Set column widths
    if col_widths:
        for i, width in enumerate(col_widths):
            col_letter = chr(65 + i) if i < 26 else "A" + chr(65 + i - 26)
            ws.column_dimensions[col_letter].width = width
    else:
        for i, col in enumerate(df.columns):
            col_letter = chr(65 + i) if i < 26 else "A" + chr(65 + i - 26)
            ws.column_dimensions[col_letter].width = max(len(str(col)) + 2, 10)


# ---------------------------------------------------------------------------
# GENOME VISUALIZATION
# ---------------------------------------------------------------------------

def get_replicon_order(all_chroms):
    """
    Return replicons in display order: chromosome first, megaplasmid second.
    Handles whatever accession format is in the data.
    """
    order = []
    # Prefer known accessions first
    for acc in ["NC_008313.2", "NC_005241.1"]:
        if acc in all_chroms:
            order.append(acc)
    # Add any remaining chroms not already included
    for c in sorted(all_chroms):
        if c not in order:
            order.append(c)
    return order


def get_replicon_size(chrom):
    """Return known size or estimate from data."""
    return REPLICON_SIZES.get(chrom, None)


def get_replicon_label(chrom):
    return REPLICON_LABELS.get(chrom, chrom)


def plot_genome_insertions(sites_dict, title, outpath, show_strand=True):
    """
    sites_dict: {sample_name: annotated_df}
    Produces a multi-panel figure: one row per replicon.
    Each panel shows insertion positions as vertical lines colored by sample.
    Top track: + strand, bottom track: - strand (mirrored).
    """
    # Determine replicons present
    all_chroms = set()
    for df in sites_dict.values():
        all_chroms.update(df["chrom"].unique())
    replicons = get_replicon_order(all_chroms)

    n_replicons = len(replicons)
    fig_height = 3.5 * n_replicons + 1.5
    fig, axes = plt.subplots(n_replicons, 1,
                             figsize=(18, fig_height),
                             squeeze=False)

    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.01)

    samples = list(sites_dict.keys())

    for row_idx, chrom in enumerate(replicons):
        ax = axes[row_idx, 0]
        rep_size = get_replicon_size(chrom)
        rep_label = get_replicon_label(chrom)

        # Determine axis range
        max_pos = rep_size if rep_size else 0
        for df in sites_dict.values():
            sub = df[df["chrom"] == chrom]
            if not sub.empty:
                max_pos = max(max_pos, sub["ins_start"].max())

        for s_idx, sample in enumerate(samples):
            df = sites_dict[sample]
            sub = df[df["chrom"] == chrom].copy()
            if sub.empty:
                continue

            color = COLORS.get(sample, f"C{s_idx}")
            alpha = 0.6 if len(samples) > 1 else 0.7

            # Separate strands
            plus  = sub[sub["site_strand"] == "+"]
            minus = sub[sub["site_strand"] == "-"]

            # Plot plus strand upward, minus strand downward
            # Use dedup_read_count for height (log-scaled) if available,
            # otherwise uniform height
            def _plot_strand(ax, pos_df, direction, color, alpha, label):
                if pos_df.empty:
                    return
                positions = pos_df["ins_start"].values

                if "dedup_read_count" in pos_df.columns:
                    counts = pos_df["dedup_read_count"].values.astype(float)
                    heights = np.log1p(counts) * direction
                    # Normalize so tallest bar = 1 (or -1)
                    max_h = np.abs(heights).max()
                    if max_h > 0:
                        heights = heights / max_h * 0.9
                else:
                    heights = np.full(len(positions), direction * 0.9)

                ax.vlines(positions, 0, heights,
                          color=color, alpha=alpha, linewidth=0.4,
                          label=label if direction > 0 else None)

            _plot_strand(ax, plus,  1, color, alpha,
                         f"{sample} (+ strand)")
            _plot_strand(ax, minus, -1, color, alpha * 0.7,
                         f"{sample} (- strand)")

        # Chromosome backbone
        ax.axhline(0, color="black", linewidth=1.2, zorder=5)

        # Formatting
        ax.set_xlim(0, max_pos)
        ax.set_ylim(-1.15, 1.15)
        ax.set_ylabel(rep_label, fontsize=10, fontweight="bold", labelpad=8)
        ax.set_yticks([1, 0, -1])
        ax.set_yticklabels(["+", "0", "−"], fontsize=9)

        # X axis in Mb or kb
        if max_pos > 1e6:
            ax.xaxis.set_major_formatter(
                ticker.FuncFormatter(lambda x, _: f"{x/1e6:.1f} Mb"))
        else:
            ax.xaxis.set_major_formatter(
                ticker.FuncFormatter(lambda x, _: f"{x/1e3:.0f} kb"))

        ax.set_xlabel("Genomic position", fontsize=9)

        # Replicon size annotation
        if rep_size:
            ax.text(0.99, 0.95, f"{rep_size/1e6:.3f} Mb",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=8, color="gray")

        # Strand labels
        ax.text(-0.01, 0.75, "+ strand", transform=ax.transAxes,
                ha="right", va="center", fontsize=8, color="gray")
        ax.text(-0.01, 0.25, "− strand", transform=ax.transAxes,
                ha="right", va="center", fontsize=8, color="gray")

        sns.despine(ax=ax, left=False)

    # Legend
    handles = []
    for s_idx, sample in enumerate(samples):
        color = COLORS.get(sample, f"C{s_idx}")
        handles.append(mpatches.Patch(color=color, label=sample))
    handles.append(mpatches.Patch(color="white", label="Height = log(dedup reads)"))
    fig.legend(handles=handles, loc="upper right",
               bbox_to_anchor=(1.0, 1.0), fontsize=10, framealpha=0.9)

    plt.tight_layout()
    fig.savefig(outpath, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Figure: {outpath}")


def plot_insertion_density(sites_dict, title, outpath, bin_size=10000):
    """
    Sliding-window insertion density plot.
    One panel per replicon, one line per sample.
    Y-axis: unique insertion sites per bin.
    """
    all_chroms = set()
    for df in sites_dict.values():
        all_chroms.update(df["chrom"].unique())
    replicons = get_replicon_order(all_chroms)

    n = len(replicons)
    fig, axes = plt.subplots(n, 1, figsize=(18, 3 * n + 1), squeeze=False)
    fig.suptitle(f"{title} — Insertion Density ({bin_size//1000} kb bins)",
                 fontsize=13, fontweight="bold", y=1.01)

    for row_idx, chrom in enumerate(replicons):
        ax = axes[row_idx, 0]
        rep_size = get_replicon_size(chrom)
        max_pos = rep_size if rep_size else 0
        for df in sites_dict.values():
            sub = df[df["chrom"] == chrom]
            if not sub.empty:
                max_pos = max(max_pos, sub["ins_start"].max())

        bins = np.arange(0, max_pos + bin_size, bin_size)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        for s_idx, (sample, df) in enumerate(sites_dict.items()):
            sub = df[df["chrom"] == chrom]
            if sub.empty:
                continue
            counts, _ = np.histogram(sub["ins_start"].values, bins=bins)
            color = COLORS.get(sample, f"C{s_idx}")
            ax.fill_between(bin_centers, counts, alpha=0.4, color=color)
            ax.plot(bin_centers, counts, color=color, linewidth=1.2, label=sample)

        ax.set_xlim(0, max_pos)
        ax.set_ylabel("Sites per bin", fontsize=9)
        ax.set_title(get_replicon_label(chrom), fontsize=10, fontweight="bold", pad=4)

        if max_pos > 1e6:
            ax.xaxis.set_major_formatter(
                ticker.FuncFormatter(lambda x, _: f"{x/1e6:.1f} Mb"))
        else:
            ax.xaxis.set_major_formatter(
                ticker.FuncFormatter(lambda x, _: f"{x/1e3:.0f} kb"))

        ax.set_xlabel("Genomic position", fontsize=9)
        ax.legend(fontsize=9, loc="upper right")
        sns.despine(ax=ax)

    plt.tight_layout()
    fig.savefig(outpath, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Density figure: {outpath}")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build Excel workbook and genome-wide insertion plots from pipeline output"
    )
    parser.add_argument("--run-dir", required=True,
                        help="Path to the pipeline run directory (contains per-sample subdirs)")
    parser.add_argument("--samples", nargs="+", required=True,
                        help="Sample names (must match subdirectory names)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for Excel and figures")
    parser.add_argument("--bin-size", type=int, default=10000,
                        help="Bin size in bp for density plot [default: 10000]")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"\nLoading data for samples: {', '.join(args.samples)}")

    # Load annotated sites for all samples
    annotated = {}
    for s in args.samples:
        print(f"  Loading {s}...")
        annotated[s] = load_annotated_sites(args.run_dir, s)
        n = len(annotated[s])
        chroms = annotated[s]["chrom"].unique().tolist()
        print(f"    {n} sites across {len(chroms)} replicons: {chroms}")

    # ------------------------------------------------------------------
    # Excel workbook
    # ------------------------------------------------------------------
    print("\nBuilding Excel workbook...")
    build_excel(args.run_dir, args.samples, args.outdir)

    # ------------------------------------------------------------------
    # Per-sample insertion maps
    # ------------------------------------------------------------------
    print("\nGenerating genome-wide insertion maps...")
    for s in args.samples:
        plot_genome_insertions(
            {s: annotated[s]},
            title=f"Tn5 Insertion Sites — {s}",
            outpath=os.path.join(args.outdir, f"genome_insertions_{s}.pdf"),
        )
        plot_insertion_density(
            {s: annotated[s]},
            title=f"{s}",
            outpath=os.path.join(args.outdir, f"insertion_density_{s}.pdf"),
            bin_size=args.bin_size,
        )

    # ------------------------------------------------------------------
    # Combined figure (all samples overlaid)
    # ------------------------------------------------------------------
    if len(args.samples) > 1:
        print("\nGenerating combined insertion map...")
        plot_genome_insertions(
            annotated,
            title="Tn5 Insertion Sites — All Samples",
            outpath=os.path.join(args.outdir, "genome_insertions_combined.pdf"),
        )
        plot_insertion_density(
            annotated,
            title="All Samples",
            outpath=os.path.join(args.outdir, "insertion_density_combined.pdf"),
            bin_size=args.bin_size,
        )

    print(f"\nAll outputs written to: {args.outdir}")
    print("\nFiles:")
    for f in sorted(os.listdir(args.outdir)):
        size = os.path.getsize(os.path.join(args.outdir, f))
        print(f"  {f:50s}  {size/1024:.0f} KB")


if __name__ == "__main__":
    main()
