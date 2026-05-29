#!/usr/bin/env python3
"""
qc_report.py
------------
Generate a QC report for the Tn5 insertion mapping pipeline.

Produces:
  - SAMPLE_qc_report.pdf   : multi-page PDF with all QC plots
  - SAMPLE_qc_summary.tsv  : flat TSV of all QC metrics (for multi-sample aggregation)

Plots included:
  1. Read fate funnel        — reads at each pipeline stage (input → trimmed → mapped → dedup)
  2. Trimming outcome        — fraction passing / discarded-untrimmed / too-short
  3. Insert size distribution — histogram of inferred insert sizes from samtools stats
  4. Mapping quality dist.   — MAPQ score histogram from aligned BAM
  5. Duplication rate        — raw vs dedup read counts and duplication %
  6. Insertion site read depth — histogram of reads-per-site (raw and dedup)
  7. Cumulative site coverage — what fraction of sites are covered by ≥N reads
  8. Genomic feature distribution — intragenic vs intergenic insertion fractions

Usage (called by tn5_insertion_mapping.sh Step 7):
    python3 qc_report.py \
        --sample         SAMPLE_NAME \
        --cutadapt-json  stats/SAMPLE_cutadapt.json \
        --bowtie2-stats  stats/SAMPLE_bowtie2_stats.txt \
        --samtools-stats stats/SAMPLE_samtools_stats.txt \
        --dedup-stats    stats/SAMPLE_dedup_stats.txt \
        --raw-sites      sites/SAMPLE_raw_insertion_sites.bed \
        --dedup-sites    sites/SAMPLE_dedup_insertion_sites.bed \
        --summary-stats  stats/SAMPLE_summary_stats.tsv \
        --outdir         qc/

Requires: python3, matplotlib >= 3.5, seaborn >= 0.11, pandas >= 1.3
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# STYLE
# ---------------------------------------------------------------------------
sns.set_theme(style="whitegrid", font_scale=1.1)
PALETTE = {
    "pass":       "#2ecc71",
    "fail":       "#e74c3c",
    "warn":       "#f39c12",
    "raw":        "#3498db",
    "dedup":      "#9b59b6",
    "intragenic": "#1abc9c",
    "intergenic": "#95a5a6",
    "neutral":    "#34495e",
}

# QC thresholds — flag values outside these ranges in the summary table
THRESHOLDS = {
    "pct_reads_trimmed_passing":  (50.0, 100.0),   # % of input reads passing trim
    "pct_mapped":                 (50.0, 100.0),
    "pct_properly_paired":        (50.0, 100.0),
    "duplication_rate_pct":       (0.0,  80.0),     # flag if >80% duplication
    "median_insert_size_bp":      (50.0, 900.0),    # flag if very short/long
    "pct_intragenic":             (20.0, 100.0),    # flag if <20% intragenic (poor enrichment)
}


# ---------------------------------------------------------------------------
# PARSERS
# ---------------------------------------------------------------------------

def parse_cutadapt_json(path):
    """Parse cutadapt --json output. Handles cutadapt 3.x and 5.x schema."""
    with open(path) as fh:
        d = json.load(fh)

    rc = d.get("read_counts", {})

    # cutadapt 5.x: rc["input"] is an int, not a dict
    # cutadapt 3.x: rc["input"] is a dict with a "reads" key
    def _get_reads(val):
        if isinstance(val, dict):
            return val.get("reads", 0)
        return int(val) if val else 0

    total_reads  = _get_reads(rc.get("input",  0))
    output_reads = _get_reads(rc.get("output", 0))
    total_pairs  = total_reads  // 2
    passing      = output_reads // 2

    filt = rc.get("filtered", {})
    def _get_filtered(key):
        val = filt.get(key, 0)
        return _get_reads(val) // 2

    too_short        = _get_filtered("too_short")
    discarded_untrim = _get_filtered("discard_untrimmed")

    adapters   = d.get("adapters_read1", [])
    trim_rate  = 0.0
    if adapters:
        trim_rate = adapters[0].get("total_matches", 0) / max(total_pairs, 1)

    return {
        "input_pairs":               total_pairs,
        "passing_pairs":             passing,
        "discarded_untrimmed":       discarded_untrim,
        "discarded_too_short":       too_short,
        "pct_reads_trimmed_passing": passing / max(total_pairs, 1) * 100,
        "r1_tn5_trim_rate":          trim_rate,
    }

def parse_cutadapt_txt(path):
    """
    Fallback: parse cutadapt plain-text output when JSON is unavailable.
    Returns the same keys as parse_cutadapt_json with best-effort values.
    """
    metrics = {
        "input_pairs": 0, "passing_pairs": 0,
        "discarded_untrimmed": 0, "discarded_too_short": 0,
        "pct_reads_trimmed_passing": 0.0, "r1_tn5_trim_rate": 0.0,
    }
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            m = re.search(r"Total read pairs processed:\s+([\d,]+)", line)
            if m:
                metrics["input_pairs"] = int(m.group(1).replace(",", ""))
            m = re.search(r"Pairs written \(passing filters\):\s+([\d,]+)", line)
            if m:
                metrics["passing_pairs"] = int(m.group(1).replace(",", ""))
            m = re.search(r"Pairs discarded as untrimmed:\s+([\d,]+)", line)
            if m:
                metrics["discarded_untrimmed"] = int(m.group(1).replace(",", ""))
            m = re.search(r"Pairs that were too short:\s+([\d,]+)", line)
            if m:
                metrics["discarded_too_short"] = int(m.group(1).replace(",", ""))

    total = metrics["input_pairs"]
    if total > 0:
        metrics["pct_reads_trimmed_passing"] = metrics["passing_pairs"] / total * 100
    return metrics


def parse_bowtie2_stats(path):
    """Parse bowtie2 stderr alignment stats."""
    metrics = {}
    with open(path) as fh:
        text = fh.read()

    m = re.search(r"(\d+) reads; of these:", text)
    if m:
        metrics["total_reads_bt2"] = int(m.group(1))

    m = re.search(r"(\d+) \(.+\) aligned concordantly exactly 1 time", text)
    if m:
        metrics["concordant_1"] = int(m.group(1))

    m = re.search(r"(\d+) \(.+\) aligned concordantly >1 times", text)
    if m:
        metrics["concordant_multi"] = int(m.group(1))

    m = re.search(r"(\d+) \(.+\) aligned concordantly 0 times", text)
    if m:
        metrics["concordant_0"] = int(m.group(1))

    m = re.search(r"(\d+\.\d+)% overall alignment rate", text)
    if m:
        metrics["pct_overall_aligned"] = float(m.group(1))

    total = metrics.get("total_reads_bt2", 1)
    mapped = metrics.get("concordant_1", 0) + metrics.get("concordant_multi", 0)
    metrics["pct_mapped"] = mapped / max(total, 1) * 100
    metrics["pct_multi"]  = metrics.get("concordant_multi", 0) / max(total, 1) * 100

    return metrics


def parse_samtools_stats(path):
    """
    Parse samtools stats output.
    Extracts insert size histogram (IS lines) and summary stats (SN lines).
    """
    sn = {}
    insert_sizes = []   # list of (size, count) tuples

    with open(path) as fh:
        for line in fh:
            if line.startswith("SN\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    key = parts[1].rstrip(":").strip()
                    try:
                        val = float(parts[2])
                    except ValueError:
                        val = parts[2]
                    sn[key] = val
            elif line.startswith("IS\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    try:
                        size  = int(parts[1])
                        count = int(parts[2])
                        if count > 0:
                            insert_sizes.append((size, count))
                    except ValueError:
                        pass

    insert_df = pd.DataFrame(insert_sizes, columns=["insert_size", "count"]) \
                if insert_sizes else pd.DataFrame(columns=["insert_size", "count"])

    # Compute median insert size from histogram
    median_insert = 0.0
    if not insert_df.empty:
        expanded = np.repeat(insert_df["insert_size"].values,
                             insert_df["count"].values.astype(int))
        if len(expanded) > 0:
            median_insert = float(np.median(expanded))

    return {
        "sn":                  sn,
        "insert_df":           insert_df,
        "median_insert_size":  median_insert,
        "mean_insert_size":    float(sn.get("insert size average", 0)),
        "pct_properly_paired": float(sn.get("percentage of properly paired reads (%)", 0)),
        "reads_mapped":        int(sn.get("reads mapped", 0)),
        "reads_total":         int(sn.get("raw total sequences", 0)),
    }


def parse_dedup_stats(path):
    """Parse coordinate-pair dedup stats written by the pipeline."""
    d = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                try:
                    d[parts[0]] = float(parts[1])
                except ValueError:
                    d[parts[0]] = parts[1]
    d["duplication_rate_pct"] = float(d.get("duplication_rate", 0)) * 100
    return d


def parse_sites_bed(path):
    """
    Read merged insertion sites BED.
    Col 5 = read count per site.
    Returns a DataFrame with columns: chrom, start, end, name, count, strand.
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=["chrom","start","end","name","count","strand"])
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom","start","end","name","count","strand"])
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(int)
    return df


def parse_summary_stats(path):
    """Parse the per-sample summary_stats.tsv produced by annotate_insertions.py."""
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["metric"], df["value"]))


# ---------------------------------------------------------------------------
# PLOTTING HELPERS
# ---------------------------------------------------------------------------

def add_value_labels(ax, fmt="{:.0f}", fontsize=9, color="white", ha="center", va="center"):
    """Add text labels to the center of each bar in a bar chart."""
    for patch in ax.patches:
        h = patch.get_height()
        if h > 0:
            ax.text(
                patch.get_x() + patch.get_width() / 2,
                patch.get_y() + h / 2,
                fmt.format(h),
                ha=ha, va=va, fontsize=fontsize, color=color, fontweight="bold"
            )


def flag_metric(val, key):
    """Return color string based on whether value is within threshold."""
    if key not in THRESHOLDS:
        return PALETTE["neutral"]
    lo, hi = THRESHOLDS[key]
    if lo <= float(val) <= hi:
        return PALETTE["pass"]
    return PALETTE["fail"]


# ---------------------------------------------------------------------------
# INDIVIDUAL PLOT FUNCTIONS
# ---------------------------------------------------------------------------

def plot_read_funnel(ax, cutadapt, bowtie2, dedup):
    """Bar chart: reads at each stage of the pipeline."""
    stages = [
        ("Input\npairs",    cutadapt.get("input_pairs",   0)),
        ("After\ntrimming", cutadapt.get("passing_pairs", 0)),
        ("After\nalignment",bowtie2.get("reads_mapped",   0) // 2),
        ("After\ndedup",    int(dedup.get("unique_coordinate_pairs", 0))),
    ]
    labels, values = zip(*stages)
    colors = [PALETTE["neutral"], PALETTE["pass"], PALETTE["raw"], PALETTE["dedup"]]
    bars = ax.bar(labels, values, color=colors, edgecolor="white", linewidth=0.8)
    ax.set_title("Read Fate Funnel", fontweight="bold")
    ax.set_ylabel("Read pairs")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda x, _: f"{int(x):,}"))
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02,
                f"{val:,}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_ylim(0, max(values) * 1.15)
    sns.despine(ax=ax)


def plot_trimming_outcome(ax, cutadapt):
    """Stacked bar: disposition of input read pairs after trimming."""
    total      = cutadapt.get("input_pairs",          1)
    passing    = cutadapt.get("passing_pairs",         0)
    no_adapter = cutadapt.get("discarded_untrimmed",   0)
    too_short  = cutadapt.get("discarded_too_short",   0)
    other      = max(0, total - passing - no_adapter - too_short)

    fractions = {
        "Passed":             passing    / total * 100,
        "No Tn5 ME found":    no_adapter / total * 100,
        "Too short":          too_short  / total * 100,
        "Other":              other      / total * 100,
    }
    colors_map = {
        "Passed":           PALETTE["pass"],
        "No Tn5 ME found":  PALETTE["fail"],
        "Too short":        PALETTE["warn"],
        "Other":            PALETTE["neutral"],
    }

    bottom = 0
    for label, pct in fractions.items():
        if pct > 0:
            ax.bar(0, pct, bottom=bottom, color=colors_map[label],
                   label=f"{label} ({pct:.1f}%)", width=0.5)
            if pct > 3:
                ax.text(0, bottom + pct/2, f"{pct:.1f}%",
                        ha="center", va="center", fontsize=9,
                        color="white", fontweight="bold")
            bottom += pct

    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0, 105)
    ax.set_xticks([])
    ax.set_ylabel("% of input read pairs")
    ax.set_title("Trimming Outcome", fontweight="bold")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.8)
    sns.despine(ax=ax)


def plot_insert_size(ax, samtools_stats):
    """Histogram of insert sizes from samtools stats IS table."""
    df = samtools_stats.get("insert_df", pd.DataFrame())
    if df.empty or df["count"].sum() == 0:
        ax.text(0.5, 0.5, "No insert size data", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="gray")
        ax.set_title("Insert Size Distribution", fontweight="bold")
        return

    median = samtools_stats.get("median_insert_size", 0)
    mean   = samtools_stats.get("mean_insert_size",   0)

    # Clip to sensible range for display
    df = df[df["insert_size"] <= 1000]

    ax.fill_between(df["insert_size"], df["count"],
                    alpha=0.6, color=PALETTE["raw"], step="mid")
    ax.plot(df["insert_size"], df["count"],
            color=PALETTE["raw"], linewidth=1.2)
    ax.axvline(median, color=PALETTE["fail"],  linestyle="--",
               linewidth=1.5, label=f"Median: {median:.0f} bp")
    ax.axvline(mean,   color=PALETTE["warn"],  linestyle=":",
               linewidth=1.5, label=f"Mean: {mean:.0f} bp")
    ax.set_title("Insert Size Distribution", fontweight="bold")
    ax.set_xlabel("Insert size (bp)")
    ax.set_ylabel("Read pairs")
    ax.legend(fontsize=9)
    sns.despine(ax=ax)


def plot_mapping_stats(ax, bowtie2):
    """Horizontal bar: mapping categories."""
    total = max(bowtie2.get("total_reads_bt2", 1), 1)
    uniq  = bowtie2.get("concordant_1",     0)
    multi = bowtie2.get("concordant_multi", 0)
    unmapped = bowtie2.get("concordant_0",  0)
    other = max(0, total - uniq - multi - unmapped)

    categories = [
        ("Unique concordant",  uniq     / total * 100, PALETTE["pass"]),
        ("Multi-mapping",      multi    / total * 100, PALETTE["warn"]),
        ("Unmapped",           unmapped / total * 100, PALETTE["fail"]),
        ("Other",              other    / total * 100, PALETTE["neutral"]),
    ]

    labels = [c[0] for c in categories]
    values = [c[1] for c in categories]
    colors = [c[2] for c in categories]

    bars = ax.barh(labels, values, color=colors, edgecolor="white")
    for bar, val in zip(bars, values):
        if val > 1:
            ax.text(val + 0.5, bar.get_y() + bar.get_height()/2,
                    f"{val:.1f}%", va="center", fontsize=9)
    ax.set_xlim(0, 110)
    ax.set_xlabel("% of reads")
    ax.set_title("Alignment Categories", fontweight="bold")
    sns.despine(ax=ax)


def plot_duplication(ax, dedup):
    """Side-by-side bar: total vs unique coordinate pairs."""
    total   = int(dedup.get("total_properly_paired",   0))
    unique  = int(dedup.get("unique_coordinate_pairs", 0))
    dup_pct = dedup.get("duplication_rate_pct",        0.0)

    bars = ax.bar(["All properly\npaired", "Unique coord\npairs (dedup)"],
                  [total, unique],
                  color=[PALETTE["raw"], PALETTE["dedup"]],
                  edgecolor="white", width=0.5)
    for bar, val in zip(bars, [total, unique]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02,
                f"{val:,}", ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_title(f"Duplication Rate: {dup_pct:.1f}%", fontweight="bold",
                 color=flag_metric(dup_pct, "duplication_rate_pct"))
    ax.set_ylabel("Read pairs")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda x, _: f"{int(x):,}"))
    sns.despine(ax=ax)


def plot_site_depth_histogram(ax, raw_df, dedup_df):
    """
    Histogram of reads-per-insertion-site for both raw and dedup counts.
    X-axis capped at 99th percentile for readability.
    """
    for df, label, color in [
        (raw_df,   "Raw",   PALETTE["raw"]),
        (dedup_df, "Dedup", PALETTE["dedup"]),
    ]:
        if df.empty or "count" not in df.columns:
            continue
        counts = df["count"].values
        cap = int(np.percentile(counts, 99)) + 1
        counts_clipped = np.clip(counts, 0, cap)
        ax.hist(counts_clipped, bins=min(50, cap),
                alpha=0.55, color=color, label=label,
                edgecolor="white", linewidth=0.5)

    ax.set_xlabel("Reads per insertion site")
    ax.set_ylabel("Number of sites")
    ax.set_title("Site Read Depth Distribution", fontweight="bold")
    ax.legend(fontsize=9)
    sns.despine(ax=ax)


def plot_cumulative_coverage(ax, raw_df, dedup_df):
    """
    Cumulative fraction of sites with ≥N reads.
    Shows how much of the site library is represented at various depth thresholds.
    """
    for df, label, color in [
        (raw_df,   "Raw",   PALETTE["raw"]),
        (dedup_df, "Dedup", PALETTE["dedup"]),
    ]:
        if df.empty or "count" not in df.columns:
            continue
        counts = np.sort(df["count"].values)[::-1]
        cum_frac = np.arange(1, len(counts)+1) / len(counts)
        ax.plot(counts, cum_frac * 100, color=color, label=label, linewidth=1.8)

    ax.set_xscale("log")
    ax.set_xlabel("Minimum reads per site (log scale)")
    ax.set_ylabel("% of sites with ≥ this depth")
    ax.set_title("Cumulative Site Coverage", fontweight="bold")
    ax.axvline(1, color="gray", linestyle=":", linewidth=1)
    ax.axvline(5, color="gray", linestyle=":", linewidth=1)
    ax.text(1.1, 5, "1×", color="gray", fontsize=8)
    ax.text(5.5, 5, "5×", color="gray", fontsize=8)
    ax.legend(fontsize=9)
    sns.despine(ax=ax)


def plot_genomic_features(ax, summary):
    """Pie/bar chart of intragenic vs intergenic insertion fractions."""
    intragenic = float(summary.get("intragenic_sites",  0))
    intergenic = float(summary.get("intergenic_sites",  0))
    total      = intragenic + intergenic

    if total == 0:
        ax.text(0.5, 0.5, "No annotation data", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="gray")
        ax.set_title("Genomic Feature Distribution", fontweight="bold")
        return

    labels = ["Intragenic", "Intergenic"]
    sizes  = [intragenic, intergenic]
    colors = [PALETTE["intragenic"], PALETTE["intergenic"]]
    explode = (0.03, 0)

    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, colors=colors, explode=explode,
        autopct="%1.1f%%", startangle=90,
        wedgeprops={"edgecolor": "white", "linewidth": 1.5},
        textprops={"fontsize": 10},
    )
    for at in autotexts:
        at.set_fontweight("bold")
        at.set_color("white")
        at.set_fontsize(10)

    genes_hit   = int(float(summary.get("genes_with_insertions",  0)))
    total_genes = int(float(summary.get("total_annotated_genes",  0)))
    frac_genes  = float(summary.get("fraction_genes_hit",         0))

    ax.set_title(
        f"Genomic Feature Distribution\n"
        f"{genes_hit}/{total_genes} genes hit ({frac_genes*100:.1f}%)",
        fontweight="bold"
    )


def plot_qc_metrics_table(ax, all_metrics):
    """
    Summary metrics table on the final page.
    Values outside thresholds are highlighted in red.
    """
    ax.axis("off")

    rows = [
        ("Input read pairs",             f"{int(float(all_metrics.get('input_pairs',0))):,}",
         "pct_reads_trimmed_passing"),
        ("Passing trimming (%)",          f"{float(all_metrics.get('pct_reads_trimmed_passing',0)):.1f}%",
         "pct_reads_trimmed_passing"),
        ("Discarded (no Tn5 ME) (%)",     f"{float(all_metrics.get('discarded_untrimmed',0)) / max(float(all_metrics.get('input_pairs',1)),1)*100:.1f}%",
         None),
        ("Alignment rate (%)",            f"{float(all_metrics.get('pct_overall_aligned',0)):.1f}%",
         "pct_mapped"),
        ("Unique concordant (%)",         f"{float(all_metrics.get('pct_mapped',0)):.1f}%",
         "pct_mapped"),
        ("Multi-mapping (%)",             f"{float(all_metrics.get('pct_multi',0)):.1f}%",
         None),
        ("Properly paired (%)",           f"{float(all_metrics.get('pct_properly_paired',0)):.1f}%",
         "pct_properly_paired"),
        ("Median insert size (bp)",       f"{float(all_metrics.get('median_insert_size',0)):.0f}",
         "median_insert_size_bp"),
        ("Duplication rate (%)",          f"{float(all_metrics.get('duplication_rate_pct',0)):.1f}%",
         "duplication_rate_pct"),
        ("Total unique sites (raw)",      f"{int(float(all_metrics.get('total_unique_insertion_sites',0))):,}",
         None),
        ("Intragenic sites (%)",          f"{float(all_metrics.get('pct_intragenic',0))*100:.1f}%",
         "pct_intragenic"),
        ("Genes with insertions",         f"{int(float(all_metrics.get('genes_with_insertions',0))):,}",
         None),
        ("Fraction of genes hit (%)",     f"{float(all_metrics.get('fraction_genes_hit',0))*100:.1f}%",
         None),
    ]

    table_data  = [[r[0], r[1]] for r in rows]
    cell_colors = []
    for _, _, threshold_key in rows:
        if threshold_key and threshold_key in all_metrics:
            try:
                val = float(str(all_metrics[threshold_key]).replace("%",""))
                color = flag_metric(val, threshold_key)
                cell_colors.append(["#f8f9fa", color + "33"])   # light tint for value cell
            except (ValueError, TypeError):
                cell_colors.append(["#f8f9fa", "#f8f9fa"])
        else:
            cell_colors.append(["#f8f9fa", "#f8f9fa"])

    tbl = ax.table(
        cellText=table_data,
        colLabels=["Metric", "Value"],
        cellLoc="left",
        loc="center",
        cellColours=cell_colors,
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.2, 1.8)

    # Style header
    for j in range(2):
        tbl[0, j].set_facecolor(PALETTE["neutral"])
        tbl[0, j].set_text_props(color="white", fontweight="bold")

    ax.set_title(f"QC Summary Metrics", fontweight="bold", pad=20)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate Tn5 pipeline QC report")
    parser.add_argument("--sample",         required=True)
    parser.add_argument("--cutadapt-json",  required=True,
                        help="cutadapt --json output (or plain text if JSON unavailable)")
    parser.add_argument("--bowtie2-stats",  required=True)
    parser.add_argument("--samtools-stats", required=True)
    parser.add_argument("--dedup-stats",    required=True)
    parser.add_argument("--raw-sites",      required=True)
    parser.add_argument("--dedup-sites",    required=True)
    parser.add_argument("--summary-stats",  required=True)
    parser.add_argument("--outdir",         required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ------------------------------------------------------------------
    # Parse all inputs
    # ------------------------------------------------------------------
    print(f"  Parsing QC inputs for sample: {args.sample}")

    # cutadapt: try JSON first, fall back to plain text
    if args.cutadapt_json.endswith(".json") and os.path.exists(args.cutadapt_json):
        try:
            cutadapt = parse_cutadapt_json(args.cutadapt_json)
        except (KeyError, json.JSONDecodeError):
            # Some cutadapt versions have different JSON schemas; fall back
            txt_path = args.cutadapt_json.replace(".json", ".txt")
            cutadapt = parse_cutadapt_txt(txt_path) if os.path.exists(txt_path) else {}
    else:
        cutadapt = parse_cutadapt_txt(args.cutadapt_json)

    bowtie2  = parse_bowtie2_stats(args.bowtie2_stats)
    samstats = parse_samtools_stats(args.samtools_stats)
    dedup    = parse_dedup_stats(args.dedup_stats)
    raw_df   = parse_sites_bed(args.raw_sites)
    dedup_df = parse_sites_bed(args.dedup_sites)
    summary  = parse_summary_stats(args.summary_stats)

    # ------------------------------------------------------------------
    # Merge all metrics into one flat dict for the summary table + TSV
    # ------------------------------------------------------------------
    all_metrics = {}
    all_metrics.update(cutadapt)
    all_metrics.update(bowtie2)
    all_metrics["median_insert_size"]    = samstats.get("median_insert_size", 0)
    all_metrics["mean_insert_size"]      = samstats.get("mean_insert_size",   0)
    all_metrics["pct_properly_paired"]   = samstats.get("pct_properly_paired", 0)
    all_metrics["reads_mapped"]          = samstats.get("reads_mapped",        0)
    all_metrics.update({k: v for k, v in dedup.items()})
    all_metrics["n_raw_sites"]   = len(raw_df)
    all_metrics["n_dedup_sites"] = len(dedup_df)
    all_metrics.update({k: v for k, v in summary.items()
                        if k not in ("sample",)})
    all_metrics["pct_intragenic"] = float(summary.get("fraction_intragenic", 0))
    all_metrics["sample"] = args.sample

    # ------------------------------------------------------------------
    # Save QC summary TSV
    # ------------------------------------------------------------------
    tsv_path = os.path.join(args.outdir, f"{args.sample}_qc_summary.tsv")
    pd.DataFrame([all_metrics]).T.reset_index().rename(
        columns={"index": "metric", 0: args.sample}
    ).to_csv(tsv_path, sep="\t", index=False)
    print(f"  QC summary TSV: {tsv_path}")

    # ------------------------------------------------------------------
    # Build PDF
    # ------------------------------------------------------------------
    pdf_path = os.path.join(args.outdir, f"{args.sample}_qc_report.pdf")

    with PdfPages(pdf_path) as pdf:

        # --- Page 1: Read fate, trimming outcome, insert size, mapping ---
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f"QC Report — {args.sample}\n(Tn5 Insertion Mapping Pipeline)",
                     fontsize=14, fontweight="bold", y=1.01)
        plot_read_funnel(       axes[0, 0], cutadapt, samstats, dedup)
        plot_trimming_outcome(  axes[0, 1], cutadapt)
        plot_insert_size(       axes[1, 0], samstats)
        plot_mapping_stats(     axes[1, 1], bowtie2)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # --- Page 2: Duplication, site depth, cumulative coverage, features ---
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f"QC Report — {args.sample}  (continued)",
                     fontsize=14, fontweight="bold", y=1.01)
        plot_duplication(           axes[0, 0], dedup)
        plot_site_depth_histogram(  axes[0, 1], raw_df, dedup_df)
        plot_cumulative_coverage(   axes[1, 0], raw_df, dedup_df)
        plot_genomic_features(      axes[1, 1], summary)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # --- Page 3: Summary metrics table ---
        fig, ax = plt.subplots(figsize=(10, 8))
        plot_qc_metrics_table(ax, all_metrics)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # PDF metadata
        d = pdf.infodict()
        d["Title"]   = f"Tn5 Insertion Mapping QC — {args.sample}"
        d["Subject"] = "Transposon insertion library QC"

    print(f"  QC report PDF: {pdf_path}")
    print(f"  Pages: 3  (read funnel + trimming + insert size + mapping |"
          f" dedup + site depth + cumulative coverage + features | summary table)")


if __name__ == "__main__":
    main()
