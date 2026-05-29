#!/usr/bin/env python3
"""
annotate_insertions.py
----------------------
Annotate Tn5 insertion sites against a GFF3 gene annotation.

Inputs:
  --raw-sites   : BED6 of insertion sites from ALL mapped reads (before dedup)
                  Score column (col 5) = raw read count supporting the site
  --dedup-sites : BED6 of insertion sites from coordinate-pair-dedup reads
                  Score column (col 5) = dedup read count supporting the site
  --gff         : GFF3 annotation (C. necator H16 or any standard GFF3)
  --sample      : Sample name (used in output header)

Outputs:
  --outfile      : Per-site annotated TSV with both raw and dedup counts
  --per-gene-out : Per-gene aggregated count table
  --stats-out    : Summary statistics TSV

Usage (called by tn5_insertion_mapping.sh):
    python3 annotate_insertions.py \
        --raw-sites   sample_raw_insertion_sites.bed \
        --dedup-sites sample_dedup_insertion_sites.bed \
        --gff         cnecator_H16.gff \
        --sample      SAMPLE \
        --outfile     sample_insertions_annotated.tsv \
        --per-gene-out sample_per_gene_counts.tsv \
        --stats-out    sample_summary_stats.tsv

Requires: pandas >= 1.3, bedtools on PATH
"""

import argparse
import os
import re
import subprocess
import sys

import pandas as pd


# ---------------------------------------------------------------------------
# GFF3 PARSER
# ---------------------------------------------------------------------------

def parse_gff3_genes(gff_path):
    """
    Parse GFF3 and return a DataFrame of gene-level rows with columns:
      chrom, start (0-based), end, strand, gene_id, gene_name, locus_tag
    Only 'gene' feature type rows are used for overlap assignment.
    """
    genes = []

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _, ftype, start, end, _, strand, _, attrs = parts

            if ftype != "gene":
                continue

            start = int(start) - 1  # GFF3 is 1-based inclusive → 0-based BED
            end   = int(end)

            attr_dict = {}
            for kv in attrs.split(";"):
                kv = kv.strip()
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_dict[k.strip()] = v.strip()

            genes.append({
                "chrom":     chrom,
                "start":     start,
                "end":       end,
                "strand":    strand,
                "gene_id":   attr_dict.get("ID", ""),
                "gene_name": attr_dict.get("Name", attr_dict.get("gene", "")),
                "locus_tag": attr_dict.get("locus_tag", ""),
            })

    if not genes:
        sys.exit("ERROR: No 'gene' features found in GFF3. Check file format.")

    return pd.DataFrame(genes)


# ---------------------------------------------------------------------------
# BED READER
# ---------------------------------------------------------------------------

def read_sites_bed(bed_path, count_col_name):
    """
    Read a BED6 file of insertion sites.
    Column 5 (score) holds the read count supporting each merged site.
    Returns DataFrame with: chrom, start, end, strand, <count_col_name>
    """
    df = pd.read_csv(
        bed_path, sep="\t", header=None,
        names=["chrom", "ins_start", "ins_end", "name", count_col_name, "strand"]
    )
    return df[["chrom", "ins_start", "ins_end", "strand", count_col_name]]


# ---------------------------------------------------------------------------
# BEDTOOLS INTERSECT (subprocess)
# ---------------------------------------------------------------------------

def intersect_sites_with_genes(sites_df, genes_df, tmp_prefix):
    """
    Intersect insertion sites BED with gene features using bedtools.
    Returns the sites DataFrame with gene annotation columns added.
    Intergenic sites receive gene_id='intergenic', all gene fields null/empty.
    """
    sites_bed = f"{tmp_prefix}_sites.bed"
    genes_bed = f"{tmp_prefix}_genes.bed"

    # Write sites BED (BED6 format expected by bedtools)
    sites_df.assign(name="INS", score=0)[
        ["chrom", "ins_start", "ins_end", "name", "score", "strand"]
    ].to_csv(sites_bed, sep="\t", index=False, header=False)

    # Write genes BED6: chrom, start, end, gene_id, gene_name, strand
    genes_df[["chrom", "start", "end", "gene_id", "gene_name", "strand"]].to_csv(
        genes_bed, sep="\t", index=False, header=False
    )

    cmd = [
        "bedtools", "intersect",
        "-a", sites_bed,
        "-b", genes_bed,
        "-wa", "-wb", "-loj"   # left outer join: keep all sites
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    os.remove(sites_bed)
    os.remove(genes_bed)

    rows = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        f = line.split("\t")
        # f[0:6]  = site BED6  (chrom, start, end, name, score, strand)
        # f[6:12] = gene BED6  (chrom, start, end, gene_id, gene_name, gene_strand)
        #           or '.' if no overlap
        site_chrom  = f[0]
        site_start  = int(f[1])
        site_end    = int(f[2])
        site_strand = f[5]

        if f[6] == ".":
            gene_id    = "intergenic"
            gene_name  = ""
            gene_start = None
            gene_end   = None
            gene_strand = "."
        else:
            gene_id    = f[9]
            gene_name  = f[10]
            gene_start = int(f[7])
            gene_end   = int(f[8])
            gene_strand = f[11]

        rows.append({
            "chrom":       site_chrom,
            "ins_start":   site_start,
            "ins_end":     site_end,
            "site_strand": site_strand,
            "gene_id":     gene_id,
            "gene_name":   gene_name,
            "gene_start":  gene_start,
            "gene_end":    gene_end,
            "gene_strand": gene_strand,
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# RELATIVE POSITION WITHIN GENE
# ---------------------------------------------------------------------------

def relative_position(row):
    """
    Fractional position of insertion within gene body (0.0 = TSS, 1.0 = stop),
    in the direction of transcription. Returns None for intergenic sites.
    """
    if row["gene_id"] == "intergenic" or pd.isna(row.get("gene_start")):
        return None
    gene_len = row["gene_end"] - row["gene_start"]
    if gene_len == 0:
        return None
    ins = row["ins_start"]
    if row["gene_strand"] == "+":
        frac = (ins - row["gene_start"]) / gene_len
    else:
        frac = (row["gene_end"] - ins) / gene_len
    return round(max(0.0, min(1.0, frac)), 4)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Annotate Tn5 insertion sites with raw and dedup counts")
    parser.add_argument("--raw-sites",    required=True,  help="BED6 of raw (all-reads) insertion sites")
    parser.add_argument("--dedup-sites",  required=True,  help="BED6 of coordinate-pair-dedup insertion sites")
    parser.add_argument("--gff",          required=True,  help="GFF3 annotation file")
    parser.add_argument("--sample",       required=True,  help="Sample name (for output labeling)")
    parser.add_argument("--outfile",      required=True,  help="Output: per-site annotated TSV")
    parser.add_argument("--per-gene-out", required=True,  help="Output: per-gene count table TSV")
    parser.add_argument("--stats-out",    required=True,  help="Output: summary statistics TSV")
    args = parser.parse_args()

    tmp_prefix = args.outfile.replace(".tsv", "_tmp")

    # ------------------------------------------------------------------
    # Load inputs
    # ------------------------------------------------------------------
    print(f"  Parsing GFF3: {args.gff}")
    genes_df = parse_gff3_genes(args.gff)
    n_genes  = len(genes_df)
    print(f"  Loaded {n_genes} gene features.")

    print(f"  Reading raw sites:   {args.raw_sites}")
    raw_df   = read_sites_bed(args.raw_sites,   "raw_read_count")

    print(f"  Reading dedup sites: {args.dedup_sites}")
    dedup_df = read_sites_bed(args.dedup_sites, "dedup_read_count")

    # ------------------------------------------------------------------
    # Annotate raw sites (primary intersection)
    # ------------------------------------------------------------------
    if raw_df.empty:
        print("  WARNING: No insertion sites found. Writing empty output files.")
        pd.DataFrame(columns=["chrom","ins_start","ins_end","site_strand",
            "raw_read_count","dedup_read_count","feature_type","gene_id",
            "gene_name","locus_tag","gene_start","gene_end","gene_strand",
            "relative_pos_in_gene"]).to_csv(args.outfile, sep="\t", index=False)
        pd.DataFrame().to_csv(args.per_gene_out, sep="\t", index=False)
        pd.DataFrame([{"metric":"total_unique_insertion_sites","value":0}]).to_csv(
            args.stats_out, sep="\t", index=False)
        print("  Check trimming and alignment steps.")
        return

    print("  Intersecting raw sites with gene features...")
    annotated = intersect_sites_with_genes(raw_df, genes_df, tmp_prefix + "_raw")

    # Merge raw counts back by position (sites_df may have been split by bedtools
    # if a site overlaps multiple gene features — keep the first/best hit)
    raw_count_map = raw_df.set_index(["chrom", "ins_start"])["raw_read_count"].to_dict()
    annotated["raw_read_count"] = annotated.apply(
        lambda r: raw_count_map.get((r["chrom"], r["ins_start"]), 0), axis=1
    )

    # Drop duplicate site entries keeping intragenic over intergenic when both exist
    annotated["_is_genic"] = annotated["gene_id"] != "intergenic"
    annotated = (
        annotated
        .sort_values(["chrom", "ins_start", "_is_genic"], ascending=[True, True, False])
        .drop_duplicates(subset=["chrom", "ins_start"], keep="first")
        .drop(columns=["_is_genic"])
    )

    # ------------------------------------------------------------------
    # Merge dedup counts by coordinate
    # ------------------------------------------------------------------
    dedup_count_map = dedup_df.set_index(["chrom", "ins_start"])["dedup_read_count"].to_dict()
    annotated["dedup_read_count"] = annotated.apply(
        lambda r: dedup_count_map.get((r["chrom"], r["ins_start"]), 0), axis=1
    )

    # ------------------------------------------------------------------
    # Add locus_tag and relative position
    # ------------------------------------------------------------------
    locus_map = genes_df.set_index("gene_id")["locus_tag"].to_dict()
    annotated["locus_tag"] = annotated["gene_id"].map(locus_map).fillna("")

    annotated["relative_pos_in_gene"] = annotated.apply(relative_position, axis=1)

    annotated["feature_type"] = annotated["gene_id"].apply(
        lambda g: "intragenic" if g != "intergenic" else "intergenic"
    )

    # ------------------------------------------------------------------
    # Write annotated sites TSV
    # ------------------------------------------------------------------
    out_cols = [
        "chrom", "ins_start", "ins_end", "site_strand",
        "raw_read_count", "dedup_read_count",
        "feature_type", "gene_id", "gene_name", "locus_tag",
        "gene_start", "gene_end", "gene_strand",
        "relative_pos_in_gene"
    ]
    annotated[out_cols].sort_values(["chrom", "ins_start"]).to_csv(
        args.outfile, sep="\t", index=False
    )
    print(f"  Annotated sites written: {args.outfile}")

    # ------------------------------------------------------------------
    # Per-gene count table
    # ------------------------------------------------------------------
    intragenic = annotated[annotated["feature_type"] == "intragenic"].copy()

    per_gene = (
        intragenic
        .groupby(["gene_id", "gene_name", "locus_tag",
                  "chrom", "gene_start", "gene_end", "gene_strand"])
        .agg(
            n_unique_sites_raw      = ("ins_start",        "nunique"),
            total_raw_reads         = ("raw_read_count",   "sum"),
            n_unique_sites_dedup    = ("ins_start",        lambda x:
                                       (intragenic.loc[x.index, "dedup_read_count"] > 0).sum()),
            total_dedup_reads       = ("dedup_read_count", "sum"),
            mean_relative_pos       = ("relative_pos_in_gene", "mean"),
        )
        .reset_index()
        .sort_values("total_raw_reads", ascending=False)
    )

    per_gene["gene_length_bp"] = per_gene["gene_end"] - per_gene["gene_start"]
    per_gene["unique_sites_per_kb_raw"] = (
        per_gene["n_unique_sites_raw"] / (per_gene["gene_length_bp"] / 1000)
    ).round(3)
    per_gene["unique_sites_per_kb_dedup"] = (
        per_gene["n_unique_sites_dedup"] / (per_gene["gene_length_bp"] / 1000)
    ).round(3)
    per_gene["mean_relative_pos"] = per_gene["mean_relative_pos"].round(4)

    per_gene.to_csv(args.per_gene_out, sep="\t", index=False)
    print(f"  Per-gene counts written: {args.per_gene_out}")

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------
    n_sites_total      = len(annotated)
    n_sites_intragenic = (annotated["feature_type"] == "intragenic").sum()
    n_sites_intergenic = (annotated["feature_type"] == "intergenic").sum()
    n_genes_hit        = per_gene["gene_id"].nunique()

    total_raw_reads    = int(annotated["raw_read_count"].sum())
    total_dedup_reads  = int(annotated["dedup_read_count"].sum())
    duplication_pct    = (1 - total_dedup_reads / total_raw_reads) * 100 if total_raw_reads > 0 else 0

    stats = {
        "sample":                           args.sample,
        "total_raw_reads_at_sites":         total_raw_reads,
        "total_dedup_reads_at_sites":       total_dedup_reads,
        "duplication_pct_at_sites":         round(duplication_pct, 2),
        "total_unique_insertion_sites":     n_sites_total,
        "intragenic_sites":                 int(n_sites_intragenic),
        "intergenic_sites":                 int(n_sites_intergenic),
        "fraction_intragenic":              round(n_sites_intragenic / max(n_sites_total, 1), 4),
        "genes_with_insertions":            int(n_genes_hit),
        "total_annotated_genes":            n_genes,
        "fraction_genes_hit":               round(n_genes_hit / max(n_genes, 1), 4),
        "median_raw_reads_per_site":        float(annotated["raw_read_count"].median()),
        "max_raw_reads_per_site":           int(annotated["raw_read_count"].max()),
        "median_dedup_reads_per_site":      float(annotated["dedup_read_count"].median()),
        "max_dedup_reads_per_site":         int(annotated["dedup_read_count"].max()),
    }

    stats_df = pd.DataFrame([stats]).T.reset_index()
    stats_df.columns = ["metric", "value"]
    stats_df.to_csv(args.stats_out, sep="\t", index=False)

    print("")
    print("  === Summary Statistics ===")
    for _, row in stats_df.iterrows():
        print(f"  {row['metric']:45s}: {row['value']}")

    print(f"\n  Summary stats written: {args.stats_out}")


if __name__ == "__main__":
    main()
