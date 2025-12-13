#!/usr/bin/env python3
import pandas as pd
import numpy as np

# -----------------------------
# Helper: normalize REF/ALT for indels
# -----------------------------

def normalize_ref_alt(df, ref_col="ref", alt_col="alt"):
    """
    Apply the rule:
      - If REF and ALT are both non-empty/non-dot and have different lengths,
        drop the first character from each.
      - Otherwise leave them as is.

    This is done in-place on the given columns and also returns the df.
    """
    ref = df[ref_col].astype(str)
    alt = df[alt_col].astype(str)

    # Identify rows where we want to trim
    mask = (
        (ref != "") & (ref != ".") &
        (alt != "") & (alt != ".") &
        (ref.str.len() != alt.str.len())
    )

    df.loc[mask, ref_col] = ref[mask].str.slice(1)
    df.loc[mask, alt_col] = alt[mask].str.slice(1)

    return df

# -----------------------------
# 1. Load per-variant correctness
# -----------------------------

bial_path = "bial_per_site.tsv"    # baseline run (e.g. og/bi)
super_path = "super_per_site.tsv"  # supersite run
supersite_annot_path = "supersite_annotation.tsv"

bial = pd.read_csv(bial_path, sep="\t")
sup  = pd.read_csv(super_path, sep="\t")

required_cols = {"chrom", "pos", "ref", "alt", "sample", "is_correct"}
for name, df in [("bial_per_site.tsv", bial), ("super_per_site.tsv", sup)]:
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")

# Make sure types are consistent
for df in (bial, sup):
    df["chrom"] = df["chrom"].astype(str)
    df["pos"]   = df["pos"].astype(int)
    df["ref"]   = df["ref"].fillna("").astype(str)
    df["alt"]   = df["alt"].fillna("").astype(str)
    df["sample"]= df["sample"].astype(str)
    df["is_correct"] = df["is_correct"].astype(int)

# Normalize REF/ALT in both runs
bial = normalize_ref_alt(bial, "ref", "alt")
sup  = normalize_ref_alt(sup,  "ref", "alt")

# Rename correctness columns
bial = bial.rename(columns={"is_correct": "B_correct"})
sup  = sup.rename(columns={"is_correct": "S_correct"})

# -------------------------------------------
# 2. Inner join on same variant & sample
# -------------------------------------------

merge_keys = ["chrom", "pos", "ref", "alt", "sample"]

df = bial.merge(
    sup[merge_keys + ["S_correct"]],
    on=merge_keys,
    how="inner",
)

df["B_correct"] = df["B_correct"].astype(int)
df["S_correct"] = df["S_correct"].astype(int)

# -----------------------------------
# 3. Classify B/S correctness combos (per sample)
# -----------------------------------

def classify_row(row):
    b = row["B_correct"]
    s = row["S_correct"]
    if b == 1 and s == 1:
        return "both_correct"
    elif b == 1 and s == 0:
        return "regression"  # correct in baseline, wrong in supersite
    elif b == 0 and s == 1:
        return "supersite_improved"
    elif b == 0 and s == 0:
        return "both_wrong"
    else:
        return "unknown"

df["combo"] = df.apply(classify_row, axis=1)

# --------------------------------
# 4. Load & normalize supersite annotations
# --------------------------------

annot = pd.read_csv(supersite_annot_path, sep="\t")

# Expect a column that has the combined variant id like chrom_pos_ref_alt
if "variant_id" not in annot.columns:
    raise ValueError("supersite_annotation.tsv must have a 'variant_id' column of form chrom_pos_ref_alt")

# Parse variant_id -> chrom,pos,ref,alt
# Example variant_id: chr22_19057961_A_G
split_cols = annot["variant_id"].astype(str).str.split("_", n=3, expand=True)
if split_cols.shape[1] != 4:
    raise ValueError("Expected variant_id of the form chrom_pos_ref_alt (3 underscores).")

annot["chrom"] = split_cols[0].astype(str)
annot["pos"]   = split_cols[1].astype(int)
annot["ref"]   = split_cols[2].fillna("").astype(str)
annot["alt"]   = split_cols[3].fillna("").astype(str)

# Normalize ref/alt in annotations with the same rule
annot = normalize_ref_alt(annot, "ref", "alt")

# Keep only relevant columns
keep_cols = ["chrom", "pos", "ref", "alt", "variant_id"]
for col in ["supersite_id", "supersite_role"]:
    if col in annot.columns:
        keep_cols.append(col)

annot = annot[keep_cols].drop_duplicates()

# Merge by chrom,pos,ref,alt
df = df.merge(
    annot,
    on=["chrom", "pos", "ref", "alt"],
    how="left",
)

if "supersite_role" in df.columns:
    df["supersite_role"] = df["supersite_role"].fillna("none")
else:
    df["supersite_role"] = "none"

if "supersite_id" in df.columns:
    df["supersite_id"] = df["supersite_id"].fillna("none")

# --------------------------------------------
# 5. Global per-sample stats (as before)
# --------------------------------------------

total_rows = len(df)
print(f"Total sample-variant rows (shared between runs): {total_rows:,}\n")

combo_counts = df["combo"].value_counts().sort_index()
print("Global counts by combo (per sample-variant):")
for combo, count in combo_counts.items():
    frac = count / total_rows * 100.0
    print(f"  {combo:18s} : {count:10d} ({frac:6.3f}%)")
print()

bial_error_rate = (df["B_correct"] == 0).mean() * 100.0
super_error_rate = (df["S_correct"] == 0).mean() * 100.0
print(f"Baseline error rate (per sample-variant)      : {bial_error_rate:6.3f}%")
print(f"Supersite run error rate (per sample-variant) : {super_error_rate:6.3f}%\n")

print("Per-sample counts by supersite_role × combo:")
role_combo = (
    df
    .groupby(["supersite_role", "combo"])
    .size()
    .unstack(fill_value=0)
)
role_combo["TOTAL"] = role_combo.sum(axis=1)
role_combo_frac = role_combo.div(role_combo["TOTAL"], axis=0) * 100.0

print(role_combo)
print("\nRow-wise % fractions (per sample-variant):")
print(role_combo_frac.round(3))

# --------------------------------------------
# 6. Per-variant statistics (chrom,pos,ref,alt define a variant)
# --------------------------------------------

variant_keys = ["chrom", "pos", "ref", "alt", "supersite_role", "supersite_id", "variant_id"]

# Ensure all exist (some may not be in annot; fill if missing)
for col in ["supersite_id", "variant_id"]:
    if col not in df.columns:
        df[col] = "none"

grouped = df.groupby(variant_keys, dropna=False)

variant_stats = grouped.agg(
    n_samples=("sample", "nunique"),
    B_correct_mean=("B_correct", "mean"),
    S_correct_mean=("S_correct", "mean"),
    B_error_rate=("B_correct", lambda x: 1.0 - float(np.mean(x))),
    S_error_rate=("S_correct", lambda x: 1.0 - float(np.mean(x))),
    n_regression=("combo", lambda x: (x == "regression").sum()),
    n_supersite_improved=("combo", lambda x: (x == "supersite_improved").sum()),
    n_both_wrong=("combo", lambda x: (x == "both_wrong").sum()),
    n_both_correct=("combo", lambda x: (x == "both_correct").sum()),
)

variant_stats = variant_stats.reset_index()

# Quick sanity summary
n_variants = len(variant_stats)
print(f"\nTotal unique variants (chrom,pos,ref,alt): {n_variants:,}")
worse = (variant_stats["S_error_rate"] > variant_stats["B_error_rate"]).sum()
better = (variant_stats["S_error_rate"] < variant_stats["B_error_rate"]).sum()
same   = (variant_stats["S_error_rate"] == variant_stats["B_error_rate"]).sum()
print(f"Variants where supersite is worse : {worse:,}")
print(f"Variants where supersite is better: {better:,}")
print(f"Variants where equal              : {same:,}")

# --------------------------------------------
# 7. Write outputs
# --------------------------------------------

out_per_row = "bial_vs_supersite_per_row_annotated.tsv"
out_per_var = "bial_vs_supersite_per_variant_stats.tsv"

df.to_csv(out_per_row, sep="\t", index=False)
variant_stats.to_csv(out_per_var, sep="\t", index=False)

print(f"\nWrote annotated per sample-variant table to: {out_per_row}")
print(f"Wrote per-variant statistics to:           {out_per_var}")
