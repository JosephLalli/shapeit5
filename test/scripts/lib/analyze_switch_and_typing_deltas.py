#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Switch files: space-delimited, no header
# Columns: var_id pos n_switch_err n_biallelic SER
COLS_SWITCH = ["id", "pos", "nse", "nb", "SER"]

# Typing files: space-delimited, no header
# Columns: id pos errors ground_truth_genotypes gt_err_rate
COLS_TYPING = ["id", "pos", "errors", "gt", "gt_err_rate"]


def guess_typing_path(switch_path: str) -> str:
    """Replace 'switch' with 'typing' once in the filename."""
    if "switch" in switch_path:
        return switch_path.replace("switch", "typing", 1)
    # fallback: if the pattern isn't present, just complain and append
    base, ext = os.path.splitext(switch_path)
    print(f"[WARN] No 'switch' in {switch_path}; guessing {base}.typing{ext}")
    return f"{base}.typing{ext}"

def annotate_supersites(df: pd.DataFrame) -> pd.DataFrame:
    """Annotate supersites in a switch DataFrame."""
    # Count variants per position
    df["var_count"] = df.reset_index().groupby("pos")["index"].transform("count")
    
    # Assign supersite index
    df["ss_idx"] = df.groupby("pos").cumcount() + 1
    
    # Determine sup_status
    def determine_status(row):
        if row["var_count"] == 1:
            return "non"
        elif row["ss_idx"] == 1:
            return "anchor"
        else:
            return "sibling"
    
    df["sup_status"] = df.apply(determine_status, axis=1)
    
    return df


def load_switch_supersites(path: str) -> pd.DataFrame:
    """
    Load supersites switch file and auto-detect supersites:
      - supersite = position with >1 variant
      - anchor    = first variant at that position (file order)
      - sibling   = other variants at that position
      - non       = positions with only one variant
    Returns columns:
      id, pos, SER_super, sup_status, var_count, ss_idx
    """
    df = pd.read_csv(path, sep=" ", header=None,
                     names=COLS_SWITCH, compression="infer")

    # # Compute SER if missing or all zero
    # if df["SER"].isna().all() or (df["SER"] == 0).all():
    #     df["SER"] = np.where(
    #         df["nb"] > 0,
    #         100.0 * df["nse"] / df["nb"],
    #         np.nan
    #     )

    # Preserve file order
    df["row_order"] = np.arange(len(df))

    df = annotate_supersites(df)

    out = df[["id", "pos", "nse", "nb", "SER", "sup_status", "var_count", "ss_idx"]].copy()
    out = out.rename(columns={"SER": "SER_super",
                              "nse": "nse_super",
                              "nb": "nb_super"
                           })
    return out


def load_switch_general(path: str, label: str) -> pd.DataFrame:
    """Load a generic switch file -> id, pos, SER_<label>."""
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=COLS_SWITCH, compression="infer")

    # Compute SER if needed
    # if df["SER"].isna().all() or (df["SER"] == 0).all():
    #     df["SER"] = np.where(
    #         df["nb"] > 0,
    #         100.0 * df["nse"] / df["nb"],
    #         np.nan
    #     )

    # df = df[["id", "pos", "SER"]].copy()
    df = df.rename(columns={
        "nse": f"nse_{label}",
        "nb": f"nb_{label}",
        "pos": f"pos_{label}",
        "SER": f"SER_{label}",
    })
    return df


def load_typing(path: str, label: str) -> pd.DataFrame:
    """Load a typing file -> id, pos_<label>, GTERR_<label>."""
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=COLS_TYPING, compression="infer")
    df = df[["id", "pos", "gt_err_rate"]].copy()
    df = df.rename(columns={
        "pos": f"pos_{label}",
        "gt_err_rate": f"GTERR_{label}",
    })
    return df


def resolve_position(df: pd.DataFrame) -> pd.Series:
    """Choose a single 'pos' column, preferring supersites, then main, then og."""
    pos = pd.Series(np.nan, index=df.index, dtype=float)
    for col in ["pos", "pos_main", "pos_og"]:
        if col in df.columns:
            pos = pos.fillna(df[col].astype(float))
    return pos


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--supersites", required=True,
                   help="supersites switch file (space-delimited, no header)")
    p.add_argument("--main", required=True,
                   help="main_algo switch file (space-delimited, no header)")
    p.add_argument("--og", required=True,
                   help="original switch file (space-delimited, no header)")
    p.add_argument("--input-positions", required=False,
                   help="optional file with raw variant positions for determining supersite locations (one per line)")
    p.add_argument("--out-tsv", default="variant_error_comparison.tsv")
    p.add_argument("--out-prefix", default="plots")
    args = p.parse_args()

    # --- Load switch data ---
    df_super = load_switch_supersites(args.supersites).rename(columns={'pos_super': 'pos'})
    df_main  = load_switch_general(args.main, "main").drop('pos_main', axis=1)
    df_og    = load_switch_general(args.og,   "og").drop('pos_og', axis=1)

    # --- Load typing data for each run (auto path) ---
    super_typ_path = guess_typing_path(args.supersites)
    main_typ_path  = guess_typing_path(args.main)
    og_typ_path    = guess_typing_path(args.og)

    print(f"[INFO] Supersite typing: {super_typ_path}")
    print(f"[INFO] Main typing:      {main_typ_path}")
    print(f"[INFO] OG typing:        {og_typ_path}")

    # df_super_typ = load_typing(super_typ_path, "super").drop('pos_super', axis=1)
    # df_main_typ  = load_typing(main_typ_path,  "main").drop('pos_main', axis=1)
    # df_og_typ    = load_typing(og_typ_path,    "og").drop('pos_og', axis=1)


    # --- Merge everything on variant ID ---
    # Start from supersite switch (as it carries sup_status)
    df = df_super.merge(df_main, on="id", how="outer") \
                 .merge(df_og,   on="id", how="outer") #\
                #  .merge(df_super_typ, on="id", how="outer") \
                #  .merge(df_main_typ,  on="id", how="outer") \
                #  .merge(df_og_typ,    on="id", how="outer")
    # --- If input positions provided, use to identify supersites ---
    if args.input_positions:
        pos=pd.read_csv(args.input_positions, sep="\t", header=None, names=["pos"])
        pos = pos.loc[pos.pos.isin(df.pos)]
        pos=annotate_supersites(pos.reset_index(drop=True))
        
        # Merge sup_status back into df_super
        if 'sup_status' in df.columns:
            df = df.drop(columns=["sup_status", "var_count"])
        
        df = df.merge(
                pos[["pos", "sup_status", "var_count", "ss_idx"]], on=['pos', "ss_idx"], how='left')

    # --- Switch error deltas ---
    df["d_SER_super_main"] = df["SER_super"] - df["SER_main"]
    df["d_SER_super_og"]   = df["SER_super"] - df["SER_og"]
    df["d_SER_main_og"]    = df["SER_main"]  - df["SER_og"]

    # # --- Genotype error deltas ---
    # df["d_GT_super_main"] = df["GTERR_super"] - df["GTERR_main"]
    # df["d_GT_super_og"]   = df["GTERR_super"] - df["GTERR_og"]
    # df["d_GT_main_og"]    = df["GTERR_main"]  - df["GTERR_og"]

    # --- Write out full table ---
    # Keep everything; still space-delimited, no header
    df.to_csv(args.out_tsv, sep=" ", index=False, header=True)
    print(f"[INFO] Wrote {args.out_tsv} with {len(df):,} rows")

    # ==================== PLOTS ====================

    # Helper to safely get mask
    def not_nan(*cols):
        m = pd.Series(True, index=df.index)
        for c in cols:
            m &= df[c].notna()
        return m

    # 1) SER_super vs SER_main scatter
    plt.figure(figsize=(6, 6))
    m = not_nan("SER_super", "SER_main")
    plt.scatter(df.loc[m, "SER_main"], df.loc[m, "SER_super"], s=5, alpha=0.4)
    lim = float(np.nanmax([df["SER_main"].max(), df["SER_super"].max()]))
    plt.plot([0, lim], [0, lim], "--")
    plt.xlabel("SER_main (%)")
    plt.ylabel("SER_super (%)")
    plt.title("Switch error: supersites vs main")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.SER_super_vs_main.png", dpi=200)
    plt.close()

    # 2) Histogram of ΔSER (super - main)
    plt.figure(figsize=(6, 4))
    plt.hist(df["d_SER_super_main"].dropna(), bins=80)
    plt.axvline(0, linestyle="--")
    plt.xlabel("ΔSER (super − main)")
    plt.ylabel("Count")
    plt.title("Distribution of Δ switch error (supersites − main)")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.d_SER_super_main.hist.png", dpi=200)
    plt.close()

    # # 3) Histogram of ΔGTERR (super - main)
    # plt.figure(figsize=(6, 4))
    # plt.hist(df["d_GT_super_main"].dropna(), bins=80)
    # plt.axvline(0, linestyle="--")
    # plt.xlabel("ΔGTERR (super − main)")
    # plt.ylabel("Count")
    # plt.title("Distribution of Δ genotype error (supersites − main)")
    # plt.tight_layout()
    # plt.savefig(f"{args.out_prefix}.d_GT_super_main.hist.png", dpi=200)
    # plt.close()

    # 4) ΔSER vs position, colored by sup_status
    plt.figure(figsize=(10, 4))
    if "sup_status" in df.columns:
        m = not_nan("d_SER_super_main", "pos")
        for status, sub in df.loc[m].groupby("sup_status"):
            plt.scatter(sub["pos"], sub["d_SER_super_main"],
                        s=6, alpha=0.6, label=status)
        plt.legend(title="sup_status")
    else:
        plt.scatter(df["pos"], df["d_SER_super_main"], s=3, alpha=0.4)
    plt.axhline(0, linestyle="--")
    plt.xlabel("Position (bp)")
    plt.ylabel("ΔSER (super − main)")
    plt.title("Δ switch error along chromosome")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.d_SER_super_main.by_pos.png", dpi=200)
    plt.close()

    # # 5) ΔGTERR vs position, colored by sup_status
    # plt.figure(figsize=(10, 4))
    # if "sup_status" in df.columns:
    #     m = not_nan("d_GT_super_main", "pos")
    #     for status, sub in df.loc[m].groupby("sup_status"):
    #         plt.scatter(sub["pos"], sub["d_GT_super_main"],
    #                     s=6, alpha=0.6, label=status)
    #     plt.legend(title="sup_status")
    # else:
    #     plt.scatter(df["pos"], df["d_GT_super_main"], s=3, alpha=0.4)
    # plt.axhline(0, linestyle="--")
    # plt.xlabel("Position (bp)")
    # plt.ylabel("ΔGTERR (super − main)")
    # plt.title("Δ genotype error along chromosome")
    # plt.tight_layout()
    # plt.savefig(f"{args.out_prefix}.d_GT_super_main.by_pos.png", dpi=200)
    # plt.close()

    # 6) ΔSER vs ΔGTERR scatter (are sites where supersites hurt genotyping also where phasing gets worse?)
    # plt.figure(figsize=(6, 6))
    # m = not_nan("d_SER_super_main", "d_GT_super_main")
    # if "sup_status" in df.columns:
    #     for status, sub in df.loc[m].groupby("sup_status"):
    #         plt.scatter(sub["d_GT_super_main"], sub["d_SER_super_main"],
    #                     s=6, alpha=0.6, label=status)
    #     plt.legend(title="sup_status")
    # else:
    #     plt.scatter(df.loc[m, "d_GT_super_main"],
    #                 df.loc[m, "d_SER_super_main"],
    #                 s=6, alpha=0.6)
    # plt.axhline(0, linestyle="--")
    # plt.axvline(0, linestyle="--")
    # plt.xlabel("ΔGTERR (super − main)")
    # plt.ylabel("ΔSER (super − main)")
    # plt.title("Δ switch error vs Δ genotype error")
    # plt.tight_layout()
    # plt.savefig(f"{args.out_prefix}.d_SER_vs_d_GT_super_main.png", dpi=200)
    # plt.close()


if __name__ == "__main__":
    main()
