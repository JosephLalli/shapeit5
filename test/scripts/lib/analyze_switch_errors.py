#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

COLS_SWITCH = ["id","pos","nse","nb","SER"]
COLS_TYPING = ["id","pos","errors","gt","gt_err_rate"]
COLS_SUPER  = ["id","pos","ss_idx","is_anchor","var_count","member_rank"]

def load_switch(path, label):
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=COLS_SWITCH, compression="infer")
    # Compute SER if needed
    if df["SER"].isna().all():
        df["SER"] = np.where(
            df["nb"] > 0,
            100.0 * df["nse"] / df["nb"],
            np.nan
        )
    df = df[["id","pos","SER"]].copy()
    df = df.rename(columns={"SER": f"SER_{label}"})
    return df

def load_typing(path):
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=COLS_TYPING, compression="infer")
    return df[["id","pos","gt_err_rate"]]

def load_supersite_meta(path):
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=COLS_SUPER, compression="infer")
    # derive sup_status
    df["sup_status"] = "non"
    df.loc[df["ss_idx"] >= 0, "sup_status"] = "sibling"
    df.loc[(df["ss_idx"] >= 0) & (df["is_anchor"] == 1), "sup_status"] = "anchor"
    return df

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--supersites", required=True)
    p.add_argument("--main",       required=True)
    p.add_argument("--og",         required=True)
    p.add_argument("--typing",     required=False)
    p.add_argument("--supersite-meta", required=False,
                   help="space-delimited, no header: id pos ss_idx is_anchor var_count member_rank")
    p.add_argument("--out-tsv",    default="variant_error_comparison.tsv")
    p.add_argument("--out-prefix", default="plots")
    args = p.parse_args()

    df_super = load_switch(args.supersites, "super")
    df_main  = load_switch(args.main,       "main")
    df_og    = load_switch(args.og,         "og")

    # Merge switch stats
    df = df_super.merge(df_main, on=["id","pos"], how="outer") \
                 .merge(df_og,   on=["id","pos"], how="outer")

    # Deltas
    df["d_super_main"] = df["SER_super"] - df["SER_main"]
    df["d_super_og"]   = df["SER_super"] - df["SER_og"]
    df["d_main_og"]    = df["SER_main"]  - df["SER_og"]

    # Typing
    if args.typing:
        df = df.merge(load_typing(args.typing), on=["id","pos"], how="left")

    # Supersite metadata
    if args.supersite_meta:
        df_meta = load_supersite_meta(args.supersite_meta)
        df = df.merge(df_meta, on=["id","pos"], how="left")

    # Save (space-delimited, no header)
    df.to_csv(args.out_tsv, sep=" ", index=False, header=False)
    print(f"Wrote {args.out_tsv} ({len(df)} rows)")

    # ---- PLOTS ----

    # 1) Scatter: supersites vs main
    plt.figure(figsize=(6,6))
    m = df["SER_super"].notna() & df["SER_main"].notna()
    plt.scatter(df.loc[m,"SER_main"], df.loc[m,"SER_super"], s=5, alpha=0.4)
    lim = float(np.nanmax([df["SER_main"].max(), df["SER_super"].max()]))
    plt.plot([0,lim],[0,lim],"--")
    plt.xlabel("SER_main (%)")
    plt.ylabel("SER_super (%)")
    plt.title("supersites vs main switch error")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.super_vs_main.png", dpi=200)
    plt.close()

    # 2) ΔSER histogram
    plt.figure(figsize=(6,4))
    plt.hist(df["d_super_main"].dropna(), bins=80)
    plt.axvline(0, linestyle="--")
    plt.xlabel("ΔSER (super − main)")
    plt.ylabel("count")
    plt.title("ΔSER distribution")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.delta_hist.png", dpi=200)
    plt.close()

    # 3) ΔSER vs position (colored by sup_status if present)
    plt.figure(figsize=(10,4))
    if "sup_status" in df.columns:
        m = df["d_super_main"].notna()
        for status, sub in df.loc[m].groupby("sup_status"):
            plt.scatter(sub["pos"], sub["d_super_main"],
                        s=5, alpha=0.5, label=status)
        plt.legend(title="sup_status")
    else:
        plt.scatter(df["pos"], df["d_super_main"], s=3, alpha=0.4)
    plt.axhline(0, linestyle="--")
    plt.xlabel("position")
    plt.ylabel("ΔSER (super − main)")
    plt.title("ΔSER along chromosome")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.delta_by_pos.png", dpi=200)
    plt.close()

    # 4) ΔSER vs genotype error
    if "gt_err_rate" in df.columns:
        mask = df["gt_err_rate"].notna() & df["d_super_main"].notna()
        plt.figure(figsize=(6,4))
        if "sup_status" in df.columns:
            for status, sub in df.loc[mask].groupby("sup_status"):
                plt.scatter(sub["gt_err_rate"], sub["d_super_main"],
                            s=5, alpha=0.5, label=status)
            plt.legend(title="sup_status")
        else:
            plt.scatter(df.loc[mask,"gt_err_rate"],
                        df.loc[mask,"d_super_main"],
                        s=5, alpha=0.4)
        plt.axhline(0, linestyle="--")
        plt.xlabel("genotype error rate (%)")
        plt.ylabel("ΔSER (super − main)")
        plt.title("ΔSER vs genotype error")
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}.delta_vs_gt_err.png", dpi=200)
        plt.close()

if __name__ == "__main__":
    main()
