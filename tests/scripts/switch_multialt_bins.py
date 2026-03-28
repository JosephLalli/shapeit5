#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(info_field):
    info = {}
    for item in info_field.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def parse_alt_mafs(info, n_alts):
    if n_alts <= 0:
        return []
    ac_list = []
    an = None
    if "AN" in info:
        try:
            an = int(info["AN"])
        except ValueError:
            an = None
    if "AC" in info:
        ac_list = [int(x) for x in info["AC"].split(",") if x != "."]
    elif "AF" in info:
        if an is None:
            return []
        af_list = [float(x) for x in info["AF"].split(",") if x != "."]
        ac_list = [int(round(af * an)) for af in af_list]
    if an is None or an <= 0:
        return []
    if len(ac_list) < n_alts:
        ac_list.extend([0] * (n_alts - len(ac_list)))
    mafs = []
    for ac in ac_list[:n_alts]:
        af = float(ac) / float(an)
        maf = min(af, 1.0 - af)
        mafs.append(maf)
    return mafs


def load_switch_variants(path):
    by_pos = {}
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                pos = int(parts[1])
                errors = int(parts[2])
                checked = int(parts[3])
            except ValueError:
                continue
            if pos not in by_pos:
                by_pos[pos] = [0, 0]
            by_pos[pos][0] += errors
            by_pos[pos][1] += checked
    return by_pos


def build_bins(edges):
    bins = []
    for i in range(len(edges) - 1):
        bins.append((edges[i], edges[i + 1]))
    bins.append((edges[-1], None))
    return bins


def bin_label(low, high):
    if high is None:
        return f">={low}"
    return f"{low}-{high}"


def pick_maf(mafs, mode):
    if not mafs:
        return None
    if mode == "min":
        return min(mafs)
    if mode == "max":
        return max(mafs)
    return sum(mafs) / float(len(mafs))


def main():
    parser = argparse.ArgumentParser(
        description="Bin switch error rates by biallelic vs multiallelic and per-ALT MAF."
    )
    parser.add_argument("--switch-variant", required=True,
                        help="switch output .variant.switch.txt[.gz]")
    parser.add_argument("--frequency-vcf", required=True,
                        help="VCF/VCF.gz with AC/AN (or AF) for MAF bins")
    parser.add_argument("--out", required=True, help="Output CSV path")
    parser.add_argument("--alt-maf-mode", choices=["min", "max", "mean"], default="min",
                        help="How to summarize per-ALT MAF for multiallelic sites")
    parser.add_argument("--bins", default="0,0.001,0.01,0.05,0.5",
                        help="Comma-separated MAF bin edges (ascending)")
    args = parser.parse_args()

    edges = [float(x) for x in args.bins.split(",") if x != ""]
    if len(edges) < 2:
        raise SystemExit("Need at least two bin edges")
    bins = build_bins(edges)

    switch_by_pos = load_switch_variants(args.switch_variant)
    stats = defaultdict(lambda: {"sites": 0, "errors": 0, "checked": 0})
    skipped = 0
    missing_switch = 0

    with open_text(args.frequency_vcf) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            alt_field = parts[4]
            if alt_field == ".":
                continue
            alts = alt_field.split(",")
            n_alts = len(alts)
            info = parse_info(parts[7])
            mafs = parse_alt_mafs(info, n_alts)
            maf_value = pick_maf(mafs, args.alt_maf_mode)
            if maf_value is None:
                skipped += 1
                continue
            if pos not in switch_by_pos:
                missing_switch += 1
                continue
            errors, checked = switch_by_pos[pos]
            if checked == 0:
                continue
            category = "biallelic" if n_alts == 1 else "multiallelic"
            label = None
            for low, high in bins:
                if high is None:
                    if maf_value >= low:
                        label = bin_label(low, high)
                        break
                else:
                    if low <= maf_value < high:
                        label = bin_label(low, high)
                        break
            if label is None:
                skipped += 1
                continue
            key = (category, label)
            stats[key]["sites"] += 1
            stats[key]["errors"] += errors
            stats[key]["checked"] += checked

    with open(args.out, "w") as out:
        out.write("category,maf_bin,sites,errors,checked,error_rate\n")
        for key in sorted(stats.keys()):
            category, label = key
            row = stats[key]
            rate = (row["errors"] / row["checked"]) if row["checked"] else 0.0
            out.write(f"{category},{label},{row['sites']},{row['errors']},{row['checked']},{rate:.6f}\n")

    if skipped or missing_switch:
        print(f"Skipped variants (missing MAF bin): {skipped}")
        print(f"Skipped variants (missing switch data): {missing_switch}")


if __name__ == "__main__":
    main()
