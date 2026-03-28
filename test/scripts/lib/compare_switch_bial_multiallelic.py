#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_variant_switch(path):
    data = defaultdict(lambda: [0, 0])
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
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
            data[pos][0] += errors
            data[pos][1] += checked
    return data


def read_sample_switch(path):
    data = {}
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                errors = int(parts[1])
                checked = int(parts[2])
            except ValueError:
                continue
            data[parts[0]] = (errors, checked)
    return data


def read_site_log(path):
    if not path:
        return None
    with open_text(path) as handle:
        header = None
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                header = line.lstrip("#").split("\t")
                break
        if header is None:
            return None
        try:
            pos_idx = header.index("truth_pos")
            ref_idx = header.index("truth_ref")
            alt_idx = header.index("truth_alt")
        except ValueError:
            return None
        stats = defaultdict(lambda: {"refs": set(), "alts": set(), "rows": 0})
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= max(pos_idx, ref_idx, alt_idx):
                continue
            try:
                pos = int(parts[pos_idx])
            except ValueError:
                continue
            ref = parts[ref_idx]
            alt = parts[alt_idx]
            entry = stats[pos]
            entry["refs"].add(ref)
            entry["alts"].add(alt)
            entry["rows"] += 1
    return stats


def summarize_totals(label, data):
    total_errors = sum(v[0] for v in data.values())
    total_checked = sum(v[1] for v in data.values())
    positions = len(data)
    checked_positions = sum(1 for v in data.values() if v[1] > 0)
    ser = (total_errors * 100.0 / total_checked) if total_checked else float("nan")
    print(f"{label}: positions={positions} checked_positions={checked_positions} "
          f"errors={total_errors} checked={total_checked} SER={ser:.5f}")


def write_diffs(path, rows, header):
    if not path:
        return
    with open(path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(x) for x in row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Compare switch results between multiallelic and biallelic representations."
    )
    parser.add_argument("--multi-variant", required=True, help="multiallelic variant.switch.txt.gz")
    parser.add_argument("--bial-variant", required=True, help="biallelic variant.switch.txt.gz")
    parser.add_argument("--multi-sample", required=True, help="multiallelic sample.switch.txt.gz")
    parser.add_argument("--bial-sample", required=True, help="biallelic sample.switch.txt.gz")
    parser.add_argument("--multi-site-log", help="multiallelic site log (optional)")
    parser.add_argument("--bial-site-log", help="biallelic site log (optional)")
    parser.add_argument("--top", type=int, default=10, help="show top N diffs (default: 10)")
    parser.add_argument("--out-prefix", help="write TSV outputs with this prefix")
    args = parser.parse_args()

    multi_var = read_variant_switch(args.multi_variant)
    bial_var = read_variant_switch(args.bial_variant)
    multi_sample = read_sample_switch(args.multi_sample)
    bial_sample = read_sample_switch(args.bial_sample)

    print("== Per-position comparison (biallelic collapsed by position) ==")
    summarize_totals("multiallelic", multi_var)
    summarize_totals("biallelic", bial_var)

    all_pos = sorted(set(multi_var) | set(bial_var))
    diff_pos = []
    for pos in all_pos:
        me, mc = multi_var.get(pos, (0, 0))
        be, bc = bial_var.get(pos, (0, 0))
        if (me, mc) != (be, bc):
            diff_pos.append((pos, me, mc, be, bc, be - me, bc - mc))

    print(f"positions with differing (errors,checked): {len(diff_pos)}")

    diff_pos_sorted = sorted(
        diff_pos,
        key=lambda r: (abs(r[5]), abs(r[6]), r[0]),
        reverse=True,
    )
    for row in diff_pos_sorted[: max(args.top, 0)]:
        pos, me, mc, be, bc, de, dc = row
        print(f"  pos {pos}: multi=({me},{mc}) bial=({be},{bc}) delta=({de},{dc})")

    print("\n== Per-sample comparison ==")
    summarize_totals("multiallelic", {k: list(v) for k, v in multi_sample.items()})
    summarize_totals("biallelic", {k: list(v) for k, v in bial_sample.items()})

    all_samples = sorted(set(multi_sample) | set(bial_sample))
    diff_samples = []
    for sample in all_samples:
        me, mc = multi_sample.get(sample, (0, 0))
        be, bc = bial_sample.get(sample, (0, 0))
        if (me, mc) != (be, bc):
            diff_samples.append((sample, me, mc, be, bc, be - me, bc - mc))

    print(f"samples with differing (errors,checked): {len(diff_samples)}")

    diff_samples_sorted = sorted(
        diff_samples,
        key=lambda r: (abs(r[5]), abs(r[6]), r[0]),
        reverse=True,
    )
    for row in diff_samples_sorted[: max(args.top, 0)]:
        sample, me, mc, be, bc, de, dc = row
        print(f"  sample {sample}: multi=({me},{mc}) bial=({be},{bc}) delta=({de},{dc})")

    if args.out_prefix:
        pos_out = args.out_prefix + ".position_diffs.tsv"
        sample_out = args.out_prefix + ".sample_diffs.tsv"
        write_diffs(
            pos_out,
            diff_pos_sorted,
            ["pos", "multi_errors", "multi_checked", "bial_errors", "bial_checked", "delta_errors", "delta_checked"],
        )
        write_diffs(
            sample_out,
            diff_samples_sorted,
            ["sample", "multi_errors", "multi_checked", "bial_errors", "bial_checked", "delta_errors", "delta_checked"],
        )
        print(f"\nWrote: {pos_out}")
        print(f"Wrote: {sample_out}")

    if args.multi_site_log:
        stats = read_site_log(args.multi_site_log)
        if stats is None:
            print("\nSite log summary (multiallelic): unavailable (missing headers or parse error)")
        else:
            dup = sum(1 for v in stats.values() if v["rows"] > 1)
            ref_conflict = sum(1 for v in stats.values() if len(v["refs"]) > 1)
            alt_conflict = sum(1 for v in stats.values() if len(v["alts"]) > 1)
            print("\nSite log summary (multiallelic):")
            print(f"  positions={len(stats)} duplicate_positions={dup} ref_conflicts={ref_conflict} alt_conflicts={alt_conflict}")

    if args.bial_site_log:
        stats = read_site_log(args.bial_site_log)
        if stats is None:
            print("\nSite log summary (biallelic): unavailable (missing headers or parse error)")
        else:
            ref_conflict = sum(1 for v in stats.values() if len(v["refs"]) > 1)
            alt_multi = sum(1 for v in stats.values() if len(v["alts"]) > 1)
            print("\nSite log summary (biallelic):")
            print(f"  positions={len(stats)} positions_with_multiple_alts={alt_multi} ref_conflicts={ref_conflict}")


if __name__ == "__main__":
    main()
