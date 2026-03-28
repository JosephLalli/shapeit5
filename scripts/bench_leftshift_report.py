#!/usr/bin/env python3
import argparse
import gzip
import re
from pathlib import Path


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_variant_metrics(path):
    data = {}
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            vid = parts[0]
            try:
                pos = int(parts[1])
                err = int(parts[2])
                total = int(parts[3])
                rate = float(parts[4])
            except ValueError:
                continue
            data[(pos, vid)] = (err, total, rate)
    return data


def compare_typing(vt_a, vt_b, out_path):
    a = read_variant_metrics(vt_a)
    b = read_variant_metrics(vt_b)
    common = set(a.keys()) & set(b.keys())
    rows = []
    for key in common:
        a_err, a_tot, a_rate = a[key]
        b_err, b_tot, b_rate = b[key]
        delta_rate = b_rate - a_rate
        delta_err = b_err - a_err
        rows.append((key[0], key[1], a_err, a_tot, a_rate,
                     b_err, b_tot, b_rate, delta_rate, delta_err))
    rows.sort(key=lambda r: abs(r[8]), reverse=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as out:
        out.write(
            "pos\tid\t"
            "err_biallelic\ttotal_biallelic\trate_biallelic\t"
            "err_multiallelic_split\ttotal_multiallelic_split\trate_multiallelic_split\t"
            "delta_rate\tdelta_err\n"
        )
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    better = sum(1 for r in rows if r[9] < 0)
    worse = sum(1 for r in rows if r[9] > 0)
    same = len(rows) - better - worse
    summary = {
        "common_sites": len(rows),
        "only_biallelic": len(a) - len(common),
        "only_multiallelic_split": len(b) - len(common),
        "multiallelic_split_better": better,
        "multiallelic_split_worse": worse,
        "multiallelic_split_equal": same,
    }
    return summary


def compute_confusion_and_corr(a, b):
    common = set(a.keys()) & set(b.keys())
    if not common:
        return None, None, 0
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    rates_a = []
    rates_b = []
    for key in common:
        a_err, a_tot, a_rate = a[key]
        b_err, b_tot, b_rate = b[key]
        a_bad = a_err > 0
        b_bad = b_err > 0
        if a_bad and b_bad:
            tp += 1
        elif a_bad and not b_bad:
            fn += 1
        elif not a_bad and b_bad:
            fp += 1
        else:
            tn += 1
        rates_a.append(a_rate)
        rates_b.append(b_rate)
    corr = correlation(rates_a, rates_b)
    return (tn, fp, fn, tp), corr, len(common)


def correlation(xs, ys):
    if len(xs) != len(ys) or not xs:
        return None
    mean_x = sum(xs) / float(len(xs))
    mean_y = sum(ys) / float(len(ys))
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    if var_x <= 0.0 or var_y <= 0.0:
        return None
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    return cov / ((var_x ** 0.5) * (var_y ** 0.5))


def grab(pattern, text, default="NA"):
    match = re.search(pattern, text)
    return match.group(1) if match else default


def read_switch_log(path):
    if not path.exists():
        return None
    text = path.read_text(errors="ignore")
    return {
        "overlap_samples": grab(r"#Overlapping samples\s*=\s*(\d+)", text),
        "overlap_variants": grab(r"#Overlapping variants\s*=\s*(\d+)", text),
        "genotyping_errors": grab(r"#Genotyping errors\s*=\s*(\d+)", text),
        "overall_SER": grab(r"#Overall switch error rate\s*=\s*([0-9.]+)", text),
        "pure_SER": grab(r"#Overall pure switch error rate\s*=\s*([0-9.]+)", text),
    }


def parse_summary_kv(path):
    if not path.exists():
        return None
    data = {}
    with path.open("r") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) != 2:
                continue
            data[parts[0]] = parts[1]
    return data


def format_table(headers, rows):
    cols = list(zip(*([headers] + rows)))
    widths = [max(len(str(x)) for x in col) for col in cols]

    def fmt_row(row):
        return "  ".join(str(val).ljust(widths[i]) for i, val in enumerate(row))

    lines = [fmt_row(headers)]
    lines.extend(fmt_row(row) for row in rows)
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize leftshift bench outputs (switch logs + typing deltas)."
    )
    parser.add_argument("--out-dir", default="tmp_multi/bench_leftshift",
                        help="Benchmark output directory")
    parser.add_argument("--vt-biallelic",
                        help="Variant typing output for biallelic estimate")
    parser.add_argument("--vt-multiallelic-split",
                        help="Variant typing output for multiallelic split estimate")
    parser.add_argument("--delta-out",
                        help="Output path for variant typing delta TSV")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    vt_biallelic = Path(args.vt_biallelic) if args.vt_biallelic else (
        out_dir / "eval.1kgp.multi.unphased.split_bial.phased.switch.variant.typing.txt.gz"
    )
    vt_multiallelic = Path(args.vt_multiallelic_split) if args.vt_multiallelic_split else (
        out_dir / "eval.1kgp.multi.noshift.split_bial.switch.variant.typing.txt.gz"
    )
    delta_out = Path(args.delta_out) if args.delta_out else (
        out_dir / "compare.variant_typing.delta.tsv"
    )

    rows = [
        ("multi", "multi truth", "1KGP multi phased",
         out_dir / "eval.1kgp.multi.noshift.switch.log"),
        ("post-phase split multi", "bi truth",
         "1KGP multi phased (split->bi)",
         out_dir / "eval.1kgp.multi.noshift.split_bial.switch.log"),
        ("split->bi (from multi)", "bi truth",
         "1KGP split biallelic phased",
         out_dir / "eval.1kgp.multi.unphased.split_bial.phased.switch.log"),
    ]

    parsed_rows = []
    for label, truth, est, log_path in rows:
        stats = read_switch_log(log_path)
        if stats is None:
            parsed_rows.append((label, truth, est, "NA", "NA", "NA", "NA", "NA"))
            continue
        parsed_rows.append((
            label, truth, est,
            stats["overlap_samples"],
            stats["overlap_variants"],
            stats["genotyping_errors"],
            stats["overall_SER"],
            stats["pure_SER"],
        ))

    headers = [
        "comparison",
        "verification",
        "estimation",
        "overlap_samples",
        "overlap_variants",
        "genotyping_errors",
        "overall_SER",
        "pure_SER",
    ]

    print(format_table(headers, parsed_rows))

    if vt_biallelic.exists() and vt_multiallelic.exists():
        summary = compare_typing(vt_biallelic, vt_multiallelic, delta_out)
        print("")
        print("Variant typing delta summary")
        for key in [
            "common_sites",
            "only_biallelic",
            "only_multiallelic_split",
            "multiallelic_split_better",
            "multiallelic_split_worse",
            "multiallelic_split_equal",
        ]:
            print(f"{key}={summary[key]}")
        print(f"wrote {delta_out}")

        a = read_variant_metrics(vt_biallelic)
        b = read_variant_metrics(vt_multiallelic)
        conf, corr, n_common = compute_confusion_and_corr(a, b)
        print("")
        print("Genotype mismatch (typing) comparison")
        if conf is None:
            print("No common variants for confusion/correlation")
        else:
            tn, fp, fn, tp = conf
            print(f"common_variants={n_common}")
            print("confusion_matrix (rows=biallelic, cols=multiallelic_split)")
            print(f"tn={tn} fp={fp} fn={fn} tp={tp}")
            if corr is None:
                print("correlation_matrix=NA")
            else:
                print("correlation_matrix")
                print(f"1.0  {corr:.6f}")
                print(f"{corr:.6f}  1.0")

    sw_biallelic = vt_biallelic.parent / "eval.1kgp.multi.unphased.split_bial.phased.switch.variant.switch.txt.gz"
    sw_multiallelic = vt_multiallelic.parent / "eval.1kgp.multi.noshift.split_bial.switch.variant.switch.txt.gz"
    if sw_biallelic.exists() and sw_multiallelic.exists():
        a = read_variant_metrics(sw_biallelic)
        b = read_variant_metrics(sw_multiallelic)
        conf, corr, n_common = compute_confusion_and_corr(a, b)
        print("")
        print("Switch error comparison")
        if conf is None:
            print("No common variants for confusion/correlation")
        else:
            tn, fp, fn, tp = conf
            print(f"common_variants={n_common}")
            print("confusion_matrix (rows=biallelic, cols=multiallelic_split)")
            print(f"tn={tn} fp={fp} fn={fn} tp={tp}")
            if corr is None:
                print("correlation_matrix=NA")
            else:
                print("correlation_matrix")
                print(f"1.0  {corr:.6f}")
                print(f"{corr:.6f}  1.0")

    gt_summary = out_dir / "compare.gt_mismatch.biallelic_vs_multiallelic_split.summary.txt"
    gt_data = parse_summary_kv(gt_summary)
    if gt_data:
        print("")
        print("Genotype mismatch summary")
        for key in [
            "variants_common",
            "variants_only_a",
            "variants_only_b",
            "genotypes_compared",
            "genotypes_mismatched",
            "genotype_mismatch_rate",
        ]:
            if key in gt_data:
                print(f"{key}={gt_data[key]}")


if __name__ == "__main__":
    main()
