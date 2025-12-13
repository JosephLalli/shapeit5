#!/usr/bin/env bash
set -euo pipefail

# Aggregate change in number of switch errors (nse) between supersites and main
# across multiple *.ser_deltas.tsv files (different seeds).
#
# For each variant present in ALL input files, this script computes:
#   - nse_super and nse_main per seed
#   - delta_nse_super_main = nse_super - nse_main per seed
#   - mean and std of delta across seeds
#   - counts of seeds with positive / negative / zero delta
#
# Usage:
#   scripts/triage_nse_multi.sh ser_deltas.seed1.tsv ser_deltas.seed2.tsv [...]
#
# The order of files defines the seed index in the output.

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <ser_deltas.tsv> <ser_deltas2.tsv> [ser_deltas3.tsv ...]" >&2
  exit 1
fi

files=("$@")

awk '
  BEGIN {
    file_idx = 0;
  }

  # Skip header in each file, track seed index
  FNR == 1 {
    file_idx++;
    next;
  }

  {
    id  = $1;
    pos = $2;
    sup = $13;  # sup_status
    vc  = $14;  # var_count

    nse_super = $3 + 0;  # nse_super
    nse_main  = $7 + 0;  # nse_main
    d         = nse_super - nse_main;

    key = id;
    meta_pos[key] = pos;
    meta_sup[key] = sup;
    meta_vc[key]  = vc;

    nseS[key, file_idx] = nse_super;
    nseM[key, file_idx] = nse_main;
    dval[key, file_idx] = d;
  }

  END {
    # Header
    printf("id\tpos\tsup_status\tvar_count");
    for (i = 1; i <= file_idx; i++) {
      printf("\tnse_super_seed%d", i);
    }
    for (i = 1; i <= file_idx; i++) {
      printf("\tnse_main_seed%d", i);
    }
    for (i = 1; i <= file_idx; i++) {
      printf("\tdelta_nse_super_main_seed%d", i);
    }
    printf("\tmean_delta_nse_super_main\tstd_delta_nse_super_main\tn_seeds_pos\tn_seeds_neg\tn_seeds_zero\n");

    # Iterate variants
    for (k in meta_pos) {
      ok = 1;
      # Require presence in all seeds
      for (i = 1; i <= file_idx; i++) {
        if (!((k, i) in dval)) {
          ok = 0;
          break;
        }
      }
      if (!ok) continue;

      # Collect deltas per seed
      sum = 0.0;
      n   = 0;
      for (i = 1; i <= file_idx; i++) {
        dtmp[i] = dval[k, i] + 0.0;
        sum += dtmp[i];
        n++;
      }
      if (n == 0) continue;

      mean = sum / n;

      sq = 0.0;
      n_pos = n_neg = n_zero = 0;
      for (i = 1; i <= file_idx; i++) {
        d = dtmp[i];
        diff = d - mean;
        sq += diff * diff;
        if (d > 0) n_pos++;
        else if (d < 0) n_neg++;
        else n_zero++;
      }
      std = (n > 1 ? sqrt(sq / n) : 0.0);

      # Row output
      printf("%s\t%s\t%s\t%s", k, meta_pos[k], meta_sup[k], meta_vc[k]);
      for (i = 1; i <= file_idx; i++) {
        printf("\t%d", nseS[k, i]);
      }
      for (i = 1; i <= file_idx; i++) {
        printf("\t%d", nseM[k, i]);
      }
      for (i = 1; i <= file_idx; i++) {
        printf("\t%d", dtmp[i]);
      }
      printf("\t%.6g\t%.6g\t%d\t%d\t%d\n", mean, std, n_pos, n_neg, n_zero);
    }
  }
' "${files[@]}"

