#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <ser_deltas.tsv> <ser_deltas2.tsv> [ser_deltas3.tsv ...]" >&2
  exit 1
fi

files=("$@")

awk '
  BEGIN {
    file_idx = 0;
  }
  FNR == 1 {
    file_idx++;
    next;
  }
  {
    id = $1;
    pos = $2;
    sup = $13;
    vc  = $14;
    d   = $15;
    key = id;
    meta_pos[key] = pos;
    meta_sup[key] = sup;
    meta_vc[key]  = vc;
    dval[key, file_idx] = d;
  }
  END {
    printf("id\tpos\tsup_status\tvar_count");
    for (i = 1; i <= file_idx; i++) {
      printf("\td_SER_super_main_seed%d", i);
    }
    printf("\tmean_d_SER_super_main\tstd_d_SER_super_main\tn_seeds_pos\tn_seeds_neg\tn_seeds_zero\n");

    for (k in meta_pos) {
      ok = 1;
      n  = 0;
      sum = 0.0;

      for (i = 1; i <= file_idx; i++) {
        if ((k, i) in dval) {
          dtmp[i] = dval[k, i] + 0.0;
        } else {
          ok = 0;
          break;
        }
      }
      if (!ok) continue;

      for (i = 1; i <= file_idx; i++) {
        d = dtmp[i];
        sum += d;
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

      printf("%s\t%s\t%s\t%s", k, meta_pos[k], meta_sup[k], meta_vc[k]);
      for (i = 1; i <= file_idx; i++) {
        printf("\t%.6g", dtmp[i]);
      }
      printf("\t%.6g\t%.6g\t%d\t%d\t%d\n", mean, std, n_pos, n_neg, n_zero);
    }
  }
' "${files[@]}"

