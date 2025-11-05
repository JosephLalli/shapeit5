#!/bin/bash

set -euo pipefail

_prepend_unique_path() {
  local var_name=$1
  local dir_path=$2
  if [[ -d "$dir_path" ]]; then
    local current=""
    if [[ -n ${!var_name+x} ]]; then
      current="${!var_name}"
    fi
    case ":${current}:" in
      *":${dir_path}:"*) ;;
      *)
        if [[ -n "$current" ]]; then
          local new_value="${dir_path}:${current}"
          printf -v "$var_name" '%s' "$new_value"
        else
          printf -v "$var_name" '%s' "$dir_path"
        fi
        export "$var_name"
        ;;
    esac
  fi
}

# Ensure shared libraries from user and linuxbrew prefixes are visible
_prepend_unique_path LD_LIBRARY_PATH "$HOME/usr/local/lib"
_prepend_unique_path LD_LIBRARY_PATH "$HOME/.linuxbrew/lib"

normalize_bcf() {
  local bcf_path=$1
  if [[ ! -f $bcf_path ]]; then
    echo "normalize_bcf: missing file $bcf_path" >&2
    return 1
  fi
  SSH_AUTH_SOCK= bcftools view -Ov "$bcf_path" |
    grep -Ev '^##(fileDate|bcftools_view(Command|Version))='
}

# Test result reporting functions
test_pass() {
  local test_name="${1:-unknown}"
  local duration="${2:-0.000}"
  printf "%s: PASS (%.3fs)\n" "$test_name" "$duration"
}

test_fail() {
  local test_name="${1:-unknown}"
  local duration="${2:-0.000}"
  local reason="${3:-unknown failure}"
  printf "%s: FAIL (%.3fs) - %s\n" "$test_name" "$duration" "$reason"
}

# Enhanced assertion function that doesn't crash
assert_same_variants() {
  local actual_bcf=$1
  local expected_vcf=$2
  local test_name="${3:-variant_comparison}"
  local start_time="${4:-$(date +%s.%N)}"
  
  if [[ ! -f $expected_vcf ]]; then
    local duration
    duration=$(echo "$(date +%s.%N) - $start_time" | bc -l 2>/dev/null || echo "0.000")
    test_fail "$test_name" "$duration" "missing expected VCF $expected_vcf"
    return 1
  fi
  
  if [[ ! -f $actual_bcf ]]; then
    local duration
    duration=$(echo "$(date +%s.%N) - $start_time" | bc -l 2>/dev/null || echo "0.000")
    test_fail "$test_name" "$duration" "missing actual BCF $actual_bcf"
    return 1
  fi
  
  local tmp_dir
  tmp_dir=$(mktemp -d)
  trap 'rm -rf "$tmp_dir"' RETURN
  local normalized_actual=$tmp_dir/actual.vcf
  
  # Capture normalization errors
  local norm_error_file=$tmp_dir/norm_error.log
  if ! normalize_bcf "$actual_bcf" >"$normalized_actual" 2>"$norm_error_file"; then
    local duration
    duration=$(echo "$(date +%s.%N) - $start_time" | bc -l 2>/dev/null || echo "0.000")
    test_fail "$test_name" "$duration" "failed to normalize BCF: $(cat "$norm_error_file")"
    return 1
  fi
  
  # Capture diff errors
  local diff_error_file=$tmp_dir/diff_error.log
  if ! diff -u "$expected_vcf" "$normalized_actual" >"$tmp_dir/diff_output.txt" 2>"$diff_error_file"; then
    local duration
    duration=$(echo "$(date +%s.%N) - $start_time" | bc -l 2>/dev/null || echo "0.000")
    local diff_lines
    diff_lines=$(wc -l < "$tmp_dir/diff_output.txt" || echo "0")
    test_fail "$test_name" "$duration" "variants differ ($diff_lines diff lines)"
    echo "Diff output (first 20 lines):" >&2
    head -20 "$tmp_dir/diff_output.txt" >&2
    return 1
  fi
  
  local duration
  duration=$(echo "$(date +%s.%N) - $start_time" | bc -l 2>/dev/null || echo "0.000")
  test_pass "$test_name" "$duration"
  return 0
}
