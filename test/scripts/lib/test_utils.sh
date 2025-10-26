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

assert_same_variants() {
  local actual_bcf=$1
  local expected_vcf=$2
  if [[ ! -f $expected_vcf ]]; then
    echo "assert_same_variants: missing expected VCF $expected_vcf" >&2
    return 1
  fi
  local tmp_dir
  tmp_dir=$(mktemp -d)
  trap 'rm -rf "$tmp_dir"' RETURN
  local normalized_actual=$tmp_dir/actual.vcf
  normalize_bcf "$actual_bcf" >"$normalized_actual"
  diff -u "$expected_vcf" "$normalized_actual"
}

# New headerless VCF.gz + MD5 validation functions

extract_variants() {
  local bcf_path=$1
  local output_path=$2
  if [[ ! -f $bcf_path ]]; then
    echo "extract_variants: missing file $bcf_path" >&2
    return 1
  fi
  # Extract headerless variants and compress
  bcftools view -H -Oz "$bcf_path" -o "$output_path"
}

generate_headerless_vcf_gz() {
  local bcf_path=$1
  local output_base=$2
  if [[ ! -f $bcf_path ]]; then
    echo "generate_headerless_vcf_gz: missing file $bcf_path" >&2
    return 1
  fi
  
  local vcf_gz="${output_base}.vcf.gz"
  local md5_file="${output_base}.md5"
  
  # Extract headerless variants and compress
  extract_variants "$bcf_path" "$vcf_gz"
  
  # Generate MD5 checksum
  md5sum "$vcf_gz" | cut -d' ' -f1 > "$md5_file"
  
  echo "Generated: $vcf_gz and $md5_file"
}

assert_same_md5() {
  local actual_bcf=$1
  local expected_base=$2
  
  local expected_vcf_gz="$SCRIPT_DIR/expected/${expected_base}.vcf.gz"
  local expected_md5="$SCRIPT_DIR/expected/${expected_base}.md5"
  
  # Optional: regenerate expected when UPDATE_EXPECTED=1
  if [[ ${UPDATE_EXPECTED:-0} == 1 ]]; then
    mkdir -p "$(dirname "$expected_vcf_gz")"
    rm -f "$expected_vcf_gz" "$expected_md5"
    local output_base="$SCRIPT_DIR/expected/${expected_base}"
    generate_headerless_vcf_gz "$actual_bcf" "$output_base"
  fi
  
  if [[ ! -f $expected_vcf_gz ]]; then
    echo "assert_same_md5: missing expected VCF.gz $expected_vcf_gz" >&2
    echo "(set UPDATE_EXPECTED=1 to generate)" >&2
    return 1
  fi
  
  if [[ ! -f $expected_md5 ]]; then
    echo "assert_same_md5: missing expected MD5 $expected_md5" >&2
    echo "(set UPDATE_EXPECTED=1 to generate)" >&2
    return 1
  fi
  
  local tmp_dir
  tmp_dir=$(mktemp -d)
  trap 'rm -rf "$tmp_dir"' RETURN
  
  local actual_vcf_gz="$tmp_dir/actual.vcf.gz"
  
  # Extract variants from actual BCF
  extract_variants "$actual_bcf" "$actual_vcf_gz"
  
  # Calculate MD5 of actual output
  local actual_md5
  actual_md5=$(md5sum "$actual_vcf_gz" | cut -d' ' -f1)
  
  # Read expected MD5
  local expected_md5_value
  expected_md5_value=$(cat "$expected_md5")
  
  # Compare MD5 checksums
  if [[ "$actual_md5" == "$expected_md5_value" ]]; then
    echo "✓ MD5 validation passed: $actual_md5"
    return 0
  else
    echo "✗ MD5 validation failed:" >&2
    echo "  Expected: $expected_md5_value" >&2
    echo "  Actual:   $actual_md5" >&2
    echo "  Files: $expected_vcf_gz vs $actual_vcf_gz" >&2
    return 1
  fi
}

# Assert zero one-allele violations by parsing tools/check_oneallele output
assert_no_oneallele_violations() {
  local bcf_path=$1
  local tool_path="${TEST_DIR:-$(dirname "${BASH_SOURCE[0]}")/..}/../tools/check_oneallele"
  if [[ ! -x "$tool_path" ]]; then
    echo "assert_no_oneallele_violations: missing tool $tool_path" >&2
    return 1
  fi
  local out
  out=$("$tool_path" "$bcf_path")
  echo "$out"
  local hv
  hv=$(awk -F'\t' '{for(i=1;i<=NF;i++){if($i~ /^hap_violations=/){split($i,a,"=");print a[2]}}}' <<<"$out")
  if [[ -z "$hv" ]]; then
    echo "assert_no_oneallele_violations: could not parse output" >&2
    return 1
  fi
  if [[ "$hv" -eq 0 ]]; then
    echo "✓ One-allele violations: 0"
    return 0
  else
    echo "✗ One-allele violations: $hv" >&2
    return 1
  fi
}
