#!/usr/bin/env bash
# Find the most recent commit (within a window) where the supersite expansion
# regression test still passes. This script walks commits chronologically,
# builds the test binary, runs it with deterministic sorting, and records the
# last passing SHA. It restores your original branch/commit on exit.

set -euo pipefail

WINDOW_SPEC=${WINDOW_SPEC:-"36 hours ago"}        # Override to widen/narrow search
TEST_CMD=${TEST_CMD:-"SHAPEIT5_DETERMINISTIC_SORT=1 tests/bin/test_supersite_expansion_epochs"}
TEST_TARGET=${TEST_TARGET:-"test_supersite_expansion_epochs"}
LD_LIBRARY_PATH_OVERRIDE=${LD_LIBRARY_PATH_OVERRIDE:-"$HOME/.linuxbrew/lib:/usr/local/lib:$LD_LIBRARY_PATH"}

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"

if ! git diff --quiet || ! git diff --quiet --cached; then
  echo "Working tree is dirty; stash or commit your changes first." >&2
  exit 1
fi

start_ref=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || true)
[ "$start_ref" = "HEAD" ] && start_ref=$(git rev-parse HEAD)

commits=($(git rev-list --reverse --since="$WINDOW_SPEC" HEAD))
if [ ${#commits[@]} -eq 0 ]; then
  echo "No commits found in window '$WINDOW_SPEC'." >&2
  exit 1
fi

timestamp=$(date +"%Y%m%d-%H%M%S")
log_dir="$repo_root/logs/supersite-bisect"
mkdir -p "$log_dir"
report="$log_dir/run-$timestamp.log"

last_good=""
echo "Searching commits since $WINDOW_SPEC" | tee "$report"
for sha in "${commits[@]}"; do
  echo "------------------------------------------------------------" | tee -a "$report"
  echo "Checking $sha" | tee -a "$report"
  git checkout -q "$sha"

  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH_OVERRIDE"
  build_log="$log_dir/build-$sha.log"
  run_log="$log_dir/run-$sha.log"
  rm -f "$build_log" "$run_log"

  status="FAIL"
  if make -C tests "$TEST_TARGET" >"$build_log" 2>&1; then
    if eval "$TEST_CMD" >"$run_log" 2>&1; then
      status="PASS"
      last_good="$sha"
    fi
  fi

  echo "Result: $status" | tee -a "$report"
  echo "  build log: $build_log" | tee -a "$report"
  echo "  run log  : $run_log" | tee -a "$report"
done

echo "------------------------------------------------------------" | tee -a "$report"
echo "Last passing commit: ${last_good:-NONE}" | tee -a "$report"

git checkout -q "$start_ref"
