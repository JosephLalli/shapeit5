#!/bin/bash
set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

# Regenerate expected MD5s for the curated suite
UPDATE_EXPECTED=1 "${SCRIPT_DIR}/run.sh" --update

