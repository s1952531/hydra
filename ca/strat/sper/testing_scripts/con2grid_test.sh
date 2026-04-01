#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CAPS_DIR="$(cd "$SCRIPT_DIR/../caps" && pwd)"

echo "=========================================="
echo "Con2Grid Baseline Test"
echo "=========================================="
echo "Working directory: $CAPS_DIR"

cd "$CAPS_DIR"

if [[ ! -f c2g_inputs.dat ]]; then
    echo "ERROR: c2g_inputs.dat not found in $CAPS_DIR"
    exit 1
fi

if [[ ! -f c2g_outputs.dat ]]; then
    echo "ERROR: c2g_outputs.dat not found in $CAPS_DIR"
    exit 1
fi

echo "Building test_con2grid..."
make test_con2grid

echo "Running test_con2grid..."
./test_con2grid
