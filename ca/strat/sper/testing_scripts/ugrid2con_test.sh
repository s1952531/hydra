#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CAPS_DIR="$(cd "$SCRIPT_DIR/../caps" && pwd)"

echo "=========================================="
echo "UGrid2Con Baseline Test"
echo "=========================================="
echo "Working directory: $CAPS_DIR"

cd "$CAPS_DIR"

if [[ ! -f ug2c_inputs.dat ]]; then
  echo "ERROR: ug2c_inputs.dat not found in $CAPS_DIR"
  exit 1
fi

if [[ ! -f ug2c_outputs.dat ]]; then
  echo "ERROR: ug2c_outputs.dat not found in $CAPS_DIR"
  exit 1
fi

echo "Building test_ugrid2con..."
make test_ugrid2con

echo "Running test_ugrid2con..."
./test_ugrid2con
