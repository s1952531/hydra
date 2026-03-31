#!/bin/bash
# con2grid_test.sh - Build and run con2grid baseline verification
# 
# Usage: ./con2grid_test.sh
# 
# This script compiles the con2grid test driver and runs baseline verification.
# Baseline files (c2g_inputs.dat, c2g_outputs.dat) must exist in the caps directory.

set -e  # Exit on any error

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Navigate to caps directory (one level up)
CAPS_DIR="$SCRIPT_DIR/../caps"

echo "=========================================="
echo "Con2Grid Baseline Test"
echo "=========================================="
echo ""
echo "Working directory: $CAPS_DIR"
cd "$CAPS_DIR"

echo "Building test driver..."
make clean test_con2grid

echo ""
echo "Checking for baseline files..."
if [ ! -f c2g_inputs.dat ]; then
    echo "ERROR: c2g_inputs.dat not found in $CAPS_DIR"
    echo "  → Run 'caps' simulation with log_con2grid=.true. to generate baselines"
    exit 1
fi

if [ ! -f c2g_outputs.dat ]; then
    echo "ERROR: c2g_outputs.dat not found in $CAPS_DIR"
    exit 1
fi

echo "Found baseline files. Running verification..."
echo ""
./test_con2grid

exit_code=$?
echo ""
if [ $exit_code -eq 0 ]; then
    echo "=========================================="
    echo "Test driver completed successfully"
    echo "=========================================="
else
    echo "=========================================="
    echo "Test driver failed (exit code: $exit_code)"
    echo "=========================================="
fi

exit $exit_code
