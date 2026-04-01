#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CAPS_DIR="$(cd "$SCRIPT_DIR/../caps" && pwd)"

echo "=========================================="
echo "UGrid2Con Baseline Test"
echo "=========================================="
echo "Working directory: $CAPS_DIR"

cd "$CAPS_DIR"

PARAM_FILE="$CAPS_DIR/parameters.f90"
PARAM_BACKUP=""

restore_parameters() {
  if [[ -n "$PARAM_BACKUP" && -f "$PARAM_BACKUP" ]]; then
    mv "$PARAM_BACKUP" "$PARAM_FILE"
  fi
}

trap restore_parameters EXIT

if [[ ! -f ug2c_inputs.dat ]]; then
  echo "ERROR: ug2c_inputs.dat not found in $CAPS_DIR"
  exit 1
fi

if [[ ! -f ug2c_outputs.dat ]]; then
  echo "ERROR: ug2c_outputs.dat not found in $CAPS_DIR"
  exit 1
fi

if grep -Eq '\b(N_X|L_X|Y_MIN|Y_MAX|T_SIM|T_GSAVE|T_CSAVE|N_CONTB|N_CONTZ|U_REF|N_NU|PRE_DISS)\b' "$PARAM_FILE"; then
  PARAM_BACKUP="$(mktemp "$CAPS_DIR/parameters.f90.bak.XXXXXX")"
  cp "$PARAM_FILE" "$PARAM_BACKUP"

  cat > "$PARAM_FILE" <<'EOF'
module parameters

! This module contains all the modifiable parameters for
! the suite of caps f90 files.

 !Domain grid dimensions:
integer,parameter:: nx=512,ny=64

 !Domain width in x and limits in y:
double precision,parameter:: ellx=51200.d0
double precision,parameter:: ymin=0.d0,ymax=6400.d0

 !Simulation duration and data save interval:
double precision,parameter:: tsim=900.d0,tgsave=9.d0,tcsave=90.d0

 !Number of contours used for representing buoyancy and vorticity:
integer,parameter:: ncontb=100,ncontz=20

 !Reference translational velocity (added to u):
double precision,parameter:: uref=0.d0

 !(Hyper-)viscosity parameters (for residual vorticity only):
integer,parameter:: nnu=3
double precision,parameter:: prediss=10.d0

end module
EOF

  echo "Applied temporary concrete test parameters in $PARAM_FILE"
fi

echo "Building test_ugrid2con..."
make test_ugrid2con

echo "Running test_ugrid2con..."
./test_ugrid2con
