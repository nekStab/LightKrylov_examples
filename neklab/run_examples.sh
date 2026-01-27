#!/usr/bin/env bash

set -euo pipefail

### -------------------------------
### Helper functions
### -------------------------------
error_exit() {
  echo -e "Error: $1" >&2
  exit 1
}

is_positive_integer() {
  [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

usage() {
  echo "Usage: $0 <case>"
  echo
  echo "Valid cases:"
  echo "  Newton-GMRES"
  echo "  eigs"
  echo
  echo "Example:"
  echo "  $0 Newton-GMRES"
}

AVAILABLE_CORES=$(getconf _NPROCESSORS_ONLN)

### -------------------------------
### Argument handling
### -------------------------------
OPTIONS=("Newton-GMRES" "eigs")

if [[ $# -gt 1 ]]; then
  usage && error_exit "Wrong argument."
fi

if [[ $# -eq 0 ]]; then
  echo "No example case provided. Please choose one of the following:"
  for opt in "${OPTIONS[@]}"; do
    echo "  - $opt"
  done
  exit 0
fi

CASE="$1"
if [[ "$CASE" != "Newton-GMRES" && "$CASE" != "eigs" ]]; then
  usage && error_exit "Invalid example '$CASE'."
fi

### -------------------------------
### Setup checks and cleaning
### -------------------------------
#
echo -e "\nRunning LightKrylov/neklab example $CASE\n"
#
HOME_DIR=$(pwd)
# Remove existing symlink to the example (if any)
if [[ -L "logfile" ]] || [[ -L "lightkrylov.log" ]]; then
	echo "Cleaning up log files."
	rm -vf "logfile" "lightkrylov.log"
fi
#
if [[ ! -d Nek5000 ]]; then
  echo "Nek5000 directory not found. Running Nek5000_setup.sh..."
  ./Nek5000_setup.sh
else
  echo "Nek5000 directory exists. Skipping setup."
fi

if [[ ! -d LightKrylov ]]; then
  echo "LightKrylov directory not found. Running LightKrylov_setup.sh..."
  ./LightKrylov_setup.sh
else
  echo "LightKrylov directory exists. Skipping setup."
fi

### -------------------------------
### Case execution
### -------------------------------
cd "examples/$CASE" || error_exit "Cannot enter examples/$CASE"

# Link 1cyl.re2 only if it does not already exist
if [[ -e 1cyl.re2 ]]; then
  echo "1cyl.re2 already exists — skipping symlink creation."
else
  ln -vs ../1cyl.re2 1cyl.re2
fi

# Run genmap only if 1cyl.ma2 does not already exist
if [[ -e 1cyl.ma2 ]]; then
  echo "1cyl.ma2 already exists — skipping genmap."
else
  echo "Running genmap..."
  echo -e "1cyl\n0.0001" | genmap
fi

### -------------------------------
### Number of cores
### -------------------------------
read -r -p "Number of cores to use for nek5000 [12]: " NPROC
NPROC=${NPROC:-12}

is_positive_integer "$NPROC" || error_exit "Number of cores must be a positive integer"

if (( NPROC > 12 )); then
  echo "Requested cores > 12, setting to 12"
  NPROC=12
fi

if (( NPROC > AVAILABLE_CORES )); then
  error_exit "Requested $NPROC cores, but only $AVAILABLE_CORES available"
fi

### -------------------------------
### Update SIZE
### -------------------------------
CURRENT_LPMIN=$(grep -oE 'lpmin=[0-9]+' SIZE | head -n1 | cut -d= -f2)
if (( CURRENT_LPMIN != NPROC )); then
  echo "Updating lpmin from $CURRENT_LPMIN to $NPROC in SIZE"
  sed -i -E "s/(lpmin=)[0-9]+/\1${NPROC}/" SIZE
else
  echo "lpmin already set to $NPROC — no update needed"
fi

### -------------------------------
### Build
### -------------------------------
echo "Running makeneklab..."
makeneklab

if ! grep -q "Compilation successful" build.log && ! grep -q "Nothing to be done" build.log; then
  error_exit "Compilation failed (see build.log)"
fi

### -------------------------------
### Run mode
### -------------------------------
read -r -p "Run in foreground (1) or background (0)? [0]: " RUNMODE
RUNMODE=${RUNMODE:-0}

LPMIN=$(grep -oP 'lpmin=\K[0-9]+' SIZE)

# Validate RUNMODE (1 = foreground, 0 = background)
if [[ "$RUNMODE" != "0" && "$RUNMODE" != "1" ]]; then
  error_exit "Error: RUNMODE must be 1 (foreground) or 0 (background)"
fi

echo "Starting run..."
if (( LPMIN == 1 )); then
  if [[ "$RUNMODE" == "1" ]]; then
    nek 1cyl
  else
    echo "Launching nek in serial mode in the background."
    echo -e "\nRuntime info:\n\tlogfile:         nek5000\n\tlightkrylov.log: LightKrylov\n"
    echo -e "Note: Depending on the speed of the nek5000 setup, it may take a few seconds for the files to appear."
    nekb 1cyl
  fi
else
  if [[ "$RUNMODE" == "1" ]]; then
    nekmpi 1cyl "$NPROC"
  else
    echo "Launching nek in parallel in the background."
    echo -e "\nRuntime info:\n\tlogfile:         nek5000\n\tlightkrylov.log: LightKrylov\n"
    echo -e "Note: Depending on the speed of the nek5000 setup, it may take a few seconds for the files to appear."
    nekbmpi 1cyl "$NPROC"
  fi
fi

cd "$HOME_DIR" || exit 1
ln -vs "examples/$CASE/logfile"
ln -vs "examples/$CASE/lightkrylov.log"
