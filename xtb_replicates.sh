#!/bin/bash
set -euo pipefail

INPUT=""
#SCRIPT="${SCRIPT:-$HOME/xtb_md/rjm_xtb_md.py}"  # usage on HPC
SCRIPT="${SCRIPT:-./rjm_xtb_md.py}" # usage locally
NREP=10
SEED=12345
RUNNAME="job"
OUTNAME="pre_trajectory.xyz"

PY_ARGS=()

usage() {
    cat <<EOF
Usage: $0 -i input.xyz [wrapper options] [python options]

Wrapper options:
  -i FILE         Input xyz
  -script FILE    Path to rjm_xtb_md.py
  -nrep INT       Number of replicates
  -seed INT       Master seed
  -runname STR    Base run name
  -outname FILE   Initial xyz filename
  -h              Show help

All other options are passed to rjm_xtb_md.py. These are detailed below:

  -n              Number of molecules; type=int; default=2
  -prefix         Base name for generated files; default=None
  -mode           Choose from box or cluster; default = cluster
  -box            Box dimensions, for mode = box; type=float; default=50.0
  -dist           Distance between molecules; type=float; default=5.0
  -m              Method to use, default="GFN2-XTB"
  -t              Temperature (K); type=float, default=300.0
  -ts             Timestep (fs); type=float; default=1.0
  -df             Frequency to write to trajectory file (steps); type=int, default=50
  -ns             Number of steps; type=int, default=2000
  -shake          SHAKE algorith mode; type=int,; default=2
  -sccacc         Accuracy; type=float; default=2.0
  -velo           Do not generate random velocities
  -nprocs         Number of cpu cores to use; type=int, default=4

EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT="$2"; shift 2 ;;
        -script) SCRIPT="$2"; shift 2 ;;
        -nrep) NREP="$2"; shift 2 ;;
        -seed) SEED="$2"; shift 2 ;;
        -runname) RUNNAME="$2"; shift 2 ;;
        -outname) OUTNAME="$2"; shift 2 ;;
        -h) usage; exit 0 ;;
        *)
            PY_ARGS+=("$1")
            shift
            if [[ $# -gt 0 && "$1" != -* ]]; then
                PY_ARGS+=("$1")
                shift
            fi
            ;;
    esac
done

[[ -n "$INPUT" ]] || { echo "Error: input xyz required"; exit 1; }

INPUT=$(realpath "$INPUT")
SCRIPT=$(realpath "$SCRIPT")

for i in $(seq 0 $((NREP-1))); do
    repdir=$(printf "rep_%03d" "$i")
    mkdir -p "$repdir"

    (
        cd "$repdir"
        python "$SCRIPT" \
            -i "$INPUT" \
            -seed "$SEED" \
            -rep "$i" \
            -o "$OUTNAME" \
            -prefix "$RUNNAME" \
            "${PY_ARGS[@]}"
    )
done

