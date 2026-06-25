#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_cliper_chunks.sh \
    --cliper_obj cliper_obj.rds \
    --genes gene_list.txt \
    --gr_anno gr_anno.rds \
    --outdir cliper_chunks \
    --n_chunks 20 \
    --max_jobs 4

Optional Run_CLIPER arguments:
  --flank 500000
  --p1 0.8
  --K 5
  --n_iter 10000
  --burn_in 5000
  --alpha_conc 5
  --tau2 5000
  --rho0 0.5
  --posterior_b_cutoff 0.1
  --pip_cutoff 0.8
  --seed 2001
  --no_scale
  --scale_y
  --add_b_mu

Outputs:
  OUTDIR/cliper_chunk_XXX_of_YYY.rds
  OUTDIR/cliper_merged.rds
  OUTDIR/summary_all.csv
  OUTDIR/cliper_summary.csv
  OUTDIR/cliper_select.csv
  OUTDIR/summary_info.csv
EOF
}

CLIPER_OBJ=""
GENES=""
GR_ANNO=""
OUTDIR="cliper_chunks"
N_CHUNKS=""
MAX_JOBS=""

EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cliper_obj) CLIPER_OBJ="$2"; shift 2 ;;
    --genes) GENES="$2"; shift 2 ;;
    --gr_anno) GR_ANNO="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --n_chunks) N_CHUNKS="$2"; shift 2 ;;
    --max_jobs) MAX_JOBS="$2"; shift 2 ;;
    --flank|--p1|--K|--n_iter|--burn_in|--alpha_conc|--tau2|--rho0|--posterior_b_cutoff|--pip_cutoff|--seed)
      EXTRA_ARGS+=("$1" "$2")
      shift 2
      ;;
    --no_scale|--scale_y|--add_b_mu)
      EXTRA_ARGS+=("$1")
      shift 1
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$CLIPER_OBJ" || -z "$GENES" || -z "$GR_ANNO" || -z "$N_CHUNKS" || -z "$MAX_JOBS" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="${SCRIPT_DIR}/run_cliper_chunk.R"

mkdir -p "$OUTDIR/logs"

echo "[CLIPER] cliper_obj: $CLIPER_OBJ"
echo "[CLIPER] genes:      $GENES"
echo "[CLIPER] gr_anno:    $GR_ANNO"
echo "[CLIPER] outdir:     $OUTDIR"
echo "[CLIPER] n_chunks:   $N_CHUNKS"
echo "[CLIPER] max_jobs:   $MAX_JOBS"

for CHUNK_ID in $(seq 1 "$N_CHUNKS"); do
  LOG="${OUTDIR}/logs/chunk_${CHUNK_ID}_of_${N_CHUNKS}.log"
  echo "[CLIPER] launching chunk ${CHUNK_ID}/${N_CHUNKS}"

  Rscript "$RUNNER" run \
    --cliper_obj "$CLIPER_OBJ" \
    --genes "$GENES" \
    --gr_anno "$GR_ANNO" \
    --outdir "$OUTDIR" \
    --chunk_id "$CHUNK_ID" \
    --n_chunks "$N_CHUNKS" \
    "${EXTRA_ARGS[@]}" \
    > "$LOG" 2>&1 &

  while [[ "$(jobs -pr | wc -l)" -ge "$MAX_JOBS" ]]; do
    sleep 10
  done
done

echo "[CLIPER] waiting for chunks"
wait

echo "[CLIPER] merging chunks"
Rscript "$RUNNER" merge \
  --outdir "$OUTDIR" \
  --merged_rds "${OUTDIR}/cliper_merged.rds" \
  --summary_all_csv "${OUTDIR}/summary_all.csv" \
  --summary_csv "${OUTDIR}/cliper_summary.csv" \
  --select_csv "${OUTDIR}/cliper_select.csv" \
  --info_csv "${OUTDIR}/summary_info.csv"

echo "[CLIPER] done"
