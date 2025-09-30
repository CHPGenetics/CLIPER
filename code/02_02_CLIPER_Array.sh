#!/bin/bash
#SBATCH --account=rduerr
#SBATCH --job-name=CLIPER
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --output=/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Code/02_CLIPER/Slurm_Out/%j.out
#SBATCH --error=/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Code/02_CLIPER/Slurm_Out/%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=32000
#SBATCH --time=48:00:00
#SBATCH --array=0-249

set -euo pipefail

module purge || true
module load r/4.5.0

GENES_FILE="${GENES_FILE:-/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Data/genes_2000.txt}"
CELLTYPES_FILE="${CELLTYPES_FILE:-/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Data/celltype.txt}"
RSCRIPT="${RSCRIPT:-/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Code/02_CLIPER/02_01_PCGS_Command.R}"

GENES_PER_TASK="${GENES_PER_TASK:-8}"

task_id=${SLURM_ARRAY_TASK_ID:-0}
N_GENES=$(wc -l < "$GENES_FILE")
start_line=$(( task_id * GENES_PER_TASK + 1 ))
end_line=$(( start_line + GENES_PER_TASK - 1 ))
if (( end_line > N_GENES )); then end_line=$N_GENES; fi

echo "Processing gene lines ${start_line}-${end_line} out of ${N_GENES}"

for gene_line in $(seq "$start_line" "$end_line"); do
  gene=$(sed -n "${gene_line}p" "$GENES_FILE" | tr -d '\r')
  if [[ -z "${gene}" ]]; then
    echo "Empty gene at line ${gene_line}, skip."
    continue
  fi

  for CELLTYPE_LINE in $(seq 1 5); do
    celltype=$(sed -n "${CELLTYPE_LINE}p" "$CELLTYPES_FILE" | tr -d '\r')
    if [[ -z "${celltype}" ]]; then
      echo "Empty celltype at line ${CELLTYPE_LINE}, skip."
      continue
    fi

    OUTDIR1="/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Result/01_CLIPER/02_Output/wo_b_mu0/No_Center/${celltype}/"
    mkdir -p "$OUTDIR1"
    OUTDIR2="/ix1/wchen/Shiyue/Projects/2025_01_CE_Reg_P2G/02_P2G/01_DOGMA/Result/01_CLIPER/02_Output/wo_b_mu0/Center/${celltype}/"
    mkdir -p "$OUTDIR2"

    for stim in 'Act_IL1B_IL23'; do
      for cell in 'meta'; do
        for atac in 'peak'; do
          for rna_norm in 'delta' 'deseq'; do
            for atac_norm in 'delta' 'deseq' 'gcfq'; do

              batch="${cell}_rna_${rna_norm}_${atac}_${atac_norm}"
              echo "[`date '+%F %T'`] gene_line=${gene_line} gene='${gene}' stim=${stim} ct_line=${CELLTYPE_LINE} ct='${celltype}' batch=${batch}"

              Rscript "$RSCRIPT" "$gene" "$celltype" 5 10000 5000 5000 0.5 0.5 FALSE "$OUTDIR1" "$batch" FALSE
              Rscript "$RSCRIPT" "$gene" "$celltype" 5 10000 5000 5000 0.5 0.5 FALSE "$OUTDIR2" "$batch" TRUE

            done
          done
        done
      done
    done
  done
done

echo "Task ${task_id} done."
