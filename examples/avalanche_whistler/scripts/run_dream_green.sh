#!/bin/bash
#SBATCH --job-name=dream-green
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/dream_green_%j.log
#SBATCH --error=logs/dream_green_%j.log

mkdir -p logs

# === Edit the paths here ===
DREAM_OUTPUT='/data/zhzhou/DREAM/examples/avalanche_whistler/figures_wk40A7_0.1_1.67_t2.5_avalanche_s1_i0.2/quasilinear_whistler_output.h5'
GREEN_FILE='/data/zhzhou/QUADRE_series/SOFT_workspace/SOFT2/examples/GreenDream/data/green.mat'
RESULTS_DIR='/data/zhzhou/DREAM/examples/avalanche_whistler/figures_wk40A7_0.1_1.67_t2.5_avalanche_s1_i0.2/green_results'
SCRIPTS_DIR='/data/zhzhou/DREAM/examples/avalanche_whistler/scripts'

# === Parameters ===
RADIAL_INDEX=0
DT=0.01
T_MIN=0.5
T_MAX=1.5

cd ${SCRIPTS_DIR}

python3 dream_green_processor.py \
    ${DREAM_OUTPUT} \
    ${GREEN_FILE} \
    ${RESULTS_DIR} \
    --radial-index ${RADIAL_INDEX} \
    --dt ${DT} \
    --t-min ${T_MIN} --t-max ${T_MAX}

python3 dream_se_power.py \
    ${RESULTS_DIR} \
    --t-min ${T_MIN} --t-max ${T_MAX}
