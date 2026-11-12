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
DREAM_OUTPUT='/data/zhzhou/DREAM/examples/avalanche_whistler/figures_wk40_A4_t3_s2_i0.2_T1000/quasilinear_whistler_output.h5'
GREEN_FILE='/data/zhzhou/QUADRE_series/SOFT_workspace/SOFT2/examples/GreenDream/data/green.mat'
RESULTS_DIR='/data/zhzhou/DREAM/examples/avalanche_whistler/figures_wk40_A4_t3_s2_i0.2_T1000/green_results'
SCRIPTS_DIR='/data/zhzhou/DREAM/examples/avalanche_whistler/scripts'

# === Parameters ===
RADIAL_INDEX=0
DT=0.005
T_MIN=1.0
T_MAX=3.0
# Wave injection ON windows (start_inject=2.0, cycle=0.2 s, 50% duty)
WAVE_INTERVALS='2.0,2.1;2.2,2.3;2.4,2.5;2.6,2.7;2.8,2.9'

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
    --t-min ${T_MIN} --t-max ${T_MAX} \
    --wave-intervals "${WAVE_INTERVALS}"
