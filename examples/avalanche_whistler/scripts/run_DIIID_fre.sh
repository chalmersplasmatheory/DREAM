#!/bin/bash
#SBATCH --job-name=DREAM_whistler_sim   # 作业名称，便于在队列中识别
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
set -e

cd /data/zhzhou/DREAM/examples/avalanche_whistler/scripts

# ── Run DREAM simulation (choose one) ────────────────────────────────────
# Dreicer + avalanche kinetic simulation (with f_re grid)
# python generate_with_fre.py --a 0.4 --R 1.67 --E 0.05 --n 5e18 --T 2165 \
#   --Np-hot 100 --Np-re 200 --Nxi 40 --tMax 2.5 --Nt 2500 --output ../outputs/dreicer_with_fre_output.h5 --source on 


# ── Select output file to analyze ────────────────────────────────────────
PLOT_DIR="../figures_w_0.4_1.67_E0.05_t2.5_no_avalanche_s1_i0.2"
DATA_DIR="${PLOT_DIR}/quasilinear_whistler_output.h5"
# Create plot directory if it doesn't exist
if [ ! -d "$PLOT_DIR" ]; then
    mkdir -p "$PLOT_DIR"
    echo "Created directory: $PLOT_DIR"
fi
# Quasilinear whistler simulation
python generate_with_fre_whistler.py \
    --amplitude 6e2 --a 0.4 --R 1.67 \
    --tMax 2.5 --Nt 5000 --source off \
    --output $DATA_DIR \
    --start-inject-time 1 --inject-cycle-duration 0.2 --ramp-time 0.02 




# 基本分析（逃逸电子密度和增长率）
python analyze_dreicer.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 绘制分布函数
python plot_dreicer_dist.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 绘制时间演化
python plot_f_time_evolution.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR --time_series 1.0,1.5

# 绘制动量空间分布 (p_perp vs p_parallel) - f_hot
python plot_p_perp_parallel.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 绘制动量空间分布 (p_perp vs p_parallel) - f_re（逃逸电子）
python plot_fre_p_perp_parallel.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 综合可视化
python visual.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

python plot_logf_px_contour.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR 