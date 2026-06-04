set -e

cd /data/zhzhou/DREAM/examples/avalanche_whistler/scripts

# ── Run DREAM simulation (choose one) ────────────────────────────────────
# Dreicer + avalanche kinetic simulation (with f_re grid)
# python generate_with_fre.py --a 0.4 --R 1.67 --E 0.05 --n 5e18 --T 2165 \
#   --Np-hot 100 --Np-re 200 --Nxi 40 --tMax 2.5 --Nt 2500 --output ../outputs/dreicer_with_fre_output.h5 --source on 

# Quasilinear whistler simulation
python generate_with_fre_whistler.py \
    --amplitude 6e2 --a 0.4 --R 1.67 \
    --tMax 3.0 --Nt 6000 --source off \
    --output ../outputs/quasilinear_whistler_output.h5 \
    --start-inject-time 2 --inject-cycle-duration 5 --ramp-time 0.02

# ── Select output file to analyze ────────────────────────────────────────
# DATA_DIR="../outputs/dreicer_with_fre_output.h5"
DATA_DIR="../outputs/quasilinear_whistler_output.h5"

# ── Plot directory (derived from simulation parameters) ──────────────────
# PLOT_DIR="../figures_0.01_1.67_E0.05_t2.5_no_avalanche"
PLOT_DIR="../figures_w_0.4_1.67_E0.05_t3_no_avalanche_s2_i5"

# Create plot directory if it doesn't exist
if [ ! -d "$PLOT_DIR" ]; then
    mkdir -p "$PLOT_DIR"
    echo "Created directory: $PLOT_DIR"
fi

# 基本分析（逃逸电子密度和增长率）
python analyze_dreicer.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 绘制分布函数
python plot_dreicer_dist.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 绘制时间演化
python plot_f_time_evolution.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR --time_series 2.0,2.5

# 绘制动量空间分布 (p_perp vs p_parallel) - f_hot
python plot_p_perp_parallel.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 绘制动量空间分布 (p_perp vs p_parallel) - f_re（逃逸电子）
python plot_fre_p_perp_parallel.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 综合可视化
python visual.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

python plot_logf_px_contour.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR 