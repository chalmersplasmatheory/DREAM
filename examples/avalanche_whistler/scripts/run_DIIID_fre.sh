set -e

cd /data/zhzhou/DREAM/examples/avalanche_whistler/scripts
# python generate_with_fre.py
python generate_with_fre_whistler.py --amplitude 1e3 --a 0.4 --R 1.67

# Set output file path
DATA_DIR="../outputs/quasilinear_whistler_output.h5"
# DATA_DIR="../outputs/dreicer_with_fre_output.h5"

PLOT_DIR="../figures"

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
python plot_f_time_evolution.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 绘制动量空间分布 (p_perp vs p_parallel) - f_hot
python plot_p_perp_parallel.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 绘制动量空间分布 (p_perp vs p_parallel) - f_re（逃逸电子）
python plot_fre_p_perp_parallel.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

# 综合可视化
python visual.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR

python plot_logf_px_contour.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR