set -e

cd /data/zhzhou/DREAM/examples/Dreicer/scripts
python generate_with_fre.py
cd /data/zhzhou/DREAM/examples/Dreicer/scripts

# 基本分析（逃逸电子密度和增长率）
python analyze_dreicer.py

# 绘制分布函数
python plot_dreicer_dist.py

# 绘制时间演化
python plot_f_time_evolution.py

# 绘制动量空间分布 (p_perp vs p_parallel) - f_hot
python plot_p_perp_parallel.py

# 绘制动量空间分布 (p_perp vs p_parallel) - f_re（逃逸电子）
python plot_fre_p_perp_parallel.py

# 综合可视化
python visual.py

python plot_logf_px_contour.py