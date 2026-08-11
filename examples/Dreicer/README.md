# Dreicer 示例目录说明

此目录包含了用于验证和分析Dreicer逃逸电子生成机制的示例和相关文件。

## 目录结构

```
Dreicer/
├── README.md              # 说明文件
├── scripts/               # Python脚本和设置文件
│   ├── generateD.py       # 基本的Dreicer生成脚本
│   ├── generate_with_fre.py # 启用逃逸电子网格的生成脚本
│   ├── analyze_dreicer.py  # 分析Dreicer生成率
│   ├── plot_dreicer_dist.py # 绘制分布函数
│   ├── plot_f_time_evolution.py # 绘制时间演化
│   ├── plot_p_perp_parallel.py # 绘制p_perp-p_parallel图(f_hot)
│   ├── plot_fre_p_perp_parallel.py # 绘制p_perp-p_parallel图(f_re)
│   └── *.py               # 其他辅助脚本
├── outputs/               # 模拟输出数据
│   ├── dreicer_output.h5  # 基本Dreicer模拟结果
│   ├── dreicer_with_fre_output.h5 # 启用逃逸电子网格的模拟结果
│   └── *.h5               # 其他HDF5输出文件
├── figures/               # 生成的图表和可视化结果
│   ├── dreicer_analysis.png # 基本分析图
│   ├── *.png              # 其他PNG图像文件
│   └── *.png
└── analysis_data/         # 分析过程中产生的额外数据（暂为空）
```

## 使用说明

1. 运行基本Dreicer模拟：
   ```
   cd scripts
   python generateD.py
   ```

2. **运行启用逃逸电子网格的模拟**：
   ```
   cd scripts
   python generate_with_fre.py
   ```

3. 分析结果：
   ```
   cd scripts
   python analyze_dreicer.py
   ```

4. 生成可视化图表：
   ```
   cd scripts
   python plot_dreicer_dist.py
   python plot_f_time_evolution.py
   python plot_p_perp_parallel.py
   python plot_fre_p_perp_parallel.py
   ```

## 文件命名约定

- 以`dreicer_`开头的文件属于基础Dreicer分析
- 以`dreicer_fre_`开头的文件涉及逃逸电子分布(f_re)的分析
- 设置文件以`_settings.h5`结尾
- 输出文件以`_output.h5`结尾
- 图像文件以`.png`结尾