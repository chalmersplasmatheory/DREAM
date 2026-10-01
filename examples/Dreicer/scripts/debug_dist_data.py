#!/usr/bin/env python3
#
# Debug distribution function data
#

import numpy as np
import sys

sys.path.append('../../py')

from DREAM import *
from DREAM.DREAMOutput import DREAMOutput

# 加载数据
do = DREAMOutput('../outputs/dreicer_with_fre_output.h5')

# 获取分布函数数据
f_hot = do.eqsys.f_hot[:]  # 所有时间步
p = do.grid.hottail.p[:]
xi = do.grid.hottail.xi[:]

print("Data shapes:")
print(f"f_hot shape: {f_hot.shape}")
print(f"p shape: {p.shape}")
print(f"xi shape: {xi.shape}")
print(f"Time points: {do.grid.t[:].shape}")

print("\nTime points (in microseconds):")
for i, t in enumerate(do.grid.t[:]):
    print(f"  t[{i}] = {t*1e6:.3f} μs")

print("\nDistribution function values at different times (first radial point, first xi, all p):")
for i in range(0, f_hot.shape[0], max(1, f_hot.shape[0]//5)):  # Show at most 5 time points
    vals = f_hot[i, 0, 0, :]  # time, radius, xi, p
    print(f"  t[{i}] (t={do.grid.t[i]*1e6:.3f} μs): min={np.min(vals):.2e}, max={np.max(vals):.2e}, mean={np.mean(vals):.2e}")

print("\nDistribution function values at different times (first radial point, last xi, all p):")
for i in range(0, f_hot.shape[0], max(1, f_hot.shape[0]//5)):  # Show at most 5 time points
    vals = f_hot[i, 0, -1, :]  # time, radius, xi, p
    print(f"  t[{i}] (t={do.grid.t[i]*1e6:.3f} μs): min={np.min(vals):.2e}, max={np.max(vals):.2e}, mean={np.mean(vals):.2e}")

print("\nDistribution function values at different xi (at last time point, first radial point, all p):")
for i in range(0, f_hot.shape[2], max(1, f_hot.shape[2]//5)):  # Show at most 5 xi points
    vals = f_hot[-1, 0, i, :]  # time, radius, xi, p
    print(f"  xi[{i}] (xi={xi[i]:.3f}): min={np.min(vals):.2e}, max={np.max(vals):.2e}, mean={np.mean(vals):.2e}")

# Check if there are significant changes in the distribution
print("\nChecking for significant changes in distribution over time:")
initial_dist = f_hot[0, 0, -1, :]  # t=0, first radius, last xi
final_dist = f_hot[-1, 0, -1, :]  # t=final, first radius, last xi

diff = np.abs(final_dist - initial_dist)
print(f"Max absolute difference: {np.max(diff):.2e}")
print(f"Mean absolute difference: {np.mean(diff):.2e}")

# Close the file
do.close()