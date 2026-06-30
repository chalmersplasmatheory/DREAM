#!/usr/bin/env python3

import sys
import h5py

sys.path.append('/data/zhzhou/DREAM/py/')
from DREAM import DREAMOutput

def print_h5_tree(grp, prefix=''):
    """
    递归打印HDF5文件的树状结构
    """
    for key in grp.keys():
        item = grp[key]
        if isinstance(item, h5py.Group):
            print(f"{prefix}{key}/ (Group)")
            print_h5_tree(item, prefix + '  ')
        else:
            print(f"{prefix}{key} (Dataset)")

def main():
    # 尝试加载现有的Dreicer输出文件
    filename = '../outputs/dreicer_with_fre_output.h5'
    
    try:
        # 打开HDF5文件并打印其结构
        with h5py.File(filename, 'r') as f:
            print(f"HDF5 file structure of '{filename}':")
            print("=" * 50)
            print_h5_tree(f)
        
        print("\n" + "=" * 50)
        print("Loading with DREAMOutput to check for tauEETh...")
        
        # 使用DREAMOutput加载并检查other.quantity的内容
        do = DREAMOutput(filename)
        
        # 检查是否有tauEETh
        if hasattr(do.other, 'fluid') and hasattr(do.other.fluid, 'tauEETh'):
            print("\nFound tauEETh in other.fluid!")
            tauEETh = do.other.fluid.tauEETh[:]
            print(f"Shape: {tauEETh.shape}")
            print(f"Values: {tauEETh}")
        else:
            print("\ntauEETh not found in other.fluid")
            
        # 列出所有other.fluid中的量
        if hasattr(do.other, 'fluid'):
            print("\nAvailable quantities in other.fluid:")
            for name in sorted(do.other.fluid.getQuantityNames()):
                print(f"  - {name}")
                
    except FileNotFoundError:
        print(f"File '{filename}' not found. Please run the simulation first.")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()