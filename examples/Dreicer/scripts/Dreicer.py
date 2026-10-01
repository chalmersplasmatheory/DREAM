#!/usr/bin/env python3

import numpy as np
import sys
from pathlib import Path

# 添加DREAM路径
try:
    import DREAM
except ModuleNotFoundError:
    sys.path.append('../../py')
    import DREAM

from DREAM import DREAMSettings
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as RunawayElectrons
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions

def create_dreicer_validation_settings(E_over_Ec = 2, temperature=100, density=5e19, simulation_time=1e-4):
    """
    创建一个专门用于验证Dreicer增长率的设置
    关闭雪崩效应，只保留Dreicer过程
    
    参数:
    E_over_Ec: 电场与临界电场的比值
    temperature: 电子温度 (eV)
    density: 电子密度 (m^-3)
    simulation_time: 模拟时间 (s)
    """
    
    # 初始化设置
    ds = DREAMSettings()
    
    # 设置基本等离子体参数
    T = temperature  # Temperature (eV)
    n = density      # Density (m^-3)
    
    # 计算临界电场 (使用简单估算)
    # E_c ≈ (n/1e19) * (Te/100)^2 * 6 V/m
    Ec = (n/1e19) * (T/100)**2 * 6
    E = E_over_Ec * Ec
    
    print(f"Setting up simulation with E/Ec = {E_over_Ec}")
    print(f"Electric field: {E:.2f} V/m, Critical field: {Ec:.2f} V/m")
    print(f"Temperature: {T} eV, Density: {n:.2e} m^-3")
    
    # 设置电场
    ds.eqsys.E_field.setPrescribedData(E)
    
    # 设置温度
    ds.eqsys.T_cold.setPrescribedData(T)
    
    # 添加离子 (完全电离氘)
    ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)
    
    # 设置径向网格
    ds.radialgrid.setB0(5)           # Tesla
    ds.radialgrid.setMinorRadius(0.1) # 米
    ds.radialgrid.setNr(1)           # 只使用单个径向点
    ds.radialgrid.setWallRadius(0.15) # 米
    
    # 关闭雪崩效应 - 只关注Dreicer过程
    ds.eqsys.n_re.setAvalanche(RunawayElectrons.AVALANCHE_MODE_NEGLECT)
    
    # 启用Dreicer过程 (使用神经网络模型，这是最精确的)
    ds.eqsys.n_re.setDreicer(RunawayElectrons.DREICER_RATE_NEURAL_NETWORK)
    
    # 设置runaway electron grid
    ds.runawaygrid.setEnabled(True)
    ds.runawaygrid.setNp(10)    # 设置动量网格点数
    ds.runawaygrid.setNxi(5)    # 设置xi网格点数
    ds.runawaygrid.setPmax(1.0) # 设置最大动量
    
    # 设置热电子分布
    ds.hottailgrid.setNxi(5)   # 减少网格点数以加快计算
    ds.hottailgrid.setNp(10)
    ds.hottailgrid.setPmax(0.5)
    ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
    
    # 时间步长设置
    ds.timestep.setTmax(simulation_time)
    ds.timestep.setNt(100)
    
    # 求解器设置
    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setVerbose(False)
    
    # 输出设置
    ds.output.setTiming(True, False)
    
    return ds, E, Ec

def run_dreicer_validation(Efactor_scan = [0.5, 1.0, 1.5, 2.0, 3.0]):
    """
    运行一系列Dreicer验证模拟
    """
    print("Running Dreicer validation...")
    
    results = []
    
    for Efact in Efactor_scan:
        print(f"\n--- Running simulation with E/Ec = {Efact} ---")
        
        # 创建设置
        ds, E, Ec = create_dreicer_validation_settings(E_over_Ec=Efact)
        
        # 保存设置
        filename = f'dreicer_validation_Efac{Efact}.h5'
        settings_file = f'dreicer_settings_Efac{Efact}.h5'
        ds.save(settings_file)
        
        # 运行模拟
        try:
            do = DREAM.runiface(ds, filename, quiet=False)
            
            # 提取结果
            t = do.grid.t[:]
            nre = do.eqsys.n_re.get()[1:, 0]  # 排除初始时刻
            t = t[1:]  # 排除初始时刻
            
            # 计算增长率 (使用最后几个时间步的平均值)
            gamma_numerical = np.diff(np.log(nre)) / np.diff(t)
            gamma_avg = np.mean(gamma_numerical[-20:])  # 平均最后20个时间步的增长率
            
            results.append({
                'E_factor': Efact,
                'E_field': E,
                'E_critical': Ec,
                'gamma_numerical': gamma_avg,
                'time': t,
                'n_re': nre,
                'gamma_profile': gamma_numerical
            })
            
            print(f"E/Ec = {Efact}: Numerical growth rate = {gamma_avg:.4e} s^-1")
            
        except Exception as e:
            print(f"ERROR: Failed to run simulation with E/Ec = {Efact}")
            print(f"Error: {e}")
            continue
    
    return results

def compare_with_chen_formula(results):
    """
    与Chen教科书中的经典Dreicer公式进行比较
    """
    print("\n" + "="*60)
    print("Comparing with Chen's classic Dreicer formula:")
    print("="*60)
    
    # Chen教科书中的Dreicer公式 (近似形式)
    # γ_D ≈ ν_e * exp(-E_c/E)
    # 其中 ν_e 是电子碰撞频率
    
    for res in results:
        E_fact = res['E_factor']
        E = res['E_field']
        Ec = res['E_critical']
        gamma_num = res['gamma_numerical']
        
        # 简化的理论估计
        # 使用ν_e ≈ 1e12 * n_e(1e19) / T_e(100eV)^(3/2) s^-1
        T_norm = 100  # 参考温度 (eV)
        n_norm = 1e19 # 参考密度 (m^-3)
        nu_e = 1e12 * (res['E_critical']/6/Ec) * (T_norm/100)**(-1.5) 
        
        # 经典Dreicer指数因子
        chen_gamma = nu_e * np.exp(-Ec/E)
        
        ratio = gamma_num / chen_gamma
        
        print(f"E/Ec = {E_fact:4.1f} | Numerical: {gamma_num:.2e} | Chen formula: {chen_gamma:.2e} | Ratio: {ratio:.2f}")

def main():
    """
    主函数
    """
    print("DREAM Dreicer Validation Script")
    print("===============================")
    
    # 运行验证
    results = run_dreicer_validation()
    
    if len(results) > 0:
        # 与经典公式比较
        compare_with_chen_formula(results)
        
        print("\nValidation completed.")
        print("Results saved to dreicer_validation_Efac*.h5 files")
    else:
        print("ERROR: No simulations completed successfully.")

if __name__ == "__main__":
    main()