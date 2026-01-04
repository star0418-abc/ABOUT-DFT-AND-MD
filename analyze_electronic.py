#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_electronic.py - VASP 电子性质后处理

功能：
  - 功函数 (wf): 从 LOCPOT/OUTCAR 计算功函数，绘制真空电势剖面
  - DOS (dos): 解析 DOSCAR，绘制 total DOS

用法：
  python3 analyze_electronic.py --calcdir calc_wf/wf_static --mode wf
  python3 analyze_electronic.py --calcdir calc_dos/dos_nscf --mode dos

依赖：
  pip install numpy matplotlib

作者：STAR0418-ABC
"""

import argparse
import os
import sys
import re
from typing import Tuple, Optional
import numpy as np

# 检查 matplotlib
try:
    import matplotlib
    matplotlib.use('Agg')  # 无头模式
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("[WARN] matplotlib 未安装，无法生成图片")
    print("[INFO] 安装: pip install matplotlib")


def parse_efermi(outcar_path: str) -> float:
    """
    从 OUTCAR 解析费米能级 E-fermi
    
    搜索行: "E-fermi :   X.XXXXX"
    
    返回: 费米能级 (eV)
    """
    if not os.path.isfile(outcar_path):
        raise FileNotFoundError(f"OUTCAR 不存在: {outcar_path}")
    
    efermi = None
    pattern = re.compile(r'E-fermi\s*:\s*([-+]?\d+\.?\d*)')
    
    with open(outcar_path, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                efermi = float(match.group(1))
    
    if efermi is None:
        raise ValueError("未在 OUTCAR 中找到 E-fermi")
    
    return efermi


def parse_locpot(locpot_path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    解析 LOCPOT 文件
    
    LOCPOT 格式:
      - 前几行是 POSCAR 头（系统名、缩放因子、晶格矢量、元素、原子数）
      - 然后是坐标（Direct 或 Cartesian）
      - 空行
      - 网格大小: nx ny nz
      - 电势数据（按 x-fast, y-medium, z-slow 顺序）
    
    返回:
        cell_z: z 方向晶胞长度
        z_coords: z 坐标数组
        potential_3d: 3D 电势数组
        planar_avg: z 方向平面平均电势
    """
    if not os.path.isfile(locpot_path):
        raise FileNotFoundError(f"LOCPOT 不存在: {locpot_path}")
    
    with open(locpot_path, 'r') as f:
        lines = f.readlines()
    
    # 解析晶格矢量
    scale = float(lines[1].strip())
    lattice = np.zeros((3, 3))
    for i in range(3):
        lattice[i] = [float(x) for x in lines[2 + i].split()]
    lattice *= scale
    
    cell_z = lattice[2, 2]
    
    # 找到网格大小行
    # 跳过 POSCAR 部分（元素名、原子数、坐标）
    idx = 5
    
    # 跳过元素名行
    if not lines[idx].strip()[0].isdigit():
        idx += 1
    
    # 原子数行
    natoms = sum(int(x) for x in lines[idx].split())
    idx += 1
    
    # 坐标类型
    if lines[idx].strip()[0] in 'SsDd':
        idx += 1
    if lines[idx].strip()[0] in 'CcDd':
        idx += 1
    
    # 跳过坐标行
    idx += natoms
    
    # 跳过空行
    while idx < len(lines) and lines[idx].strip() == '':
        idx += 1
    
    # 读取网格大小
    grid_line = lines[idx].strip().split()
    nx, ny, nz = int(grid_line[0]), int(grid_line[1]), int(grid_line[2])
    idx += 1
    
    # 读取电势数据
    potential_data = []
    for i in range(idx, len(lines)):
        line = lines[i].strip()
        if line:
            potential_data.extend([float(x) for x in line.split()])
    
    potential_data = np.array(potential_data[:nx * ny * nz])
    
    # 重塑为 3D 数组 (注意 VASP 的顺序是 x-fast)
    potential_3d = potential_data.reshape((nz, ny, nx))
    
    # 计算 z 方向平面平均
    planar_avg = np.mean(potential_3d, axis=(1, 2))
    
    # z 坐标
    z_coords = np.linspace(0, cell_z, nz, endpoint=False)
    
    return cell_z, z_coords, potential_3d, planar_avg


def estimate_vacuum_potential(z_coords: np.ndarray, planar_avg: np.ndarray,
                              fraction: float = 0.15, both_sides: bool = False) -> float:
    """
    估计真空电势 V_vac
    
    对于居中的 slab，真空区域在 z 的两端
    
    参数:
        z_coords: z 坐标
        planar_avg: 平面平均电势
        fraction: 取 z 末端的比例（默认 15%）
        both_sides: 是否取两端平均
    
    返回:
        真空电势 V_vac (eV)
    """
    nz = len(planar_avg)
    n_sample = max(1, int(nz * fraction))
    
    if both_sides:
        # 取两端的平均
        left_avg = np.mean(planar_avg[:n_sample])
        right_avg = np.mean(planar_avg[-n_sample:])
        v_vac = (left_avg + right_avg) / 2
    else:
        # 只取末端
        v_vac = np.mean(planar_avg[-n_sample:])
    
    return v_vac


def analyze_work_function(calcdir: str, both_sides: bool = False, 
                          vac_fraction: float = 0.15):
    """
    分析功函数
    
    功函数 Φ = V_vac - E_F
    
    其中:
      - V_vac: 真空电势（从 LOCPOT 的平面平均电势在真空区域的值）
      - E_F: 费米能级（从 OUTCAR 读取）
    """
    print("\n>>> 分析功函数...")
    
    outcar_path = os.path.join(calcdir, 'OUTCAR')
    locpot_path = os.path.join(calcdir, 'LOCPOT')
    
    # 解析 E-fermi
    print(f"    读取 OUTCAR: {outcar_path}")
    efermi = parse_efermi(outcar_path)
    print(f"    E-fermi = {efermi:.4f} eV")
    
    # 解析 LOCPOT
    print(f"    读取 LOCPOT: {locpot_path}")
    cell_z, z_coords, potential_3d, planar_avg = parse_locpot(locpot_path)
    print(f"    Cell z = {cell_z:.2f} Å")
    print(f"    Grid nz = {len(z_coords)}")
    
    # 估计真空电势
    v_vac = estimate_vacuum_potential(z_coords, planar_avg, vac_fraction, both_sides)
    print(f"    V_vac = {v_vac:.4f} eV (取末端 {vac_fraction*100:.0f}%)")
    
    # 计算功函数
    phi = v_vac - efermi
    print(f"    Φ = V_vac - E_F = {v_vac:.4f} - {efermi:.4f} = {phi:.4f} eV")
    
    # 保存数据
    data_path = os.path.join(calcdir, 'vacuum_potential_z.dat')
    with open(data_path, 'w') as f:
        f.write("# z (Angstrom)    V(z) (eV)\n")
        f.write(f"# E-fermi = {efermi:.6f} eV\n")
        f.write(f"# V_vac = {v_vac:.6f} eV\n")
        f.write(f"# Work function Phi = {phi:.6f} eV\n")
        for z, v in zip(z_coords, planar_avg):
            f.write(f"{z:12.6f}  {v:12.6f}\n")
    print(f"\n>>> 数据已保存: {data_path}")
    
    # 绘图
    if HAS_MPL:
        fig_path = os.path.join(calcdir, 'wf_profile.png')
        
        fig, ax = plt.subplots(figsize=(8, 5))
        
        ax.plot(z_coords, planar_avg, linewidth=1.5)
        ax.axhline(y=efermi, linestyle='--', linewidth=1, label=f'E_F = {efermi:.2f} eV')
        ax.axhline(y=v_vac, linestyle=':', linewidth=1, label=f'V_vac = {v_vac:.2f} eV')
        
        ax.set_xlabel('z (Å)')
        ax.set_ylabel('Planar Average Potential (eV)')
        ax.set_title(f'Work Function Profile (Φ = {phi:.2f} eV)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(fig_path, dpi=150)
        plt.close()
        
        print(f">>> 图已保存: {fig_path}")
    
    # 打印摘要
    print("\n" + "=" * 50)
    print("功函数计算结果")
    print("=" * 50)
    print(f"E_F    = {efermi:.4f} eV")
    print(f"V_vac  = {v_vac:.4f} eV")
    print(f"Φ      = {phi:.4f} eV")
    print("=" * 50)
    
    return phi


def parse_doscar(doscar_path: str) -> Tuple[np.ndarray, np.ndarray, float, int]:
    """
    解析 DOSCAR 文件
    
    DOSCAR 格式:
      - 第 1 行: 原子数等信息
      - 第 2-5 行: 其他信息
      - 第 6 行: Emax, Emin, NEDOS, E_F, 1.0
      - 第 7 行开始: 能量, DOS_up, (DOS_down), integrated...
    
    返回:
        energy: 能量数组 (相对 E_F)
        dos: DOS 数组
        efermi: 费米能级
        nedos: 采样点数
    """
    if not os.path.isfile(doscar_path):
        raise FileNotFoundError(f"DOSCAR 不存在: {doscar_path}")
    
    with open(doscar_path, 'r') as f:
        lines = f.readlines()
    
    # 第 6 行包含 NEDOS 和 E_F
    header = lines[5].split()
    emax = float(header[0])
    emin = float(header[1])
    nedos = int(header[2])
    efermi = float(header[3])
    
    # 读取 total DOS（从第 7 行开始，共 NEDOS 行）
    energy = []
    dos = []
    
    for i in range(6, 6 + nedos):
        parts = lines[i].split()
        energy.append(float(parts[0]))
        dos.append(float(parts[1]))  # Total DOS (up)
        # 如果是自旋极化，parts[2] 是 DOS_down
    
    energy = np.array(energy) - efermi  # 相对费米能级
    dos = np.array(dos)
    
    return energy, dos, efermi, nedos


def analyze_dos(calcdir: str):
    """
    分析 DOS
    
    解析 DOSCAR，绘制 total DOS
    """
    print("\n>>> 分析 DOS...")
    
    doscar_path = os.path.join(calcdir, 'DOSCAR')
    
    if not os.path.isfile(doscar_path):
        print(f"[WARN] DOSCAR 不存在: {doscar_path}")
        print("\n推荐使用以下工具进行 DOS 分析:")
        print("  - sumo: pip install sumo; sumo-dosplot")
        print("  - p4vasp: GUI 工具")
        print("  - VASPKIT: vaspkit -task 111")
        return
    
    # 解析 DOSCAR
    print(f"    读取 DOSCAR: {doscar_path}")
    energy, dos, efermi, nedos = parse_doscar(doscar_path)
    print(f"    E-fermi = {efermi:.4f} eV")
    print(f"    NEDOS = {nedos}")
    print(f"    能量范围: {energy.min():.2f} ~ {energy.max():.2f} eV (相对 E_F)")
    
    # 保存数据
    csv_path = os.path.join(calcdir, 'dos_total.csv')
    with open(csv_path, 'w') as f:
        f.write("# Energy (eV, relative to E_F), DOS (states/eV)\n")
        f.write(f"# E_F = {efermi:.6f} eV\n")
        for e, d in zip(energy, dos):
            f.write(f"{e:12.6f},{d:12.6f}\n")
    print(f"\n>>> 数据已保存: {csv_path}")
    
    # 绘图
    if HAS_MPL:
        fig_path = os.path.join(calcdir, 'dos.png')
        
        fig, ax = plt.subplots(figsize=(8, 5))
        
        ax.plot(energy, dos, linewidth=1)
        ax.axvline(x=0, linestyle='--', linewidth=0.8, alpha=0.7, label='E_F')
        ax.fill_between(energy, dos, where=(energy <= 0), alpha=0.3)
        
        ax.set_xlabel('Energy - E_F (eV)')
        ax.set_ylabel('DOS (states/eV)')
        ax.set_title('Total Density of States')
        ax.set_xlim(energy.min(), energy.max())
        ax.set_ylim(0, None)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(fig_path, dpi=150)
        plt.close()
        
        print(f">>> 图已保存: {fig_path}")
    
    # 打印摘要
    print("\n" + "=" * 50)
    print("DOS 分析结果")
    print("=" * 50)
    print(f"E_F = {efermi:.4f} eV")
    print(f"NEDOS = {nedos}")
    print(f"DOS @ E_F = {dos[np.abs(energy).argmin()]:.4f} states/eV")
    print("=" * 50)
    
    print("\n推荐进一步分析:")
    print("  - PDOS: 使用 sumo-dosplot 或 p4vasp")
    print("  - 带结构: 使用 sumo-bandplot")


def main():
    parser = argparse.ArgumentParser(
        description="VASP 电子性质后处理（功函数/DOS）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 功函数后处理
  python3 analyze_electronic.py --calcdir calc_wf/wf_static --mode wf

  # DOS 后处理
  python3 analyze_electronic.py --calcdir calc_dos/dos_nscf --mode dos

  # 功函数，取两端真空平均
  python3 analyze_electronic.py --calcdir calc_wf/wf_static --mode wf --both_sides
        """
    )
    
    parser.add_argument("--calcdir", required=True,
                        help="VASP 计算目录（包含 OUTCAR, LOCPOT 或 DOSCAR）")
    parser.add_argument("--mode", required=True, choices=['wf', 'dos'],
                        help="分析模式: wf=功函数, dos=DOS")
    parser.add_argument("--both_sides", action="store_true",
                        help="功函数: 取两端真空区域平均 (默认只取末端)")
    parser.add_argument("--vac_fraction", type=float, default=0.15,
                        help="功函数: 真空区域取样比例 (默认: 0.15)")
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("analyze_electronic.py - VASP 电子性质后处理")
    print("=" * 70)
    print(f"计算目录: {args.calcdir}")
    print(f"分析模式: {args.mode}")
    
    # 检查目录
    if not os.path.isdir(args.calcdir):
        print(f"[ERROR] 目录不存在: {args.calcdir}")
        sys.exit(1)
    
    # 执行分析
    if args.mode == 'wf':
        analyze_work_function(args.calcdir, args.both_sides, args.vac_fraction)
        
    elif args.mode == 'dos':
        analyze_dos(args.calcdir)


if __name__ == "__main__":
    main()

