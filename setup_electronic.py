#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
setup_electronic.py - 生成 VASP 电子性质计算输入文件

功能：
  - 功函数 (wf): 生成 slab 静态计算输入，输出 LOCPOT
  - DOS/PDOS (dos): 生成两步法输入（SCF + NSCF）

用法：
  python3 setup_electronic.py --src CONTCAR --mode wf --vacuum 20
  python3 setup_electronic.py --src CONTCAR --mode dos --two_step

依赖：
  pip install ase

作者：STAR0418-ABC
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, Any, Optional, Tuple

# 检查 ASE
try:
    from ase.io import read, write
    from ase.calculators.vasp import Vasp
    HAS_ASE = True
except ImportError:
    HAS_ASE = False
    print("[ERROR] 需要 ASE 库: pip install ase")
    sys.exit(1)


def check_vasp_pp_path() -> Tuple[bool, str]:
    """
    检查 VASP_PP_PATH 环境变量
    
    返回: (是否存在, 路径)
    """
    pp_path = os.environ.get('VASP_PP_PATH', '')
    if pp_path and os.path.isdir(pp_path):
        return True, pp_path
    return False, pp_path


def parse_kpts(kpts_str: str) -> Tuple[int, int, int]:
    """解析 K 点字符串 '8 8 1' -> (8, 8, 1)"""
    parts = kpts_str.strip().split()
    if len(parts) != 3:
        raise ValueError(f"K 点格式错误: {kpts_str}，应为 'k1 k2 k3'")
    return tuple(int(x) for x in parts)


def detect_structure_type(atoms) -> Tuple[str, str]:
    """
    自动检测结构类型
    
    返回: (type, recommended_kpts)
        type: 'cluster' / 'slab' / 'bulk'
        recommended_kpts: 推荐的 K 点字符串
    """
    if not atoms.pbc.any():
        # 无 PBC，一定是 cluster
        return 'cluster', '1 1 1'
    
    cell = atoms.get_cell()
    positions = atoms.get_positions()
    
    # 计算每个方向的尺寸和原子分布
    cell_lengths = cell.lengths()
    
    # 检查每个方向的原子分布
    vacuum_directions = []
    for axis in range(3):
        coords = positions[:, axis]
        span = coords.max() - coords.min()
        cell_len = cell_lengths[axis]
        
        # 如果原子跨度远小于 cell 长度，可能有真空
        if cell_len > 0 and span / cell_len < 0.6:
            vacuum_directions.append(axis)
    
    if len(vacuum_directions) == 3:
        # 所有方向都有真空 -> cluster
        return 'cluster', '1 1 1'
    elif len(vacuum_directions) == 1:
        # 一个方向有真空 -> slab
        axis = vacuum_directions[0]
        if axis == 2:
            return 'slab', '8 8 1'
        elif axis == 1:
            return 'slab', '8 1 8'
        else:
            return 'slab', '1 8 8'
    elif len(vacuum_directions) == 2:
        # 两个方向有真空 -> wire/1D
        # 简化处理为 cluster
        return 'cluster', '1 1 1'
    else:
        # 无明显真空 -> bulk
        return 'bulk', '12 12 12'


def get_common_incar_settings(encut: float, ediff: float, ncore: Optional[int]) -> Dict[str, Any]:
    """
    获取通用 INCAR 设置
    
    这些是电子性质计算的基本参数
    """
    settings = {
        # 精度控制
        'prec': 'Accurate',      # 高精度模式
        'encut': encut,          # 截断能 (eV)
        'ediff': ediff,          # 电子步收敛判据
        'nelm': 200,             # 最大电子步数
        'algo': 'Normal',        # 电子步算法
        
        # 输出控制
        'lwave': False,          # 不写 WAVECAR（节省空间）
    }
    
    if ncore is not None:
        settings['ncore'] = ncore
    
    return settings


def setup_work_function(atoms, outdir: str, vacuum: float, kpts: Tuple[int, int, int],
                        gamma: bool, encut: float, ediff: float, ncore: Optional[int],
                        xc: str) -> str:
    """
    设置功函数计算
    
    功函数 Φ = V_vac - E_F
    需要输出 LOCPOT（静电势）用于提取真空电势 V_vac
    
    参数:
        atoms: ASE Atoms 对象
        outdir: 输出目录
        vacuum: 真空层厚度 (Å)
        kpts: K 点网格
        gamma: 是否 Gamma-centered
        encut: 截断能
        ediff: 收敛判据
        ncore: 并行参数
        xc: 交换关联泛函
    
    返回:
        计算目录路径
    """
    calcdir = os.path.join(outdir, 'wf_static')
    os.makedirs(calcdir, exist_ok=True)
    
    # 检查是否是 slab 结构
    cell = atoms.get_cell()
    z_length = cell[2, 2]
    positions = atoms.get_positions()
    z_min, z_max = positions[:, 2].min(), positions[:, 2].max()
    slab_thickness = z_max - z_min
    
    print(f"\n>>> 结构分析:")
    print(f"    原始 cell z 方向: {z_length:.2f} Å")
    print(f"    原子 z 范围: {z_min:.2f} ~ {z_max:.2f} Å (厚度 {slab_thickness:.2f} Å)")
    
    # 添加真空并居中
    # atoms.center(vacuum=vacuum, axis=2) 会在 z 方向两侧各加 vacuum/2 的真空
    atoms.center(vacuum=vacuum, axis=2)
    
    new_cell = atoms.get_cell()
    new_z = new_cell[2, 2]
    print(f"    添加真空后 z 方向: {new_z:.2f} Å (真空层 ~{vacuum:.1f} Å)")
    
    if slab_thickness > z_length * 0.8:
        print("[WARN] 原结构可能是 bulk 而非 slab，功函数计算需要 slab + 真空结构")
    
    # 获取通用设置
    incar_settings = get_common_incar_settings(encut, ediff, ncore)
    
    # 功函数特定设置
    incar_settings.update({
        # 静态计算
        'ibrion': -1,            # 不进行离子弛豫
        'nsw': 0,                # 离子步数为 0
        
        # 功函数关键参数
        'lvhar': True,           # LVHAR=.TRUE. 输出 LOCPOT（局域静电势）
                                 # LOCPOT 包含 Hartree 势 + 交换关联势
                                 # 用于提取真空电势 V_vac
        
        'ldipol': True,          # LDIPOL=.TRUE. 启用偶极修正
                                 # 对于不对称 slab（如表面吸附）消除周期性镜像的人工偶极
        
        'idipol': 3,             # IDIPOL=3 在 z 方向应用偶极修正
                                 # 1/2/3 分别对应 x/y/z 方向
        
        # 展宽设置
        'ismear': 0,             # ISMEAR=0 Gaussian 展宽
                                 # 适用于 slab/分子/金属表面
        'sigma': 0.05,           # SIGMA=0.05 eV 展宽宽度
        
        # 输出
        'lcharg': False,         # 不需要 CHGCAR
    })
    
    # 写入 VASP 输入
    has_pp, pp_path = check_vasp_pp_path()
    
    try:
        calc = Vasp(
            directory=calcdir,
            xc=xc,
            kpts=kpts,
            gamma=gamma,
            **incar_settings
        )
        calc.write_input(atoms)
        print(f"\n>>> 已写入 VASP 输入到: {calcdir}")
        
    except Exception as e:
        # 如果 ASE 写入失败（可能是 POTCAR 问题），手动写入
        print(f"[WARN] ASE 写入出错: {e}")
        print("[INFO] 尝试手动写入 POSCAR/INCAR/KPOINTS...")
        
        write(os.path.join(calcdir, 'POSCAR'), atoms, format='vasp')
        _write_incar_manual(calcdir, incar_settings)
        _write_kpoints_manual(calcdir, kpts, gamma)
    
    # 检查 POTCAR
    if not has_pp:
        print(f"\n[WARN] VASP_PP_PATH 未设置或不存在: {pp_path}")
        print("[INFO] 请手动准备 POTCAR 文件")
        print("[INFO] 设置方法: export VASP_PP_PATH=/path/to/potentials")
    
    return calcdir


def setup_dos(atoms, outdir: str, kpts: Tuple[int, int, int], gamma: bool,
              encut: float, ediff: float, ncore: Optional[int], xc: str,
              ismear: int, two_step: bool) -> str:
    """
    设置 DOS 计算
    
    两步法：
      1. SCF 自洽计算 → 产生 CHGCAR
      2. NSCF 非自洽计算 (ICHARG=11) → 产生 DOSCAR
    
    参数:
        atoms: ASE Atoms 对象
        outdir: 输出目录
        kpts: K 点网格
        gamma: 是否 Gamma-centered
        encut: 截断能
        ediff: 收敛判据
        ncore: 并行参数
        xc: 交换关联泛函
        ismear: 展宽方法
        two_step: 是否生成两步目录
    
    返回:
        计算目录路径
    """
    # 获取通用设置
    base_settings = get_common_incar_settings(encut, ediff, ncore)
    
    if two_step:
        # ============ 步骤 1: SCF 自洽 ============
        scf_dir = os.path.join(outdir, 'dos_scf')
        os.makedirs(scf_dir, exist_ok=True)
        
        scf_settings = base_settings.copy()
        scf_settings.update({
            'ibrion': -1,
            'nsw': 0,
            'ismear': 0,             # SCF 用 Gaussian 更稳定
            'sigma': 0.05,
            'lcharg': True,          # LCHARG=.TRUE. 输出 CHGCAR
                                     # CHGCAR 包含自洽电荷密度，供 NSCF 使用
        })
        
        _write_vasp_input(atoms, scf_dir, scf_settings, kpts, gamma, xc)
        print(f"\n>>> 已写入 SCF 输入到: {scf_dir}")
        
        # ============ 步骤 2: NSCF DOS ============
        nscf_dir = os.path.join(outdir, 'dos_nscf')
        os.makedirs(nscf_dir, exist_ok=True)
        
        nscf_settings = base_settings.copy()
        nscf_settings.update({
            'ibrion': -1,
            'nsw': 0,
            
            'icharg': 11,            # ICHARG=11 从 CHGCAR 读取电荷密度
                                     # 不进行自洽，只计算本征值和 DOS
                                     # 需要将 SCF 的 CHGCAR 拷贝到此目录
            
            'lorbit': 11,            # LORBIT=11 输出投影 DOS (PDOS)
                                     # 将 DOS 分解到各原子和轨道 (s/p/d)
                                     # 输出 DOSCAR 和 PROCAR
            
            'nedos': 3000,           # NEDOS=3000 DOS 采样点数
                                     # 更多点 = 更平滑的 DOS 曲线
            
            # 展宽设置
            # ISMEAR=-5: 四面体法，适合半导体/绝缘体，DOS 更准确
            # ISMEAR=0: Gaussian 展宽，适合金属
            'ismear': ismear,
            'sigma': 0.05 if ismear != -5 else 0.01,
            
            'lcharg': False,
        })
        
        _write_vasp_input(atoms, nscf_dir, nscf_settings, kpts, gamma, xc)
        print(f">>> 已写入 NSCF 输入到: {nscf_dir}")
        
        return nscf_dir
    
    else:
        # 单步 DOS（不推荐，但有时可用）
        dos_dir = os.path.join(outdir, 'dos_single')
        os.makedirs(dos_dir, exist_ok=True)
        
        settings = base_settings.copy()
        settings.update({
            'ibrion': -1,
            'nsw': 0,
            'lorbit': 11,
            'nedos': 3000,
            'ismear': ismear,
            'sigma': 0.05 if ismear != -5 else 0.01,
            'lcharg': False,
        })
        
        _write_vasp_input(atoms, dos_dir, settings, kpts, gamma, xc)
        print(f"\n>>> 已写入 DOS 输入到: {dos_dir}")
        
        return dos_dir


def _write_vasp_input(atoms, calcdir: str, settings: Dict, kpts: Tuple[int, int, int],
                      gamma: bool, xc: str):
    """写入 VASP 输入文件"""
    has_pp, pp_path = check_vasp_pp_path()
    
    try:
        calc = Vasp(
            directory=calcdir,
            xc=xc,
            kpts=kpts,
            gamma=gamma,
            **settings
        )
        calc.write_input(atoms)
        
    except Exception as e:
        print(f"[WARN] ASE 写入出错: {e}")
        print("[INFO] 手动写入 POSCAR/INCAR/KPOINTS...")
        
        write(os.path.join(calcdir, 'POSCAR'), atoms, format='vasp')
        _write_incar_manual(calcdir, settings)
        _write_kpoints_manual(calcdir, kpts, gamma)
    
    if not has_pp:
        potcar_path = os.path.join(calcdir, 'POTCAR')
        if not os.path.exists(potcar_path):
            print(f"[WARN] POTCAR 未生成，请手动准备: {potcar_path}")


def _write_incar_manual(calcdir: str, settings: Dict):
    """手动写入 INCAR"""
    incar_path = os.path.join(calcdir, 'INCAR')
    
    with open(incar_path, 'w') as f:
        f.write("# INCAR generated by setup_electronic.py\n\n")
        
        for key, value in settings.items():
            if isinstance(value, bool):
                val_str = '.TRUE.' if value else '.FALSE.'
            elif isinstance(value, float):
                if abs(value) < 1e-4:
                    val_str = f"{value:.0E}"
                else:
                    val_str = str(value)
            else:
                val_str = str(value)
            
            f.write(f"{key.upper()} = {val_str}\n")


def _write_kpoints_manual(calcdir: str, kpts: Tuple[int, int, int], gamma: bool):
    """手动写入 KPOINTS"""
    kpoints_path = os.path.join(calcdir, 'KPOINTS')
    
    with open(kpoints_path, 'w') as f:
        f.write("Automatic mesh\n")
        f.write("0\n")
        f.write("Gamma\n" if gamma else "Monkhorst-Pack\n")
        f.write(f"{kpts[0]} {kpts[1]} {kpts[2]}\n")
        f.write("0 0 0\n")


def print_next_steps(mode: str, calcdir: str, outdir: str, two_step: bool):
    """打印后续步骤提示"""
    print("\n" + "=" * 70)
    print("下一步操作")
    print("=" * 70)
    
    if mode == 'wf':
        print(f"""
1. 检查 POTCAR 是否存在:
   ls {calcdir}/POTCAR

2. 运行功函数计算:
   cd {calcdir}
   NP=16 EXE=vasp_std run_vasp.sh

3. 计算完成后进行后处理:
   python3 analyze_electronic.py --calcdir {calcdir} --mode wf

输出文件:
   - OUTCAR: 包含 E-fermi
   - LOCPOT: 静电势数据
""")
    
    elif mode == 'dos':
        if two_step:
            scf_dir = os.path.join(outdir, 'dos_scf')
            nscf_dir = os.path.join(outdir, 'dos_nscf')
            print(f"""
=== 两步法 DOS 计算 ===

步骤 1: 运行 SCF 自洽计算
   cd {scf_dir}
   NP=16 EXE=vasp_std run_vasp.sh

步骤 2: 拷贝 CHGCAR 到 NSCF 目录
   cp {scf_dir}/CHGCAR {nscf_dir}/

步骤 3: 运行 NSCF DOS 计算
   cd {nscf_dir}
   NP=16 EXE=vasp_std run_vasp.sh

步骤 4: 后处理
   python3 analyze_electronic.py --calcdir {nscf_dir} --mode dos

输出文件:
   - DOSCAR: DOS 数据
   - PROCAR: 投影信息 (如果 LORBIT=11)
""")
        else:
            print(f"""
运行 DOS 计算:
   cd {calcdir}
   NP=16 EXE=vasp_std run_vasp.sh

后处理:
   python3 analyze_electronic.py --calcdir {calcdir} --mode dos
""")
    
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="生成 VASP 电子性质计算输入（功函数/DOS）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 功函数计算
  python3 setup_electronic.py --src CONTCAR --mode wf --vacuum 20 --ncore 8

  # DOS 计算（两步法）
  python3 setup_electronic.py --src CONTCAR --mode dos --two_step

  # DOS 计算（金属体系）
  python3 setup_electronic.py --src POSCAR --mode dos --ismear_dos 0
        """
    )
    
    parser.add_argument("--src", required=True,
                        help="输入结构文件路径 (POSCAR/CONTCAR/cif/xyz)")
    parser.add_argument("--mode", required=True, choices=['wf', 'dos'],
                        help="计算模式: wf=功函数, dos=DOS")
    parser.add_argument("--outdir", default="calc_electronic",
                        help="输出目录 (默认: calc_electronic)")
    
    # 功函数参数
    parser.add_argument("--vacuum", type=float, default=20.0,
                        help="真空层厚度 Å (默认: 20, 仅 wf 模式)")
    
    # 通用计算参数
    parser.add_argument("--ncore", type=int, default=None,
                        help="NCORE 并行参数 (可选)")
    parser.add_argument("--encut", type=float, default=500.0,
                        help="截断能 eV (默认: 500)")
    parser.add_argument("--ediff", type=float, default=1e-6,
                        help="电子步收敛判据 (默认: 1e-6)")
    parser.add_argument("--xc", default="PBE",
                        help="交换关联泛函 (默认: PBE)")
    
    # K 点参数
    parser.add_argument("--kpts_wf", default="8 8 1",
                        help="功函数 K 点 (默认: '8 8 1')")
    parser.add_argument("--kpts_dos", default="12 12 1",
                        help="DOS K 点 (默认: '12 12 1')")
    parser.add_argument("--gamma", type=bool, default=True,
                        help="是否 Gamma-centered (默认: True)")
    
    # DOS 参数
    parser.add_argument("--ismear_dos", type=int, default=-5,
                        help="DOS 的 ISMEAR (默认: -5 四面体法; 金属用 0)")
    parser.add_argument("--two_step", action="store_true", default=True,
                        help="DOS 使用两步法 (默认: True)")
    parser.add_argument("--no_two_step", action="store_false", dest="two_step",
                        help="DOS 使用单步法")
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("setup_electronic.py - VASP 电子性质输入生成器")
    print("=" * 70)
    print(f"结构文件: {args.src}")
    print(f"计算模式: {args.mode}")
    print(f"输出目录: {args.outdir}")
    
    # 检查输入文件
    if not os.path.isfile(args.src):
        print(f"[ERROR] 输入文件不存在: {args.src}")
        sys.exit(1)
    
    # 读取结构
    print(f"\n>>> 读取结构: {args.src}")
    try:
        atoms = read(args.src)
        print(f"    原子数: {len(atoms)}")
        print(f"    元素: {set(atoms.get_chemical_symbols())}")
    except Exception as e:
        print(f"[ERROR] 读取结构失败: {e}")
        sys.exit(1)
    
    # 自动检测结构类型
    print("\n>>> 结构类型检测...")
    struct_type, recommended_kpts = detect_structure_type(atoms)
    print(f"    检测到类型: {struct_type}")
    print(f"    推荐 K 点: {recommended_kpts}")
    
    if struct_type == 'cluster':
        print("\n[INFO] 检测到 cluster 结构，将使用 Gamma-only K 点 (1 1 1)")
        if args.mode == 'wf':
            args.kpts_wf = "1 1 1"
        else:
            args.kpts_dos = "1 1 1"
    elif struct_type == 'bulk' and args.mode == 'wf':
        print("\n[WARN] 检测到 bulk 结构，功函数计算需要 slab + 真空")
        print("[INFO] 请确保结构是正确的 slab 模型")
    
    # 检查 VASP_PP_PATH
    has_pp, pp_path = check_vasp_pp_path()
    if has_pp:
        print(f"\n>>> VASP_PP_PATH: {pp_path}")
    else:
        print(f"\n[WARN] VASP_PP_PATH 未设置或无效")
        print("[INFO] 将写入 POSCAR/INCAR/KPOINTS，但 POTCAR 需手动准备")
    
    # 创建输出目录
    os.makedirs(args.outdir, exist_ok=True)
    
    # 执行设置
    if args.mode == 'wf':
        kpts = parse_kpts(args.kpts_wf)
        calcdir = setup_work_function(
            atoms, args.outdir, args.vacuum, kpts, args.gamma,
            args.encut, args.ediff, args.ncore, args.xc
        )
        print_next_steps('wf', calcdir, args.outdir, False)
        
    elif args.mode == 'dos':
        kpts = parse_kpts(args.kpts_dos)
        calcdir = setup_dos(
            atoms, args.outdir, kpts, args.gamma,
            args.encut, args.ediff, args.ncore, args.xc,
            args.ismear_dos, args.two_step
        )
        print_next_steps('dos', calcdir, args.outdir, args.two_step)


if __name__ == "__main__":
    main()

