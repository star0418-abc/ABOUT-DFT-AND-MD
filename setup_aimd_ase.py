#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
setup_aimd_ase.py - 从大体系结构切割 AIMD cluster 并生成 VASP 输入

用法:
  python3 setup_aimd_ase.py --src equilibrated.pdb --center_atom Li --radius 8 --temp 350 --outdir aimd_Li8A
  python3 setup_aimd_ase.py --src equilibrated.gro --center_atom 0 --radius 10 --selection molecule

功能:
  - 从大体系（GROMACS/Packmol 输出）切割局部量子区域
  - 支持 PBC 最小镜像距离
  - 支持 sphere/molecule 两种选择模式
  - 生成 VASP AIMD 输入文件（POSCAR/INCAR/KPOINTS/POTCAR）

依赖:
  pip install ase numpy

作者: STAR0418-ABC
"""

from __future__ import annotations

import argparse
import os
import sys
import shutil
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Set
from datetime import datetime

import numpy as np

# 检查 ASE
try:
    from ase import Atoms
    from ase.io import read, write
    from ase.geometry import get_distances
    HAS_ASE = True
except ImportError:
    HAS_ASE = False
    print("[ERROR] 需要 ASE 库: pip install ase")
    sys.exit(1)


# ============================================================================
# 常量与配置
# ============================================================================

# 离子电荷估计表（粗略，仅用于警告）
CHARGE_MAP: Dict[str, int] = {
    # 阳离子
    'Li': +1, 'Na': +1, 'K': +1, 'Rb': +1, 'Cs': +1,
    'Mg': +2, 'Ca': +2, 'Sr': +2, 'Ba': +2,
    'Zn': +2, 'Cu': +2, 'Fe': +2, 'Co': +2, 'Ni': +2,
    'Al': +3, 'Fe3': +3,
    # 阴离子
    'F': -1, 'Cl': -1, 'Br': -1, 'I': -1,
    'O2': -2, 'S2': -2,
    # 有机元素（通常中性）
    'C': 0, 'H': 0, 'O': 0, 'N': 0, 'S': 0, 'P': 0, 'Si': 0, 'B': 0,
}


# ============================================================================
# 辅助函数
# ============================================================================

def parse_center_atom(
    center_str: str,
    atoms: Atoms,
    one_based: bool = False
) -> int:
    """
    解析中心原子参数
    
    支持:
      - 整数索引（0-based 或 1-based）
      - 元素符号（如 'Li'，返回第一个匹配的原子索引）
    
    参数:
        center_str: 用户输入的中心原子字符串
        atoms: ASE Atoms 对象
        one_based: 是否按 1-based 解释数字索引
    
    返回:
        中心原子索引（0-based）
    
    异常:
        ValueError: 找不到指定的原子
    """
    n_atoms = len(atoms)
    symbols = atoms.get_chemical_symbols()
    
    # 尝试解析为整数
    try:
        idx = int(center_str)
        
        # 处理 1-based
        if one_based:
            idx -= 1
        
        # 智能检测：如果索引等于原子数，可能是 1-based
        if idx == n_atoms:
            print(f"[WARN] 索引 {center_str} 等于原子总数 {n_atoms}，可能是 1-based")
            print(f"[INFO] 尝试按 1-based 解释为 {idx - 1}")
            idx -= 1
        
        if idx < 0 or idx >= n_atoms:
            raise ValueError(f"原子索引 {idx} 超出范围 [0, {n_atoms - 1}]")
        
        return idx
        
    except ValueError:
        pass
    
    # 作为元素符号解析
    element = center_str.strip()
    for i, sym in enumerate(symbols):
        if sym == element:
            print(f"[INFO] 找到第一个 {element} 原子，索引 {i}")
            return i
    
    raise ValueError(f"未找到元素 '{element}'，可用元素: {set(symbols)}")


def select_indices_within_radius_mic(
    atoms: Atoms,
    center_idx: int,
    radius: float,
    use_mic: bool = True
) -> np.ndarray:
    """
    选择中心原子周围指定半径内的所有原子索引
    
    使用最小镜像距离 (MIC) 处理周期性边界条件
    
    参数:
        atoms: ASE Atoms 对象
        center_idx: 中心原子索引
        radius: 切割半径 (Å)
        use_mic: 是否使用最小镜像距离
    
    返回:
        选中的原子索引数组
    """
    n_atoms = len(atoms)
    positions = atoms.get_positions()
    center_pos = positions[center_idx]
    
    if use_mic and atoms.pbc.any() and atoms.cell.any():
        # 使用 ASE 的 get_distances 计算最小镜像距离
        try:
            # 计算中心原子到所有原子的距离
            _, distances = get_distances(
                [center_pos],
                positions,
                cell=atoms.cell,
                pbc=atoms.pbc
            )
            distances = distances[0]  # shape: (n_atoms,)
        except Exception as e:
            print(f"[WARN] MIC 距离计算失败: {e}")
            print("[INFO] 回退到普通欧氏距离")
            distances = np.linalg.norm(positions - center_pos, axis=1)
    else:
        # 普通欧氏距离
        distances = np.linalg.norm(positions - center_pos, axis=1)
    
    # 选择半径内的原子（包括中心原子）
    indices = np.where(distances <= radius)[0]
    
    return indices


def get_molecule_indices(atoms: Atoms) -> Optional[List[Set[int]]]:
    """
    获取分子分组信息
    
    尝试从 atoms.arrays 中读取 residue/molecule 信息
    
    返回:
        分子列表，每个分子是一个原子索引集合；如果无法获取则返回 None
    """
    # 尝试从 arrays 获取分子信息
    arrays = atoms.arrays if hasattr(atoms, 'arrays') else {}
    
    # 常见的分子标识字段
    mol_keys = ['residuenumbers', 'molid', 'mol', 'residue', 'resid']
    
    mol_ids = None
    for key in mol_keys:
        if key in arrays:
            mol_ids = arrays[key]
            print(f"[INFO] 使用 '{key}' 字段进行分子分组")
            break
    
    if mol_ids is None:
        return None
    
    # 构建分子列表
    unique_mols = np.unique(mol_ids)
    molecules = []
    for mol_id in unique_mols:
        indices = set(np.where(mol_ids == mol_id)[0])
        molecules.append(indices)
    
    return molecules


def expand_to_molecules_if_possible(
    atoms: Atoms,
    selected_indices: np.ndarray,
    selection_mode: str
) -> np.ndarray:
    """
    根据选择模式扩展原子选择
    
    参数:
        atoms: ASE Atoms 对象
        selected_indices: 初步选中的原子索引
        selection_mode: 'sphere' 或 'molecule'
    
    返回:
        最终选中的原子索引
    """
    if selection_mode == 'sphere':
        return selected_indices
    
    # molecule 模式：扩展到完整分子
    molecules = get_molecule_indices(atoms)
    
    if molecules is None:
        print("[WARN] 无法获取分子信息，回退到 sphere 模式")
        print("[INFO] 建议使用 .pdb 格式或确保结构文件包含分子信息")
        return selected_indices
    
    selected_set = set(selected_indices)
    expanded_set = set()
    
    for mol_indices in molecules:
        # 如果分子中任何原子被选中，则整个分子都被选中
        if mol_indices & selected_set:
            expanded_set.update(mol_indices)
    
    expanded_indices = np.array(sorted(expanded_set))
    
    n_original = len(selected_indices)
    n_expanded = len(expanded_indices)
    if n_expanded > n_original:
        print(f"[INFO] molecule 模式: {n_original} -> {n_expanded} 原子 (扩展了 {n_expanded - n_original} 个)")
    
    return expanded_indices


def estimate_charge(atoms: Atoms) -> Tuple[int, Dict[str, int]]:
    """
    估算体系总电荷（粗略）
    
    基于元素符号和预定义的离子电荷表
    
    返回:
        (总电荷, 元素计数字典)
    """
    symbols = atoms.get_chemical_symbols()
    element_counts: Dict[str, int] = {}
    total_charge = 0
    
    for sym in symbols:
        element_counts[sym] = element_counts.get(sym, 0) + 1
        charge = CHARGE_MAP.get(sym, 0)
        total_charge += charge
    
    return total_charge, element_counts


def create_vacuum_box(cluster: Atoms, vacuum: float) -> Atoms:
    """
    将 cluster 放入带真空层的盒子中
    
    参数:
        cluster: 提取的原子簇
        vacuum: 真空层厚度 (Å)
    
    返回:
        处于大盒子中居中的 cluster
    """
    positions = cluster.get_positions()
    
    # 计算 bounding box
    min_pos = positions.min(axis=0)
    max_pos = positions.max(axis=0)
    size = max_pos - min_pos
    
    # 新盒子尺寸 = bounding box + 2 * vacuum
    cell = size + 2 * vacuum
    
    # 确保最小尺寸
    cell = np.maximum(cell, vacuum * 2)
    
    cluster.set_cell(cell)
    cluster.set_pbc([True, True, True])
    cluster.center()
    
    return cluster


def write_vasp_inputs(
    cluster: Atoms,
    outdir: str,
    temp: float,
    steps: int,
    potim: float,
    kpoints: Tuple[int, int, int],
    ncore: Optional[int],
    encut: float = 400.0
) -> bool:
    """
    生成 VASP AIMD 输入文件
    
    参数:
        cluster: 提取的原子簇
        outdir: 输出目录
        temp: AIMD 温度 (K)
        steps: AIMD 步数
        potim: 时间步长 (fs)
        kpoints: K 点网格
        ncore: NCORE 参数
        encut: 截断能 (eV)
    
    返回:
        是否成功生成 POTCAR
    """
    os.makedirs(outdir, exist_ok=True)
    
    # 检查 VASP_PP_PATH
    pp_path = os.environ.get('VASP_PP_PATH', '')
    has_pp = bool(pp_path and os.path.isdir(pp_path))
    
    potcar_generated = False
    
    if has_pp:
        # 使用 ASE Vasp calculator
        try:
            from ase.calculators.vasp import Vasp
            
            calc = Vasp(
                directory=outdir,
                xc='PBE',
                encut=encut,
                prec='Normal',
                algo='Normal',
                ediff=1e-5,
                nelm=200,
                ismear=0,
                sigma=0.05,
                ibrion=0,          # MD 模式
                mdalgo=2,          # Nosé-Hoover 恒温器
                isym=0,            # 禁用对称性（AIMD 必须）
                tebeg=temp,        # 初始温度
                teend=temp,        # 终止温度
                potim=potim,       # 时间步长
                nsw=steps,         # MD 步数
                lreal='Auto',      # 实空间投影
                lwave=False,       # 不写 WAVECAR
                lcharg=False,      # 不写 CHGCAR
                kpts=kpoints,
                gamma=True,
                ncore=ncore if ncore else 4,
            )
            calc.write_input(cluster)
            potcar_generated = True
            
        except Exception as e:
            print(f"[WARN] ASE Vasp 写入失败: {e}")
            has_pp = False
    
    if not has_pp:
        # 手动写入 POSCAR
        write(os.path.join(outdir, 'POSCAR'), cluster, format='vasp')
        
        # 手动写入 INCAR
        _write_incar_manual(outdir, temp, steps, potim, encut, ncore)
        
        # 手动写入 KPOINTS
        _write_kpoints_manual(outdir, kpoints)
        
        if not pp_path:
            print("\n[WARN] VASP_PP_PATH 未设置，未生成 POTCAR")
        else:
            print(f"\n[WARN] VASP_PP_PATH 无效: {pp_path}")
        print("[INFO] 请手动准备 POTCAR 文件")
        print("[INFO] 设置方法: export VASP_PP_PATH=/path/to/potentials")
    
    # 写入 cluster_visual.xyz
    xyz_path = os.path.join(outdir, 'cluster_visual.xyz')
    write(xyz_path, cluster, format='xyz')
    
    return potcar_generated


def _write_incar_manual(
    outdir: str,
    temp: float,
    steps: int,
    potim: float,
    encut: float,
    ncore: Optional[int]
) -> None:
    """手动写入 INCAR 文件"""
    incar_path = os.path.join(outdir, 'INCAR')
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    with open(incar_path, 'w') as f:
        f.write(f"# INCAR for AIMD - Generated by setup_aimd_ase.py\n")
        f.write(f"# {timestamp}\n")
        f.write(f"# Temperature: {temp} K, Steps: {steps}, POTIM: {potim} fs\n\n")
        
        f.write("# ============ 基础参数 ============\n")
        f.write("PREC = Normal\n")
        f.write(f"ENCUT = {encut}\n")
        f.write("ALGO = Normal\n")
        f.write("EDIFF = 1E-5    # MD 可用较宽松收敛\n")
        f.write("NELM = 200\n")
        f.write("LREAL = Auto    # 实空间投影，大体系推荐\n\n")
        
        f.write("# ============ 展宽 ============\n")
        f.write("ISMEAR = 0      # Gaussian 展宽（MD 推荐，避免 -5）\n")
        f.write("SIGMA = 0.05\n\n")
        
        f.write("# ============ 分子动力学 ============\n")
        f.write("IBRION = 0      # MD 模式（不是优化）\n")
        f.write("MDALGO = 2      # Nosé-Hoover 恒温器\n")
        f.write("ISYM = 0        # 禁用对称性（AIMD 必须！）\n")
        f.write(f"TEBEG = {temp}  # 初始温度 (K)\n")
        f.write(f"TEEND = {temp}  # 终止温度 (K)\n")
        f.write(f"POTIM = {potim}  # 时间步长 (fs)，含 H 建议 0.5-1.0\n")
        f.write(f"NSW = {steps}    # MD 步数\n")
        f.write("SMASS = -3      # Nosé-Hoover 质量参数（-3=自动）\n\n")
        
        f.write("# ============ 输出控制 ============\n")
        f.write("LWAVE = .FALSE.  # 不写 WAVECAR（省空间）\n")
        f.write("LCHARG = .FALSE. # 不写 CHGCAR\n\n")
        
        f.write("# ============ 并行参数 ============\n")
        if ncore:
            f.write(f"NCORE = {ncore}\n")
        else:
            f.write("# NCORE = 4     # 根据核数设置\n")


def _write_kpoints_manual(outdir: str, kpoints: Tuple[int, int, int]) -> None:
    """手动写入 KPOINTS 文件"""
    kpoints_path = os.path.join(outdir, 'KPOINTS')
    
    with open(kpoints_path, 'w') as f:
        f.write("Automatic mesh\n")
        f.write("0\n")
        f.write("Gamma\n")
        f.write(f"{kpoints[0]} {kpoints[1]} {kpoints[2]}\n")
        f.write("0 0 0\n")


def write_selected_indices(
    outdir: str,
    center_idx: int,
    selected_indices: np.ndarray
) -> None:
    """写入选中的原子索引信息"""
    path = os.path.join(outdir, 'selected_indices.txt')
    
    with open(path, 'w') as f:
        f.write(f"# Center atom index (0-based): {center_idx}\n")
        f.write(f"# Total selected atoms: {len(selected_indices)}\n")
        f.write(f"# Indices (0-based):\n")
        for idx in selected_indices:
            f.write(f"{idx}\n")


# ============================================================================
# 主函数
# ============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="从大体系结构切割 AIMD cluster 并生成 VASP 输入",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 以 Li 为中心，半径 8 Å
  python3 setup_aimd_ase.py --src equilibrated.pdb --center_atom Li --radius 8 --temp 350 --outdir aimd_Li8A

  # 使用原子索引（0-based）
  python3 setup_aimd_ase.py --src equilibrated.gro --center_atom 0 --radius 10

  # 保留完整分子
  python3 setup_aimd_ase.py --src system.pdb --center_atom Li --radius 8 --selection molecule

注意:
  - .gro 文件可能需要先用 gmx trjconv 转换为 .pdb/.xyz
  - 含 H 原子时建议 --potim 0.5-1.0
  - molecule 模式需要结构文件包含分子信息
        """
    )
    
    # 必需参数
    parser.add_argument("--src", required=True,
                        help="输入结构文件路径 (.pdb/.gro/.xyz/.cif)")
    parser.add_argument("--center_atom", required=True,
                        help="中心原子：整数索引（0-based）或元素符号（如 Li）")
    
    # 切割参数
    parser.add_argument("--radius", type=float, default=8.0,
                        help="切割半径 Å (默认: 8.0)")
    parser.add_argument("--selection", choices=['sphere', 'molecule'], default='sphere',
                        help="选择模式: sphere=纯半径, molecule=保留完整分子 (默认: sphere)")
    parser.add_argument("--one_based", action="store_true",
                        help="原子索引按 1-based 解释")
    
    # 盒子参数
    parser.add_argument("--vacuum", type=float, default=20.0,
                        help="真空层厚度 Å (默认: 20)")
    
    # AIMD 参数
    parser.add_argument("--temp", type=float, default=350.0,
                        help="AIMD 温度 K (默认: 350)")
    parser.add_argument("--steps", type=int, default=2000,
                        help="AIMD 步数 (默认: 2000)")
    parser.add_argument("--potim", type=float, default=2.0,
                        help="时间步长 fs (默认: 2.0; 含 H 建议 0.5-1.0)")
    parser.add_argument("--kpoints", default="1 1 1",
                        help="K 点网格 (默认: '1 1 1' Gamma-only)")
    parser.add_argument("--ncore", type=int, default=None,
                        help="NCORE 并行参数 (可选)")
    parser.add_argument("--encut", type=float, default=400.0,
                        help="截断能 eV (默认: 400)")
    
    # 输出参数
    parser.add_argument("--outdir", default="aimd_cluster",
                        help="输出目录 (默认: aimd_cluster)")
    parser.add_argument("--overwrite", action="store_true",
                        help="覆盖已存在的输出目录")
    
    # 其他
    parser.add_argument("--seed", type=int, default=None,
                        help="随机种子（保留参数，暂未使用）")
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("setup_aimd_ase.py - AIMD Cluster 切割与输入生成")
    print("=" * 70)
    print(f"输入文件: {args.src}")
    print(f"中心原子: {args.center_atom}")
    print(f"切割半径: {args.radius} Å")
    print(f"选择模式: {args.selection}")
    print(f"真空层: {args.vacuum} Å")
    print(f"温度: {args.temp} K")
    print(f"步数: {args.steps}")
    print(f"时间步长: {args.potim} fs")
    print(f"输出目录: {args.outdir}")
    print("=" * 70)
    
    # 检查输入文件
    if not os.path.isfile(args.src):
        print(f"[ERROR] 输入文件不存在: {args.src}")
        sys.exit(1)
    
    # 检查输出目录
    if os.path.exists(args.outdir):
        if args.overwrite:
            backup = f"{args.outdir}.bak.{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            print(f"[INFO] 备份已存在目录: {args.outdir} -> {backup}")
            shutil.move(args.outdir, backup)
        else:
            print(f"[ERROR] 输出目录已存在: {args.outdir}")
            print("[INFO] 使用 --overwrite 覆盖或指定其他目录")
            sys.exit(1)
    
    # 读取结构
    print("\n>>> 读取结构...")
    try:
        atoms = read(args.src)
        print(f"    原子数: {len(atoms)}")
        print(f"    元素: {set(atoms.get_chemical_symbols())}")
        
        if atoms.pbc.any():
            print(f"    PBC: {atoms.pbc}")
            print(f"    Cell: {atoms.cell.lengths()}")
            atoms.wrap()
            print("    [INFO] 已执行 wrap() 将原子包回主胞")
        else:
            print("    PBC: 无周期性边界")
            
    except Exception as e:
        print(f"[ERROR] 读取结构失败: {e}")
        if args.src.endswith('.gro'):
            print("[INFO] .gro 文件读取问题，建议使用:")
            print("       gmx trjconv -f input.gro -o output.pdb")
        sys.exit(1)
    
    # 解析中心原子
    print("\n>>> 解析中心原子...")
    try:
        center_idx = parse_center_atom(args.center_atom, atoms, args.one_based)
        center_sym = atoms.get_chemical_symbols()[center_idx]
        center_pos = atoms.get_positions()[center_idx]
        print(f"    中心原子索引: {center_idx} (0-based)")
        print(f"    中心原子元素: {center_sym}")
        print(f"    中心原子位置: {center_pos}")
    except ValueError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
    
    # 选择半径内的原子
    print("\n>>> 选择半径内的原子...")
    use_mic = atoms.pbc.any() and atoms.cell.any()
    if use_mic:
        print("    [INFO] 使用最小镜像距离 (MIC)")
    
    selected_indices = select_indices_within_radius_mic(
        atoms, center_idx, args.radius, use_mic
    )
    print(f"    初步选中: {len(selected_indices)} 原子")
    
    # 扩展到完整分子（如果需要）
    if args.selection == 'molecule':
        print("\n>>> 扩展到完整分子...")
        selected_indices = expand_to_molecules_if_possible(
            atoms, selected_indices, args.selection
        )
    
    print(f"    最终选中: {len(selected_indices)} 原子")
    
    # 创建 cluster
    cluster = atoms[selected_indices]
    
    # 检查是否包含氢
    if 'H' in cluster.get_chemical_symbols():
        print(f"\n[WARN] cluster 包含 H 原子，建议 POTIM=0.5-1.0 (当前: {args.potim})")
    
    # 估算电荷
    print("\n>>> 估算电荷...")
    total_charge, element_counts = estimate_charge(cluster)
    print(f"    元素组成: {element_counts}")
    print(f"    估算总电荷: {total_charge:+d}")
    
    if abs(total_charge) >= 1:
        print(f"\n[WARN] cluster 可能带电 ({total_charge:+d})！")
        print("[INFO] 可能需要调整 VASP NELECT 或重新选择中心/半径")
    
    # 创建真空盒
    print("\n>>> 创建真空盒...")
    cluster = create_vacuum_box(cluster, args.vacuum)
    cell = cluster.get_cell().lengths()
    print(f"    盒子尺寸: {cell[0]:.1f} x {cell[1]:.1f} x {cell[2]:.1f} Å")
    
    # 解析 K 点
    kpts_parts = args.kpoints.strip().split()
    if len(kpts_parts) != 3:
        print(f"[ERROR] K 点格式错误: {args.kpoints}，应为 'k1 k2 k3'")
        sys.exit(1)
    kpoints = tuple(int(x) for x in kpts_parts)
    
    # 生成 VASP 输入
    print("\n>>> 生成 VASP 输入文件...")
    potcar_ok = write_vasp_inputs(
        cluster, args.outdir,
        args.temp, args.steps, args.potim, kpoints,
        args.ncore, args.encut
    )
    
    # 保存选中索引
    write_selected_indices(args.outdir, center_idx, selected_indices)
    
    # 输出摘要
    print("\n" + "=" * 70)
    print("切割完成！")
    print("=" * 70)
    print(f"提取原子数: {len(cluster)}")
    print(f"中心原子: {center_sym} (索引 {center_idx})")
    print(f"切割半径: {args.radius} Å")
    print(f"选择模式: {args.selection}")
    print(f"估算电荷: {total_charge:+d}")
    print(f"输出目录: {args.outdir}")
    print("=" * 70)
    
    print("\n输出文件:")
    print(f"  - {args.outdir}/POSCAR")
    print(f"  - {args.outdir}/INCAR")
    print(f"  - {args.outdir}/KPOINTS")
    if potcar_ok:
        print(f"  - {args.outdir}/POTCAR")
    print(f"  - {args.outdir}/cluster_visual.xyz")
    print(f"  - {args.outdir}/selected_indices.txt")
    
    print("\n下一步:")
    if not potcar_ok:
        print("  1. 准备 POTCAR 文件:")
        print(f"     export VASP_PP_PATH=/path/to/potentials")
        print(f"     cat POT_*/POTCAR > {args.outdir}/POTCAR")
        print("  2. 运行 AIMD:")
    else:
        print("  1. 运行 AIMD:")
    
    print(f"     cd {args.outdir}")
    print(f"     NP=16 EXE=vasp_std run_vasp.sh")
    print("  3. 监控: aimd_watch.sh")
    print("  4. 可视化: 使用 OVITO/VMD 打开 cluster_visual.xyz")


if __name__ == "__main__":
    main()

