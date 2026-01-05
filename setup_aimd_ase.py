#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
setup_aimd_ase.py v2.3.1 - 从大体系结构切割 AIMD 子体系并生成 VASP 输入

================================================================================
CHANGELOG v2.3.1 (Production Hardening):
================================================================================
1. [PERF] Clash detection speedup: use get_all_distances(mic=True) once,
   instead of O(N^2) individual safe_get_distances() calls
2. [FIX] Multi-frame file: force read(filepath, index=0) to ensure first frame
3. [FIX] LANGEVIN_GAMMA order: parse POSCAR line 6 to match exact element order
   - New function: parse_poscar_element_order()
   - Guarantees NTYP order consistency with POSCAR
4. [FIX] Relax script: backup INCAR to INCAR.aimd before overwriting
5. [FIX] Cut-bond graph: use has_valid_cell() for cell/pbc consistency
6. [DOC] Clarified clash threshold logic comments

================================================================================
CHANGELOG v2.3:
================================================================================
1. [FIX] get_distances() return order guard (ASE version compatibility)
2. [FIX] Safer cell checks: atoms.pbc.any() AND cell.volume > 1e-8
3. [FIX] LANGEVIN_GAMMA format: per-element-type (NTYP) instead of per-atom
4. [FIX] Density sanity checks: warn if <0.5 or >3.0 g/cm³
5. [FIX] Neutralization verification: re-check charge after adding counterions
6. [FIX] Multi-frame file handling: safely take first frame from trajectories
7. [FIX] Strict --kpoints validation: must be exactly 3 integers
8. [NEW] Clash detection: detect atomic overlaps after density compression
9. [NEW] Relaxation guidance: RELAX_GUIDE.txt + optional INCAR.relax

关键改进 (v2.3 → v2.3.1):
  - ✅ 碰撞检测加速 ~100x (N=400: ~8万次调用 → 1次矩阵运算)
  - ✅ 多帧文件严格取第一帧
  - ✅ LANGEVIN_GAMMA 顺序与 POSCAR 元素行严格一致
  - ✅ 弛豫脚本自动备份 INCAR → INCAR.aimd
  - ✅ 切断键构建统一使用 has_valid_cell() 判定

用法:
  # bulk 模式（默认，按原体系密度定盒子）
  python3 setup_aimd_ase.py --src eq.pdb --center_atom Li --radius 8 --mode bulk

  # bulk 模式，指定目标密度（带弛豫输入）
  python3 setup_aimd_ase.py --src eq.pdb --center_atom Li --radius 8 \\
      --density_g_cm3 1.2 --write_relax_inputs

  # cluster 模式（真空簇，需显式指定）
  python3 setup_aimd_ase.py --src eq.pdb --center_atom Li --radius 8 \\
      --mode cluster --vacuum 20

  # 使用 bond_hops 避免切断聚合物链
  python3 setup_aimd_ase.py --src eq.pdb --center_atom Li --radius 8 --bond_hops 3

依赖:
  pip install ase numpy

作者: STAR0418-ABC
版本: v2.3
"""

from __future__ import annotations

import argparse
import os
import sys
import shutil
import json
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Set, Any
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

# 导入本地工具模块
try:
    from utils.connectivity import (
        build_bond_graph, find_connected_components,
        expand_by_bond_hops, detect_cut_bonds, write_cut_bonds_report
    )
    from utils.charges import (
        ELEMENT_CHARGE_MAP, RESIDUE_CHARGE_MAP,
        estimate_charge_by_residue, estimate_charge_by_element,
        find_counterion_residues, load_charge_map_file,
        CATION_RESIDUES, ANION_RESIDUES
    )
    from utils.units import (
        compute_mass, compute_density, volume_from_density,
        scale_cell_to_volume, ATOMIC_MASSES
    )
    HAS_UTILS = True
except ImportError:
    HAS_UTILS = False
    print("[WARN] utils 模块未找到，使用内置简化版本")

VERSION = "v2.3.1"

# ==============================================================================
# 共价半径表（用于碰撞检测）
# ==============================================================================
COVALENT_RADII = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'S': 1.05, 'P': 1.07,
    'Li': 1.28, 'Na': 1.66, 'K': 2.03, 'Mg': 1.41, 'Ca': 1.76, 'Zn': 1.22, 'Al': 1.21,
    'Cl': 1.02, 'Br': 1.20, 'I': 1.39, 'Si': 1.11, 'B': 0.84,
}

# ==============================================================================
# 内置简化版工具（当 utils 不可用时）
# ==============================================================================

if not HAS_UTILS:
    ELEMENT_CHARGE_MAP = {
        'Li': +1, 'Na': +1, 'K': +1, 'Mg': +2, 'Ca': +2, 'Zn': +2, 'Al': +3,
        'F': -1, 'Cl': -1, 'Br': -1, 'I': -1,
        'C': 0, 'H': 0, 'O': 0, 'N': 0, 'S': 0, 'P': 0,
    }
    
    ATOMIC_MASSES = {
        'H': 1.008, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'F': 19.00, 'S': 32.07,
        'Li': 6.941, 'Na': 22.99, 'K': 39.10, 'Mg': 24.31, 'Ca': 40.08,
        'Zn': 65.38, 'Al': 26.98, 'Cl': 35.45, 'Br': 79.90,
    }
    
    def compute_mass(symbols):
        return sum(ATOMIC_MASSES.get(s, 10.0) for s in symbols)
    
    def compute_density(symbols, volume_A3):
        mass_g = compute_mass(symbols) / 6.022e23
        volume_cm3 = volume_A3 * 1e-24
        return mass_g / volume_cm3
    
    def volume_from_density(symbols, density):
        mass_g = compute_mass(symbols) / 6.022e23
        volume_cm3 = mass_g / density
        return volume_cm3 / 1e-24
    
    def scale_cell_to_volume(cell, target_vol, mode='scale_proportional'):
        if mode == 'cubic':
            L = target_vol ** (1/3)
            return np.diag([L, L, L])
        current_vol = abs(np.linalg.det(cell))
        if current_vol < 1e-10:
            L = target_vol ** (1/3)
            return np.diag([L, L, L])
        scale = (target_vol / current_vol) ** (1/3)
        return cell * scale
    
    def build_bond_graph(positions, symbols, cell=None, pbc=None, **kwargs):
        return {i: set() for i in range(len(positions))}
    
    def detect_cut_bonds(graph, selected):
        return []
    
    def expand_by_bond_hops(graph, seed, max_hops, max_size):
        return seed


# ==============================================================================
# 辅助函数：安全检查与兼容性
# ==============================================================================

def has_valid_cell(atoms: Atoms) -> bool:
    """
    安全检查原子对象是否有有效的周期性晶胞
    
    替代不安全的 atoms.cell.any() 检查
    """
    if not atoms.pbc.any():
        return False
    try:
        vol = atoms.cell.volume
        return vol > 1e-8
    except:
        return False


def safe_get_distances(
    center_pos: np.ndarray,
    positions: np.ndarray,
    cell: Optional[np.ndarray] = None,
    pbc: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    ASE get_distances() 的安全包装，处理返回值顺序歧义
    
    不同 ASE 版本可能返回 (vectors, distances) 或 (distances, vectors)
    通过形状检测来确定正确顺序:
    - vectors: shape (1, n_atoms, 3) 或 (n_atoms, 3)，最后一维 = 3
    - distances: shape (1, n_atoms) 或 (n_atoms,)，是标量数组
    
    Returns:
        mic_vectors: (n_atoms, 3) MIC 位移向量
        distances: (n_atoms,) 距离数组
    """
    result = get_distances([center_pos], positions, cell=cell, pbc=pbc)
    
    ret0, ret1 = result[0], result[1]
    
    # 检测哪个是向量（最后一维是 3）
    shape0 = np.array(ret0).shape
    shape1 = np.array(ret1).shape
    
    # vectors 的最后一维应该是 3
    if len(shape0) >= 2 and shape0[-1] == 3:
        vectors = np.array(ret0)
        distances = np.array(ret1)
    elif len(shape1) >= 2 and shape1[-1] == 3:
        vectors = np.array(ret1)
        distances = np.array(ret0)
    else:
        # 回退：假设 (vectors, distances) 顺序
        vectors = np.array(ret0)
        distances = np.array(ret1)
    
    # 归一化形状
    if vectors.ndim == 3:
        vectors = vectors[0]  # (1, n, 3) -> (n, 3)
    if distances.ndim == 2:
        distances = distances[0]  # (1, n) -> (n,)
    
    return vectors, distances


def validate_kpoints(kpts_str: str) -> Tuple[int, int, int]:
    """验证 K 点字符串，必须是 3 个整数"""
    parts = kpts_str.split()
    if len(parts) != 3:
        raise ValueError(f"--kpoints 必须是 3 个整数 (如 '1 1 1')，收到: '{kpts_str}'")
    try:
        kpts = tuple(int(x) for x in parts)
    except ValueError:
        raise ValueError(f"--kpoints 必须是整数，收到: '{kpts_str}'")
    if any(k < 1 for k in kpts):
        raise ValueError(f"K 点值必须 >= 1，收到: {kpts}")
    return kpts


def safe_read_structure(filepath: str) -> Atoms:
    """
    安全读取结构文件，强制取第一帧
    
    v2.3.1: 使用 index=0 确保读取第一帧，避免不同 reader 行为差异
    """
    try:
        # 优先使用 index=0 强制取第一帧
        result = read(filepath, index=0)
    except TypeError:
        # 某些格式可能不支持 index 参数
        result = read(filepath)
    
    # 双重保险：如果仍然返回 list，取第一个
    if isinstance(result, list):
        if len(result) == 0:
            raise ValueError(f"文件 {filepath} 不包含任何结构")
        print(f"[INFO] 文件包含 {len(result)} 帧，使用第一帧")
        return result[0]
    
    return result


def parse_poscar_element_order(poscar_path: str) -> List[str]:
    """
    解析 POSCAR 第 6 行获取元素顺序
    
    v2.3.1: 确保 LANGEVIN_GAMMA 与 POSCAR 元素顺序严格一致
    
    POSCAR 格式:
        Line 1: Comment
        Line 2: Scale
        Line 3-5: Lattice vectors
        Line 6: Element symbols (e.g., "C H O Li F S N")
        Line 7: Atom counts
        ...
    
    Returns:
        elements: 元素符号列表，按 POSCAR 顺序
    """
    with open(poscar_path, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 7:
        raise ValueError(f"POSCAR 格式错误: {poscar_path}")
    
    # 第 6 行是元素符号（0-indexed line 5）
    element_line = lines[5].strip()
    elements = element_line.split()
    
    # 验证：应该都是元素符号，不是数字
    for elem in elements:
        if elem.isdigit():
            raise ValueError(f"POSCAR 第 6 行应为元素符号，但发现数字: {element_line}")
    
    return elements


# ==============================================================================
# 核心函数
# ==============================================================================

def parse_center_atom(center_str: str, atoms: Atoms, one_based: bool = False) -> int:
    """解析中心原子参数"""
    n_atoms = len(atoms)
    symbols = atoms.get_chemical_symbols()
    
    try:
        idx = int(center_str)
        if one_based:
            idx -= 1
        if idx == n_atoms:
            print(f"[WARN] 索引 {center_str} = 原子总数，按 1-based 解释为 {idx - 1}")
            idx -= 1
        if idx < 0 or idx >= n_atoms:
            raise ValueError(f"原子索引 {idx} 超出范围 [0, {n_atoms - 1}]")
        return idx
    except ValueError:
        pass
    
    element = center_str.strip()
    for i, sym in enumerate(symbols):
        if sym == element:
            print(f"[INFO] 找到第一个 {element} 原子，索引 {i}")
            return i
    
    raise ValueError(f"未找到元素 '{element}'，可用: {set(symbols)}")


def select_indices_with_mic_vectors(
    atoms: Atoms,
    center_idx: int,
    radius: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    选择中心原子周围指定半径内的所有原子，返回 MIC 位移向量
    
    Returns:
        indices: 选中的原子索引
        distances: 到中心的距离
        mic_vectors: MIC 位移向量 (center -> atom)
    """
    positions = atoms.get_positions()
    center_pos = positions[center_idx]
    
    # 使用安全的晶胞检查
    use_mic = has_valid_cell(atoms)
    
    if use_mic:
        try:
            # 使用安全包装处理 ASE 版本兼容性
            mic_vectors, distances = safe_get_distances(
                center_pos, positions,
                cell=atoms.cell, pbc=atoms.pbc
            )
        except Exception as e:
            print(f"[WARN] MIC 计算失败: {e}，回退到普通距离")
            mic_vectors = positions - center_pos
            distances = np.linalg.norm(mic_vectors, axis=1)
    else:
        mic_vectors = positions - center_pos
        distances = np.linalg.norm(mic_vectors, axis=1)
    
    indices = np.where(distances <= radius)[0]
    return indices, distances, mic_vectors


def reimage_atoms_around_center(
    atoms: Atoms,
    selected_indices: np.ndarray,
    center_idx: int,
    mic_vectors: np.ndarray
) -> Atoms:
    """
    将选中的原子重新成像到中心附近（空间连贯）
    
    使用 MIC 位移向量重建坐标: r_i = r_center + mic_vector(center→i)
    """
    center_pos = atoms.get_positions()[center_idx]
    symbols = atoms.get_chemical_symbols()
    
    # 构建新坐标
    new_positions = []
    new_symbols = []
    
    for idx in selected_indices:
        new_pos = center_pos + mic_vectors[idx]
        new_positions.append(new_pos)
        new_symbols.append(symbols[idx])
    
    # 创建新 Atoms 对象
    cluster = Atoms(
        symbols=new_symbols,
        positions=new_positions
    )
    
    # 复制数组属性（如 residuenumbers）
    if hasattr(atoms, 'arrays'):
        for key, arr in atoms.arrays.items():
            if key in ['positions', 'numbers']:
                continue
            if len(arr) == len(atoms):
                try:
                    cluster.arrays[key] = arr[selected_indices]
                except:
                    pass
    
    return cluster


def get_residue_info(atoms: Atoms) -> Tuple[Optional[List[str]], Optional[List[int]]]:
    """获取残基信息"""
    arrays = atoms.arrays if hasattr(atoms, 'arrays') else {}
    
    # 残基名称
    resnames = None
    for key in ['residuenames', 'resname', 'resnames']:
        if key in arrays:
            resnames = [str(r) for r in arrays[key]]
            break
    
    # 残基编号
    resids = None
    for key in ['residuenumbers', 'resid', 'resids', 'molid']:
        if key in arrays:
            resids = list(arrays[key])
            break
    
    return resnames, resids


def expand_selection_with_molecules(
    atoms: Atoms,
    selected_indices: np.ndarray,
    selection_mode: str,
    bond_hops: int,
    max_atoms: int,
    allow_exceed: bool
) -> Tuple[np.ndarray, bool, int]:
    """
    扩展选择到完整分子/键跳
    
    Returns:
        final_indices: 最终选中的索引
        was_truncated: 是否被截断
        n_cut_bonds: 切断键数量
    """
    if selection_mode == 'sphere' and bond_hops <= 0:
        return selected_indices, False, 0
    
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    # 使用安全的晶胞检查
    valid_cell = has_valid_cell(atoms)
    cell = atoms.cell if valid_cell else None
    pbc = atoms.pbc if valid_cell else None
    
    # 构建键图
    if HAS_UTILS:
        graph = build_bond_graph(positions, symbols, cell, pbc)
    else:
        graph = {i: set() for i in range(len(atoms))}
    
    selected_set = set(selected_indices)
    
    # 根据模式扩展
    if selection_mode == 'molecule':
        # 尝试基于残基信息扩展
        resnames, resids = get_residue_info(atoms)
        
        if resnames is not None and resids is not None:
            # 按残基扩展
            touched_residues = set()
            for idx in selected_indices:
                touched_residues.add((resnames[idx], resids[idx]))
            
            expanded_set = set()
            for i in range(len(atoms)):
                if (resnames[i], resids[i]) in touched_residues:
                    expanded_set.add(i)
            
            selected_set = expanded_set
        else:
            # 回退到连通分量
            print("[INFO] 无残基信息，使用键图连通性扩展")
            if HAS_UTILS:
                components = find_connected_components(graph)
                touched_components = []
                for comp in components:
                    if comp & set(selected_indices):
                        touched_components.append(comp)
                
                for comp in touched_components:
                    if len(comp) <= max_atoms or allow_exceed:
                        selected_set.update(comp)
    
    # 键跳扩展
    if bond_hops > 0 and HAS_UTILS:
        print(f"[INFO] 执行 {bond_hops} 步键跳扩展...")
        selected_set = expand_by_bond_hops(graph, selected_set, bond_hops, max_atoms)
    
    # 检查原子数限制
    was_truncated = False
    if len(selected_set) > max_atoms and not allow_exceed:
        print(f"[WARN] 扩展后原子数 ({len(selected_set)}) 超过 max_atoms ({max_atoms})")
        print("[INFO] 回退到原始 sphere 选择")
        selected_set = set(selected_indices)
        was_truncated = True
    
    final_indices = np.array(sorted(selected_set))
    
    # 检测切断键
    if HAS_UTILS:
        cut_bonds = detect_cut_bonds(graph, set(final_indices))
        n_cut_bonds = len(cut_bonds)
    else:
        n_cut_bonds = 0
    
    return final_indices, was_truncated, n_cut_bonds


def create_density_based_bulk_box(
    cluster: Atoms,
    original_density: Optional[float],
    target_density: Optional[float],
    cell_shape: str = 'scale_parent',
    parent_cell: Optional[np.ndarray] = None,
    min_cell_length: float = 10.0
) -> Tuple[Atoms, float, float]:
    """
    创建基于密度的周期性盒子
    
    公式: V_target = M_sub / ρ_target
    
    Args:
        cluster: 子体系原子
        original_density: 原体系密度 (g/cm³)
        target_density: 目标密度 (g/cm³)，None 则使用 original_density
        cell_shape: 'scale_parent' 或 'cubic'
        parent_cell: 父体系晶胞（用于 scale_parent）
        min_cell_length: 最小盒子边长 (Å)
    
    Returns:
        cluster: 带新晶胞的子体系
        target_density: 使用的目标密度
        achieved_density: 实际达到的密度
    """
    symbols = cluster.get_chemical_symbols()
    
    # 确定目标密度
    if target_density is None:
        if original_density is not None:
            target_density = original_density
        else:
            # 默认凝胶电解质密度
            target_density = 1.2
            print(f"[WARN] 无法确定原体系密度，使用默认值 {target_density} g/cm³")
    
    # 计算目标体积
    target_volume = volume_from_density(symbols, target_density)
    
    # 确保最小盒子尺寸
    min_volume = min_cell_length ** 3
    if target_volume < min_volume:
        print(f"[WARN] 目标体积 ({target_volume:.1f} Å³) 太小，调整到最小 ({min_volume:.1f} Å³)")
        target_volume = min_volume
    
    # 生成新晶胞
    if cell_shape == 'scale_parent' and parent_cell is not None:
        new_cell = scale_cell_to_volume(parent_cell, target_volume, 'scale_proportional')
    else:
        new_cell = scale_cell_to_volume(np.eye(3) * 10, target_volume, 'cubic')
    
    # 设置晶胞并居中
    cluster.set_cell(new_cell)
    cluster.set_pbc([True, True, True])
    cluster.center()
    
    # 计算实际达到的密度
    actual_volume = abs(np.linalg.det(new_cell))
    achieved_density = compute_density(symbols, actual_volume)
    
    return cluster, target_density, achieved_density


def create_vacuum_box(cluster: Atoms, vacuum: float) -> Atoms:
    """创建真空盒子（cluster 模式）"""
    positions = cluster.get_positions()
    
    min_pos = positions.min(axis=0)
    max_pos = positions.max(axis=0)
    size = max_pos - min_pos
    
    cell = size + 2 * vacuum
    cell = np.maximum(cell, vacuum * 2)
    
    cluster.set_cell(cell)
    cluster.set_pbc([True, True, True])
    cluster.center()
    
    return cluster


def estimate_charge_comprehensive(
    cluster: Atoms,
    charge_map_file: Optional[str] = None
) -> Tuple[int, Dict[str, int], bool, List[str]]:
    """
    综合电荷估计：优先残基，回退元素
    
    Returns:
        total_charge: 总电荷
        counts: 组分计数
        is_reliable: 是否可靠
        warnings: 警告信息
    """
    symbols = cluster.get_chemical_symbols()
    resnames, resids = get_residue_info(cluster)
    
    # 加载自定义电荷表
    custom_res_map = None
    custom_elem_map = None
    if charge_map_file and HAS_UTILS:
        try:
            custom_res_map, custom_elem_map = load_charge_map_file(charge_map_file)
            print(f"[INFO] 已加载自定义电荷表: {charge_map_file}")
        except Exception as e:
            print(f"[WARN] 加载电荷表失败: {e}")
    
    # 尝试残基估计
    if resnames is not None and resids is not None and HAS_UTILS:
        selected_set = set(range(len(cluster)))
        charge, counts, reliable, warnings = estimate_charge_by_residue(
            resnames, resids, selected_set,
            custom_res_map, symbols, custom_elem_map
        )
        return charge, counts, reliable, warnings
    
    # 回退到元素估计
    charge, counts, reliable = estimate_charge_by_element(
        symbols, custom_elem_map
    )
    warnings = ["使用元素级电荷估计（低可信度）"]
    return charge, counts, reliable, warnings


def neutralize_by_residue(
    atoms: Atoms,
    selected_indices: np.ndarray,
    center_idx: int,
    current_charge: int,
    target_charge: int,
    all_distances: np.ndarray,
    mic_vectors: np.ndarray
) -> Tuple[np.ndarray, int, List[str]]:
    """
    按残基/分子中和电荷
    
    Returns:
        new_indices: 添加反离子后的索引
        added_count: 添加的原子数
        info: 信息列表
    """
    if current_charge == target_charge:
        return selected_indices, 0, []
    
    resnames, resids = get_residue_info(atoms)
    positions = atoms.get_positions()
    center_pos = positions[center_idx]
    # 使用安全的晶胞检查
    cell = atoms.cell if has_valid_cell(atoms) else None
    
    info = []
    
    if resnames is not None and resids is not None and HAS_UTILS:
        # 确保选中集合不包含残基的“半截”
        selected_set = set(selected_indices)
        residue_groups: Dict[Tuple[str, int], Set[int]] = {}
        for idx, (resname, resid) in enumerate(zip(resnames, resids)):
            residue_groups.setdefault((resname, resid), set()).add(idx)

        expanded_selected = set()
        for (resname, resid), indices in residue_groups.items():
            if indices & selected_set:
                expanded_selected.update(indices)

        if expanded_selected != selected_set:
            selected_set = expanded_selected
            selected_indices = np.array(sorted(selected_set))
            info.append("[INFO] 已将选中集合扩展为完整残基，避免部分残基参与中和")

        # 按残基中和
        counterion_groups, remaining = find_counterion_residues(
            resnames, resids, selected_set,
            current_charge, target_charge,
            positions, center_pos, cell
        )
        
        added_indices = []
        for group in counterion_groups:
            added_indices.extend(group)
            # 获取残基名称
            if group:
                idx = list(group)[0]
                rname = resnames[idx] if idx < len(resnames) else "?"
                info.append(f"添加残基 {rname} ({len(group)} 原子)")
        
        if added_indices:
            new_indices = np.concatenate([selected_indices, np.array(added_indices)])
            new_indices = np.unique(new_indices)
            return new_indices, len(added_indices), info
    
    # 回退：按元素（单原子）中和
    info.append("[WARN] 无残基信息，按元素原子中和（可能不准确）")
    
    symbols = atoms.get_chemical_symbols()
    selected_set = set(selected_indices)
    
    needed_charge = current_charge - target_charge
    if needed_charge > 0:
        target_elements = {'F', 'Cl', 'Br', 'I'}
    else:
        target_elements = {'Li', 'Na', 'K', 'Mg', 'Ca', 'Zn'}
    
    candidates = []
    for i, sym in enumerate(symbols):
        if i not in selected_set and sym in target_elements:
            charge = abs(ELEMENT_CHARGE_MAP.get(sym, 0))
            if charge > 0:
                candidates.append((i, all_distances[i], sym, charge))
    
    candidates.sort(key=lambda x: x[1])
    
    added = []
    remaining = abs(needed_charge)
    for idx, dist, sym, charge in candidates:
        if remaining <= 0:
            break
        added.append(idx)
        remaining -= charge
        info.append(f"添加 {sym} (索引 {idx}, 距离 {dist:.2f} Å)")
    
    if added:
        new_indices = np.concatenate([selected_indices, np.array(added)])
        new_indices = np.unique(new_indices)
        return new_indices, len(added), info
    
    return selected_indices, 0, info


def detect_atomic_clashes(
    cluster: Atoms,
    scale: float = 0.75,
    global_min_h: float = 0.8,
    global_min_other: float = 1.2
) -> Dict[str, Any]:
    """
    检测原子碰撞/重叠（v2.3.1 加速版）
    
    当密度压缩后，边界原子可能重叠，导致 AIMD 不稳定。
    
    v2.3.1: 使用 get_all_distances(mic=True) 一次性计算距离矩阵，
            避免 O(N^2) 次 safe_get_distances() 调用
            N=400 时加速约 100x
    
    Args:
        cluster: 原子对象（已设置 PBC）
        scale: 共价半径比例阈值 (默认 0.75)
        global_min_h: 涉及 H 的全局最小距离阈值 (Å)
        global_min_other: 其他原子对的全局最小距离阈值 (Å)
    
    阈值逻辑说明:
        对每对原子 (i, j)，计算两种阈值:
        1. threshold_cov = (r_cov_i + r_cov_j) * scale  # 基于共价半径
        2. threshold_global = 0.8 Å (涉及H) 或 1.2 Å (其他)
        取 min(threshold_cov, threshold_global) 作为最终阈值。
        
        这是"宽松"策略：只捕捉严重重叠（d < 较小阈值），
        避免误报正常的短距离（如 H-O 氢键 ~1.8Å）。
    
    Returns:
        clash_info: 包含 d_min, clash_count, clash_pairs, has_clashes
    """
    n_atoms = len(cluster)
    symbols = cluster.get_chemical_symbols()
    
    # v2.3.1: 一次性计算全距离矩阵（关键加速）
    use_mic = has_valid_cell(cluster)
    dist_mat = cluster.get_all_distances(mic=use_mic)  # (N, N) 对称矩阵
    
    clash_pairs = []
    d_min = float('inf')
    d_min_pair = None
    
    # 遍历上三角矩阵（i < j）
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            d_ij = dist_mat[i, j]
            
            # 更新最小距离
            if d_ij < d_min:
                d_min = d_ij
                d_min_pair = (i, j)
            
            # 检测碰撞
            sym_i, sym_j = symbols[i], symbols[j]
            
            # 方法1：基于共价半径比例
            r_cov_i = COVALENT_RADII.get(sym_i, 1.0)
            r_cov_j = COVALENT_RADII.get(sym_j, 1.0)
            threshold_cov = (r_cov_i + r_cov_j) * scale
            
            # 方法2：全局硬阈值
            if 'H' in (sym_i, sym_j):
                threshold_global = global_min_h
            else:
                threshold_global = global_min_other
            
            # 取较小值（宽松策略：只抓严重重叠，避免误报）
            threshold = min(threshold_cov, threshold_global)
            
            if d_ij < threshold:
                clash_pairs.append({
                    'i': int(i),
                    'j': int(j),
                    'sym_i': sym_i,
                    'sym_j': sym_j,
                    'distance': float(d_ij),
                    'threshold': float(threshold)
                })
    
    # 按距离排序
    clash_pairs.sort(key=lambda x: x['distance'])
    
    return {
        'd_min': float(d_min) if d_min < float('inf') else None,
        'd_min_pair': d_min_pair,
        'clash_count': len(clash_pairs),
        'clash_pairs': clash_pairs[:20],  # 只保留前 20 个
        'has_clashes': len(clash_pairs) > 0
    }


def write_relax_guide(
    outdir: str,
    clash_info: Dict[str, Any],
    has_h: bool,
    gamma_1ps: float,
    volume_compression_ratio: float
) -> None:
    """
    写入弛豫/预平衡指导文件
    
    当检测到碰撞或大体积压缩时，AIMD 前需要弛豫
    """
    guide_path = os.path.join(outdir, 'RELAX_GUIDE.txt')
    
    with open(guide_path, 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\n")
        f.write("AIMD 前处理指南 - 生成自 setup_aimd_ase.py v2.3\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("【重要警告】\n")
        f.write("-" * 70 + "\n")
        
        if clash_info['has_clashes']:
            f.write(f"⚠️  检测到 {clash_info['clash_count']} 对原子碰撞/重叠！\n")
            f.write(f"⚠️  最小原子间距: {clash_info['d_min']:.3f} Å\n")
            f.write("⚠️  直接运行 AIMD 可能导致:\n")
            f.write("    - 巨大的初始力和能量\n")
            f.write("    - SCF 不收敛\n")
            f.write("    - 模拟爆炸/崩溃\n\n")
        
        if volume_compression_ratio > 0.3:
            f.write(f"⚠️  体积压缩比例: {volume_compression_ratio*100:.1f}%\n")
            f.write("⚠️  大幅压缩可能导致原子重叠\n\n")
        
        f.write("\n【推荐工作流】\n")
        f.write("-" * 70 + "\n")
        f.write("步骤 1: 离子弛豫 (必须！)\n")
        f.write("        - 使用 INCAR.relax (如已生成)\n")
        f.write("        - IBRION = 2 (共轭梯度)\n")
        f.write("        - NSW = 200-500\n")
        f.write("        - ISIF = 2 (固定晶胞)\n")
        f.write(f"        - POTIM = {0.5 if has_h else 1.0} (含 H 用更小步长)\n")
        f.write("        - EDIFFG = -0.02 到 -0.05\n\n")
        
        f.write("步骤 2: NVT 预平衡\n")
        f.write("        - 从弛豫后的 CONTCAR 启动\n")
        f.write("        - 较小 POTIM (0.5-1.0 fs)\n")
        f.write("        - 较大 LANGEVIN_GAMMA (10-20) 快速控温\n")
        f.write("        - NSW = 2000-5000\n\n")
        
        f.write("步骤 3: 生产 AIMD\n")
        f.write("        - 从预平衡 CONTCAR 启动\n")
        f.write("        - 正常 POTIM (1.0-2.0 fs)\n")
        f.write("        - 较小 LANGEVIN_GAMMA (1-5) 减少对动力学干扰\n\n")
        
        if gamma_1ps >= 10:
            f.write(f"【注意】当前 LANGEVIN_GAMMA = {gamma_1ps}\n")
            f.write("        适合平衡段，生产段建议降至 1-5\n\n")
        
        f.write("\n【碰撞详情】\n")
        f.write("-" * 70 + "\n")
        if clash_info['has_clashes']:
            f.write("最差的原子对 (前 20 个):\n")
            for i, cp in enumerate(clash_info['clash_pairs'][:20]):
                f.write(f"  {i+1:2d}. 原子 {cp['i']:4d} ({cp['sym_i']}) - "
                       f"原子 {cp['j']:4d} ({cp['sym_j']}): "
                       f"{cp['distance']:.3f} Å (阈值 {cp['threshold']:.3f} Å)\n")
        else:
            f.write("未检测到严重碰撞（d_min 可能仍较小，建议弛豫）\n")
        
        f.write("\n" + "=" * 70 + "\n")


def write_incar_relax(
    outdir: str,
    encut: float,
    ncore: Optional[int],
    has_h: bool,
    n_element_types: int
) -> None:
    """写入弛豫用 INCAR"""
    incar_path = os.path.join(outdir, 'INCAR.relax')
    
    potim = 0.5 if has_h else 1.0
    
    with open(incar_path, 'w', encoding='utf-8') as f:
        f.write(f"# INCAR for relaxation - Generated by setup_aimd_ase.py {VERSION}\n")
        f.write("# 用于密度压缩后的离子弛豫（AIMD 前必须步骤）\n\n")
        
        f.write("# ============ 基础参数 ============\n")
        f.write("PREC = Accurate\n")
        f.write(f"ENCUT = {encut}\n")
        f.write("ALGO = Normal\n")
        f.write("EDIFF = 1E-6\n")
        f.write("NELM = 200\n")
        f.write("LREAL = Auto\n\n")
        
        f.write("# ============ 展宽 ============\n")
        f.write("ISMEAR = 0\n")
        f.write("SIGMA = 0.05\n\n")
        
        f.write("# ============ 离子弛豫 ============\n")
        f.write("IBRION = 2       # 共轭梯度\n")
        f.write("ISIF = 2         # 固定晶胞，弛豫离子\n")
        f.write("NSW = 300        # 最大步数\n")
        f.write(f"POTIM = {potim}        # 步长 (含 H 用更小值)\n")
        f.write("EDIFFG = -0.03   # 力收敛标准 (eV/Å)\n")
        f.write("ISYM = 0         # 关闭对称性\n\n")
        
        f.write("# ============ 输出控制 ============\n")
        f.write("LWAVE = .FALSE.\n")
        f.write("LCHARG = .TRUE.  # 保留 CHGCAR 供后续计算\n")
        
        if ncore:
            f.write(f"\nNCORE = {ncore}\n")


def write_relax_script(outdir: str) -> None:
    """
    写入弛豫-AIMD 运行脚本
    
    v2.3.1: 自动备份 INCAR 到 INCAR.aimd
    """
    script_path = os.path.join(outdir, 'run_relax_then_aimd.sh')
    
    # v2.3.1: 如果 INCAR 存在，先备份为 INCAR.aimd
    incar_path = os.path.join(outdir, 'INCAR')
    incar_aimd_path = os.path.join(outdir, 'INCAR.aimd')
    if os.path.exists(incar_path) and not os.path.exists(incar_aimd_path):
        shutil.copy(incar_path, incar_aimd_path)
        print(f"    [INFO] 已备份 INCAR → INCAR.aimd")
    
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write("#!/bin/bash\n")
        f.write("# ================================================================\n")
        f.write("# 弛豫 + AIMD 分步运行脚本\n")
        f.write(f"# 生成自 setup_aimd_ase.py {VERSION}\n")
        f.write("# ================================================================\n\n")
        
        f.write("set -e  # 出错即停\n\n")
        
        f.write("# ================== Step 1: 离子弛豫 ==================\n")
        f.write("echo '>>> Step 1: 离子弛豫'\n")
        f.write("if [ -f INCAR.relax ]; then\n")
        f.write("    # 备份 AIMD INCAR（如果尚未备份）\n")
        f.write("    if [ -f INCAR ] && [ ! -f INCAR.aimd ]; then\n")
        f.write("        cp INCAR INCAR.aimd\n")
        f.write("        echo '[INFO] 已备份 INCAR → INCAR.aimd'\n")
        f.write("    fi\n")
        f.write("    cp INCAR.relax INCAR\n")
        f.write("    echo '[INFO] 已设置 INCAR.relax → INCAR'\n")
        f.write("    echo '[INFO] 请运行 VASP 进行弛豫:'\n")
        f.write("    echo '        mpirun -np 8 vasp_std > vasp_relax.out 2>&1'\n")
        f.write("else\n")
        f.write("    echo '[ERROR] 未找到 INCAR.relax'\n")
        f.write("    exit 1\n")
        f.write("fi\n\n")
        
        f.write("# ================== Step 2: 弛豫后准备 ==================\n")
        f.write("# 弛豫完成后运行以下命令:\n")
        f.write("#   cp CONTCAR POSCAR\n")
        f.write("#   cp INCAR.aimd INCAR\n")
        f.write("#   # 然后运行 AIMD\n")
        f.write("#   mpirun -np 8 vasp_std > vasp_aimd.out 2>&1\n\n")
        
        f.write("echo ''\n")
        f.write("echo '================================================================'\n")
        f.write("echo '弛豫完成后的命令:'\n")
        f.write("echo '  cp CONTCAR POSCAR'\n")
        f.write("echo '  cp INCAR.aimd INCAR'\n")
        f.write("echo '  # 运行 AIMD'\n")
        f.write("echo '================================================================'\n")
    
    os.chmod(script_path, 0o755)


def write_model_meta(
    outdir: str,
    args: argparse.Namespace,
    n_atoms: int,
    center_idx: int,
    center_symbol: str,
    total_charge: int,
    element_counts: Dict[str, int],
    was_neutralized: bool,
    neutralized_count: int,
    selection_truncated: bool,
    cell_lengths: np.ndarray,
    density_original: Optional[float],
    density_target: Optional[float],
    density_achieved: float,
    n_cut_bonds: int,
    charge_warnings: List[str],
    clash_info: Optional[Dict[str, Any]] = None,
    neutralize_verified: bool = True
) -> None:
    """写入模型元数据"""
    meta = {
        "generator": f"setup_aimd_ase.py {VERSION}",
        "timestamp": datetime.now().isoformat(),
        "model_type": args.mode,
        "source_file": args.src,
        "center_atom": {
            "index": int(center_idx),
            "element": center_symbol
        },
        "radius_angstrom": args.radius,
        "selection_mode": args.selection,
        "bond_hops": args.bond_hops,
        "n_atoms": n_atoms,
        "estimated_charge": total_charge,
        "element_counts": element_counts,
        "cell_angstrom": [float(x) for x in cell_lengths],
        "pbc": [True, True, True],
        "has_vacuum": args.mode == "cluster",
        "vacuum_angstrom": args.vacuum if args.mode == "cluster" else None,
        "density": {
            "original_g_cm3": density_original,
            "target_g_cm3": density_target,
            "achieved_g_cm3": density_achieved
        },
        "neutralization": {
            "enabled": args.neutralize != "none",
            "method": args.neutralize,
            "target_charge": args.target_charge,
            "atoms_added": neutralized_count,
            "was_neutralized": was_neutralized,
            "verified": neutralize_verified
        },
        "cut_bonds": {
            "count": n_cut_bonds,
            "report_file": "cut_bonds_report.txt" if n_cut_bonds > 0 else None
        },
        "aimd_params": {
            "temperature_K": args.temp,
            "steps": args.steps,
            "potim_fs": args.potim,
            "thermostat": args.thermostat,
            "gamma_1ps": args.gamma_1ps
        },
        "warnings": []
    }
    
    # 添加碰撞检测信息
    if clash_info is not None:
        meta["clash_check"] = {
            "d_min_angstrom": clash_info['d_min'],
            "clash_count": clash_info['clash_count'],
            "has_clashes": clash_info['has_clashes'],
            "worst_pairs": clash_info['clash_pairs'][:10],
            "relaxation_recommended": clash_info['has_clashes']
        }
    
    # 添加警告
    if abs(total_charge) >= 1:
        meta["warnings"].append(f"体系可能带电 ({total_charge:+d})")
    if selection_truncated:
        meta["warnings"].append("molecule/bond_hops 模式被截断为 sphere")
    if n_cut_bonds > 0:
        meta["warnings"].append(f"检测到 {n_cut_bonds} 个切断的化学键")
    if args.mode == "cluster":
        meta["warnings"].append("CLUSTER 模式：有限簇模型，不能用于 bulk 输运性质！")
    if density_target and density_achieved:
        error = abs(density_achieved - density_target) / density_target * 100
        if error > 5:
            meta["warnings"].append(f"密度偏差 {error:.1f}%")
    if clash_info and clash_info['has_clashes']:
        meta["warnings"].append(f"检测到 {clash_info['clash_count']} 对原子碰撞！需要弛豫！")
    if not neutralize_verified and args.neutralize != 'none':
        meta["warnings"].append("中和后电荷验证失败，可能未完全中和")
    meta["warnings"].extend(charge_warnings)
    
    json_path = os.path.join(outdir, "model_meta.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)


def write_vasp_inputs(
    cluster: Atoms,
    outdir: str,
    temp: float,
    steps: int,
    potim: float,
    kpoints: Tuple[int, int, int],
    ncore: Optional[int],
    encut: float,
    is_bulk: bool,
    thermostat: str,
    gamma_1ps: float
) -> Tuple[bool, int]:
    """
    生成 VASP AIMD 输入文件
    
    v2.3.1: LANGEVIN_GAMMA 顺序从 POSCAR 元素行解析，确保一致性
    
    Returns:
        potcar_generated: 是否生成了 POTCAR
        n_element_types: 元素种类数 (NTYP)
    """
    os.makedirs(outdir, exist_ok=True)
    
    pp_path = os.environ.get('VASP_PP_PATH', '')
    has_pp = bool(pp_path and os.path.isdir(pp_path))
    potcar_generated = False
    
    # 检查含 H
    symbols = cluster.get_chemical_symbols()
    has_h = 'H' in symbols
    if has_h and potim > 1.5:
        print(f"[WARN] 体系含 H，POTIM={potim} fs 可能过大，建议 0.5-1.0 fs")
    
    poscar_path = os.path.join(outdir, 'POSCAR')
    
    if has_pp:
        try:
            from ase.calculators.vasp import Vasp
            
            mdalgo = 3 if thermostat == 'langevin' else 2
            
            calc_params = {
                'directory': outdir,
                'xc': 'PBE',
                'encut': encut,
                'prec': 'Normal',
                'algo': 'VeryFast',
                'ediff': 1e-5,
                'nelm': 200,
                'ismear': 0,
                'sigma': 0.05,
                'ibrion': 0,
                'mdalgo': mdalgo,
                'isym': 0,
                'tebeg': temp,
                'teend': temp,
                'potim': potim,
                'nsw': steps,
                'lreal': 'Auto',
                'lwave': False,
                'lcharg': False,
                'kpts': kpoints,
                'gamma': True,
            }
            
            if ncore:
                calc_params['ncore'] = ncore
            
            calc = Vasp(**calc_params)
            calc.write_input(cluster)
            potcar_generated = True
            
            # v2.3.1: 从 POSCAR 解析元素顺序，确保 LANGEVIN_GAMMA 一致
            if thermostat == 'langevin':
                poscar_elements = parse_poscar_element_order(poscar_path)
                n_element_types = len(poscar_elements)
                _append_langevin_to_incar(outdir, gamma_1ps, n_element_types, poscar_elements)
            else:
                poscar_elements = parse_poscar_element_order(poscar_path)
                n_element_types = len(poscar_elements)
            
        except Exception as e:
            print(f"[WARN] ASE Vasp 写入失败: {e}")
            has_pp = False
    
    if not has_pp:
        # 先写 POSCAR
        write(poscar_path, cluster, format='vasp')
        
        # v2.3.1: 从 POSCAR 解析元素顺序
        poscar_elements = parse_poscar_element_order(poscar_path)
        n_element_types = len(poscar_elements)
        
        _write_incar_manual(outdir, temp, steps, potim, encut, ncore, is_bulk, 
                          thermostat, gamma_1ps, n_element_types, poscar_elements)
        _write_kpoints_manual(outdir, kpoints)
        
        if not pp_path:
            print("\n[WARN] VASP_PP_PATH 未设置，未生成 POTCAR")
        print("[INFO] 请手动准备 POTCAR")
    
    write(os.path.join(outdir, 'cluster_visual.xyz'), cluster, format='xyz')
    
    return potcar_generated, n_element_types


def _append_langevin_to_incar(
    outdir: str, 
    gamma_1ps: float, 
    n_element_types: int,
    poscar_elements: Optional[List[str]] = None
):
    """
    追加 Langevin 参数到 INCAR
    
    v2.3.1: 元素顺序从 POSCAR 第 6 行解析，确保与 VASP 一致
    """
    incar_path = os.path.join(outdir, 'INCAR')
    
    element_str = " ".join(poscar_elements) if poscar_elements else f"(NTYP={n_element_types})"
    
    with open(incar_path, 'a') as f:
        f.write("\n# Langevin thermostat parameters\n")
        f.write(f"# LANGEVIN_GAMMA: per-element-type, order follows POSCAR: {element_str}\n")
        f.write(f"LANGEVIN_GAMMA = {' '.join([str(gamma_1ps)] * n_element_types)}\n")


def _write_incar_manual(
    outdir: str,
    temp: float,
    steps: int,
    potim: float,
    encut: float,
    ncore: Optional[int],
    is_bulk: bool,
    thermostat: str,
    gamma_1ps: float,
    n_element_types: int,
    poscar_elements: Optional[List[str]] = None
) -> None:
    """
    手动写入 INCAR
    
    v2.3.1: LANGEVIN_GAMMA 顺序与 POSCAR 元素行严格一致
    """
    incar_path = os.path.join(outdir, 'INCAR')
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    mode_str = "BULK（周期性子胞）" if is_bulk else "CLUSTER（真空簇）"
    
    # MDALGO: 2=Nosé-Hoover, 3=Langevin
    mdalgo = 3 if thermostat == 'langevin' else 2
    mdalgo_comment = "Langevin 恒温器" if mdalgo == 3 else "Nosé-Hoover 恒温器"
    
    element_str = " ".join(poscar_elements) if poscar_elements else f"NTYP={n_element_types}"
    
    with open(incar_path, 'w') as f:
        f.write(f"# INCAR for AIMD - Generated by setup_aimd_ase.py {VERSION}\n")
        f.write(f"# {timestamp}\n")
        f.write(f"# Mode: {mode_str}\n")
        f.write(f"# Thermostat: {thermostat} (MDALGO={mdalgo})\n")
        f.write(f"# Temperature: {temp} K, Steps: {steps}, POTIM: {potim} fs\n")
        
        if not is_bulk:
            f.write("#\n")
            f.write("# ⚠️ WARNING: CLUSTER 模式 - 有限簇模型！\n")
            f.write("# ⚠️ 扩散系数等输运性质不能与 bulk 结果比较！\n")
            f.write("#\n")
        
        f.write("\n# ============ 基础参数 ============\n")
        f.write("PREC = Normal\n")
        f.write(f"ENCUT = {encut}\n")
        f.write("ALGO = VeryFast\n")
        f.write("EDIFF = 1E-5\n")
        f.write("NELM = 200\n")
        f.write("LREAL = Auto\n\n")
        
        f.write("# ============ 展宽 ============\n")
        f.write("ISMEAR = 0\n")
        f.write("SIGMA = 0.05\n\n")
        
        f.write("# ============ 分子动力学 ============\n")
        f.write("IBRION = 0      # MD 模式\n")
        f.write(f"MDALGO = {mdalgo}      # {mdalgo_comment}\n")
        f.write("ISYM = 0        # AIMD 必须关闭对称性！\n")
        f.write(f"TEBEG = {temp}\n")
        f.write(f"TEEND = {temp}\n")
        f.write(f"POTIM = {potim}  # fs，含 H 建议 0.5-1.0\n")
        f.write(f"NSW = {steps}\n")
        
        if thermostat == 'langevin':
            f.write(f"\n# Langevin 恒温器摩擦系数 (1/ps)\n")
            f.write(f"# 元素顺序与 POSCAR 第 6 行一致: {element_str}\n")
            if gamma_1ps >= 10:
                f.write(f"# ⚠️ gamma={gamma_1ps} 较大，适合平衡段；生产段建议 1-5\n")
            f.write(f"LANGEVIN_GAMMA = {' '.join([str(gamma_1ps)] * n_element_types)}\n")
        else:
            f.write("SMASS = -3      # Nosé-Hoover\n")
        
        f.write("\n# ============ 输出控制 ============\n")
        f.write("LWAVE = .FALSE.\n")
        f.write("LCHARG = .FALSE.\n")
        
        if ncore:
            f.write(f"\nNCORE = {ncore}\n")


def _write_kpoints_manual(outdir: str, kpoints: Tuple[int, int, int]):
    """手动写入 KPOINTS"""
    with open(os.path.join(outdir, 'KPOINTS'), 'w') as f:
        f.write("Automatic mesh\n0\nGamma\n")
        f.write(f"{kpoints[0]} {kpoints[1]} {kpoints[2]}\n0 0 0\n")


def write_selected_indices(outdir: str, center_idx: int, selected_indices: np.ndarray):
    """写入选中索引"""
    with open(os.path.join(outdir, 'selected_indices.txt'), 'w') as f:
        f.write(f"# Center atom index (0-based): {center_idx}\n")
        f.write(f"# Total selected atoms: {len(selected_indices)}\n")
        for idx in selected_indices:
            f.write(f"{idx}\n")


# ==============================================================================
# 主函数
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=f"从大体系切割 AIMD 子体系 (setup_aimd_ase.py {VERSION})",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # bulk 模式（默认，按原体系密度定盒子）
  python3 setup_aimd_ase.py --src eq.pdb --center_atom Li --radius 8 --mode bulk

  # bulk 模式，指定目标密度
  python3 setup_aimd_ase.py --src eq.pdb --center_atom Li --radius 8 --density_g_cm3 1.2

  # cluster 模式（真空簇）
  python3 setup_aimd_ase.py --src eq.pdb --center_atom Li --radius 8 --mode cluster --vacuum 20

  # 使用 bond_hops 避免切断聚合物链
  python3 setup_aimd_ase.py --src eq.pdb --center_atom Li --radius 8 --bond_hops 3

注意:
  - 默认 bulk 模式按原体系密度计算盒子体积（避免"低压气相"陷阱）
  - cluster 模式必须显式指定 --mode cluster
  - 切断化学键会生成警告报告
        """
    )
    
    # 必需参数
    parser.add_argument("--src", required=True, help="输入结构文件")
    parser.add_argument("--center_atom", required=True, help="中心原子索引或元素")
    
    # 模式
    parser.add_argument("--mode", choices=['bulk', 'cluster'], default='bulk',
                        help="bulk=周期性(默认), cluster=真空簇")
    
    # 切割参数
    parser.add_argument("--radius", type=float, default=8.0, help="切割半径 Å")
    parser.add_argument("--selection", choices=['sphere', 'molecule'], default='sphere',
                        help="选择模式")
    parser.add_argument("--bond_hops", type=int, default=0,
                        help="键跳扩展步数（避免切断聚合物链）")
    parser.add_argument("--one_based", action="store_true", help="1-based 索引")
    
    # 盒子/密度参数
    parser.add_argument("--vacuum", type=float, default=20.0,
                        help="真空层 Å（仅 cluster）")
    parser.add_argument("--density_g_cm3", type=float, default=None,
                        help="目标密度 g/cm³（bulk 模式）")
    parser.add_argument("--cell_shape", choices=['scale_parent', 'cubic'], default='scale_parent',
                        help="盒子形状")
    
    # 电荷
    parser.add_argument("--neutralize", choices=['none', 'nearest_counterions'], default='none',
                        help="电荷中和")
    parser.add_argument("--target_charge", type=int, default=0, help="目标电荷")
    parser.add_argument("--charge_map_file", type=str, default=None,
                        help="自定义电荷映射文件")
    
    # 限制
    parser.add_argument("--max_atoms", type=int, default=400, help="最大原子数")
    parser.add_argument("--allow_exceed_max_atoms", action="store_true")
    
    # AIMD 参数
    parser.add_argument("--temp", type=float, default=350.0, help="温度 K")
    parser.add_argument("--steps", type=int, default=2000, help="步数")
    parser.add_argument("--potim", type=float, default=1.0, help="POTIM fs")
    parser.add_argument("--thermostat", choices=['langevin', 'nose'], default='langevin',
                        help="恒温器")
    parser.add_argument("--gamma_1ps", type=float, default=10.0,
                        help="Langevin gamma (1/ps)")
    parser.add_argument("--kpoints", default="1 1 1", help="K 点")
    parser.add_argument("--ncore", type=int, default=None)
    parser.add_argument("--encut", type=float, default=400.0)
    
    # 输出
    parser.add_argument("--outdir", default="aimd_sub")
    parser.add_argument("--overwrite", action="store_true")
    
    # v2.3 新增: 碰撞检测与弛豫输入
    parser.add_argument("--write_relax_inputs", action="store_true",
                        help="生成弛豫辅助文件 (INCAR.relax, RELAX_GUIDE.txt)")
    parser.add_argument("--clash_threshold_scale", type=float, default=0.75,
                        help="碰撞检测阈值比例 (默认 0.75)")
    parser.add_argument("--force_density", action="store_true",
                        help="强制使用可疑密度（跳过警告）")
    
    args = parser.parse_args()
    
    # 验证
    if args.mode == 'bulk' and args.vacuum != 20.0:
        print("[ERROR] --vacuum 仅在 cluster 模式有效")
        sys.exit(1)
    
    if args.mode == 'cluster' and args.density_g_cm3 is not None:
        print("[ERROR] --density_g_cm3 仅在 bulk 模式有效")
        sys.exit(1)
    
    print("=" * 70)
    print(f"setup_aimd_ase.py {VERSION} - AIMD 子体系切割")
    print("=" * 70)
    
    mode_desc = "BULK（周期性子胞，按密度定盒子）" if args.mode == "bulk" else "CLUSTER（真空簇）"
    print(f"模式: {mode_desc}")
    
    if args.mode == "cluster":
        print("")
        print("!" * 70)
        print("!!! 警告: CLUSTER 模式 - 有限真空簇！")
        print("!!! 扩散系数等输运性质不能与 bulk 结果比较！")
        print("!" * 70)
    
    print(f"输入: {args.src}")
    print(f"中心: {args.center_atom}, 半径: {args.radius} Å")
    print(f"选择: {args.selection}, bond_hops: {args.bond_hops}")
    if args.density_g_cm3:
        print(f"目标密度: {args.density_g_cm3} g/cm³")
    print("=" * 70)
    
    # 检查文件
    if not os.path.isfile(args.src):
        print(f"[ERROR] 文件不存在: {args.src}")
        sys.exit(1)
    
    # 检查输出目录
    if os.path.exists(args.outdir):
        if args.overwrite:
            backup = f"{args.outdir}.bak.{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            shutil.move(args.outdir, backup)
        else:
            print(f"[ERROR] 目录已存在: {args.outdir}，使用 --overwrite")
            sys.exit(1)
    
    # 验证 K 点
    try:
        kpts = validate_kpoints(args.kpoints)
    except ValueError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
    
    # 读取结构
    print("\n>>> 读取结构...")
    try:
        atoms = safe_read_structure(args.src)
        print(f"    原子数: {len(atoms)}")
        print(f"    元素: {set(atoms.get_chemical_symbols())}")
        
        original_density = None
        parent_cell = None
        original_volume = None
        
        # 使用安全的晶胞检查
        if has_valid_cell(atoms):
            atoms.wrap()
            parent_cell = atoms.cell.array.copy()
            original_volume = abs(np.linalg.det(parent_cell))
            original_density = compute_density(atoms.get_chemical_symbols(), original_volume)
            print(f"    原体系密度: {original_density:.4f} g/cm³")
            
            # v2.3: 密度健全性检查
            if args.mode == 'bulk' and args.density_g_cm3 is None:
                if original_density < 0.5 or original_density > 3.0:
                    print("")
                    print("!" * 70)
                    print(f"!!! 警告: 原体系密度 {original_density:.4f} g/cm³ 异常！")
                    if original_density < 0.5:
                        print("!!! 密度过低，可能是气相或未正确设置晶胞")
                    else:
                        print("!!! 密度过高，可能是晶胞设置错误")
                    print("!!! 强烈建议使用 --density_g_cm3 显式指定目标密度")
                    print("!" * 70)
                    if not args.force_density:
                        print("[ERROR] 使用 --force_density 强制继续，或指定 --density_g_cm3")
                        sys.exit(1)
        else:
            print("    [WARN] 无 PBC/cell")
            if args.mode == 'bulk' and args.density_g_cm3 is None:
                print("[ERROR] bulk 模式需要 --density_g_cm3（无法从输入获取）")
                sys.exit(1)
    except Exception as e:
        print(f"[ERROR] 读取失败: {e}")
        sys.exit(1)
    
    # 解析中心
    print("\n>>> 解析中心原子...")
    try:
        center_idx = parse_center_atom(args.center_atom, atoms, args.one_based)
        center_sym = atoms.get_chemical_symbols()[center_idx]
        print(f"    中心: {center_sym} (索引 {center_idx})")
    except ValueError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
    
    # 选择原子（带 MIC 向量）
    print("\n>>> 选择原子...")
    selected_indices, all_distances, mic_vectors = select_indices_with_mic_vectors(
        atoms, center_idx, args.radius
    )
    print(f"    初步选中: {len(selected_indices)} 原子")
    
    # 扩展选择
    selected_indices, truncated, n_cut_bonds = expand_selection_with_molecules(
        atoms, selected_indices, args.selection, args.bond_hops,
        args.max_atoms, args.allow_exceed_max_atoms
    )
    
    if truncated:
        print("[WARN] 选择被截断到 sphere 模式")
    
    # 电荷中和
    neutralized = False
    neutralized_count = 0
    neutralize_info = []
    neutralize_verified = True
    
    if args.neutralize == 'nearest_counterions':
        print("\n>>> 电荷中和...")
        temp_cluster = atoms[selected_indices]
        current_charge, _, _, _ = estimate_charge_comprehensive(temp_cluster)
        print(f"    当前电荷: {current_charge:+d}, 目标: {args.target_charge:+d}")
        
        if current_charge != args.target_charge:
            selected_indices, neutralized_count, neutralize_info = neutralize_by_residue(
                atoms, selected_indices, center_idx,
                current_charge, args.target_charge,
                all_distances, mic_vectors
            )
            if neutralized_count > 0:
                neutralized = True
                print(f"    添加 {neutralized_count} 原子")
                for info in neutralize_info:
                    print(f"    {info}")
                
                # v2.3: 中和验证
                temp_cluster_verify = atoms[selected_indices]
                verify_charge, _, _, _ = estimate_charge_comprehensive(temp_cluster_verify)
                if verify_charge != args.target_charge:
                    neutralize_verified = False
                    print(f"    [WARN] 中和验证: 电荷 {verify_charge:+d} ≠ 目标 {args.target_charge:+d}")
                    print(f"    [WARN] 可能需要更多反离子或调整选择半径")
                else:
                    print(f"    [OK] 中和验证通过: 电荷 = {verify_charge:+d}")
    
    print(f"\n>>> 最终: {len(selected_indices)} 原子")
    
    # 检查限制
    if len(selected_indices) > args.max_atoms and not args.allow_exceed_max_atoms:
        print(f"[ERROR] 超过 max_atoms ({args.max_atoms})")
        sys.exit(1)
    
    # 重新成像
    print("\n>>> 重新成像原子到中心附近...")
    cluster = reimage_atoms_around_center(atoms, selected_indices, center_idx, mic_vectors)
    
    # 估算电荷
    total_charge, elem_counts, charge_reliable, charge_warnings = estimate_charge_comprehensive(
        cluster, args.charge_map_file
    )
    print(f"    估算电荷: {total_charge:+d}")
    if not charge_reliable:
        print("    [WARN] 电荷估计可能不准确")
    
    # 创建盒子
    print("\n>>> 创建盒子...")
    if args.mode == 'bulk':
        cluster, target_density, achieved_density = create_density_based_bulk_box(
            cluster, original_density, args.density_g_cm3,
            args.cell_shape, parent_cell
        )
        print(f"    目标密度: {target_density:.4f} g/cm³")
        print(f"    实际密度: {achieved_density:.4f} g/cm³")
    else:
        cluster = create_vacuum_box(cluster, args.vacuum)
        achieved_density = 0.0
        target_density = None
        print(f"    真空层: {args.vacuum} Å")
    
    cell = cluster.get_cell().lengths()
    print(f"    盒子: {cell[0]:.1f} x {cell[1]:.1f} x {cell[2]:.1f} Å")
    
    # 切断键检测
    if n_cut_bonds > 0:
        print(f"\n[WARN] 检测到 {n_cut_bonds} 个切断的化学键！")
        print("[INFO] 建议: 增大半径 / 使用 --bond_hops / molecule 模式")
        
        # 写入报告
        if HAS_UTILS:
            os.makedirs(args.outdir, exist_ok=True)
            positions = atoms.get_positions()
            symbols = atoms.get_chemical_symbols()
            # v2.3.1: 使用统一的 has_valid_cell() 判定
            valid_cell = has_valid_cell(atoms)
            cell_arr = atoms.cell if valid_cell else None
            pbc_arr = atoms.pbc if valid_cell else None
            graph = build_bond_graph(positions, symbols, cell_arr, pbc_arr)
            cut_bonds = detect_cut_bonds(graph, set(selected_indices))
            write_cut_bonds_report(
                os.path.join(args.outdir, 'cut_bonds_report.txt'),
                cut_bonds, symbols, positions, cell_arr
            )
    
    # v2.3: 碰撞检测
    print("\n>>> 碰撞检测...")
    clash_info = detect_atomic_clashes(cluster, scale=args.clash_threshold_scale)
    
    if clash_info['d_min'] is not None:
        print(f"    最小原子间距: {clash_info['d_min']:.3f} Å")
    
    # 计算体积压缩比例
    volume_compression_ratio = 0.0
    if original_volume is not None and args.mode == 'bulk':
        current_volume = abs(np.linalg.det(cluster.get_cell()))
        # 根据原子数比例估算等效原始体积
        n_ratio = len(cluster) / len(atoms)
        equivalent_original_vol = original_volume * n_ratio
        if equivalent_original_vol > 0:
            volume_compression_ratio = max(0, (equivalent_original_vol - current_volume) / equivalent_original_vol)
    
    if clash_info['has_clashes']:
        print("")
        print("!" * 70)
        print(f"!!! 检测到 {clash_info['clash_count']} 对原子碰撞/重叠！")
        print("!!! 直接运行 AIMD 可能导致模拟崩溃！")
        print("!!! 必须先进行离子弛豫！")
        print("!" * 70)
        print(f"    最差的 5 对原子:")
        for i, cp in enumerate(clash_info['clash_pairs'][:5]):
            print(f"      {i+1}. 原子 {cp['i']} ({cp['sym_i']}) - "
                  f"原子 {cp['j']} ({cp['sym_j']}): {cp['distance']:.3f} Å")
    elif volume_compression_ratio > 0.3:
        print(f"    [WARN] 体积压缩 {volume_compression_ratio*100:.1f}%，建议弛豫")
    else:
        print(f"    [OK] 未检测到严重碰撞")
    
    # 生成 VASP 输入
    print("\n>>> 生成 VASP 输入...")
    potcar_ok, n_element_types = write_vasp_inputs(
        cluster, args.outdir, args.temp, args.steps, args.potim,
        kpts, args.ncore, args.encut, args.mode == 'bulk',
        args.thermostat, args.gamma_1ps
    )
    
    # v2.3: 弛豫指导与辅助输入
    has_h = 'H' in cluster.get_chemical_symbols()
    need_relax_guide = clash_info['has_clashes'] or volume_compression_ratio > 0.3
    
    if need_relax_guide or args.write_relax_inputs:
        print("\n>>> 生成弛豫指导...")
        write_relax_guide(args.outdir, clash_info, has_h, args.gamma_1ps, volume_compression_ratio)
        print(f"    [OK] RELAX_GUIDE.txt")
    
    if args.write_relax_inputs or clash_info['has_clashes']:
        write_incar_relax(args.outdir, args.encut, args.ncore, has_h, n_element_types)
        write_relax_script(args.outdir)
        print(f"    [OK] INCAR.relax")
        print(f"    [OK] run_relax_then_aimd.sh")
    
    # 元数据
    write_model_meta(
        args.outdir, args, len(cluster), center_idx, center_sym,
        total_charge, elem_counts, neutralized, neutralized_count,
        truncated, cell, original_density, target_density, achieved_density,
        n_cut_bonds, charge_warnings, clash_info, neutralize_verified
    )
    
    write_selected_indices(args.outdir, center_idx, selected_indices)
    
    # 摘要
    print("\n" + "=" * 70)
    print("完成！")
    print("=" * 70)
    print(f"模式: {mode_desc}")
    print(f"原子数: {len(cluster)}")
    print(f"电荷: {total_charge:+d}")
    if args.mode == 'bulk':
        print(f"密度: {achieved_density:.4f} g/cm³")
    if clash_info['d_min'] is not None:
        print(f"最小原子间距: {clash_info['d_min']:.3f} Å")
    if clash_info['has_clashes']:
        print(f"碰撞原子对: {clash_info['clash_count']} 对 ⚠️")
    if n_cut_bonds > 0:
        print(f"切断键: {n_cut_bonds} 个（见 cut_bonds_report.txt）")
    print(f"输出: {args.outdir}")
    print("=" * 70)
    
    if args.mode == "bulk":
        print("\n>>> 【周期性子胞】模型")
        print("    ✓ 按原体系密度定盒子，物理合理")
        print("    ✓ 可用于局域结构分析")
        print("    ⚠️ AIMD 时间有限（ps），长程扩散用经典 MD")
    else:
        print("\n>>> 【有限真空簇】模型")
        print("    ⚠️ 不能用于 bulk 输运性质！")
        print("    ⚠️ 表面效应显著")
    
    if clash_info['has_clashes']:
        print("")
        print("!" * 70)
        print("!!! 重要: 检测到原子碰撞，必须先弛豫！")
        print("!!! 请查看 RELAX_GUIDE.txt 和 INCAR.relax")
        print("!" * 70)
    
    if args.gamma_1ps >= 10:
        print(f"\n[WARN] Langevin gamma={args.gamma_1ps} 较大")
        print("[INFO] 适合平衡段；生产段建议 gamma=1-5")


if __name__ == "__main__":
    main()
