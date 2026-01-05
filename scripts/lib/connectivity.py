# -*- coding: utf-8 -*-

# Author: Star



"""
connectivity.py - 分子连通性检查模块
=====================================

提供 MOL2 文件的连通性验证功能。

用法:
    from lib.connectivity import check_mol2_connectivity
    
    result = check_mol2_connectivity("output.mol2")
    if result['num_fragments'] != 1:
        raise RuntimeError("分子不是单一连通体")
"""

import os
from collections import defaultdict, deque
from typing import Dict, List, Tuple, Optional


class ConnectivityError(Exception):
    """连通性错误"""
    pass


def parse_mol2_structure(filepath: str) -> Tuple[List[Tuple[int, str, Tuple[float, float, float]]], List[Tuple[int, int, str]]]:
    """
    解析 MOL2 文件的原子和键信息
    
    Returns:
        (atoms, bonds):
            atoms: [(atom_id, atom_name, (x, y, z)), ...]
            bonds: [(atom1_id, atom2_id, bond_type), ...]
    """
    atoms = []
    bonds = []
    
    in_atom_section = False
    in_bond_section = False
    
    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            
            if stripped.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                in_bond_section = False
                continue
            elif stripped.startswith("@<TRIPOS>BOND"):
                in_atom_section = False
                in_bond_section = True
                continue
            elif stripped.startswith("@<TRIPOS>"):
                in_atom_section = False
                in_bond_section = False
                continue
            
            if in_atom_section and stripped:
                parts = stripped.split()
                if len(parts) >= 5:
                    atom_id = int(parts[0])
                    atom_name = parts[1]
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atoms.append((atom_id, atom_name, (x, y, z)))
            
            if in_bond_section and stripped:
                parts = stripped.split()
                if len(parts) >= 4:
                    atom1 = int(parts[1])
                    atom2 = int(parts[2])
                    bond_type = parts[3]
                    bonds.append((atom1, atom2, bond_type))
    
    return atoms, bonds


def count_connected_components(atoms: List, bonds: List) -> Tuple[int, List[List[int]]]:
    """
    使用 BFS 计算连通分量数
    
    Returns:
        (num_components, components_list)
    """
    if not atoms:
        return 0, []
    
    atom_ids = {a[0] for a in atoms}
    
    # 构建邻接表
    adj = defaultdict(set)
    for a1, a2, _ in bonds:
        adj[a1].add(a2)
        adj[a2].add(a1)
    
    visited = set()
    components = []
    
    for start in sorted(atom_ids):
        if start in visited:
            continue
        
        # BFS
        component = []
        queue = deque([start])
        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            component.append(node)
            for neighbor in adj[node]:
                if neighbor not in visited:
                    queue.append(neighbor)
        
        components.append(sorted(component))
    
    return len(components), components


def check_mol2_connectivity(filepath: str, expected_atoms: int = 0) -> Dict:
    检查 MOL2 文件的连通性
    
    Args:
        filepath: MOL2 文件路径
        expected_atoms: 期望的原子数 (0 表示不检查)
    
    Returns:
        {
            'num_atoms': int,
            'num_bonds': int,
            'num_fragments': int,
            'is_connected': bool,  # 是否单一连通
            'fragments': [[atom_ids], ...],
        }
    
    Raises:
        ConnectivityError: 文件不存在或解析失败
    """
    if not os.path.exists(filepath):
        raise ConnectivityError(f"文件不存在: {filepath}")
    
    atoms, bonds = parse_mol2_structure(filepath)
    
    if not atoms:
        raise ConnectivityError(f"文件中没有原子: {filepath}")
    
    num_components, components = count_connected_components(atoms, bonds)
    
    result = {
        'num_atoms': len(atoms),
        'num_bonds': len(bonds),
        'num_fragments': num_components,
        'is_connected': num_components == 1,
        'fragments': components,
    }
    
    if expected_atoms > 0 and len(atoms) != expected_atoms:
        result['atom_mismatch'] = True
    
    return result


def validate_single_molecule(filepath: str) -> None:
    """
    验证 MOL2 文件是单一连通分子
    
    Raises:
        ConnectivityError: 如果不是单一连通分子
    """
    result = check_mol2_connectivity(filepath)
    
    if not result['is_connected']:
        fragments_str = "; ".join(
            f"碎片{i+1}: {len(frag)} 原子" 
            for i, frag in enumerate(result['fragments'])
        )
        raise ConnectivityError(
            f"分子不是单一连通体！发现 {result['num_fragments']} 个碎片。\n"
            f"  文件: {filepath}\n"
            f"  碎片: {fragments_str}"
        )
    
    print(f"  ✓ 连通性检查通过: 单一分子, {result['num_atoms']} 原子, {result['num_bonds']} 键")


def calculate_bond_lengths(filepath: str) -> List[Dict]:
    """
    计算所有键的长度
    
    Returns:
        [{'atom1': id, 'atom2': id, 'length': float, 'type': str}, ...]
    """
    import math
    
    atoms, bonds = parse_mol2_structure(filepath)
    
    # 创建坐标字典
    coords = {a[0]: a[2] for a in atoms}
    
    bond_lengths = []
    for a1, a2, btype in bonds:
        if a1 in coords and a2 in coords:
            c1, c2 = coords[a1], coords[a2]
            length = math.sqrt(sum((x - y) ** 2 for x, y in zip(c1, c2)))
            bond_lengths.append({
                'atom1': a1,
                'atom2': a2,
                'length': length,
                'type': btype
            })
    
    return bond_lengths


def check_bond_lengths(filepath: str, max_length: float = 3.0) -> List[Dict]:
    """
    检查异常长的键（可能表示碎片）
    
    Returns:
        异常键列表
    """
    bond_lengths = calculate_bond_lengths(filepath)
    
    abnormal = [b for b in bond_lengths if b['length'] > max_length]
    return abnormal

