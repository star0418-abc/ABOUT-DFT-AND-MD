#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
recipe_to_counts.py - 配方 wt% 转换为分子/原子数量

功能：
  - 读取 recipe.yaml 配方文件
  - 将 wt% 换算为整数 entity 数量（分子数）
  - 支持两种模式：
    1) --total_mass_g: 按总质量换算
    2) --target_atoms: 按目标总原子数换算（适合 VASP AIMD）
  - 输出 counts.csv 和 counts.json

用法：
  # 模式1: 按总质量 (需指定缩放目标)
  python3 recipe_to_counts.py --total_mass_g 1.0 --scale_to_molecules 100
  python3 recipe_to_counts.py --total_mass_g 1.0 --scale_to_atoms 5000

  # 模式2: 按目标原子数
  python3 recipe_to_counts.py --target_atoms 5000

作者：STAR0418-ABC
"""

import argparse
import sys
import os
import json
import csv
from typing import Dict, List, Any, Optional, Tuple

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

# 可选 numpy
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

# 固定的 8 类顺序
CATEGORY_ORDER = [
    ("salt_solution", "盐溶液"),
    ("polymer_matrix", "聚合物基质"),
    ("crosslinker", "交联剂"),
    ("photoinitiator", "引发剂"),
    ("plasticizer_solvent", "增塑剂/溶剂"),
    ("functional_monomer", "功能单体"),
    ("stabilizer", "稳定剂"),
    ("functional_filler", "功能填料"),
]


def load_yaml(filepath: str) -> Dict[str, Any]:
    """加载 YAML 文件"""
    if not os.path.isfile(filepath):
        print(f"[ERROR] 配方文件不存在: {filepath}")
        sys.exit(1)

    if not HAS_YAML:
        print("[ERROR] 需要 PyYAML 库: pip install pyyaml")
        sys.exit(1)

    with open(filepath, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)

    return data if data else {}


def flatten_recipe(data: Dict) -> List[Dict]:
    """
    将配方数据展平为条目列表，按固定顺序

    返回: 条目列表，每个条目包含 category 字段
    """
    entries = []

    for cat_key, cat_name in CATEGORY_ORDER:
        cat_entries = data.get(cat_key, []) or []

        for entry in cat_entries:
            if not isinstance(entry, dict):
                continue

            item = entry.copy()
            item['category'] = cat_key
            item['category_cn'] = cat_name
            entries.append(item)

    return entries


def compute_moles_from_mass(entries: List[Dict], total_mass_g: float) -> List[Dict]:
    """
    根据总质量计算每个组分的物质的量 (mol)

    公式: n_i = (wt_pct/100) * total_mass_g / mw_g_mol
    """
    results = []

    for entry in entries:
        result = entry.copy()
        wt_pct = entry.get('wt_pct', 0)
        mw = entry.get('mw_g_mol', None)

        if mw and mw > 0:
            mass_g = (wt_pct / 100.0) * total_mass_g
            moles = mass_g / mw
            result['mass_g'] = mass_g
            result['moles'] = moles
        else:
            result['mass_g'] = None
            result['moles'] = None
            result['skip_reason'] = "缺少 mw_g_mol"

        results.append(result)

    return results


def scale_to_molecules(entries: List[Dict], target_molecules: int) -> List[Dict]:
    """
    按目标分子总数缩放

    所有有 moles 的条目按比例缩放到整数分子数
    """
    # 计算总 moles
    total_moles = sum(e.get('moles', 0) or 0 for e in entries)

    if total_moles <= 0:
        print("[ERROR] 无法计算总物质的量，请检查 mw_g_mol 是否完整")
        sys.exit(1)

    results = []
    for entry in entries:
        result = entry.copy()
        moles = entry.get('moles', None)

        if moles is not None and moles > 0:
            # 按比例缩放
            fraction = moles / total_moles
            scaled = fraction * target_molecules
            result['scaled_count'] = max(1, round(scaled))  # 至少 1 个
        else:
            result['scaled_count'] = None

        results.append(result)

    return results


def scale_to_atoms(entries: List[Dict], target_atoms: int) -> List[Dict]:
    """
    按目标原子总数缩放

    需要 atoms_per_entity 字段
    """
    # 计算加权原子数（moles * atoms_per_entity）
    total_weighted = 0.0
    valid_entries = []

    for entry in entries:
        moles = entry.get('moles', None)
        atoms = entry.get('atoms_per_entity', None) or entry.get('atoms_per_molecule', None)

        if moles and atoms and moles > 0 and atoms > 0:
            total_weighted += moles * atoms
            valid_entries.append((entry, moles, atoms))

    if total_weighted <= 0:
        print("[ERROR] 无法计算加权原子数，请检查 atoms_per_entity 是否完整")
        sys.exit(1)

    # 缩放因子
    scale_factor = target_atoms / total_weighted

    results = []
    for entry in entries:
        result = entry.copy()
        moles = entry.get('moles', None)
        atoms = entry.get('atoms_per_entity', None) or entry.get('atoms_per_molecule', None)

        if moles and atoms and moles > 0 and atoms > 0:
            # 缩放后的分子数
            scaled_moles = moles * scale_factor
            scaled_count = max(1, round(scaled_moles))
            result['scaled_count'] = scaled_count
            result['scaled_atoms'] = scaled_count * atoms
        else:
            result['scaled_count'] = None
            result['scaled_atoms'] = None
            if 'skip_reason' not in result:
                result['skip_reason'] = "缺少 atoms_per_entity"

        results.append(result)

    return results


def compute_by_target_atoms(entries: List[Dict], target_atoms: int) -> List[Dict]:
    """
    模式2: 直接按目标原子数计算

    不需要 total_mass_g，直接从 wt_pct 比例计算
    需要所有条目提供 mw_g_mol 和 atoms_per_entity
    """
    # 步骤1: 假设总质量为 1g，计算 moles
    results = []
    total_weighted = 0.0

    for entry in entries:
        result = entry.copy()
        wt_pct = entry.get('wt_pct', 0)
        mw = entry.get('mw_g_mol', None)
        atoms = entry.get('atoms_per_entity', None) or entry.get('atoms_per_molecule', None)

        if mw and mw > 0 and atoms and atoms > 0:
            # 假设总质量 1g
            mass_g = wt_pct / 100.0
            moles = mass_g / mw
            weighted = moles * atoms
            total_weighted += weighted

            result['moles'] = moles
            result['weighted_atoms'] = weighted
        else:
            result['moles'] = None
            result['weighted_atoms'] = None
            skip_reasons = []
            if not mw:
                skip_reasons.append("mw_g_mol")
            if not atoms:
                skip_reasons.append("atoms_per_entity")
            result['skip_reason'] = f"缺少 {', '.join(skip_reasons)}"

        results.append(result)

    if total_weighted <= 0:
        print("[ERROR] 无法计算加权原子数，请检查数据完整性")
        sys.exit(1)

    # 步骤2: 按比例缩放到目标原子数
    scale_factor = target_atoms / total_weighted

    for result in results:
        moles = result.get('moles', None)
        atoms = result.get('atoms_per_entity', None) or result.get('atoms_per_molecule', None)

        if moles and atoms:
            scaled_moles = moles * scale_factor
            scaled_count = max(1, round(scaled_moles))
            result['scaled_count'] = scaled_count
            result['scaled_atoms'] = scaled_count * atoms
        else:
            result['scaled_count'] = None
            result['scaled_atoms'] = None

    return results


def write_csv(results: List[Dict], filepath: str):
    """输出 CSV 文件"""
    fieldnames = [
        'category', 'kind', 'name', 'wt_pct',
        'mw_g_mol', 'moles', 'scaled_count',
        'atoms_per_entity', 'scaled_atoms', 'structure_file', 'skip_reason'
    ]

    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()

        for r in results:
            row = {k: r.get(k, '') for k in fieldnames}
            # 格式化数值
            if row['moles'] and isinstance(row['moles'], float):
                row['moles'] = f"{row['moles']:.6e}"
            writer.writerow(row)


def write_json(results: List[Dict], filepath: str):
    """输出 JSON 文件"""
    # 清理结果，只保留关键字段
    output = []
    for r in results:
        item = {
            'category': r.get('category', ''),
            'kind': r.get('kind', ''),
            'name': r.get('name', ''),
            'wt_pct': r.get('wt_pct', 0),
            'mw_g_mol': r.get('mw_g_mol'),
            'moles': r.get('moles'),
            'scaled_count': r.get('scaled_count'),
            'atoms_per_entity': r.get('atoms_per_entity') or r.get('atoms_per_molecule'),
            'scaled_atoms': r.get('scaled_atoms'),
            'structure_file': r.get('structure_file'),
        }
        if 'skip_reason' in r:
            item['skip_reason'] = r['skip_reason']
        output.append(item)

    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(output, f, ensure_ascii=False, indent=2)


def print_summary(results: List[Dict], mode: str):
    """打印摘要"""
    print("\n" + "=" * 90)
    print(f"换算结果摘要 (模式: {mode})")
    print("=" * 90)

    # 按类别分组
    current_cat = None
    total_count = 0
    total_atoms = 0
    skipped = []

    print(f"\n{'类别':<20} {'名称':<30} {'wt%':>8} {'count':>10} {'atoms':>10}")
    print("-" * 90)

    for r in results:
        cat = r.get('category', '')
        if cat != current_cat:
            current_cat = cat
            cat_cn = r.get('category_cn', cat)
            print(f"\n>>> {cat} ({cat_cn})")

        name = r.get('name', 'N/A')
        if len(name) > 28:
            name = name[:25] + "..."

        wt = r.get('wt_pct', 0)
        count = r.get('scaled_count', None)
        atoms = r.get('scaled_atoms', None)

        count_str = str(count) if count is not None else "---"
        atoms_str = str(atoms) if atoms is not None else "---"

        print(f"    {name:<30} {wt:>8.2f} {count_str:>10} {atoms_str:>10}")

        if count is not None:
            total_count += count
        if atoms is not None:
            total_atoms += atoms

        if 'skip_reason' in r:
            skipped.append((r.get('name', 'N/A'), r['skip_reason']))

    print("\n" + "-" * 90)
    print(f"{'总计:':<50} {' ':>8} {total_count:>10} {total_atoms:>10}")

    if skipped:
        print("\n[INFO] 以下条目未参与 counts 换算:")
        for name, reason in skipped:
            print(f"  - {name}: {reason}")

    print("=" * 90)


def main():
    parser = argparse.ArgumentParser(
        description="配方 wt% 转换为分子/原子数量",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 模式1: 按总质量换算
  python3 recipe_to_counts.py --total_mass_g 1.0 --scale_to_molecules 100
  python3 recipe_to_counts.py --total_mass_g 1.0 --scale_to_atoms 5000

  # 模式2: 按目标原子数（适合 VASP AIMD）
  python3 recipe_to_counts.py --target_atoms 5000

  # 指定配方文件和输出
  python3 recipe_to_counts.py --recipe my_recipe.yaml --target_atoms 3000 --output my_counts
        """
    )

    parser.add_argument("--recipe", default="recipe.yaml",
                        help="配方文件路径 (默认: recipe.yaml)")
    parser.add_argument("--output", default="counts",
                        help="输出文件前缀 (默认: counts -> counts.csv, counts.json)")

    # 模式1: 按总质量
    parser.add_argument("--total_mass_g", type=float,
                        help="模式1: 总质量 (g)")
    parser.add_argument("--scale_to_molecules", type=int,
                        help="模式1: 缩放到目标分子总数")
    parser.add_argument("--scale_to_atoms", type=int,
                        help="模式1: 缩放到目标原子总数")

    # 模式2: 按目标原子数
    parser.add_argument("--target_atoms", type=int,
                        help="模式2: 目标总原子数（直接换算）")

    args = parser.parse_args()

    # 验证参数
    mode = None
    if args.target_atoms:
        mode = "target_atoms"
    elif args.total_mass_g:
        if args.scale_to_molecules:
            mode = "total_mass_molecules"
        elif args.scale_to_atoms:
            mode = "total_mass_atoms"
        else:
            print("[ERROR] --total_mass_g 需要配合 --scale_to_molecules 或 --scale_to_atoms")
            sys.exit(1)
    else:
        print("[ERROR] 请指定换算模式:")
        print("  模式1: --total_mass_g M --scale_to_molecules N")
        print("  模式1: --total_mass_g M --scale_to_atoms N")
        print("  模式2: --target_atoms N")
        sys.exit(1)

    print("=" * 90)
    print("recipe_to_counts.py - 配方换算工具")
    print("=" * 90)
    print(f"配方文件: {args.recipe}")
    print(f"换算模式: {mode}")

    if mode == "target_atoms":
        print(f"目标原子数: {args.target_atoms}")
    elif mode == "total_mass_molecules":
        print(f"总质量: {args.total_mass_g} g")
        print(f"目标分子数: {args.scale_to_molecules}")
    elif mode == "total_mass_atoms":
        print(f"总质量: {args.total_mass_g} g")
        print(f"目标原子数: {args.scale_to_atoms}")

    # 加载配方
    data = load_yaml(args.recipe)

    # 展平
    entries = flatten_recipe(data)

    if not entries:
        print("[ERROR] 配方为空")
        sys.exit(1)

    print(f"\n>>> 读取 {len(entries)} 个条目")

    # 执行换算
    if mode == "target_atoms":
        results = compute_by_target_atoms(entries, args.target_atoms)

    elif mode == "total_mass_molecules":
        results = compute_moles_from_mass(entries, args.total_mass_g)
        results = scale_to_molecules(results, args.scale_to_molecules)

    elif mode == "total_mass_atoms":
        results = compute_moles_from_mass(entries, args.total_mass_g)
        results = scale_to_atoms(results, args.scale_to_atoms)

    # 输出
    csv_path = f"{args.output}.csv"
    json_path = f"{args.output}.json"

    print(f"\n>>> 输出文件:")
    write_csv(results, csv_path)
    print(f"    [OK] {csv_path}")

    write_json(results, json_path)
    print(f"    [OK] {json_path}")

    # 打印摘要
    print_summary(results, mode)

    # 计算实际总原子数
    actual_atoms = sum(r.get('scaled_atoms', 0) or 0 for r in results)
    print(f"\n[INFO] 实际总原子数: {actual_atoms}")

    if mode == "target_atoms":
        diff = abs(actual_atoms - args.target_atoms)
        print(f"[INFO] 与目标差异: {diff} 原子 ({diff/args.target_atoms*100:.2f}%)")


if __name__ == "__main__":
    main()

