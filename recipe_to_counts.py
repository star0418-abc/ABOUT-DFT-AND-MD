#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
recipe_to_counts.py - 配方 wt% 转换为分子/原子数量（含误差分析）

功能：
  - 读取 recipe.yaml 配方文件
  - 将 wt% 换算为整数 entity 数量
  - 使用 Largest Remainder 方法最小化误差
  - 输出目标/实际 wt% 对比与误差报告
  - 支持 min_count 强制最小数量

用法：
  python3 recipe_to_counts.py --target_atoms 200
  python3 recipe_to_counts.py --total_mass_g 1.0 --scale_to_atoms 5000

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
    """将配方数据展平为条目列表"""
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


def largest_remainder_allocation(
    entries: List[Dict],
    target_total: int,
    count_field: str = 'float_count',
    min_count_field: str = 'min_count'
) -> List[Dict]:
    """
    Largest Remainder 方法分配整数
    
    1. 先取 floor
    2. 强制满足 min_count
    3. 按小数部分从大到小补齐剩余配额
    """
    results = []
    
    # 第一步：计算浮点 count 并取 floor
    total_floor = 0
    for entry in entries:
        result = entry.copy()
        float_val = entry.get(count_field, 0) or 0
        min_val = entry.get(min_count_field, 0) or 0
        
        floor_val = int(float_val)
        # 强制满足 min_count
        floor_val = max(floor_val, min_val)
        
        result['floor_count'] = floor_val
        result['remainder'] = float_val - floor_val
        result['min_count'] = min_val
        total_floor += floor_val
        results.append(result)
    
    # 第二步：分配剩余配额
    remaining = target_total - total_floor
    
    if remaining > 0:
        # 按 remainder 从大到小排序
        sorted_indices = sorted(
            range(len(results)),
            key=lambda i: results[i]['remainder'],
            reverse=True
        )
        
        for i in sorted_indices:
            if remaining <= 0:
                break
            results[i]['floor_count'] += 1
            remaining -= 1
    
    # 设置最终 count
    for result in results:
        result['scaled_count'] = result['floor_count']
    
    return results


def compute_by_target_atoms(entries: List[Dict], target_atoms: int) -> List[Dict]:
    """按目标原子数计算（使用 Largest Remainder）"""
    
    # 计算加权原子数
    total_weighted = 0.0
    valid_entries = []
    
    for entry in entries:
        wt_pct = entry.get('wt_pct', 0)
        mw = entry.get('mw_g_mol', None)
        atoms = entry.get('atoms_per_entity', None) or entry.get('atoms_per_molecule', None)
        
        if mw and mw > 0 and atoms and atoms > 0:
            mass_g = wt_pct / 100.0
            moles = mass_g / mw
            weighted = moles * atoms
            total_weighted += weighted
            entry['moles'] = moles
            entry['weighted_atoms'] = weighted
            valid_entries.append(entry)
        else:
            entry['moles'] = None
            entry['weighted_atoms'] = None
            skip_reasons = []
            if not mw:
                skip_reasons.append("mw_g_mol")
            if not atoms:
                skip_reasons.append("atoms_per_entity")
            entry['skip_reason'] = f"缺少 {', '.join(skip_reasons)}"
    
    if total_weighted <= 0:
        print("[ERROR] 无法计算加权原子数")
        sys.exit(1)
    
    # 计算浮点 count
    scale_factor = target_atoms / total_weighted
    target_count = 0
    
    for entry in entries:
        moles = entry.get('moles', None)
        atoms = entry.get('atoms_per_entity', None) or entry.get('atoms_per_molecule', None)
        
        if moles and atoms:
            float_count = moles * scale_factor
            entry['float_count'] = float_count
            target_count += int(float_count)  # 估算目标分子数
        else:
            entry['float_count'] = 0
    
    # 使用 Largest Remainder 分配
    results = largest_remainder_allocation(entries, target_count + len(valid_entries))
    
    # 计算实际原子数和 wt%
    total_mass = 0.0
    for result in results:
        count = result.get('scaled_count', 0) or 0
        atoms = result.get('atoms_per_entity', None) or result.get('atoms_per_molecule', None)
        mw = result.get('mw_g_mol', None)
        
        if count > 0 and atoms:
            result['scaled_atoms'] = count * atoms
        else:
            result['scaled_atoms'] = None
        
        if count > 0 and mw:
            result['actual_mass'] = count * mw / 6.022e23 * 1e15  # 相对值
            total_mass += result['actual_mass']
        else:
            result['actual_mass'] = 0
    
    # 计算实际 wt% 和误差
    for result in results:
        actual_mass = result.get('actual_mass', 0)
        target_wt = result.get('wt_pct', 0)
        
        if total_mass > 0:
            actual_wt = actual_mass / total_mass * 100
        else:
            actual_wt = 0
        
        result['actual_wt_pct'] = actual_wt
        result['wt_error_abs'] = abs(actual_wt - target_wt)
        result['wt_error_rel'] = result['wt_error_abs'] / target_wt * 100 if target_wt > 0 else 0
    
    return results


def write_csv(results: List[Dict], filepath: str):
    """输出 CSV 文件"""
    fieldnames = [
        'category', 'kind', 'name', 'wt_pct', 'actual_wt_pct',
        'wt_error_abs', 'wt_error_rel',
        'mw_g_mol', 'scaled_count', 'min_count',
        'atoms_per_entity', 'scaled_atoms', 'structure_file', 'skip_reason'
    ]

    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()

        for r in results:
            row = {k: r.get(k, '') for k in fieldnames}
            # 格式化数值
            for key in ['actual_wt_pct', 'wt_error_abs']:
                if row[key] and isinstance(row[key], float):
                    row[key] = f"{row[key]:.2f}"
            if row['wt_error_rel'] and isinstance(row['wt_error_rel'], float):
                row['wt_error_rel'] = f"{row['wt_error_rel']:.1f}%"
            writer.writerow(row)


def write_json(results: List[Dict], filepath: str):
    """输出 JSON 文件"""
    output = []
    for r in results:
        item = {
            'category': r.get('category', ''),
            'kind': r.get('kind', ''),
            'name': r.get('name', ''),
            'wt_pct_target': r.get('wt_pct', 0),
            'wt_pct_actual': round(r.get('actual_wt_pct', 0), 2),
            'wt_error_rel_pct': round(r.get('wt_error_rel', 0), 1),
            'mw_g_mol': r.get('mw_g_mol'),
            'scaled_count': r.get('scaled_count'),
            'min_count': r.get('min_count', 0),
            'atoms_per_entity': r.get('atoms_per_entity') or r.get('atoms_per_molecule'),
            'scaled_atoms': r.get('scaled_atoms'),
            'structure_file': r.get('structure_file'),
        }
        if 'skip_reason' in r:
            item['skip_reason'] = r['skip_reason']
        output.append(item)

    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(output, f, ensure_ascii=False, indent=2)


def write_report(results: List[Dict], filepath: str, target_atoms: int):
    """输出误差报告"""
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write("=" * 90 + "\n")
        f.write("配方换算误差报告\n")
        f.write("=" * 90 + "\n\n")
        
        f.write(f"目标原子数: {target_atoms}\n")
        
        actual_atoms = sum(r.get('scaled_atoms', 0) or 0 for r in results)
        f.write(f"实际原子数: {actual_atoms}\n")
        f.write(f"差异: {abs(actual_atoms - target_atoms)} ({abs(actual_atoms - target_atoms)/target_atoms*100:.1f}%)\n\n")
        
        f.write("-" * 90 + "\n")
        f.write(f"{'名称':<35} {'目标wt%':>10} {'实际wt%':>10} {'误差':>10} {'count':>8}\n")
        f.write("-" * 90 + "\n")
        
        warnings = []
        for r in results:
            name = r.get('name', 'N/A')
            if len(name) > 33:
                name = name[:30] + "..."
            
            target = r.get('wt_pct', 0)
            actual = r.get('actual_wt_pct', 0)
            error = r.get('wt_error_rel', 0)
            count = r.get('scaled_count', 0) or 0
            
            marker = "⚠️ " if error > 50 else "   "
            f.write(f"{marker}{name:<32} {target:>10.2f} {actual:>10.2f} {error:>9.1f}% {count:>8}\n")
            
            if error > 50:
                warnings.append((name, error))
        
        f.write("-" * 90 + "\n\n")
        
        if warnings:
            f.write("⚠️ 警告: 以下组分相对误差 > 50%:\n")
            for name, error in warnings:
                f.write(f"   - {name}: {error:.1f}%\n")
            f.write("\n建议: 增大 target_atoms 以减小凑整误差\n")
        else:
            f.write("✓ 所有组分相对误差 < 50%\n")
        
        f.write("\n" + "=" * 90 + "\n")


def print_summary(results: List[Dict], target_atoms: int):
    """打印摘要"""
    print("\n" + "=" * 90)
    print("换算结果摘要（含误差分析）")
    print("=" * 90)

    print(f"\n{'类别':<20} {'名称':<30} {'目标%':>8} {'实际%':>8} {'误差':>8} {'count':>8}")
    print("-" * 90)

    current_cat = None
    total_count = 0
    total_atoms = 0
    warnings = []

    for r in results:
        cat = r.get('category', '')
        if cat != current_cat:
            current_cat = cat
            cat_cn = r.get('category_cn', cat)
            print(f"\n>>> {cat} ({cat_cn})")

        name = r.get('name', 'N/A')
        if len(name) > 28:
            name = name[:25] + "..."

        target = r.get('wt_pct', 0)
        actual = r.get('actual_wt_pct', 0)
        error = r.get('wt_error_rel', 0)
        count = r.get('scaled_count', None)
        atoms = r.get('scaled_atoms', None)

        count_str = str(count) if count is not None else "---"
        error_str = f"{error:.1f}%" if error > 0 else "---"
        
        marker = "⚠️" if error > 50 else "  "
        print(f"{marker}  {name:<28} {target:>8.2f} {actual:>8.2f} {error_str:>8} {count_str:>8}")

        if count is not None:
            total_count += count
        if atoms is not None:
            total_atoms += atoms
        
        if error > 50:
            warnings.append((r.get('name', 'N/A'), error))

    print("\n" + "-" * 90)
    print(f"总计: {total_count} 分子, {total_atoms} 原子 (目标: {target_atoms})")
    
    if warnings:
        print("\n⚠️ 警告: 以下组分相对误差 > 50%:")
        for name, error in warnings:
            name_short = name[:40] + "..." if len(name) > 40 else name
            print(f"   - {name_short}: {error:.1f}%")
        print("\n   建议: 增大 --target_atoms 以减小凑整误差")
    
    print("=" * 90)


def main():
    parser = argparse.ArgumentParser(
        description="配方 wt% 转换为分子/原子数量（含误差分析）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python3 recipe_to_counts.py --target_atoms 200
  python3 recipe_to_counts.py --target_atoms 500 --output counts_500

输出:
  - counts.csv: 详细数据
  - counts.json: JSON 格式
  - counts_report.txt: 误差报告
        """
    )

    parser.add_argument("--recipe", default="recipe.yaml",
                        help="配方文件路径 (默认: recipe.yaml)")
    parser.add_argument("--output", default="counts",
                        help="输出文件前缀 (默认: counts)")
    parser.add_argument("--target_atoms", type=int, required=True,
                        help="目标总原子数")

    args = parser.parse_args()

    print("=" * 90)
    print("recipe_to_counts.py - 配方换算工具 (v2.0)")
    print("=" * 90)
    print(f"配方文件: {args.recipe}")
    print(f"目标原子数: {args.target_atoms}")

    # 加载配方
    data = load_yaml(args.recipe)
    entries = flatten_recipe(data)

    if not entries:
        print("[ERROR] 配方为空")
        sys.exit(1)

    print(f"\n>>> 读取 {len(entries)} 个条目")

    # 执行换算
    results = compute_by_target_atoms(entries, args.target_atoms)

    # 输出文件
    csv_path = f"{args.output}.csv"
    json_path = f"{args.output}.json"
    report_path = f"{args.output}_report.txt"

    print(f"\n>>> 输出文件:")
    write_csv(results, csv_path)
    print(f"    [OK] {csv_path}")

    write_json(results, json_path)
    print(f"    [OK] {json_path}")

    write_report(results, report_path, args.target_atoms)
    print(f"    [OK] {report_path}")

    # 打印摘要
    print_summary(results, args.target_atoms)

    # 实际原子数
    actual_atoms = sum(r.get('scaled_atoms', 0) or 0 for r in results)
    diff = abs(actual_atoms - args.target_atoms)
    print(f"\n[INFO] 实际总原子数: {actual_atoms}")
    print(f"[INFO] 与目标差异: {diff} 原子 ({diff/args.target_atoms*100:.1f}%)")


if __name__ == "__main__":
    main()
