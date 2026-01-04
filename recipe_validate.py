#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
recipe_validate.py - 凝胶电解质配方验证工具

功能：
  - 读取 recipe.yaml 配方文件
  - 校验所有 wt_pct 非负
  - 校验 8 类总和 = 100%（容差 1e-3）
  - 校验 simulation 段（温度、时间步、步数等）
  - 按固定顺序打印标准化摘要

用法：
  python3 recipe_validate.py [--recipe recipe.yaml]

作者：STAR0418-ABC
"""

import argparse
import sys
import os
from typing import Dict, List, Any, Tuple

# 尝试导入 yaml
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    print("[WARN] PyYAML 未安装")
    print("[INFO] 请运行: pip install pyyaml")


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

# 容差
TOLERANCE = 1e-3

# 温度范围（摄氏度）
TEMP_MIN = -50
TEMP_MAX = 300


def load_yaml(filepath: str) -> Dict[str, Any]:
    """加载 YAML 文件"""
    if not os.path.isfile(filepath):
        print(f"[ERROR] 配方文件不存在: {filepath}")
        sys.exit(1)

    if not HAS_YAML:
        print("[ERROR] 需要 PyYAML 库来解析 YAML 文件")
        print("[INFO] 请运行: pip install pyyaml")
        sys.exit(1)

    with open(filepath, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)

    if data is None:
        data = {}

    return data


def validate_entry(entry: Dict, category: str, idx: int) -> Tuple[bool, List[str]]:
    """
    验证单个组分条目

    返回: (是否有效, 错误列表)
    """
    errors = []

    # 必需字段检查
    if 'name' not in entry:
        errors.append(f"  [{category}][{idx}] 缺少 'name' 字段")

    if 'wt_pct' not in entry:
        errors.append(f"  [{category}][{idx}] 缺少 'wt_pct' 字段")
    else:
        wt = entry['wt_pct']
        if not isinstance(wt, (int, float)):
            errors.append(f"  [{category}][{idx}] 'wt_pct' 必须是数值，当前: {type(wt)}")
        elif wt < 0:
            errors.append(f"  [{category}][{idx}] 'wt_pct' 不能为负: {wt}")

    if 'kind' not in entry:
        errors.append(f"  [{category}][{idx}] 缺少 'kind' 字段")

    # 检查 name 是否包含中文
    name = entry.get('name', '')
    has_chinese = any('\u4e00' <= c <= '\u9fff' for c in name)
    if not has_chinese:
        errors.append(f"  [{category}][{idx}] 'name' 应包含中文全称: {name}")

    return len(errors) == 0, errors


def validate_simulation(sim: Dict) -> Tuple[bool, List[str], List[str]]:
    """
    验证 simulation 段

    返回: (是否有效, 错误列表, 警告列表)
    """
    errors = []
    warnings = []

    if sim is None:
        return True, [], ["[INFO] 未定义 simulation 段，跳过模拟条件验证"]

    mode = sim.get('mode', 'static')

    if mode == 'aimd':
        # 温度检查
        if 'temperature_C' not in sim:
            errors.append("[simulation] AIMD 模式必须指定 temperature_C")
        else:
            temp_c = sim['temperature_C']
            if not isinstance(temp_c, (int, float)):
                errors.append(f"[simulation] temperature_C 必须是数值，当前: {type(temp_c)}")
            elif temp_c < TEMP_MIN or temp_c > TEMP_MAX:
                warnings.append(f"[simulation] temperature_C={temp_c}°C 超出建议范围 [{TEMP_MIN}, {TEMP_MAX}]")

        # 时间步检查
        dt = sim.get('dt_fs', None)
        if dt is None:
            errors.append("[simulation] AIMD 模式必须指定 dt_fs (时间步长)")
        elif not isinstance(dt, (int, float)) or dt <= 0:
            errors.append(f"[simulation] dt_fs 必须是正数，当前: {dt}")

        # 步数检查
        nsteps = sim.get('nsteps', None)
        if nsteps is None:
            errors.append("[simulation] AIMD 模式必须指定 nsteps (总步数)")
        elif not isinstance(nsteps, int) or nsteps <= 0:
            errors.append(f"[simulation] nsteps 必须是正整数，当前: {nsteps}")

        # 系综检查
        ensemble = sim.get('ensemble', 'nvt')
        if ensemble not in ['nvt', 'nve']:
            warnings.append(f"[simulation] ensemble='{ensemble}' 非标准值，建议使用 nvt 或 nve")

        # 恒温器检查
        thermostat = sim.get('thermostat', 'langevin')
        if thermostat not in ['langevin', 'nose_hoover']:
            warnings.append(f"[simulation] thermostat='{thermostat}' 非标准值，建议使用 langevin 或 nose_hoover")

        # gamma 检查
        gamma = sim.get('gamma_1ps', 10.0)
        if not isinstance(gamma, (int, float)) or gamma <= 0:
            warnings.append(f"[simulation] gamma_1ps 应为正数，当前: {gamma}")

        # nelm 检查
        nelm = sim.get('nelm', 100)
        if not isinstance(nelm, int) or nelm <= 0:
            warnings.append(f"[simulation] nelm 应为正整数，当前: {nelm}")

        # ediff 检查
        ediff = sim.get('ediff', 1e-5)
        if not isinstance(ediff, (int, float)) or ediff <= 0:
            warnings.append(f"[simulation] ediff 应为正数，当前: {ediff}")

    return len(errors) == 0, errors, warnings


def validate_recipe(data: Dict) -> Tuple[bool, float, List[str]]:
    """
    验证整个配方的组分部分

    返回: (是否有效, 总 wt_pct, 错误列表)
    """
    all_errors = []
    total_wt = 0.0

    for cat_key, cat_name in CATEGORY_ORDER:
        entries = data.get(cat_key, [])

        # 处理 None 或非列表
        if entries is None:
            entries = []
        if not isinstance(entries, list):
            all_errors.append(f"[{cat_key}] 应为列表类型，当前: {type(entries)}")
            continue

        for idx, entry in enumerate(entries):
            if not isinstance(entry, dict):
                all_errors.append(f"[{cat_key}][{idx}] 条目应为字典类型")
                continue

            valid, errors = validate_entry(entry, cat_key, idx)
            all_errors.extend(errors)

            wt = entry.get('wt_pct', 0)
            if isinstance(wt, (int, float)) and wt >= 0:
                total_wt += wt

    return len(all_errors) == 0, total_wt, all_errors


def print_summary(data: Dict):
    """按固定顺序打印配方摘要"""
    print("\n" + "=" * 78)
    print("配方摘要 (Recipe Summary)")
    print("=" * 78)

    total_wt = 0.0
    total_entries = 0

    for cat_key, cat_name in CATEGORY_ORDER:
        entries = data.get(cat_key, []) or []

        print(f"\n>>> {cat_key} ({cat_name})")
        print("-" * 78)

        if not entries:
            print("    (空)")
            continue

        print(f"    {'名称':<40} {'wt%':>8} {'种类':<12} {'MW':>10}")
        print("    " + "-" * 74)

        cat_wt = 0.0
        for entry in entries:
            name = entry.get('name', 'N/A')
            wt = entry.get('wt_pct', 0)
            kind = entry.get('kind', 'N/A')
            mw = entry.get('mw_g_mol', None)

            # 截断过长的名称
            if len(name) > 38:
                name_display = name[:35] + "..."
            else:
                name_display = name

            mw_str = f"{mw:.2f}" if mw else "N/A"

            print(f"    {name_display:<40} {wt:>8.2f} {kind:<12} {mw_str:>10}")
            cat_wt += wt
            total_entries += 1

        print(f"    {'小计:':<40} {cat_wt:>8.2f}")
        total_wt += cat_wt

    print("\n" + "=" * 78)
    print(f"总条目数: {total_entries}")
    print(f"总 wt%: {total_wt:.4f}")

    # 检查总和
    diff = abs(total_wt - 100.0)
    if diff <= TOLERANCE:
        print(f"[OK] 总和校验通过 (|{total_wt:.4f} - 100| = {diff:.6f} <= {TOLERANCE})")
    else:
        print(f"[ERROR] 总和校验失败: {total_wt:.4f} != 100 (差值: {diff:.4f})")

    print("=" * 78)


def print_simulation_summary(sim: Dict):
    """打印模拟条件摘要"""
    if sim is None:
        return

    print("\n" + "=" * 78)
    print("模拟条件摘要 (Simulation Settings)")
    print("=" * 78)

    mode = sim.get('mode', 'static')
    print(f"模式: {mode}")

    if mode == 'aimd':
        temp_c = sim.get('temperature_C', 'N/A')
        if isinstance(temp_c, (int, float)):
            temp_k = temp_c + 273.15
            print(f"温度: {temp_c} °C = {temp_k:.2f} K")
        else:
            print(f"温度: {temp_c} °C")

        dt = sim.get('dt_fs', 'N/A')
        print(f"时间步长: {dt} fs (POTIM)")

        nsteps = sim.get('nsteps', 'N/A')
        print(f"总步数: {nsteps} (NSW)")

        if isinstance(dt, (int, float)) and isinstance(nsteps, int):
            total_time_ps = dt * nsteps / 1000.0
            print(f"总模拟时间: {total_time_ps:.2f} ps")

        ensemble = sim.get('ensemble', 'nvt')
        print(f"系综: {ensemble.upper()}")

        thermostat = sim.get('thermostat', 'langevin')
        print(f"恒温器: {thermostat}")

        if thermostat == 'langevin':
            gamma = sim.get('gamma_1ps', 10.0)
            print(f"摩擦系数: {gamma} 1/ps (LANGEVIN_GAMMA)")
        elif thermostat == 'nose_hoover':
            smass = sim.get('smass', -3)
            print(f"质量参数: SMASS = {smass}")

        nelm = sim.get('nelm', 100)
        ediff = sim.get('ediff', 1e-5)
        print(f"电子步: NELM={nelm}, EDIFF={ediff:.0e}")

        encut = sim.get('encut', None)
        if encut:
            print(f"截断能: ENCUT = {encut} eV")
        else:
            print("截断能: 未指定 (需在 INCAR.base 或手动设置)")

    print("=" * 78)


def main():
    parser = argparse.ArgumentParser(
        description="凝胶电解质配方验证工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
    python3 recipe_validate.py
    python3 recipe_validate.py --recipe my_recipe.yaml
        """
    )
    parser.add_argument("--recipe", default="recipe.yaml",
                        help="配方文件路径 (默认: recipe.yaml)")

    args = parser.parse_args()

    print("=" * 78)
    print("recipe_validate.py - 配方验证工具")
    print("=" * 78)
    print(f"配方文件: {args.recipe}")

    # 加载配方
    data = load_yaml(args.recipe)

    # 验证组分
    print("\n>>> 验证组分...")
    valid_recipe, total_wt, recipe_errors = validate_recipe(data)

    if recipe_errors:
        print("\n[ERROR] 组分验证发现问题:")
        for err in recipe_errors:
            print(err)

    # 验证模拟条件
    print("\n>>> 验证模拟条件...")
    sim = data.get('simulation', None)
    valid_sim, sim_errors, sim_warnings = validate_simulation(sim)

    if sim_errors:
        print("\n[ERROR] 模拟条件验证发现问题:")
        for err in sim_errors:
            print(err)

    if sim_warnings:
        print("\n[WARN] 模拟条件警告:")
        for warn in sim_warnings:
            print(warn)

    # 打印摘要
    print_summary(data)
    print_simulation_summary(sim)

    # 返回状态
    has_error = False

    if not valid_recipe:
        print("\n[FAIL] 组分验证失败，请修正上述错误。")
        has_error = True

    diff = abs(total_wt - 100.0)
    if diff > TOLERANCE:
        print(f"\n[FAIL] 总和不等于 100%: {total_wt:.4f}%")
        has_error = True

    if not valid_sim:
        print("\n[FAIL] 模拟条件验证失败，请修正上述错误。")
        has_error = True

    if has_error:
        sys.exit(1)

    print("\n[PASS] 配方验证通过！")
    sys.exit(0)


if __name__ == "__main__":
    main()
