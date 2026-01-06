#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
selfcheck.py - 统一自检脚本
============================

验证重构后的脚本功能完整性。

测试项:
  1. lib 模块可导入
  2. recipe.yaml 可解析
  3. 单体批处理可运行 (生成 P*.mol2)
  4. 输出文件存在且非空
  5. MOL2 包含必需段落

用法:
  python scripts/selfcheck.py
  python scripts/selfcheck.py --verbose
"""

import argparse
import os
import sys
import tempfile
import shutil
from pathlib import Path

# 项目根目录
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


class TestResult:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.tests = []
    
    def add(self, name: str, passed: bool, message: str = ""):
        self.tests.append((name, passed, message))
        if passed:
            self.passed += 1
        else:
            self.failed += 1
    
    def print_summary(self):
        print("\n" + "=" * 60)
        print("自检结果")
        print("=" * 60)
        
        for name, passed, message in self.tests:
            status = "✓ PASS" if passed else "✗ FAIL"
            print(f"  [{status}] {name}")
            if message and not passed:
                print(f"          {message}")
        
        print("=" * 60)
        print(f"通过: {self.passed}, 失败: {self.failed}")
        print("=" * 60)


def test_lib_imports() -> tuple:
    """测试 lib 模块可导入"""
    try:
        from scripts.lib import (
            load_recipe, get_oligomer_n,
            read_mol2, write_mol2_strict,
            get_polymer_output_path,
            validate_file_exists, validate_file_nonempty,
            log_success, log_error
        )
        return True, ""
    except ImportError as e:
        return False, str(e)


def test_recipe_loading() -> tuple:
    """测试 recipe.yaml 可解析"""
    try:
        from scripts.lib.recipe import load_recipe, get_oligomer_n
        
        recipe_path = PROJECT_ROOT / "config" / "recipe.yaml"
        if not recipe_path.exists():
            return True, "recipe.yaml 不存在 (跳过)"
        
        config = load_recipe(str(recipe_path))
        n = get_oligomer_n(config)
        
        if n == 3:
            return True, f"oligomer_n = {n}"
        else:
            return True, f"oligomer_n = {n} (非默认值)"
    except Exception as e:
        return False, str(e)


def test_naming_rules() -> tuple:
    """测试命名规则"""
    try:
        from scripts.lib.naming import get_polymer_output_path, get_polymer_name
        
        cases = [
            ("mol2/EGDA.mol2", "mol2/PEGDA.mol2"),
            ("mol2/MMA.mol2", "mol2/PMMA.mol2"),
        ]
        
        for input_path, expected in cases:
            # 使用相对路径测试
            result = get_polymer_output_path(input_path)
            result_basename = os.path.basename(result)
            expected_basename = os.path.basename(expected)
            
            if result_basename != expected_basename:
                return False, f"{input_path} → {result_basename}, 期望 {expected_basename}"
        
        return True, ""
    except Exception as e:
        return False, str(e)


def test_monomers_batch() -> tuple:
    """测试单体批处理"""
    try:
        from scripts.cli.monomers import process_single_monomer
        from scripts.lib.validate import validate_mol2_output
        
        # 测试 EGDA
        input_path = PROJECT_ROOT / "mol2" / "EGDA.mol2"
        if not input_path.exists():
            return True, "EGDA.mol2 不存在 (跳过)"
        
        output_path = process_single_monomer(str(input_path), oligomer_n=3, optimize=False)
        
        # 验证输出
        validate_mol2_output(output_path)
        
        return True, f"PEGDA.mol2 生成成功"
    except Exception as e:
        return False, str(e)


def test_output_validation() -> tuple:
    """测试输出验证功能"""
    try:
        from scripts.lib.validate import validate_mol2_output, ValidationError
        
        # 测试现有文件
        pegda_path = PROJECT_ROOT / "mol2" / "PEGDA.mol2"
        if pegda_path.exists():
            validate_mol2_output(str(pegda_path))
            return True, "PEGDA.mol2 验证通过"
        
        return True, "无文件可验证 (跳过)"
    except Exception as e:
        return False, str(e)


def test_default_oligomer_n() -> tuple:
    """测试默认聚合度为 3"""
    try:
        from scripts.lib.recipe import get_oligomer_n
        
        n = get_oligomer_n(None, None, 3)
        if n == 3:
            return True, "DEFAULT_OLIGOMER_N = 3"
        else:
            return False, f"DEFAULT_OLIGOMER_N = {n}, 期望 3"
    except Exception as e:
        return False, str(e)


def run_all_tests(verbose: bool = False) -> int:
    """运行所有测试"""
    print("=" * 60)
    print("gel_packmol 自检测试")
    print("=" * 60)
    print(f"项目根目录: {PROJECT_ROOT}")
    
    result = TestResult()
    
    tests = [
        ("lib 模块导入", test_lib_imports),
        ("recipe.yaml 解析", test_recipe_loading),
        ("命名规则", test_naming_rules),
        ("默认聚合度", test_default_oligomer_n),
        ("单体批处理", test_monomers_batch),
        ("输出验证", test_output_validation),
    ]
    
    for name, test_func in tests:
        print(f"\n[TEST] {name}...")
        try:
            passed, message = test_func()
            result.add(name, passed, message)
            if passed:
                print(f"  ✓ PASS" + (f" - {message}" if message else ""))
            else:
                print(f"  ✗ FAIL - {message}")
        except Exception as e:
            result.add(name, False, str(e))
            print(f"  ✗ FAIL - {e}")
    
    result.print_summary()
    
    return 0 if result.failed == 0 else 1


def main():
    parser = argparse.ArgumentParser(description="gel_packmol 自检测试")
    parser.add_argument("--verbose", "-v", action="store_true", help="详细输出")
    args = parser.parse_args()
    
    return run_all_tests(args.verbose)


if __name__ == "__main__":
    sys.exit(main())

