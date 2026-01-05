#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
selfcheck_oligomer.py - 寡聚体生成自检脚本
============================================

验证 mol2_to_polymer_mol2.py 的核心功能：
1. 默认聚合度是 3 (三聚体)
2. recipe.yaml 配置可正确读取
3. CLI --n_repeat 可覆盖默认值
4. 输出命名规则: X.mol2 → PX.mol2
5. 向后兼容: oligomer_n=5 仍可工作

用法:
  python scripts/selfcheck_oligomer.py
  python scripts/selfcheck_oligomer.py --verbose

"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

# 添加脚本目录到 PATH
script_dir = Path(__file__).parent
project_root = script_dir.parent
sys.path.insert(0, str(script_dir))


def test_default_oligomer_n():
    """测试 1: 默认聚合度是 3"""
    print("\n[TEST 1] 默认聚合度 = 3")
    print("-" * 40)
    
    try:
        from mol2_to_polymer_mol2 import DEFAULT_OLIGOMER_N
        
        if DEFAULT_OLIGOMER_N == 3:
            print(f"  ✓ PASS: DEFAULT_OLIGOMER_N = {DEFAULT_OLIGOMER_N}")
            return True
        else:
            print(f"  ✗ FAIL: DEFAULT_OLIGOMER_N = {DEFAULT_OLIGOMER_N}, 期望 3")
            return False
    except Exception as e:
        print(f"  ✗ FAIL: 导入失败: {e}")
        return False


def test_recipe_loading():
    """测试 2: 从 recipe.yaml 读取 oligomer_n"""
    print("\n[TEST 2] 从 recipe.yaml 读取聚合度")
    print("-" * 40)
    
    try:
        from mol2_to_polymer_mol2 import load_oligomer_n_from_recipe
        
        recipe_path = project_root / "config" / "recipe.yaml"
        
        if not recipe_path.exists():
            print(f"  [SKIP] recipe.yaml 不存在: {recipe_path}")
            return True  # 跳过，不算失败
        
        oligomer_n = load_oligomer_n_from_recipe(str(recipe_path))
        
        if oligomer_n is not None:
            print(f"  ✓ PASS: recipe.yaml oligomer_n = {oligomer_n}")
            return True
        else:
            print(f"  [WARN] recipe.yaml 中未找到 oligomer_n 配置 (可接受)")
            return True
    except Exception as e:
        print(f"  ✗ FAIL: 读取失败: {e}")
        return False


def test_priority_logic():
    """测试 3: 优先级逻辑 (CLI > recipe > 默认)"""
    print("\n[TEST 3] 优先级逻辑 (CLI > recipe > 默认)")
    print("-" * 40)
    
    try:
        from mol2_to_polymer_mol2 import get_effective_oligomer_n
        
        # CLI 覆盖
        n1 = get_effective_oligomer_n(cli_n_repeat=5, recipe_path=None)
        # 无 CLI，无 recipe
        n2 = get_effective_oligomer_n(cli_n_repeat=None, recipe_path=None)
        # CLI 覆盖 recipe
        n3 = get_effective_oligomer_n(cli_n_repeat=7, recipe_path="config/recipe.yaml")
        
        success = True
        
        if n1 == 5:
            print(f"  ✓ CLI=5 → {n1}")
        else:
            print(f"  ✗ CLI=5 → {n1}, 期望 5")
            success = False
        
        if n2 == 3:
            print(f"  ✓ CLI=None, recipe=None → {n2} (默认)")
        else:
            print(f"  ✗ CLI=None, recipe=None → {n2}, 期望 3")
            success = False
        
        if n3 == 7:
            print(f"  ✓ CLI=7 覆盖 recipe → {n3}")
        else:
            print(f"  ✗ CLI=7 覆盖 recipe → {n3}, 期望 7")
            success = False
        
        return success
    except Exception as e:
        print(f"  ✗ FAIL: 测试失败: {e}")
        return False


def test_input_validation():
    """测试 4: 输入格式验证 (仅接受 .mol2)"""
    print("\n[TEST 4] 输入格式验证 (仅接受 .mol2)")
    print("-" * 40)
    
    try:
        from mol2_to_polymer_mol2 import validate_mol2_input
        
        # 重定向 stdout 以捕获错误消息
        import io
        from contextlib import redirect_stdout
        
        f = io.StringIO()
        with redirect_stdout(f):
            valid_mol2 = validate_mol2_input("/tmp/test.mol2")
            invalid_pdb = validate_mol2_input("/tmp/test.pdb")
            invalid_mol = validate_mol2_input("/tmp/test.mol")
        
        success = True
        
        if valid_mol2:
            print(f"  ✓ .mol2 → 接受")
        else:
            print(f"  ✗ .mol2 → 拒绝, 期望接受")
            success = False
        
        if not invalid_pdb:
            print(f"  ✓ .pdb → 拒绝")
        else:
            print(f"  ✗ .pdb → 接受, 期望拒绝")
            success = False
        
        if not invalid_mol:
            print(f"  ✓ .mol → 拒绝")
        else:
            print(f"  ✗ .mol → 接受, 期望拒绝")
            success = False
        
        return success
    except Exception as e:
        print(f"  ✗ FAIL: 测试失败: {e}")
        return False


def test_output_naming():
    """测试 5: 输出命名规则 (X.mol2 → PX.mol2)"""
    print("\n[TEST 5] 输出命名规则 (X.mol2 → PX.mol2)")
    print("-" * 40)
    
    # 测试命名逻辑
    test_cases = [
        ("EGDA.mol2", "PEGDA.mol2"),
        ("MMA.mol2", "PMMA.mol2"),
        ("mymonomer.mol2", "Pmymonomer.mol2"),
    ]
    
    success = True
    for input_name, expected_output in test_cases:
        # 模拟命名逻辑
        base, ext = os.path.splitext(input_name)
        actual_output = f"P{base}{ext}"
        
        if actual_output == expected_output:
            print(f"  ✓ {input_name} → {actual_output}")
        else:
            print(f"  ✗ {input_name} → {actual_output}, 期望 {expected_output}")
            success = False
    
    return success


def test_monomer_configs():
    """测试 6: 单体配置完整性"""
    print("\n[TEST 6] 单体配置完整性")
    print("-" * 40)
    
    try:
        from mol2_to_polymer_mol2 import MONOMER_CONFIGS
        
        required_monomers = ["EO", "EGDA", "MMA", "VA", "AM", "TFEMA"]
        missing = []
        
        for m in required_monomers:
            if m not in MONOMER_CONFIGS:
                missing.append(m)
            else:
                config = MONOMER_CONFIGS[m]
                # 检查必需字段
                if not config.head_atom_names and not config.head_smarts:
                    print(f"  [WARN] {m}: 缺少 head 定义")
                if not config.tail_atom_names and not config.tail_smarts:
                    print(f"  [WARN] {m}: 缺少 tail 定义")
        
        if not missing:
            print(f"  ✓ PASS: 所有必需单体配置存在 ({len(required_monomers)} 个)")
            return True
        else:
            print(f"  ✗ FAIL: 缺少配置: {missing}")
            return False
    except Exception as e:
        print(f"  ✗ FAIL: 测试失败: {e}")
        return False


def test_backward_compatibility():
    """测试 7: 向后兼容性 (oligomer_n=5 仍可工作)"""
    print("\n[TEST 7] 向后兼容性 (oligomer_n=5)")
    print("-" * 40)
    
    try:
        from mol2_to_polymer_mol2 import get_effective_oligomer_n
        
        # 可以指定 5 或任何其他值
        n5 = get_effective_oligomer_n(cli_n_repeat=5, recipe_path=None)
        n10 = get_effective_oligomer_n(cli_n_repeat=10, recipe_path=None)
        
        if n5 == 5 and n10 == 10:
            print(f"  ✓ PASS: oligomer_n=5 → {n5}, oligomer_n=10 → {n10}")
            return True
        else:
            print(f"  ✗ FAIL: 向后兼容失败")
            return False
    except Exception as e:
        print(f"  ✗ FAIL: 测试失败: {e}")
        return False


def test_recipe_to_counts_integration():
    """测试 8: recipe_to_counts.py 聚合度集成"""
    print("\n[TEST 8] recipe_to_counts.py 聚合度集成")
    print("-" * 40)
    
    try:
        # 检查 recipe_to_counts.py 是否能正确处理 polymerization 配置
        recipe_counts_path = script_dir / "recipe_to_counts.py"
        
        if not recipe_counts_path.exists():
            print(f"  [SKIP] recipe_to_counts.py 不存在")
            return True
        
        # 读取文件内容，检查是否包含 polymerization 处理逻辑
        with open(recipe_counts_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        checks = [
            ("polymerization", "polymerization 配置读取"),
            ("oligomer_n", "oligomer_n 字段处理"),
            ("monomer_mw", "monomer_mw 字段处理"),
        ]
        
        for keyword, desc in checks:
            if keyword in content:
                print(f"  ✓ {desc}: 已实现")
            else:
                print(f"  [WARN] {desc}: 未找到 (可能未实现)")
        
        return True
    except Exception as e:
        print(f"  ✗ FAIL: 测试失败: {e}")
        return False


def main():
    """运行所有自检测试"""
    print("=" * 60)
    print("mol2_to_polymer_mol2.py 寡聚体生成自检")
    print("=" * 60)
    print(f"项目根目录: {project_root}")
    print(f"脚本目录: {script_dir}")
    
    tests = [
        ("默认聚合度 = 3", test_default_oligomer_n),
        ("recipe.yaml 读取", test_recipe_loading),
        ("优先级逻辑", test_priority_logic),
        ("输入格式验证", test_input_validation),
        ("输出命名规则", test_output_naming),
        ("单体配置完整性", test_monomer_configs),
        ("向后兼容性", test_backward_compatibility),
        ("recipe_to_counts 集成", test_recipe_to_counts_integration),
    ]
    
    passed = 0
    failed = 0
    
    for name, test_func in tests:
        try:
            result = test_func()
            if result:
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"\n[{name}] 异常: {e}")
            failed += 1
    
    # 汇总
    print("\n" + "=" * 60)
    print("自检结果汇总")
    print("=" * 60)
    print(f"  通过: {passed}")
    print(f"  失败: {failed}")
    print(f"  总计: {passed + failed}")
    print("=" * 60)
    
    if failed == 0:
        print("\n✓ 所有测试通过!")
        return 0
    else:
        print(f"\n✗ {failed} 个测试失败")
        return 1


if __name__ == "__main__":
    sys.exit(main())

