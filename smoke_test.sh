#!/usr/bin/env bash
# ============================================================================
# smoke_test.sh - VASP Scripts 功能验证脚本 v2.2
# ============================================================================
# 用法: ./smoke_test.sh
#
# 验证内容:
#   1. Python 脚本语法检查（含 utils/）
#   2. --help 输出检查
#   3. recipe_validate.py 验证
#   4. make_incar_aimd.py 生成（dry-run）
#   5. recipe_to_counts.py 换算
#   6. Shell 脚本语法检查
#   7. MSD MTO 算法验证（合成数据）
#
# 注意: 不会真正运行 VASP
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "============================================"
echo "VASP Scripts Smoke Test"
echo "============================================"
echo "目录: $SCRIPT_DIR"
echo "时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================"
echo ""

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

pass_count=0
fail_count=0

pass() {
    echo -e "${GREEN}[PASS]${NC} $1"
    ((pass_count++)) || true
}

fail() {
    echo -e "${RED}[FAIL]${NC} $1"
    ((fail_count++)) || true
}

warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

# ============== 1. Python 语法检查 ==============
echo ""
echo ">>> 1. Python 语法检查..."

py_files=(
    "aimd_msd.py"
    "aimd_post.py"
    "make_incar_aimd.py"
    "recipe_validate.py"
    "recipe_to_counts.py"
    "setup_aimd_ase.py"
    "setup_electronic.py"
    "analyze_electronic.py"
    "utils/__init__.py"
    "utils/connectivity.py"
    "utils/charges.py"
    "utils/units.py"
)

for pyfile in "${py_files[@]}"; do
    if [[ -f "$pyfile" ]]; then
        if python3 -m py_compile "$pyfile" 2>/dev/null; then
            pass "$pyfile 语法正确"
        else
            fail "$pyfile 语法错误"
        fi
    else
        warn "$pyfile 不存在"
    fi
done

# ============== 2. --help 检查 ==============
echo ""
echo ">>> 2. --help 输出检查..."

help_scripts=(
    "aimd_msd.py"
    "make_incar_aimd.py"
    "recipe_validate.py"
    "recipe_to_counts.py"
    "setup_aimd_ase.py"
    "setup_electronic.py"
)

for script in "${help_scripts[@]}"; do
    if [[ -f "$script" ]]; then
        if python3 "$script" --help >/dev/null 2>&1; then
            pass "$script --help"
        else
            fail "$script --help"
        fi
    fi
done

# ============== 3. recipe_validate.py ==============
echo ""
echo ">>> 3. recipe_validate.py 测试..."

if [[ -f "examples/minimal_recipe.yaml" ]]; then
    if python3 recipe_validate.py --recipe examples/minimal_recipe.yaml >/dev/null 2>&1; then
        pass "recipe_validate.py 验证 minimal_recipe.yaml"
    else
        fail "recipe_validate.py 验证失败"
    fi
else
    warn "examples/minimal_recipe.yaml 不存在"
fi

# ============== 4. make_incar_aimd.py ==============
echo ""
echo ">>> 4. make_incar_aimd.py 测试..."

TEST_DIR=$(mktemp -d)
trap "rm -rf $TEST_DIR" EXIT

cp examples/minimal_recipe.yaml "$TEST_DIR/recipe.yaml" 2>/dev/null || true

if [[ -f "$TEST_DIR/recipe.yaml" ]]; then
    cd "$TEST_DIR"
    if python3 "$SCRIPT_DIR/make_incar_aimd.py" --recipe recipe.yaml --out INCAR.test >/dev/null 2>&1; then
        if [[ -f "INCAR.test" ]]; then
            pass "make_incar_aimd.py 生成 INCAR"
        else
            fail "INCAR 未生成"
        fi
    else
        fail "make_incar_aimd.py 执行失败"
    fi
    cd "$SCRIPT_DIR"
fi

# ============== 5. recipe_to_counts.py ==============
echo ""
echo ">>> 5. recipe_to_counts.py 测试..."

if [[ -f "examples/minimal_recipe.yaml" ]]; then
    cd "$TEST_DIR"
    if python3 "$SCRIPT_DIR/recipe_to_counts.py" --recipe "$SCRIPT_DIR/examples/minimal_recipe.yaml" --target_atoms 100 >/dev/null 2>&1; then
        if [[ -f "counts.csv" && -f "counts.json" ]]; then
            pass "recipe_to_counts.py 生成 counts 文件"
        else
            fail "counts 文件未生成"
        fi
    else
        fail "recipe_to_counts.py 执行失败"
    fi
    cd "$SCRIPT_DIR"
fi

# ============== 6. Shell 脚本语法检查 ==============
echo ""
echo ">>> 6. Shell 脚本语法检查..."

sh_files=(
    "vasp_env.sh"
    "run_vasp.sh"
    "check_vasp.sh"
    "aimd_watch.sh"
    "aimd_setup.sh"
    "clean_vasp.sh"
)

for shfile in "${sh_files[@]}"; do
    if [[ -f "$shfile" ]]; then
        if bash -n "$shfile" 2>/dev/null; then
            pass "$shfile 语法正确"
        else
            fail "$shfile 语法错误"
        fi
    else
        warn "$shfile 不存在"
    fi
done

# ============== 7. MSD MTO 算法验证 ==============
echo ""
echo ">>> 7. MSD MTO 算法验证（合成数据）..."

if [[ -f "tests/test_msd_synthetic.py" ]]; then
    if python3 tests/test_msd_synthetic.py 2>/dev/null | tail -5; then
        # 检查退出码
        if python3 tests/test_msd_synthetic.py >/dev/null 2>&1; then
            pass "MSD MTO 算法验证"
        else
            fail "MSD MTO 算法验证失败"
        fi
    else
        fail "MSD MTO 测试执行失败"
    fi
else
    warn "tests/test_msd_synthetic.py 不存在"
fi

# ============== 结果汇总 ==============
echo ""
echo "============================================"
echo "测试结果汇总"
echo "============================================"
echo -e "通过: ${GREEN}$pass_count${NC}"
echo -e "失败: ${RED}$fail_count${NC}"
echo "============================================"

if [[ $fail_count -gt 0 ]]; then
    echo -e "${RED}存在失败的测试，请检查！${NC}"
    exit 1
else
    echo -e "${GREEN}所有测试通过！${NC}"
    exit 0
fi

