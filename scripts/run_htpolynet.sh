#!/bin/bash

# Author: Star



# ==============================================================================
# run_htpolynet.sh - HTPolyNet 交联聚合物模拟驱动脚本 (v2.2)
# ==============================================================================
# 功能：运行 HTPolyNet CLI 进行交联反应模拟
#
# 特点：
#   - 完全不依赖 conda init / conda activate
#   - 直接使用系统 PATH 中的 htpolynet 和 obabel
#   - 官方示例 run.sh 中的 conda activate 会被自动替换
#
# 用法:
#   # 官方示例模式（推荐先用这个验证环境）
#   ./scripts/run_htpolynet.sh --example 0 -o outputs/smoke_full/htpolynet
#
#   # 自定义配置模式
#   ./scripts/run_htpolynet.sh -i gel.pdb -c config.yaml -o outputs/htpolynet
#
# 依赖: htpolynet, obabel, gmx
# ==============================================================================

set -euo pipefail

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'

# 获取脚本目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# 默认参数
EXAMPLE_NUM=""
INPUT_PDB=""
CONFIG_YAML=""
RECIPE_YAML=""
OUTPUT_DIR=""
SKIP_PREFLIGHT=false

# 打印函数
info()    { echo -e "${GREEN}[INFO]${NC} $1"; }
warn()    { echo -e "${YELLOW}[WARN]${NC} $1"; }
error()   { echo -e "${RED}[ERROR]${NC} $1" >&2; }
step()    { echo -e "${CYAN}[STEP]${NC} $1"; }
success() { echo -e "${MAGENTA}[SUCCESS]${NC} $1"; }

# 使用帮助
usage() {
    cat << EOF
用法: $0 [选项]

===== 两种互斥模式 =====

【模式 1: 官方示例模式（推荐先跑这个验证环境）】
  $0 --example 0 -o outputs/smoke_full/htpolynet
  
  示例列表:
    0: 0-liquid-styrene (最快，推荐)
    1: 1-polystyrene
    2: 2-polymethylstyrene
    3: 3-bisgma-styrene-thermoset
    4: 4-pacm-dgeba-epoxy-thermoset
    5: 5-dfda-fde-epoxy-thermoset

【模式 2: 自定义配置模式】
  $0 -i outputs/smoke/packmol/gel.pdb -c configs/htpolynet.yaml -o outputs/smoke/htpolynet
  
  （可选）添加 --recipe 进行 preflight 原子命名检查:
  $0 -i gel.pdb -c htpolynet.yaml --recipe configs/recipe_smoke.yaml -o outputs/htpolynet

===== 通用选项 =====
  --example <n>       运行 HTPolyNet 官方示例 (n=0-5)
  -i, --input <pdb>   输入 PDB 结构文件
  -c, --config <yaml> HTPolyNet 配置文件 (cfg.yaml)
  --recipe <yaml>     Packmol recipe.yaml 用于 preflight 原子命名检查
  -o, --output <dir>  输出目录
  --skip-preflight    跳过依赖检查
  -h, --help          显示此帮助信息
EOF
    exit 0
}

# 解析参数
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --example)
                EXAMPLE_NUM="$2"
                shift 2
                ;;
            -i|--input)
                INPUT_PDB="$2"
                shift 2
                ;;
            -c|--config)
                CONFIG_YAML="$2"
                shift 2
                ;;
            --recipe)
                RECIPE_YAML="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --skip-preflight)
                SKIP_PREFLIGHT=true
                shift
                ;;
            -h|--help)
                usage
                ;;
            *)
                warn "未知参数: $1"
                shift
                ;;
        esac
    done
}

# ============================================================================
# 初始化 Conda（仅在需要时，且不依赖 conda init）
# ============================================================================
init_conda_if_needed() {
    # 如果 htpolynet 已在 PATH 中，直接使用，无需 conda
    if command -v htpolynet &> /dev/null; then
        return 0
    fi
    
    # 尝试初始化 conda（不依赖 conda init）
    local conda_base=""
    
    # 方法 1: 使用 conda info --base
    if command -v conda &> /dev/null; then
        conda_base=$(conda info --base 2>/dev/null || echo "")
    fi
    
    # 方法 2: 常见路径
    if [[ -z "$conda_base" ]] || [[ ! -d "$conda_base" ]]; then
        for candidate in "/home/$USER/miniforge3" "/home/$USER/miniconda3" "/home/edu/soft/miniforge3" "/opt/conda"; do
            if [[ -d "$candidate" ]]; then
                conda_base="$candidate"
                break
            fi
        done
    fi
    
    if [[ -n "$conda_base" ]] && [[ -f "$conda_base/etc/profile.d/conda.sh" ]]; then
        info "初始化 Conda: $conda_base"
        # shellcheck source=/dev/null
        source "$conda_base/etc/profile.d/conda.sh"
        
        # 尝试激活可能的环境
        for env in "iongel311" "iongel" "htpolynet" "base"; do
            if conda env list 2>/dev/null | grep -q "^${env} "; then
                info "激活 Conda 环境: $env"
                conda activate "$env" 2>/dev/null || true
                break
            fi
        done
    fi
}

# ============================================================================
# PDB 原子命名 Preflight 检查
# ============================================================================
pdb_preflight_check() {
    local recipe_yaml="$1"
    
    if [[ -z "$recipe_yaml" ]] || [[ ! -f "$recipe_yaml" ]]; then
        info "未指定 --recipe，跳过 PDB 原子命名检查"
        return 0
    fi
    
    step "╔════════════════════════════════════════════════════════════╗"
    step "║        PDB 原子命名 Preflight 检查                         ║"
    step "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    local check_script="${SCRIPT_DIR}/check_pdb_atomnames.py"
    
    if [[ ! -f "$check_script" ]]; then
        warn "找不到 check_pdb_atomnames.py，跳过检查"
        return 0
    fi
    
    info "使用 recipe: $recipe_yaml"
    
    if python3 "$check_script" --recipe "$recipe_yaml"; then
        info "✓ PDB 原子命名检查通过"
        return 0
    else
        error "✗ PDB 原子命名检查失败"
        echo ""
        error "HTPolyNet 将无法识别反应位点！"
        echo "请按上述提示修复 PDB 原子命名后重试。"
        echo ""
        echo "如需跳过检查（不推荐），使用 --skip-preflight"
        return 1
    fi
}

# ============================================================================
# Preflight 依赖检查
# ============================================================================
preflight_check() {
    if [[ "$SKIP_PREFLIGHT" == "true" ]]; then
        info "跳过 preflight 检查 (--skip-preflight)"
        return 0
    fi
    
    step "╔════════════════════════════════════════════════════════════╗"
    step "║              Preflight 依赖检查                            ║"
    step "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    local critical_missing=()
    
    # 检查 htpolynet
    if command -v htpolynet &> /dev/null; then
        local version=$(htpolynet --version 2>/dev/null || echo "unknown")
        info "✓ htpolynet: $(which htpolynet)"
        info "  版本: $version"
    else
        critical_missing+=("htpolynet")
        error "✗ htpolynet 未找到"
        echo "  安装方法: pip install htpolynet"
        echo "  或: conda install -c conda-forge htpolynet"
    fi
    
    # 检查 obabel (OpenBabel) - 关键依赖
    if command -v obabel &> /dev/null; then
        local obabel_version=$(obabel -V 2>&1 | head -1)
        info "✓ obabel: $(which obabel)"
        info "  版本: $obabel_version"
    else
        critical_missing+=("obabel")
        error "✗ obabel (OpenBabel) 未找到"
        echo ""
        echo "  这是 HTPolyNet 的关键依赖！"
        echo "  安装方法: sudo apt-get install -y openbabel"
        echo ""
    fi
    
    # 检查 gmx (GROMACS)
    if command -v gmx &> /dev/null; then
        info "✓ gmx: $(which gmx)"
    else
        warn "⚠ gmx (GROMACS) 未找到"
        echo "  HTPolyNet 需要 GROMACS 运行 MD 模拟"
        echo "  安装方法: sudo apt-get install -y gromacs"
    fi
    
    echo ""
    
    # 如果有缺失的关键依赖，退出
    if [[ ${#critical_missing[@]} -gt 0 ]]; then
        error "╔════════════════════════════════════════════════════════════╗"
        error "║  缺少关键依赖: ${critical_missing[*]}"
        error "╚════════════════════════════════════════════════════════════╝"
        echo ""
        echo "下一步排障:"
        echo "  1. 安装缺失的依赖（见上方提示）"
        echo "  2. 确保 PATH 包含这些工具"
        echo "  3. 重新运行本脚本"
        exit 1
    fi
}

# ============================================================================
# 打印产物文件位置
# ============================================================================
print_output_files() {
    local output="$1"
    
    echo ""
    info "=========================================="
    info " HTPolyNet 产物文件"
    info "=========================================="
    
    # 查找 .gro 文件
    local gro_files=$(find "$output" -name "*.gro" -type f 2>/dev/null | head -10)
    if [[ -n "$gro_files" ]]; then
        info "结构文件 (.gro):"
        echo "$gro_files" | while read -r f; do
            echo "  ✓ $f"
        done
    else
        warn "未找到 .gro 文件"
    fi
    
    # 查找 topol.top
    local top_files=$(find "$output" -name "*.top" -type f 2>/dev/null | head -5)
    if [[ -n "$top_files" ]]; then
        info "拓扑文件 (.top):"
        echo "$top_files" | while read -r f; do
            echo "  ✓ $f"
        done
    else
        warn "未找到 .top 文件"
    fi
    
    # 查找 .itp 文件
    local itp_count=$(find "$output" -name "*.itp" -type f 2>/dev/null | wc -l)
    if [[ "$itp_count" -gt 0 ]]; then
        info "包含文件 (.itp): ${itp_count} 个"
    fi
    
    echo ""
}

# ============================================================================
# 打印下一步提示
# ============================================================================
print_next_steps() {
    local status="$1"
    local output="$2"
    local log_file="${output}/htpolynet.log"
    
    echo ""
    if [[ "$status" == "success" ]]; then
        success "╔════════════════════════════════════════════════════════════╗"
        success "║                 HTPolyNet 运行成功！                        ║"
        success "╚════════════════════════════════════════════════════════════╝"
        
        print_output_files "$output"
        
        info "下一步操作:"
        echo "  1. 检查产物: ls -la $output/"
        echo "  2. 运行 GROMACS: ./scripts/run_gmx.sh -i $output -o ${output%/*}/gmx"
        echo ""
        # WSL 提示
        if command -v explorer.exe &> /dev/null; then
            info "在 Windows 资源管理器中打开:"
            local open_cmd="${SCRIPT_DIR}/open_in_windows.sh"
            if [[ -x "$open_cmd" ]]; then
                echo "  $open_cmd $output"
            else
                echo "  explorer.exe \$(wslpath -w $output)"
            fi
        fi
    else
        error "╔════════════════════════════════════════════════════════════╗"
        error "║                 HTPolyNet 运行失败                          ║"
        error "╚════════════════════════════════════════════════════════════╝"
        echo ""
        warn "排障步骤:"
        echo ""
        echo "  1. 检查日志文件:"
        echo "     cat $log_file | tail -100"
        echo ""
        echo "  2. 搜索错误关键字:"
        echo "     grep -i 'error\\|failed\\|exception\\|no reaction' $log_file"
        echo ""
        echo "  3. 常见问题:"
        echo "     - 'No reaction sites found': PDB 原子命名与 cfg 不匹配"
        echo "       检查: python3 scripts/check_pdb_atomnames.py <your.pdb>"
        echo "     - 'obabel: command not found': 安装 OpenBabel"
        echo "       修复: sudo apt-get install -y openbabel"
        echo "     - 'gmx: command not found': 安装 GROMACS"
        echo "       修复: sudo apt-get install -y gromacs"
        echo ""
        echo "  4. 检查输出目录:"
        echo "     ls -la $output/"
    fi
}

# ============================================================================
# 创建不依赖 conda activate 的 run.sh 替代脚本
# ============================================================================
create_fixed_run_script() {
    local example_dir="$1"
    local fixed_script="${example_dir}/run_no_conda.sh"
    
    info "创建不依赖 conda 的执行脚本: $fixed_script"
    
    # 获取原 run.sh 的核心命令（去掉 conda activate 相关行）
    local original_run="${example_dir}/run.sh"
    
    cat > "$fixed_script" << 'SCRIPT_HEADER'
#!/bin/bash
# 自动生成的执行脚本 - 不依赖 conda activate
set -euo pipefail

# 确保系统工具在 PATH 中
export PATH="/usr/bin:/usr/local/bin:$PATH"

echo "=== 使用系统 PATH 中的工具 ==="
echo "htpolynet: $(which htpolynet 2>/dev/null || echo 'NOT FOUND')"
echo "obabel: $(which obabel 2>/dev/null || echo 'NOT FOUND')"
echo "gmx: $(which gmx 2>/dev/null || echo 'NOT FOUND')"
echo ""

SCRIPT_HEADER

    # 从原 run.sh 提取核心命令（跳过 shebang、conda 相关行）
    if [[ -f "$original_run" ]]; then
        grep -v "^#!" "$original_run" | \
        grep -v "conda activate" | \
        grep -v "conda deactivate" | \
        grep -v "^#.*conda" | \
        grep -v "^[[:space:]]*$" >> "$fixed_script" || true
    else
        # 如果没有 run.sh，使用默认命令
        cat >> "$fixed_script" << 'DEFAULT_CMD'
# 查找 cfg 文件并运行
CFG_FILE=$(find . -maxdepth 1 -name "*.yaml" -o -name "*.yml" | head -1)
if [[ -z "$CFG_FILE" ]]; then
    echo "ERROR: No cfg file found"
    exit 1
fi
echo "Running: htpolynet run $CFG_FILE"
htpolynet run "$CFG_FILE"
DEFAULT_CMD
    fi
    
    chmod +x "$fixed_script"
    echo ""
}

# ============================================================================
# 运行官方示例
# ============================================================================
run_example() {
    local example_num="$1"
    local output="$2"
    
    mkdir -p "$output"
    local log_file="${output}/htpolynet.log"
    
    step "╔════════════════════════════════════════════════════════════╗"
    step "║  运行 HTPolyNet 官方示例 (n=$example_num)                   ║"
    step "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    info "输出目录: $output"
    info "日志文件: $log_file"
    echo ""
    
    # 切换到输出目录
    pushd "$output" > /dev/null
    
    # 获取示例
    step "STEP 1: 获取示例 (htpolynet fetch-example -n $example_num)..."
    echo ""
    
    if htpolynet fetch-example -n "$example_num" 2>&1 | tee htpolynet.log; then
        info "✓ 示例获取成功"
    else
        error "获取示例失败"
        popd > /dev/null
        print_next_steps "failure" "$output"
        exit 1
    fi
    
    # 查找示例目录
    local example_dir=""
    case "$example_num" in
        0) example_dir="0-liquid-styrene" ;;
        1) example_dir="1-polystyrene" ;;
        2) example_dir="2-polymethylstyrene" ;;
        3) example_dir="3-bisgma-styrene-thermoset" ;;
        4) example_dir="4-pacm-dgeba-epoxy-thermoset" ;;
        5) example_dir="5-dfda-fde-epoxy-thermoset" ;;
    esac
    
    if [[ -z "$example_dir" ]] || [[ ! -d "$example_dir" ]]; then
        example_dir=$(find . -maxdepth 1 -type d ! -name "." | head -1 | sed 's|^\./||')
    fi
    
    if [[ -z "$example_dir" ]] || [[ ! -d "$example_dir" ]]; then
        error "找不到示例目录"
        ls -la
        popd > /dev/null
        print_next_steps "failure" "$output"
        exit 1
    fi
    
    info "示例目录: $example_dir"
    cd "$example_dir"
    
    step "STEP 2: 执行 HTPolyNet (预计需要 1-10 分钟)..."
    echo ""
    
    # 创建不依赖 conda 的执行脚本
    create_fixed_run_script "$(pwd)"
    
    local start_time=$(date +%s)
    local run_status=0
    
    # 执行不依赖 conda 的脚本
    if bash run_no_conda.sh 2>&1 | tee -a ../htpolynet.log; then
        run_status=0
    else
        run_status=$?
        warn "脚本返回非零状态: $run_status"
    fi
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    # 检查是否生成了关键产物
    local has_gro=$(find . -name "*.gro" -type f 2>/dev/null | head -1)
    local has_top=$(find . -name "*.top" -type f 2>/dev/null | head -1)
    
    popd > /dev/null
    
    if [[ -n "$has_gro" ]] || [[ -n "$has_top" ]]; then
        info "✓ HTPolyNet 完成，耗时 ${duration} 秒"
        print_next_steps "success" "$output"
    else
        if [[ $run_status -eq 0 ]]; then
            warn "HTPolyNet 完成但未找到 .gro/.top 文件"
            print_next_steps "success" "$output"
        else
            error "HTPolyNet 失败，耗时 ${duration} 秒"
            print_next_steps "failure" "$output"
            exit 1
        fi
    fi
}

# ============================================================================
# 运行自定义配置
# ============================================================================
run_custom() {
    local input="$1"
    local config="$2"
    local output="$3"
    
    mkdir -p "$output"
    local log_file="${output}/htpolynet.log"
    
    step "╔════════════════════════════════════════════════════════════╗"
    step "║  运行 HTPolyNet (自定义配置)                               ║"
    step "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    info "输入 PDB: $input"
    info "配置文件: $config"
    info "输出目录: $output"
    echo ""
    
    # 检查文件存在
    if [[ ! -f "$input" ]]; then
        error "输入文件不存在: $input"
        echo ""
        echo "排障: 请检查 Packmol 是否已成功运行并生成 PDB 文件"
        echo "      ls -la $(dirname "$input")/"
        exit 2
    fi
    
    if [[ ! -f "$config" ]]; then
        error "配置文件不存在: $config"
        echo ""
        echo "排障: 请创建 HTPolyNet 配置文件"
        echo "      参考: htpolynet fetch-example -n 0 获取示例配置"
        exit 2
    fi
    
    # 复制文件到输出目录
    cp "$input" "${output}/input.pdb"
    cp "$config" "${output}/cfg.yaml"
    
    pushd "$output" > /dev/null
    
    step "STEP 1: 验证配置 (htpolynet input-check)..."
    htpolynet input-check cfg.yaml 2>&1 | tee htpolynet.log || true
    echo ""
    
    step "STEP 2: 运行 HTPolyNet (htpolynet run)..."
    info "这可能需要较长时间..."
    echo ""
    
    local start_time=$(date +%s)
    if htpolynet run cfg.yaml 2>&1 | tee -a htpolynet.log; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        info "✓ HTPolyNet 运行完成，耗时 ${duration} 秒"
        popd > /dev/null
        print_next_steps "success" "$output"
    else
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        error "HTPolyNet 运行失败，耗时 ${duration} 秒"
        popd > /dev/null
        print_next_steps "failure" "$output"
        exit 1
    fi
}

# ============================================================================
# 主函数
# ============================================================================
main() {
    echo ""
    echo "=============================================="
    echo " HTPolyNet 交联聚合物模拟驱动脚本 (v2.2)"
    echo " 特点: 不依赖 conda init / conda activate"
    echo "=============================================="
    echo ""
    
    parse_args "$@"
    
    # 初始化 Conda（仅在需要时，不依赖 conda init）
    init_conda_if_needed
    
    # Preflight 检查
    preflight_check
    
    # 判断运行模式
    if [[ -n "$EXAMPLE_NUM" ]]; then
        # 官方示例模式
        if [[ -z "$OUTPUT_DIR" ]]; then
            OUTPUT_DIR="${PROJECT_ROOT}/outputs/smoke_full/htpolynet"
        fi
        run_example "$EXAMPLE_NUM" "$OUTPUT_DIR"
        
    elif [[ -n "$INPUT_PDB" ]] && [[ -n "$CONFIG_YAML" ]]; then
        # 自定义模式
        if [[ -z "$OUTPUT_DIR" ]]; then
            OUTPUT_DIR="${PROJECT_ROOT}/outputs/htpolynet"
        fi
        
        # PDB 原子命名检查（如果提供了 recipe）
        if [[ "$SKIP_PREFLIGHT" != "true" ]] && [[ -n "$RECIPE_YAML" ]]; then
            pdb_preflight_check "$RECIPE_YAML" || exit 1
            echo ""
        fi
        
        run_custom "$INPUT_PDB" "$CONFIG_YAML" "$OUTPUT_DIR"
        
    else
        error "请指定运行模式"
        echo ""
        echo "【推荐】先用官方示例验证环境:"
        echo "  $0 --example 0 -o outputs/smoke_full/htpolynet"
        echo ""
        echo "然后再用自定义配置:"
        echo "  $0 -i outputs/smoke/packmol/gel.pdb -c configs/htpolynet.yaml -o outputs/smoke/htpolynet"
        echo ""
        usage
    fi
}

main "$@"
