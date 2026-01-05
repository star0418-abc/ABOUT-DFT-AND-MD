#!/bin/bash

# Author: Star



# ==============================================================================
# run_packmol.sh - Packmol 凝胶电解质结构生成器
# ==============================================================================
# 功能：调用 make_packmol_from_recipe.py 生成 gel.inp，然后运行 Packmol
#
# 用法:
#   ./scripts/run_packmol.sh -c configs/recipe_smoke.yaml -o outputs/smoke
#   ./scripts/run_packmol.sh outputs/gel.inp  # 兼容旧用法
#
# 参数:
#   -c, --config   : 配方 YAML 文件路径
#   -o, --output   : 输出目录（会在其下创建 packmol/ 子目录）
#   <inp_file>     : 兼容模式：直接指定 .inp 文件运行
#
# 输出结构:
#   <output_dir>/packmol/
#     ├── gel.inp         # Packmol 输入文件
#     ├── gel.pdb         # Packmol 输出结构
#     ├── packmol.log     # Packmol 日志
#     └── summary.json    # 计算摘要
#
# 返回值: 0=成功, 1=Packmol未安装, 2=输入错误, 3=Packmol失败
# ==============================================================================

set -euo pipefail

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

# 获取脚本所在目录（项目根目录）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# 默认值
CONFIG_FILE=""
OUTPUT_DIR=""
LEGACY_INP_FILE=""

# 打印函数
info()  { echo -e "${GREEN}[INFO]${NC} $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
step()  { echo -e "${CYAN}[STEP]${NC} $1"; }

# 打印使用帮助
usage() {
    cat << EOF
用法: $0 [选项]

选项:
  -c, --config <yaml>    配方 YAML 文件路径 (如 configs/recipe_smoke.yaml)
  -o, --output <dir>     输出目录 (如 outputs/smoke)
  -h, --help             显示此帮助信息

兼容模式:
  $0 <inp_file>          直接运行已有的 .inp 文件

示例:
  $0 -c configs/recipe_smoke.yaml -o outputs/smoke
  $0 outputs/gel.inp
EOF
    exit 0
}

# 解析参数
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -c|--config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -h|--help)
                usage
                ;;
            *)
                # 兼容旧用法：直接指定 .inp 文件
                if [[ -f "$1" ]] || [[ "$1" == *.inp ]]; then
                    LEGACY_INP_FILE="$1"
                fi
                shift
                ;;
        esac
    done
}

# 检查 Packmol 是否安装
check_packmol() {
    if ! command -v packmol &> /dev/null; then
        error "Packmol 未安装或不在 PATH 中"
        echo ""
        echo "安装方法:"
        echo "  Ubuntu/Debian: sudo apt install packmol"
        echo "  Conda: conda install -c conda-forge packmol"
        exit 1
    fi
    info "找到 Packmol: $(which packmol)"
}

# 检查 Packmol 日志中的错误
check_packmol_log() {
    local log_file="$1"
    local has_error=0
    
    # 检查常见错误模式
    if grep -qi "ERROR" "$log_file" 2>/dev/null; then
        warn "日志中发现 ERROR 关键字"
        has_error=1
    fi
    
    if grep -qi "FAILED" "$log_file" 2>/dev/null; then
        warn "日志中发现 FAILED 关键字"
        has_error=1
    fi
    
    if grep -qi "Could not verify" "$log_file" 2>/dev/null; then
        warn "日志中发现 'Could not verify' - Packmol 可能未收敛"
        has_error=1
    fi
    
    # 检查 overlap（但 overlap 可能在中间出现，最终解决了）
    if grep -qi "overlap" "$log_file" 2>/dev/null; then
        # 只在最后 20 行检查 overlap
        if tail -20 "$log_file" | grep -qi "overlap"; then
            warn "日志末尾发现 overlap 警告"
        fi
    fi
    
    # 检查 STOP 退出状态（STOP 174 是正常结束）
    if grep -q "STOP 174" "$log_file" 2>/dev/null; then
        info "Packmol 正常退出 (STOP 174)"
    elif grep -q "STOP" "$log_file" 2>/dev/null; then
        local stop_code=$(grep "STOP" "$log_file" | tail -1)
        warn "Packmol 异常退出: $stop_code"
        has_error=1
    fi
    
    return $has_error
}

# 运行新模式（使用 -c/-o 参数）
run_new_mode() {
    local config="$1"
    local output_base="$2"
    
    # 创建 packmol 子目录
    local packmol_dir="${output_base}/packmol"
    mkdir -p "$packmol_dir"
    
    step "=========================================="
    step " STEP 1: 生成 Packmol 输入文件"
    step "=========================================="
    
    # 调用 Python 脚本生成 gel.inp 和 summary.json
    local python_script="${PROJECT_ROOT}/scripts/make_packmol_from_recipe.py"
    if [[ ! -f "$python_script" ]]; then
        error "找不到 Python 脚本: $python_script"
        exit 2
    fi
    
    info "配方文件: $config"
    info "输出目录: $packmol_dir"
    
    # 运行 Python 脚本
    python3 "$python_script" -c "$config" -o "$packmol_dir"
    
    local inp_file="${packmol_dir}/gel.inp"
    local log_file="${packmol_dir}/packmol.log"
    
    if [[ ! -f "$inp_file" ]]; then
        error "Python 脚本未能生成 gel.inp"
        exit 2
    fi
    
    step "=========================================="
    step " STEP 2: 运行 Packmol"
    step "=========================================="

# 运行 Packmol
    info "执行: packmol < $inp_file"
    local start_time=$(date +%s)
    
    # 切换到 packmol_dir 运行（虽然使用绝对路径，但日志更清晰）
    pushd "$packmol_dir" > /dev/null
    if packmol < gel.inp > packmol.log 2>&1; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        info "Packmol 运行完成，耗时 ${duration} 秒"
    else
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        warn "Packmol 返回非零状态，耗时 ${duration} 秒"
    fi
    popd > /dev/null
    
    step "=========================================="
    step " STEP 3: 检查结果"
    step "=========================================="
        
    # 从 inp 文件提取输出 PDB 路径
    local output_pdb=$(grep -E "^output\s+" "$inp_file" | awk '{print $2}')
    
    # 检查输出文件
    if [[ -n "$output_pdb" ]] && [[ -f "$output_pdb" ]]; then
        local file_size=$(ls -lh "$output_pdb" | awk '{print $5}')
        local atom_count=$(grep -c "^ATOM\|^HETATM" "$output_pdb" 2>/dev/null || echo "0")
        
        if [[ "$atom_count" -gt 0 ]]; then
            info "✓ 输出结构: $output_pdb"
            info "  文件大小: $file_size, 原子数: $atom_count"
        else
            error "输出文件存在但原子数为 0"
            check_packmol_log "$log_file"
            exit 3
        fi
    else
        error "未找到输出 PDB 文件: $output_pdb"
        check_packmol_log "$log_file"
        exit 3
    fi
    
    # 检查日志
    if ! check_packmol_log "$log_file"; then
        warn "日志中存在警告，请检查 $log_file"
    fi
    
    info "✓ 日志文件: $log_file"
    info "✓ 摘要文件: ${packmol_dir}/summary.json"
    
    echo ""
    info "=========================================="
    info " Packmol 完成！"
    info "=========================================="
    echo ""
    info "输出目录: $packmol_dir"
    
    # WSL 提示：用 explorer.exe 打开目录
    if command -v explorer.exe &> /dev/null; then
    echo ""
        info "在 Windows 资源管理器中打开:"
        echo "  explorer.exe $(wslpath -w "$packmol_dir" 2>/dev/null || echo "$packmol_dir")"
    fi
    
    # Density Trap 提醒
    echo ""
    warn "提醒: Packmol 初始密度可能低于目标密度（膨胀盒子策略）"
    warn "请使用 GROMACS NPT 将体系压缩到真实密度"
    warn "推荐流程: EM → NVT → Berendsen NPT (快速压缩) → PR NPT"
}

# 运行兼容模式（直接运行 .inp 文件）
run_legacy_mode() {
    local inp_file="$1"
    
    if [[ ! -f "$inp_file" ]]; then
        error "输入文件不存在: $inp_file"
        exit 2
    fi
    
    local log_file="${inp_file%.inp}.log"
    local inp_dir=$(dirname "$inp_file")
    
    step "=========================================="
    step " 运行 Packmol (兼容模式)"
    step "=========================================="
    
    info "输入文件: $inp_file"
    info "日志文件: $log_file"
    
    # 显示 seed 信息
    local seed=$(grep -E "^seed\s+" "$inp_file" 2>/dev/null | awk '{print $2}' || echo "未指定")
    info "Packmol seed: $seed"
    
    local start_time=$(date +%s)
    
    if packmol < "$inp_file" > "$log_file" 2>&1; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        info "Packmol 运行完成，耗时 ${duration} 秒"
    else
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        warn "Packmol 返回非零状态，耗时 ${duration} 秒"
    fi
    
    # 检查输出
    local output_pdb=$(grep -E "^output\s+" "$inp_file" | awk '{print $2}')
    if [[ -n "$output_pdb" ]] && [[ -f "$output_pdb" ]]; then
        local atom_count=$(grep -c "^ATOM\|^HETATM" "$output_pdb" 2>/dev/null || echo "0")
        if [[ "$atom_count" -gt 0 ]]; then
            info "✓ 输出结构: $output_pdb (原子数: $atom_count)"
        else
            error "输出文件原子数为 0"
            exit 3
        fi
    else
        error "未找到输出文件: $output_pdb"
        exit 3
    fi
    
    check_packmol_log "$log_file" || true
    info "✓ 完成！日志: $log_file"
}

# 主函数
main() {
    echo ""
    echo "=============================================="
    echo " Packmol 凝胶电解质结构生成器 (v2.1)"
    echo "=============================================="
    echo ""
    
    parse_args "$@"
    check_packmol
    
    # 判断运行模式
    if [[ -n "$CONFIG_FILE" ]] && [[ -n "$OUTPUT_DIR" ]]; then
        # 新模式：使用 -c/-o 参数
        if [[ ! -f "$CONFIG_FILE" ]]; then
            error "配方文件不存在: $CONFIG_FILE"
            exit 2
        fi
        run_new_mode "$CONFIG_FILE" "$OUTPUT_DIR"
    elif [[ -n "$LEGACY_INP_FILE" ]]; then
        # 兼容模式：直接运行 .inp 文件
        run_legacy_mode "$LEGACY_INP_FILE"
    else
        # 默认：尝试运行 outputs/gel.inp
        if [[ -f "outputs/gel.inp" ]]; then
            run_legacy_mode "outputs/gel.inp"
        else
            error "请指定参数。使用 -h 查看帮助。"
            echo ""
            usage
        fi
    fi
}

main "$@"
