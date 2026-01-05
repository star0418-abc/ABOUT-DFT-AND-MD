#!/bin/bash

# Author: Star



# ==============================================================================
# run_gmx.sh - GROMACS 分子动力学模拟驱动脚本 (v2.2)
# ==============================================================================
# 功能：运行 GROMACS 标准流程 (EM → NVT → NPT)
#
# 特点：
#   - 自动递归查找 .gro 和 .top 文件
#   - 自动检测 GPU，无 GPU 时使用 CPU-only 模式
#   - 详细的错误提示和排障指引
#
# 用法:
#   ./scripts/run_gmx.sh -i outputs/smoke_full/htpolynet -o outputs/smoke_full/gmx
#   ./scripts/run_gmx.sh -i system.gro -p topol.top -o outputs/gmx
#
# 输出结构:
#   <output_dir>/
#     ├── input.gro         # 输入结构副本
#     ├── topol.top         # 拓扑文件副本
#     ├── em.tpr, em.gro, em.log, em.edr
#     ├── nvt.tpr, nvt.gro, nvt.log, nvt.edr, nvt.xtc
#     └── npt.tpr, npt.gro, npt.log, npt.edr, npt.xtc
#
# 依赖: GROMACS (gmx)
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

# 参数默认值
INPUT_PATH=""
TOPOL_FILE=""
OUTPUT_DIR=""
MDP_DIR="${PROJECT_ROOT}/configs/mdp"
SKIP_EM=false
SKIP_NVT=false

# GPU 检测结果
HAS_GPU=false
MDRUN_GPU_FLAGS=""

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

选项:
  -i, --input <path>     输入结构文件 (.gro/.pdb) 或 HTPolyNet 输出目录
  -p, --topol <file>     拓扑文件 (默认: 自动查找)
  -o, --output <dir>     输出目录
  --mdpdir <dir>         MDP 文件目录 (默认: configs/mdp)
  --skip-em              跳过能量最小化
  --skip-nvt             跳过 NVT 平衡
  -h, --help             显示此帮助信息

示例:
  # 使用 HTPolyNet 输出目录（自动递归查找 .gro 和 topol.top）
  $0 -i outputs/smoke_full/htpolynet -o outputs/smoke_full/gmx
  
  # 使用具体文件
  $0 -i system.gro -p topol.top -o outputs/gmx

流程: EM (能量最小化) → NVT (恒温平衡) → NPT (恒压平衡)
EOF
    exit 0
}

# 解析参数
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input)
                INPUT_PATH="$2"
                shift 2
                ;;
            -p|--topol)
                TOPOL_FILE="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --mdpdir)
                MDP_DIR="$2"
                shift 2
                ;;
            --skip-em)
                SKIP_EM=true
                shift
                ;;
            --skip-nvt)
                SKIP_NVT=true
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
# GPU 检测
# ============================================================================
detect_gpu() {
    step "检测 GPU..."
    
    HAS_GPU=false
    MDRUN_GPU_FLAGS=""
    
    # 方法 1: 检查 nvidia-smi
    if command -v nvidia-smi &> /dev/null; then
        if nvidia-smi -L 2>/dev/null | grep -q "GPU"; then
            local gpu_count=$(nvidia-smi -L 2>/dev/null | wc -l)
            info "✓ 检测到 ${gpu_count} 个 NVIDIA GPU"
            HAS_GPU=true
        else
            warn "nvidia-smi 存在但无可用 GPU"
        fi
    else
        info "未检测到 nvidia-smi"
    fi
    
    # 设置 mdrun 参数
    if [[ "$HAS_GPU" == "true" ]]; then
        info "使用 GPU 加速模式"
        MDRUN_GPU_FLAGS=""  # 让 GROMACS 自动选择
    else
        warn "使用 CPU-only 模式（无 GPU 或检测失败）"
        MDRUN_GPU_FLAGS="-nb cpu -pme cpu -bonded cpu -update cpu"
    fi
    
    echo ""
}

# Preflight 检查
preflight_check() {
    step "╔════════════════════════════════════════════════════════════╗"
    step "║              Preflight 依赖检查                            ║"
    step "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    if ! command -v gmx &> /dev/null; then
        error "✗ gmx (GROMACS) 未找到"
        echo ""
        echo "安装方法:"
        echo "  Ubuntu/Debian: sudo apt-get install -y gromacs"
        echo "  Conda: conda install -c conda-forge gromacs"
        echo ""
        echo "或者加载已安装的 GROMACS:"
        echo "  source /usr/local/gromacs/bin/GMXRC"
        exit 1
    fi
    
    local version=$(gmx --version 2>&1 | grep "GROMACS version" | head -1 || echo "unknown")
    info "✓ gmx: $(which gmx)"
    info "  版本: $version"
    echo ""
    
    # GPU 检测
    detect_gpu
}

# ============================================================================
# 递归查找 HTPolyNet 输出中的结构和拓扑文件
# ============================================================================
find_htpolynet_outputs() {
    local htpolynet_dir="$1"
    
    step "在 HTPolyNet 输出目录中递归查找文件..."
    echo ""
    
    # ========== 查找 .gro 文件 ==========
    local gro_file=""
    
    # 优先级 1: 特定命名的文件
    for pattern in "system.gro" "conf.gro" "*final*.gro" "*equilibrated*.gro" "*cured*.gro"; do
        gro_file=$(find "$htpolynet_dir" -name "$pattern" -type f 2>/dev/null | head -1)
        if [[ -n "$gro_file" ]]; then
            info "  找到结构文件 (匹配 $pattern): $gro_file"
            break
        fi
    done
    
    # 优先级 2: 任意 .gro 文件
    if [[ -z "$gro_file" ]]; then
        gro_file=$(find "$htpolynet_dir" -name "*.gro" -type f 2>/dev/null | head -1)
        if [[ -n "$gro_file" ]]; then
            warn "  使用找到的第一个 .gro 文件: $gro_file"
            warn "  （未找到 system.gro/conf.gro，可能不是最终结构）"
        fi
    fi
    
    # 优先级 3: .pdb 文件
    if [[ -z "$gro_file" ]]; then
        gro_file=$(find "$htpolynet_dir" -name "*.pdb" -type f 2>/dev/null | head -1)
        if [[ -n "$gro_file" ]]; then
            warn "  使用 PDB 文件: $gro_file"
        fi
    fi
    
    # ========== 查找 .top 文件 ==========
    local top_file=""
    
    # 优先级 1: topol.top
    top_file=$(find "$htpolynet_dir" -name "topol.top" -type f 2>/dev/null | head -1)
    
    # 优先级 2: 任意 .top 文件
    if [[ -z "$top_file" ]]; then
        top_file=$(find "$htpolynet_dir" -name "*.top" -type f 2>/dev/null | head -1)
        if [[ -n "$top_file" ]]; then
            warn "  未找到 topol.top，使用: $top_file"
        fi
    fi
    
    # ========== 检查是否找到必需文件 ==========
    if [[ -z "$gro_file" ]]; then
        error "在 $htpolynet_dir 中找不到 .gro 或 .pdb 文件"
        echo ""
        echo "目录内容 (前 30 行):"
        find "$htpolynet_dir" -type f 2>/dev/null | head -30
        echo ""
        echo "排障:"
        echo "  1. HTPolyNet 可能未完成 build 阶段"
        echo "  2. 检查 htpolynet.log 是否有错误"
        echo "  3. 确认 HTPolyNet 运行完成后再运行本脚本"
        exit 2
    fi
    
    if [[ -z "$top_file" ]]; then
        error "在 $htpolynet_dir 中找不到 .top 拓扑文件"
        echo ""
        echo "目录内容 (前 30 行):"
        find "$htpolynet_dir" -type f 2>/dev/null | head -30
        echo ""
        echo "排障:"
        echo "  1. HTPolyNet 可能未跑到 build 阶段（仍在 parameterize 阶段）"
        echo "  2. 检查 htpolynet.log 看是否有错误:"
        echo "     grep -i 'error\\|failed' $htpolynet_dir/htpolynet.log"
        echo "  3. 确认 htpolynet run 完成后再运行本脚本"
        exit 2
    fi
    
    info "  拓扑文件: $top_file"
    
    # 设置全局变量
    INPUT_GRO="$gro_file"
    if [[ -z "$TOPOL_FILE" ]]; then
        TOPOL_FILE="$top_file"
    fi
}

# 确保 MDP 文件存在
ensure_mdp_files() {
    info "检查 MDP 文件..."
    
    if [[ ! -d "$MDP_DIR" ]]; then
        error "MDP 目录不存在: $MDP_DIR"
        echo ""
        echo "排障:"
        echo "  1. 创建 MDP 目录: mkdir -p $MDP_DIR"
        echo "  2. 添加 MDP 文件: em.mdp, nvt.mdp, npt.mdp"
        echo "  3. 或使用 --mdpdir 指定其他目录"
        exit 2
    fi
    
    local required_mdp=("em.mdp" "nvt.mdp" "npt.mdp")
    local missing=()
    
    for mdp in "${required_mdp[@]}"; do
        if [[ -f "${MDP_DIR}/${mdp}" ]]; then
            info "  ✓ ${mdp}"
        else
            missing+=("$mdp")
            error "  ✗ ${mdp} (缺失)"
        fi
    done
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        echo ""
        error "缺少 MDP 文件: ${missing[*]}"
        echo ""
        echo "排障:"
        echo "  1. 请在 $MDP_DIR 目录下创建以下文件:"
        for mdp in "${missing[@]}"; do
            echo "     - ${mdp}"
        done
        echo ""
        echo "  2. 可参考 GROMACS 官方教程获取示例 MDP 文件"
        exit 2
    fi
    
    info "MDP 目录: $MDP_DIR"
}

# 运行 grompp + mdrun (CPU-only 友好)
run_gmx_step() {
    local step_name="$1"
    local mdp_file="$2"
    local input_gro="$3"
    local topol="$4"
    local output_prefix="$5"
    
    step "──────────────────────────────────────────────────────────────"
    step " ${step_name}"
    step "──────────────────────────────────────────────────────────────"
    echo ""
    
    info "MDP: $mdp_file"
    info "输入: $input_gro"
    if [[ -n "$MDRUN_GPU_FLAGS" ]]; then
        info "模式: CPU-only ($MDRUN_GPU_FLAGS)"
    else
        info "模式: 自动（GPU 如可用）"
    fi
    
    local start_time=$(date +%s)
    
    # grompp
    info "运行 grompp..."
    if ! gmx grompp \
        -f "$mdp_file" \
        -c "$input_gro" \
        -p "$topol" \
        -o "${output_prefix}.tpr" \
        -maxwarn 10 \
        2>&1 | tee "${output_prefix}_grompp.log"; then
        
        error "grompp 失败"
        echo ""
        echo "排障:"
        echo "  1. 检查 grompp 日志:"
        echo "     cat ${output_prefix}_grompp.log"
        echo ""
        echo "  2. 常见错误:"
        echo "     - 'number of coordinates ... does not match': 坐标与拓扑原子数不匹配"
        echo "     - 'No default ... atomtypes': 缺少力场参数"
        echo "     - 'atom ... not found': 原子名称不匹配"
        return 1
    fi
    
    # mdrun (带 CPU-only 参数)
    info "运行 mdrun..."
    local mdrun_cmd="gmx mdrun -deffnm $output_prefix -v"
    
    # 添加 CPU-only 参数（如果需要）
    if [[ -n "$MDRUN_GPU_FLAGS" ]]; then
        mdrun_cmd="$mdrun_cmd $MDRUN_GPU_FLAGS"
    fi
    
    if ! eval "$mdrun_cmd" 2>&1 | tee "${output_prefix}.log"; then
        error "mdrun 失败"
        echo ""
        echo "排障:"
        echo "  1. 检查 mdrun 日志:"
        echo "     cat ${output_prefix}.log | tail -100"
        echo ""
        echo "  2. 搜索错误:"
        echo "     grep -i 'error\\|fatal\\|nan' ${output_prefix}.log"
        echo ""
        echo "  3. 常见问题:"
        echo "     - 'Device ID ... did not correspond': 无 GPU 但尝试使用 GPU"
        echo "       已通过 CPU-only 模式解决"
        echo "     - 'NaN': 力场参数问题或初始结构不合理"
        return 1
    fi
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    if [[ -f "${output_prefix}.gro" ]]; then
        info "✓ ${step_name} 完成，耗时 ${duration} 秒"
    else
        error "${step_name} 未生成输出结构"
        return 1
    fi
    
    echo ""
    return 0
}

# 打印成功/失败摘要
print_summary() {
    local output="$1"
    local status="$2"
    
    echo ""
    if [[ "$status" == "success" ]]; then
        success "╔════════════════════════════════════════════════════════════╗"
        success "║               GROMACS 模拟完成！                            ║"
        success "╚════════════════════════════════════════════════════════════╝"
    else
        error "╔════════════════════════════════════════════════════════════╗"
        error "║               GROMACS 模拟失败                              ║"
        error "╚════════════════════════════════════════════════════════════╝"
    fi
    
    echo ""
    info "输出目录: $output"
    echo ""
    
    # 检查关键文件
    info "成功标志 - 文件检查:"
    
    local all_present=true
    local files_check=("em.tpr" "em.gro" "nvt.tpr" "nvt.gro" "npt.tpr" "npt.gro")
    
    for file in "${files_check[@]}"; do
        if [[ -f "${output}/${file}" ]]; then
            echo -e "  ${GREEN}✓${NC} ${file}"
        else
            echo -e "  ${RED}✗${NC} ${file} (缺失)"
            all_present=false
        fi
    done
    
    echo ""
    info "成功判定 grep 命令:"
    echo "  # 检查 NPT 是否正常结束:"
    echo "  grep 'Finished mdrun' ${output}/npt.log"
    echo ""
    echo "  # 检查是否有错误:"
    echo "  grep -i 'error\\|fatal\\|nan' ${output}/*.log"
    
    echo ""
    
    # WSL 提示
    local open_script="${SCRIPT_DIR}/open_in_windows.sh"
    if command -v explorer.exe &> /dev/null; then
        info "在 Windows 资源管理器中打开:"
        if [[ -x "$open_script" ]]; then
            echo "  $open_script $output"
        else
            echo "  explorer.exe \$(wslpath -w $output)"
        fi
    fi
    
    if [[ "$status" != "success" ]]; then
        echo ""
        warn "下一步排障:"
        echo "  1. 检查最近的日志文件:"
        echo "     ls -lt ${output}/*.log | head -3"
        echo ""
        echo "  2. 查看错误详情:"
        echo "     tail -50 \$(ls -t ${output}/*.log | head -1)"
    fi
}

# 主函数
main() {
    echo ""
    echo "=============================================="
    echo " GROMACS 分子动力学模拟驱动脚本 (v2.2)"
    echo " 特点: 自动 GPU 检测，CPU-only 兼容"
    echo "=============================================="
    echo ""
    
    parse_args "$@"
    preflight_check
    
    # 检查输入
    if [[ -z "$INPUT_PATH" ]]; then
        error "请指定输入 (-i)"
        echo ""
        echo "示例:"
        echo "  $0 -i outputs/smoke_full/htpolynet -o outputs/smoke_full/gmx"
        usage
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        OUTPUT_DIR="${PROJECT_ROOT}/outputs/gmx"
    fi
    
    mkdir -p "$OUTPUT_DIR"
    
    # 判断输入类型（目录或文件）
    local INPUT_GRO=""
    
    if [[ -d "$INPUT_PATH" ]]; then
        # 输入是目录（HTPolyNet 输出）
        find_htpolynet_outputs "$INPUT_PATH"
    elif [[ -f "$INPUT_PATH" ]]; then
        # 输入是文件
        INPUT_GRO="$INPUT_PATH"
    else
        error "输入不存在: $INPUT_PATH"
        echo ""
        echo "排障:"
        echo "  1. 检查 HTPolyNet 是否已运行:"
        echo "     ls -la $(dirname "$INPUT_PATH")/"
        echo ""
        echo "  2. 先运行 HTPolyNet:"
        echo "     ./scripts/run_htpolynet.sh --example 0 -o outputs/smoke_full/htpolynet"
        exit 2
    fi
    
    # 查找拓扑文件（如果未指定）
    if [[ -z "$TOPOL_FILE" ]]; then
        local input_dir=$(dirname "$INPUT_GRO")
        
        # 递归查找 topol.top
        TOPOL_FILE=$(find "$input_dir" -name "topol.top" -type f 2>/dev/null | head -1)
        
        # 如果没找到，查找任意 .top
        if [[ -z "$TOPOL_FILE" ]]; then
            TOPOL_FILE=$(find "$input_dir" -name "*.top" -type f 2>/dev/null | head -1)
        fi
        
        # 再往上一级目录查找
        if [[ -z "$TOPOL_FILE" ]]; then
            local parent_dir=$(dirname "$input_dir")
            TOPOL_FILE=$(find "$parent_dir" -name "topol.top" -type f 2>/dev/null | head -1)
        fi
    fi
    
    if [[ -z "$TOPOL_FILE" ]] || [[ ! -f "$TOPOL_FILE" ]]; then
        error "找不到拓扑文件 (topol.top)"
        echo ""
        echo "目录内容 (前 30 行):"
        find "$(dirname "$INPUT_GRO")" -type f 2>/dev/null | head -30
        echo ""
        echo "排障:"
        echo "  1. 指定拓扑文件: -p /path/to/topol.top"
        echo ""
        echo "  2. HTPolyNet 可能未跑到 build 阶段"
        echo "     检查: ls -la $(dirname "$INPUT_GRO")/"
        echo ""
        echo "  3. 如果没有拓扑文件，可能需要用 gmx pdb2gmx 生成"
        exit 2
    fi
    
    echo ""
    info "输入结构: $INPUT_GRO"
    info "拓扑文件: $TOPOL_FILE"
    info "输出目录: $OUTPUT_DIR"
    echo ""
    
    # 确保 MDP 文件存在
    ensure_mdp_files
    echo ""
    
    # 复制输入文件到输出目录
    local work_gro="${OUTPUT_DIR}/input.gro"
    local work_top="${OUTPUT_DIR}/topol.top"
    
    # 如果输入是 PDB，转换为 GRO
    if [[ "$INPUT_GRO" == *.pdb ]]; then
        info "将 PDB 转换为 GRO..."
        gmx editconf -f "$INPUT_GRO" -o "$work_gro" 2>&1 | tee "${OUTPUT_DIR}/editconf.log"
    else
        cp "$INPUT_GRO" "$work_gro"
    fi
    
    cp "$TOPOL_FILE" "$work_top"
    
    # 复制 include 文件（.itp）
    local topol_dir=$(dirname "$TOPOL_FILE")
    local itp_count=0
    while IFS= read -r -d '' itp; do
        cp "$itp" "$OUTPUT_DIR/" 2>/dev/null || true
        ((itp_count++)) || true
    done < <(find "$topol_dir" -name "*.itp" -type f -print0 2>/dev/null)
    
    if [[ $itp_count -gt 0 ]]; then
        info "复制了 $itp_count 个 .itp 文件"
    fi
    
    # 切换到输出目录
    pushd "$OUTPUT_DIR" > /dev/null
    
    local current_gro="input.gro"
    local failed=false
    
    step "╔════════════════════════════════════════════════════════════╗"
    step "║          开始 GROMACS 模拟流程 (EM → NVT → NPT)            ║"
    step "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    # Step 1: 能量最小化
    if [[ "$SKIP_EM" != "true" ]]; then
        if run_gmx_step "能量最小化 (EM)" "${MDP_DIR}/em.mdp" "$current_gro" "topol.top" "em"; then
            current_gro="em.gro"
        else
            failed=true
        fi
    else
        warn "跳过能量最小化 (--skip-em)"
    fi
    
    # Step 2: NVT 平衡
    if [[ "$failed" != "true" ]] && [[ "$SKIP_NVT" != "true" ]]; then
        if run_gmx_step "NVT 平衡 (恒温)" "${MDP_DIR}/nvt.mdp" "$current_gro" "topol.top" "nvt"; then
            current_gro="nvt.gro"
        else
            failed=true
        fi
    elif [[ "$SKIP_NVT" == "true" ]]; then
        warn "跳过 NVT 平衡 (--skip-nvt)"
    fi
    
    # Step 3: NPT 平衡
    if [[ "$failed" != "true" ]]; then
        if run_gmx_step "NPT 平衡 (恒压)" "${MDP_DIR}/npt.mdp" "$current_gro" "topol.top" "npt"; then
            current_gro="npt.gro"
        else
            failed=true
        fi
    fi
    
    popd > /dev/null
    
    if [[ "$failed" == "true" ]]; then
        print_summary "$OUTPUT_DIR" "failure"
        exit 1
    else
        print_summary "$OUTPUT_DIR" "success"
    fi
}

main "$@"
