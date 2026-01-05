#!/bin/bash

# Author: Star



# ==============================================================================
# run_workflow.sh - 统一全流程入口脚本 (v1.0)
# ==============================================================================
# 功能: PACKMOL → HTPOLYNET → GROMACS → 配位数计算
#
# 特点:
#   - 唯一全流程入口脚本
#   - 强制使用 config/recipe.yaml 作为默认配置
#   - 所有输出集中到 outputs/<run_timestamp>/
#   - 自动清理，根目录保持干净
#
# 用法:
#   bash scripts/run_workflow.sh
#   bash scripts/run_workflow.sh -c config/my_recipe.yaml
#   bash scripts/run_workflow.sh --dry-run
#   bash scripts/run_workflow.sh --keep-big
#
# ==============================================================================

set -euo pipefail

# ==============================================================================
# 常量与配置
# ==============================================================================
VERSION="1.0.0"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
DEFAULT_RECIPE="${PROJECT_ROOT}/config/recipe.yaml"

# 运行参数
RECIPE_FILE="$DEFAULT_RECIPE"
RUN_DIR="${PROJECT_ROOT}/outputs/run_${TIMESTAMP}"
KEEP_BIG=false
DRY_RUN=false
FORCE_RERUN=false
SKIP_HTPOLYNET=false
SKIP_GROMACS=false
SKIP_CN=false
PYTHON="${PYTHON:-python3}"

# MDP 配置 (smoke 或 prod)
# 可通过环境变量 MDP_PROFILE 或 --mdp-profile 参数设置
MDP_PROFILE="${MDP_PROFILE:-smoke}"
MDP_DIR_SMOKE="${PROJECT_ROOT}/configs/mdp"
MDP_DIR_PROD="${PROJECT_ROOT}/config/mdp_prod"
MDP_DIR="$MDP_DIR_SMOKE"  # 默认使用 smoke

# 颜色
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'

# ==============================================================================
# 打印函数
# ==============================================================================
# 日志文件（初始化前使用 /dev/null）
LOG_FILE="${LOG_FILE:-/dev/null}"

info()    { echo -e "${GREEN}[INFO]${NC} $1" | tee -a "$LOG_FILE" 2>/dev/null || echo -e "${GREEN}[INFO]${NC} $1"; }
warn()    { echo -e "${YELLOW}[WARN]${NC} $1" | tee -a "$LOG_FILE" 2>/dev/null || echo -e "${YELLOW}[WARN]${NC} $1"; }
error()   { echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOG_FILE" 2>/dev/null || echo -e "${RED}[ERROR]${NC} $1" >&2; }
step()    { echo -e "${CYAN}[STEP]${NC} $1" | tee -a "$LOG_FILE" 2>/dev/null || echo -e "${CYAN}[STEP]${NC} $1"; }
success() { echo -e "${MAGENTA}[✓]${NC} $1" | tee -a "$LOG_FILE" 2>/dev/null || echo -e "${MAGENTA}[✓]${NC} $1"; }
fail()    { echo -e "${RED}[✗]${NC} $1" | tee -a "$LOG_FILE" 2>/dev/null || echo -e "${RED}[✗]${NC} $1"; }

print_banner() {
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════════╗"
    echo "║           gel_packmol 统一工作流 v${VERSION}                            ║"
    echo "║     PACKMOL → HTPOLYNET → GROMACS → Coordination Number              ║"
    echo "╚══════════════════════════════════════════════════════════════════════╝"
    echo ""
}

usage() {
    cat << EOF
用法: $0 [选项]

选项:
  -c, --config <yaml>   配方文件 (默认: config/recipe.yaml)
  -o, --output <dir>    输出目录 (默认: outputs/run_<timestamp>)
  --mdp-profile <name>  MDP 配置集 (smoke|prod，默认: smoke)
                        smoke: configs/mdp/      (快速测试)
                        prod:  config/mdp_prod/  (正式模拟)
  --force               强制重跑所有步骤
  --skip-htpolynet      跳过 HTPolyNet 阶段
  --skip-gromacs        跳过 GROMACS 阶段
  --skip-cn             跳过配位数计算
  --keep-big            保留大文件 (traj.xtc 等)
  --dry-run             模拟运行
  -h, --help            显示帮助

示例:
  $0                                    # 使用默认配置 (smoke)
  $0 --mdp-profile prod                 # 使用生产 MDP
  $0 -c config/recipe.yaml              # 显式指定配置
  $0 --skip-htpolynet                   # 跳过 HTPolyNet
  $0 --force                            # 强制重跑

环境变量:
  MDP_PROFILE=prod      # 等同于 --mdp-profile prod

MDP 配置说明:
  smoke (默认):
    - 位置: configs/mdp/
    - 用途: 流程验证、快速测试
    - 时长: 100ps ~ 1ns
  
  prod (正式):
    - 位置: config/mdp_prod/
    - 用途: 正式模拟、统计采样
    - 时长: 1ns ~ 100ns
    - 流程: em → nvt → npt_ber → npt_pr → md_prod

输出目录: outputs/run_<timestamp>/
  ├── config_used.yaml       # 使用的配方副本
  ├── logs/run.log           # 完整日志
  ├── packmol/               # Packmol 输出
  ├── htpolynet/             # HTPolyNet 输出
  ├── gromacs/               # GROMACS 输出
  └── analysis/coordination/ # 配位数结果

注意:
  - 默认配置: config/recipe.yaml (唯一权威配置)
  - 旧脚本 run_all_smoke.sh / run_full.sh 已弃用
EOF
    exit 0
}

# ==============================================================================
# 参数解析
# ==============================================================================
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -c|--config)
                RECIPE_FILE="$2"
                shift 2
                ;;
            -o|--output)
                RUN_DIR="$2"
                shift 2
                ;;
            --mdp-profile)
                MDP_PROFILE="$2"
                shift 2
                ;;
            --force)
                FORCE_RERUN=true
                shift
                ;;
            --skip-htpolynet)
                SKIP_HTPOLYNET=true
                shift
                ;;
            --skip-gromacs)
                SKIP_GROMACS=true
                shift
                ;;
            --skip-cn)
                SKIP_CN=true
                shift
                ;;
            --keep-big)
                KEEP_BIG=true
                shift
                ;;
            --dry-run)
                DRY_RUN=true
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
    
    # 设置 MDP 目录
    case "$MDP_PROFILE" in
        smoke)
            MDP_DIR="$MDP_DIR_SMOKE"
            ;;
        prod|production)
            MDP_DIR="$MDP_DIR_PROD"
            ;;
        *)
            warn "未知 MDP_PROFILE: $MDP_PROFILE，使用默认 smoke"
            MDP_PROFILE="smoke"
            MDP_DIR="$MDP_DIR_SMOKE"
            ;;
    esac
}

# ==============================================================================
# 依赖检查
# ==============================================================================
check_dependencies() {
    step "检查依赖..."
    local missing=()
    
    # Python
    if ! command -v "$PYTHON" &>/dev/null; then
        missing+=("python3: sudo apt install python3")
    fi
    
    # Packmol
    if ! command -v packmol &>/dev/null; then
        missing+=("packmol: conda install -c conda-forge packmol")
    fi
    
    # GROMACS
    if ! command -v gmx &>/dev/null; then
        missing+=("gmx (GROMACS): 请安装 GROMACS 2020+")
    fi
    
    # RDKit
    if ! "$PYTHON" -c "from rdkit import Chem" 2>/dev/null; then
        missing+=("RDKit: conda install -c conda-forge rdkit")
    fi
    
    # HTPolyNet (可选)
    if [[ "$SKIP_HTPOLYNET" != "true" ]]; then
        if ! command -v htpolynet &>/dev/null; then
            warn "HTPolyNet 未安装，将使用简化网络"
            warn "安装: pip install htpolynet"
        fi
    fi
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        error "缺少必要依赖:"
        for dep in "${missing[@]}"; do
            echo "  - $dep"
        done
        exit 1
    fi
    
    success "所有依赖已就绪"
}

# ==============================================================================
# 验证配方文件
# ==============================================================================
validate_recipe() {
    step "验证配方文件..."
    
    if [[ ! -f "$RECIPE_FILE" ]]; then
        error "配方文件不存在: $RECIPE_FILE"
        error "请确保 config/recipe.yaml 存在"
        exit 1
    fi
    
    info "使用配方: $RECIPE_FILE"
    
    # 检查必要字段
    local required_fields=("molecules" "system" "packmol")
    for field in "${required_fields[@]}"; do
        if ! grep -q "^${field}:" "$RECIPE_FILE"; then
            error "配方缺少必要字段: $field"
            exit 1
        fi
    done
    
    # 检查配方组分（至少有一个）
    local has_components=false
    for comp in "salt_solution" "polymer_matrix" "crosslinker" "photoinitiator"; do
        if grep -q "^${comp}:" "$RECIPE_FILE"; then
            has_components=true
            break
        fi
    done
    
    if [[ "$has_components" != "true" ]]; then
        warn "配方未定义任何组分 (salt_solution/polymer_matrix/crosslinker/photoinitiator)"
    fi
    
    success "配方验证通过"
}

# ==============================================================================
# 初始化运行目录
# ==============================================================================
init_run_dir() {
    step "初始化运行目录: $RUN_DIR"
    
    mkdir -p "$RUN_DIR"/{logs,packmol,htpolynet,gromacs,analysis/coordination}
    
    # 创建日志文件
    LOG_FILE="${RUN_DIR}/logs/run.log"
    touch "$LOG_FILE"
    
    # 复制配方
    cp "$RECIPE_FILE" "${RUN_DIR}/config_used.yaml"
    
    # 记录运行信息
    cat >> "$LOG_FILE" << EOF
================================================================================
gel_packmol Workflow Run
================================================================================
开始时间: $(date)
配方文件: $RECIPE_FILE
运行目录: $RUN_DIR
版本: $VERSION
================================================================================

EOF
    
    success "运行目录已初始化"
}

# ==============================================================================
# Phase 1: Packmol
# ==============================================================================
run_packmol_phase() {
    step "=== Phase 1: Packmol 分子打包 ==="
    
    local packmol_dir="${RUN_DIR}/packmol"
    local gel_pdb="${packmol_dir}/gel.pdb"
    
    # 检查是否需要跳过
    if [[ -f "$gel_pdb" ]] && [[ "$FORCE_RERUN" != "true" ]]; then
        info "Packmol 输出已存在，跳过"
        success "Phase 1 跳过（已有输出）"
        return 0
    fi
    
    cd "$PROJECT_ROOT"
    
    # 检查是否有 salt_solution 格式
    if grep -q "^salt_solution:" "$RECIPE_FILE"; then
        # 新格式：需要先换算
        info "检测到结构化配方，运行换算..."
        
        if "$PYTHON" "${SCRIPT_DIR}/recipe_to_counts.py" \
            -c "$RECIPE_FILE" \
            -o "$packmol_dir" 2>&1 | tee -a "$LOG_FILE"; then
            info "配方换算完成"
            local resolved_yaml="${packmol_dir}/recipe_resolved.yaml"
        else
            fail "配方换算失败"
            return 1
        fi
    else
        # 旧格式：直接使用
        local resolved_yaml="$RECIPE_FILE"
    fi
    
    # 运行 Packmol
    info "运行 make_packmol_from_recipe.py..."
    
    if "$PYTHON" "${SCRIPT_DIR}/make_packmol_from_recipe.py" \
        -c "$resolved_yaml" \
        -o "$packmol_dir" 2>&1 | tee -a "$LOG_FILE"; then
        
        if [[ -f "$gel_pdb" ]]; then
            local atom_count=$(grep -c "^ATOM\|^HETATM" "$gel_pdb" 2>/dev/null || echo "0")
            success "Phase 1 完成: Packmol 生成 $atom_count 原子"
        else
            fail "Packmol 未生成 gel.pdb"
            return 1
        fi
    else
        fail "Packmol 运行失败"
        return 1
    fi
}

# ==============================================================================
# Phase 2: HTPolyNet
# ==============================================================================
run_htpolynet_phase() {
    step "=== Phase 2: HTPolyNet 交联网络 ==="
    
    if [[ "$SKIP_HTPOLYNET" == "true" ]]; then
        info "跳过 HTPolyNet (--skip-htpolynet)"
        return 0
    fi
    
    local htpolynet_dir="${RUN_DIR}/htpolynet"
    local network_pdb="${htpolynet_dir}/network.pdb"
    
    # 检查是否需要跳过
    if [[ -f "$network_pdb" ]] && [[ "$FORCE_RERUN" != "true" ]]; then
        info "HTPolyNet 输出已存在，跳过"
        success "Phase 2 跳过（已有输出）"
        return 0
    fi
    
    cd "$PROJECT_ROOT"
    
    # 从 recipe.yaml 读取 HTPolyNet 参数
    local n_monomers
    local conversion
    local seed
    n_monomers=$(grep -A10 "^htpolynet:" "$RECIPE_FILE" | grep "n_monomers:" | head -1 | awk '{print $2}' || true)
    conversion=$(grep -A20 "^htpolynet:" "$RECIPE_FILE" | grep "desired_conversion:" | head -1 | awk '{print $2}' || true)
    seed=$(grep -A20 "^htpolynet:" "$RECIPE_FILE" | grep "random_seed:" | head -1 | awk '{print $2}' || true)
    
    # 使用较小的值进行快速演示
    n_monomers="${n_monomers:-50}"
    conversion="${conversion:-0.60}"
    seed="${seed:-2025}"
    
    info "HTPolyNet 参数: n_monomers=$n_monomers, conversion=$conversion, seed=$seed"
    
    # 生成 EGDA active form
    info "生成 EGDA active form..."
    
    local egda_mol
    local htpolynet_workdir
    local active_mol2="${htpolynet_dir}/EGDA_active.mol2"
    local pegda_mol2="${htpolynet_dir}/PEGDA.mol2"

    egda_mol=$(grep -A10 "^htpolynet:" "$RECIPE_FILE" | grep "source:" | head -1 | awk '{print $2}' | sed 's/"//g' || true)
    htpolynet_workdir=$(grep -A10 "^htpolynet:" "$RECIPE_FILE" | grep "workdir:" | head -1 | awk '{print $2}' | sed 's/"//g' || true)
    egda_mol="${egda_mol:-mol/EGDA.mol}"
    htpolynet_workdir="${htpolynet_workdir:-htpolynet_out/PEGDA}"
    egda_mol="${PROJECT_ROOT}/${egda_mol}"
    htpolynet_workdir="${PROJECT_ROOT}/${htpolynet_workdir}"
    
    if [[ ! -f "$egda_mol" ]]; then
        error "EGDA 单体文件不存在: $egda_mol"
        fail "Phase 2 失败"
        return 1
    fi
    
    # 使用 Python 生成 active form
    "$PYTHON" << PYTHON_SCRIPT 2>&1 | tee -a "$LOG_FILE"
import sys, os
sys.path.insert(0, '${SCRIPT_DIR}')

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule, UFFHasAllMoleculeParams
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

egda_mol_path = '${egda_mol}'
active_mol2_path = '${active_mol2}'
os.makedirs(os.path.dirname(active_mol2_path), exist_ok=True)

mol = Chem.MolFromMolFile(egda_mol_path, removeHs=False)
if mol is None:
    print('ERROR: 无法读取 EGDA.mol')
    sys.exit(1)

acrylate_pattern = Chem.MolFromSmarts('[CH2:1]=[CH:2]-[C:3](=[O:4])-[O:5]')
matches = mol.GetSubstructMatches(acrylate_pattern)

if len(matches) < 2:
    print(f'ERROR: 只找到 {len(matches)} 个丙烯酸酯端基')
    sys.exit(1)

rw_mol = Chem.RWMol(mol)
reactive_atoms = {}

for i, match in enumerate(matches[:2]):
    ch2_idx, ch_idx = match[0], match[1]
    bond = rw_mol.GetBondBetweenAtoms(ch2_idx, ch_idx)
    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
        bond.SetBondType(Chem.BondType.SINGLE)
    if i == 0:
        reactive_atoms['HA'] = ch2_idx
        reactive_atoms['TA'] = ch_idx
    else:
        reactive_atoms['HB'] = ch2_idx
        reactive_atoms['TB'] = ch_idx

mol_no_h = Chem.RemoveHs(rw_mol)
mol_with_h = Chem.AddHs(mol_no_h, addCoords=True)
if mol_with_h.GetNumConformers() == 0:
    AllChem.EmbedMolecule(mol_with_h, randomSeed=2025)
if UFFHasAllMoleculeParams(mol_with_h):
    UFFOptimizeMolecule(mol_with_h, maxIters=200)

# 写入 MOL2
conf = mol_with_h.GetConformer()
atoms = list(mol_with_h.GetAtoms())
bonds = list(mol_with_h.GetBonds())

reactive_idx_to_name = {v: k for k, v in reactive_atoms.items()}
atom_names = {}
element_count = {}

for atom in atoms:
    idx = atom.GetIdx()
    if idx in reactive_idx_to_name:
        atom_names[idx] = reactive_idx_to_name[idx]
    else:
        symbol = atom.GetSymbol()
        element_count[symbol] = element_count.get(symbol, 0) + 1
        atom_names[idx] = f'{symbol}{element_count[symbol]}'

with open(active_mol2_path, 'w') as f:
    f.write('@<TRIPOS>MOLECULE\nEGDA\n')
    f.write(f' {len(atoms)} {len(bonds)} 1 0 0\nSMALL\nNO_CHARGES\n\n')
    f.write('@<TRIPOS>ATOM\n')
    for atom in atoms:
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        symbol = atom.GetSymbol()
        sybyl = 'H' if symbol == 'H' else f'{symbol}.3'
        f.write(f'{idx+1:7d} {atom_names[idx]:<4s} {pos.x:10.4f} {pos.y:10.4f} {pos.z:10.4f} {sybyl:<6s} 1 EGD  0.0000\n')
    f.write('@<TRIPOS>BOND\n')
    for i, bond in enumerate(bonds):
        f.write(f'{i+1:6d} {bond.GetBeginAtomIdx()+1:5d} {bond.GetEndAtomIdx()+1:5d} 1\n')
    f.write('@<TRIPOS>SUBSTRUCTURE\n     1 EGD         1 RESIDUE    0 **** **** 0 ROOT\n')

print(f'Active form 已生成: {active_mol2_path}')
PYTHON_SCRIPT

    info "运行 HTPolyNet 主流程..."
    if "$PYTHON" "${SCRIPT_DIR}/build_pegda_network.py" \
        --in_mol "${egda_mol#${PROJECT_ROOT}/}" \
        --out_mol2 "${pegda_mol2#${PROJECT_ROOT}/}" \
        --workdir "${htpolynet_workdir#${PROJECT_ROOT}/}" \
        --n_monomers "$n_monomers" \
        --conversion "$conversion" \
        --seed "$seed" 2>&1 | tee -a "$LOG_FILE"; then
        info "HTPolyNet 已运行"
    else
        warn "HTPolyNet 运行失败，尝试回退到 active form"
    fi

    local htpolynet_pdb=""
    htpolynet_pdb=$(find "$htpolynet_workdir" -type f -name "*.pdb" 2>/dev/null | head -1 || true)
    if [[ -n "$htpolynet_pdb" ]]; then
        cp "$htpolynet_pdb" "$network_pdb"
        success "Phase 2 完成: 交联网络已生成"
        return 0
    fi

    if [[ -f "$pegda_mol2" ]]; then
        "$PYTHON" - << PYTHON_SCRIPT 2>&1 | tee -a "$LOG_FILE"
from rdkit import Chem
mol = Chem.MolFromMol2File("${pegda_mol2}", removeHs=False)
if mol:
    Chem.MolToPDBFile(mol, "${network_pdb}")
    print("network.pdb 已生成: ${network_pdb}")
PYTHON_SCRIPT
    fi

    if [[ -f "$network_pdb" ]]; then
        success "Phase 2 完成: 交联网络已生成"
    else
        warn "HTPolyNet 输出不可用，生成简化网络"
        if [[ -f "$active_mol2" ]]; then
            "$PYTHON" -c "
from rdkit import Chem
mol = Chem.MolFromMol2File('$active_mol2', removeHs=False)
if mol: Chem.MolToPDBFile(mol, '$network_pdb')
" 2>/dev/null || true
        fi
    fi
}

# ==============================================================================
# Phase 3: GROMACS
# ==============================================================================
run_gromacs_phase() {
    step "=== Phase 3: GROMACS 分子动力学 ==="
    
    if [[ "$SKIP_GROMACS" == "true" ]]; then
        info "跳过 GROMACS (--skip-gromacs)"
        return 0
    fi
    
    local gmx_dir="${RUN_DIR}/gromacs"
    local packmol_pdb="${RUN_DIR}/packmol/gel.pdb"
    local network_pdb="${RUN_DIR}/htpolynet/network.pdb"
    local merged_pdb="${gmx_dir}/gel_merged.pdb"
    local input_pdb=""
    local topology_input="${gmx_dir}/gel.pdb"

    mkdir -p "$gmx_dir"
    cd "$gmx_dir"

    # 选择输入结构（若两者存在则合并）
    if [[ -f "$packmol_pdb" ]] && [[ -f "$network_pdb" ]]; then
        info "检测到 Packmol + HTPolyNet，合并结构..."
        "$PYTHON" << PYTHON_SCRIPT 2>&1 | tee -a "$LOG_FILE"
import itertools
from pathlib import Path

packmol = Path("${packmol_pdb}")
network = Path("${network_pdb}")
output = Path("${merged_pdb}")

def read_atoms(path):
    lines = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                lines.append(line.rstrip("\n"))
    return lines

def parse_coords(line):
    return float(line[30:38]), float(line[38:46]), float(line[46:54])

packmol_atoms = read_atoms(packmol)
network_atoms = read_atoms(network)
atoms = packmol_atoms + network_atoms
if not atoms:
    raise SystemExit("ERROR: 未找到可合并的 ATOM/HETATM 行")

if packmol_atoms and network_atoms:
    coords = [parse_coords(line) for line in packmol_atoms]
    min_x = min(c[0] for c in coords)
    max_x = max(c[0] for c in coords)
    padding = 5.0
    dx = (max_x - min_x) + padding

    shifted = []
    for line in network_atoms:
        x, y, z = parse_coords(line)
        new_x = x + dx
        shifted.append(f"{line[:30]}{new_x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
    network_atoms = shifted

atoms = packmol_atoms + network_atoms

renumbered = []
for idx, line in enumerate(atoms, start=1):
    renumbered.append(f"{line[:6]}{idx:5d}{line[11:]}")

output.parent.mkdir(parents=True, exist_ok=True)
with output.open("w") as fh:
    for line in renumbered:
        fh.write(line + "\n")
    fh.write("END\n")

print(f"已生成合并结构: {output}")
PYTHON_SCRIPT
        input_pdb="$merged_pdb"
    elif [[ -f "$packmol_pdb" ]]; then
        input_pdb="$packmol_pdb"
    elif [[ -f "$network_pdb" ]]; then
        input_pdb="$network_pdb"
    else
        warn "没有可用的输入结构，跳过 GROMACS"
        return 0
    fi

    info "使用输入结构: $input_pdb"
    cp "$input_pdb" "$topology_input"

    info "生成 GAFF2 拓扑..."
    if ! "$PYTHON" "${SCRIPT_DIR}/generate_topology.py" \
        -i "$gmx_dir" \
        -c "$RECIPE_FILE" \
        -o "$gmx_dir" 2>&1 | tee -a "$LOG_FILE"; then
        warn "拓扑生成失败，跳过 GROMACS"
        return 0
    fi

    if [[ ! -f "${gmx_dir}/topol.top" ]]; then
        warn "未找到 topol.top，跳过 GROMACS"
        return 0
    fi

    local current_gro="${gmx_dir}/system.gro"
    if [[ ! -f "$current_gro" ]]; then
        warn "未找到 system.gro，跳过 GROMACS"
        return 0
    fi

    local run_stage
    run_stage() {
        local stage_name="$1"
        local mdp_file="$2"

        if [[ ! -f "$mdp_file" ]]; then
            warn "缺少 MDP: $mdp_file，跳过 ${stage_name}"
            return 0
        fi

        info "运行 ${stage_name}..."
        if ! gmx grompp -f "$mdp_file" -c "$current_gro" -p topol.top -o "${stage_name}.tpr" -maxwarn 10 2>&1 | tee -a "$LOG_FILE"; then
            warn "${stage_name} grompp 失败"
            return 1
        fi
        if ! gmx mdrun -v -deffnm "${stage_name}" 2>&1 | tee -a "$LOG_FILE"; then
            warn "${stage_name} mdrun 失败"
            return 1
        fi

        if [[ -f "${stage_name}.gro" ]]; then
            current_gro="${gmx_dir}/${stage_name}.gro"
        fi
        return 0
    }

    local npt_stage="npt.mdp"
    if [[ -f "${MDP_DIR}/npt_ber.mdp" ]]; then
        npt_stage="npt_ber.mdp"
    elif [[ -f "${MDP_DIR}/npt.mdp" ]]; then
        npt_stage="npt.mdp"
    fi

    run_stage "em" "${MDP_DIR}/em.mdp" || return 0
    run_stage "nvt" "${MDP_DIR}/nvt.mdp" || return 0
    run_stage "npt" "${MDP_DIR}/${npt_stage}" || return 0

    if [[ -f "${MDP_DIR}/npt_pr.mdp" ]]; then
        run_stage "npt_pr" "${MDP_DIR}/npt_pr.mdp" || return 0
    fi

    if [[ -f "${MDP_DIR}/md_prod.mdp" ]]; then
        run_stage "md_prod" "${MDP_DIR}/md_prod.mdp" || return 0
    fi

    if [[ -f "$current_gro" ]]; then
        cp "$current_gro" "${RUN_DIR}/final_structure.gro"
    fi

    success "Phase 3 完成"
}

# ==============================================================================
# Phase 4: 配位数计算
# ==============================================================================
run_coordination_phase() {
    step "=== Phase 4: 配位数(CN)计算 ==="
    
    if [[ "$SKIP_CN" == "true" ]]; then
        info "跳过配位数计算 (--skip-cn)"
        return 0
    fi
    
    local analysis_dir="${RUN_DIR}/analysis/coordination"
    local gmx_dir="${RUN_DIR}/gromacs"
    local tpr_file=""
    local traj_file=""

    cd "$analysis_dir"

    tpr_file=$(ls "${gmx_dir}"/*.tpr 2>/dev/null | tail -1 || true)
    traj_file=$(ls "${gmx_dir}"/*.xtc 2>/dev/null | tail -1 || true)

    if [[ -z "$tpr_file" || -z "$traj_file" ]]; then
        warn "缺少轨迹或 tpr，跳过配位数计算"
        return 0
    fi

    if [[ ! -f "${gmx_dir}/index.ndx" ]]; then
        warn "未找到 index.ndx，尝试自动生成..."
        if ! command -v gmx &>/dev/null; then
            warn "gmx 不可用，无法生成 index.ndx"
            return 0
        fi
        if echo -e "r LI\na O*\nq\n" | gmx make_ndx -f "$tpr_file" -o "${gmx_dir}/index.ndx" 2>&1 | tee -a "$LOG_FILE"; then
            info "已生成 index.ndx"
        else
            warn "自动生成 index.ndx 失败，请手动创建包含 LI 与 O 组的索引文件"
            return 0
        fi
    fi

    info "使用 gmx rdf 计算配位数..."
    if echo -e "LI\nO*\n" | gmx rdf -f "$traj_file" -s "$tpr_file" -n "${gmx_dir}/index.ndx" -o rdf_Li_O.xvg -cn cn_Li_O.xvg 2>&1 | tee -a "$LOG_FILE"; then
        local cn_value
        cn_value=$(awk '($1>=0.35){print $2; exit}' cn_Li_O.xvg 2>/dev/null || true)
        cat > coordination_summary.txt << EOF
# ============================================
# 配位数分析结果
# ============================================

system: $(basename "$RECIPE_FILE" .yaml)
timestamp: $(date)
run_dir: ${RUN_DIR}
recipe: ${RECIPE_FILE}

# Li-O 配位数
CN(Li-O) at r=0.35 nm: ${cn_value:-N/A}

# 说明
note: 使用 gmx rdf -cn 计算

# 文件
files:
  - rdf_Li_O.xvg
  - cn_Li_O.xvg
EOF
        success "Phase 4 完成: CN(Li-O)@0.35nm = ${cn_value:-N/A}"
    else
        warn "gmx rdf 失败，跳过配位数结果输出"
        return 0
    fi
}

# ==============================================================================
# Phase 5: 清理
# ==============================================================================
cleanup() {
    step "=== Phase 5: 清理 ==="
    
    cd "$PROJECT_ROOT"
    
    # 移动根目录临时文件
    if [[ -d "${PROJECT_ROOT}/htpolynet_out" ]]; then
        warn "移动 htpolynet_out/ 到运行目录"
        mv "${PROJECT_ROOT}/htpolynet_out" "${RUN_DIR}/_artifacts_htpolynet_out" 2>/dev/null || true
    fi
    
    # 检查散落文件
    local stray_files=$(find "${PROJECT_ROOT}" -maxdepth 1 -name "*.pdb" -o -name "*.gro" -o -name "sqm.*" 2>/dev/null | head -5)
    if [[ -n "$stray_files" ]]; then
        warn "发现散落文件（建议手动清理）:"
        echo "$stray_files" | head -5
    fi
    
    # 大文件处理
    if [[ "$KEEP_BIG" == "false" ]]; then
        info "清理大文件..."
        find "$RUN_DIR" -name "*.xtc" -delete 2>/dev/null || true
        find "$RUN_DIR" -name "*.trr" -delete 2>/dev/null || true
        find "$RUN_DIR" -name "*.edr" -size +1M -delete 2>/dev/null || true
        find "$RUN_DIR" -name "*.cpt" -delete 2>/dev/null || true
    fi
    
    success "清理完成"
}

# ==============================================================================
# 最终报告
# ==============================================================================
print_summary() {
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════════╗"
    echo "║                      Workflow 完成!                                  ║"
    echo "╚══════════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "运行目录: $RUN_DIR"
    echo "使用配方: $RECIPE_FILE"
    echo ""
    echo "关键输出:"
    echo "  配置副本:    ${RUN_DIR}/config_used.yaml"
    echo "  运行日志:    ${RUN_DIR}/logs/run.log"
    echo "  最终结构:    ${RUN_DIR}/final_structure.gro"
    echo "  配位数结果:  ${RUN_DIR}/analysis/coordination/coordination_summary.txt"
    echo ""
    
    if [[ -f "${RUN_DIR}/analysis/coordination/coordination_summary.txt" ]]; then
        echo "配位数:"
        grep "CN(Li-O)" "${RUN_DIR}/analysis/coordination/coordination_summary.txt" 2>/dev/null || true
    fi
    echo ""
    
    # 记录结束时间
    echo "结束时间: $(date)" >> "$LOG_FILE"
}

# ==============================================================================
# 主流程
# ==============================================================================
main() {
    parse_args "$@"
    
    print_banner
    
    # 检查配方
    if [[ "$RECIPE_FILE" == "$DEFAULT_RECIPE" ]]; then
        info "使用默认配方: config/recipe.yaml"
    else
        info "使用指定配方: $RECIPE_FILE"
    fi
    
    # 显示 MDP 配置
    info "MDP 配置: $MDP_PROFILE ($MDP_DIR)"
    if [[ ! -d "$MDP_DIR" ]]; then
        error "MDP 目录不存在: $MDP_DIR"
        exit 1
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        info "[DRY-RUN] 模拟运行"
        info "配方: $RECIPE_FILE"
        info "MDP: $MDP_PROFILE ($MDP_DIR)"
        info "输出: $RUN_DIR"
        exit 0
    fi
    
    check_dependencies
    validate_recipe
    init_run_dir
    
    # 运行各阶段
    run_packmol_phase || warn "Packmol 阶段有问题"
    run_htpolynet_phase || warn "HTPolyNet 阶段有问题"
    run_gromacs_phase || warn "GROMACS 阶段有问题"
    run_coordination_phase || warn "配位数计算有问题"
    cleanup
    
    print_summary
}

main "$@"
