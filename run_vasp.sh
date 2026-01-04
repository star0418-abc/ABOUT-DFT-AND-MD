#!/usr/bin/env bash
# ============================================================================
# run_vasp.sh - VASP 运行脚本（含自动备份、日志、续算支持）
# ============================================================================
# 用法: NP=16 EXE=vasp_std run_vasp.sh
#
# 环境变量:
#   NP           - MPI 进程数 (默认 8)
#   EXE          - 可执行文件 (vasp_std/vasp_gam/vasp_ncl, 默认 vasp_std)
#   OUT          - stdout 文件名 (默认 vasp.out)
#   ERR          - stderr 文件名 (默认 vasp.err)
#   RESUME       - 续算模式: 1=自动 cp CONTCAR->POSCAR (默认 0)
#   STRICT_NP    - 严格核数: 1=NP 超限时报错退出 (默认 0=自动下调)
#   RESERVE_CORES - WSL 预留核数 (默认 2)
#   MIN_FREE_GB  - 最小磁盘空间 GB (默认 20)
#   FORCE_DISK   - 忽略磁盘检查: 1=强制继续 (默认 0)
# ============================================================================
set -uo pipefail  # 不用 -e，手动捕获 mpirun 返回码

# ---------------------- 加载 VASP 环境 ----------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if ! source "$SCRIPT_DIR/vasp_env.sh"; then
    echo "[ERROR] VASP 环境加载失败"
    exit 1
fi

# ---------------------- 参数设置 ----------------------
NP="${NP:-8}"
EXE="${EXE:-vasp_std}"
OUT="${OUT:-vasp.out}"
ERR="${ERR:-vasp.err}"
RESUME="${RESUME:-0}"
STRICT_NP="${STRICT_NP:-0}"
RESERVE_CORES="${RESERVE_CORES:-2}"
MIN_FREE_GB="${MIN_FREE_GB:-20}"
FORCE_DISK="${FORCE_DISK:-0}"
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"

# ---------------------- 函数：记录到 run.log ----------------------
log_msg() {
    echo "$1" | tee -a run.log
}

# ---------------------- 文件检查 ----------------------
echo "============================================"
echo "[run_vasp] 当前目录: $(pwd)"
echo "[run_vasp] NP=$NP  EXE=$EXE  OUT=$OUT  ERR=$ERR"
echo "[run_vasp] OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "[run_vasp] 时间戳: $TIMESTAMP"
echo "============================================"

# 必需文件检查
missing=0
for f in INCAR POSCAR POTCAR; do
    if [[ ! -f "$f" ]]; then
        echo "[ERROR] 缺少必需文件: $f"
        missing=1
    else
        echo "[OK] $f 存在"
    fi
done

if [[ $missing -eq 1 ]]; then
    echo "[ABORT] 请补齐必需输入文件后重试。"
    exit 1
fi

# KPOINTS 检查（非必需但建议）
if [[ ! -f "KPOINTS" ]]; then
    echo "[WARN] KPOINTS 文件不存在，请确保 INCAR 中设置了 KSPACING"
else
    echo "[OK] KPOINTS 存在"
fi

# 检查可执行文件
if [[ ! -x "$VASP_BIN/$EXE" ]]; then
    echo "[ERROR] 可执行文件不存在或无执行权限: $VASP_BIN/$EXE"
    exit 1
fi
echo "[OK] 可执行文件: $VASP_BIN/$EXE"

# ---------------------- WSL 核数检查 ----------------------
echo ""
echo ">>> WSL 核数检查..."
TOTAL_CORES=$(nproc 2>/dev/null || echo 8)
MAX_NP=$((TOTAL_CORES - RESERVE_CORES))
if [[ $MAX_NP -lt 1 ]]; then
    MAX_NP=1
fi

echo "    总核数: $TOTAL_CORES, 预留: $RESERVE_CORES, 可用: $MAX_NP"

if [[ $NP -gt $MAX_NP ]]; then
    if [[ $STRICT_NP -eq 1 ]]; then
        echo "[ERROR] NP=$NP 超过可用核数 $MAX_NP (STRICT_NP=1)"
        echo "[INFO] 设置 STRICT_NP=0 可自动下调"
        exit 1
    else
        echo "[WARN] NP=$NP 超过可用核数，自动下调为 $MAX_NP"
        NP=$MAX_NP
    fi
fi
echo "[OK] 使用 NP=$NP"

# ---------------------- 磁盘空间检查 ----------------------
echo ""
echo ">>> 磁盘空间检查..."
FREE_KB=$(df -k . 2>/dev/null | tail -1 | awk '{print $4}')
FREE_GB=$((FREE_KB / 1024 / 1024))

echo "    可用空间: ${FREE_GB} GB, 最小要求: ${MIN_FREE_GB} GB"

if [[ $FREE_GB -lt $MIN_FREE_GB ]]; then
    if [[ $FORCE_DISK -eq 1 ]]; then
        echo "[WARN] 磁盘空间不足，但 FORCE_DISK=1，继续运行"
    else
        echo "[ERROR] 磁盘空间不足 (${FREE_GB} < ${MIN_FREE_GB} GB)"
        echo "[INFO] 设置 FORCE_DISK=1 可强制继续"
        exit 1
    fi
else
    echo "[OK] 磁盘空间充足"
fi

# ---------------------- 并行参数提示 ----------------------
echo ""
echo ">>> 并行参数检查..."

# 解析 INCAR 中的 NCORE
NCORE=$(grep -i "^[[:space:]]*NCORE" INCAR 2>/dev/null | head -1 | sed 's/.*=//;s/[^0-9]//g' || echo "")
KPAR=$(grep -i "^[[:space:]]*KPAR" INCAR 2>/dev/null | head -1 | sed 's/.*=//;s/[^0-9]//g' || echo "")

if [[ -n "$NCORE" ]]; then
    if [[ $((NP % NCORE)) -ne 0 ]]; then
        echo "[WARN] NP=$NP 不能被 NCORE=$NCORE 整除，可能影响性能"
    else
        echo "    NCORE=$NCORE, NP/NCORE=$((NP / NCORE))"
    fi
fi

if [[ -n "$KPAR" ]]; then
    if [[ $((NP % KPAR)) -ne 0 ]]; then
        echo "[WARN] NP=$NP 不能被 KPAR=$KPAR 整除，可能影响性能"
    else
        echo "    KPAR=$KPAR"
    fi
fi

# ---------------------- AIMD 续算支持 ----------------------
echo ""
echo ">>> 续算检查..."

if [[ -f CONTCAR ]]; then
    contcar_size=$(stat -c%s CONTCAR 2>/dev/null || echo 0)
    if [[ $contcar_size -gt 100 ]]; then
        if [[ $RESUME -eq 1 ]]; then
            echo "[RESUME] 检测到 CONTCAR，启用续算模式"
            echo "[RESUME] 备份 POSCAR -> POSCAR.bak.$TIMESTAMP"
            cp POSCAR "POSCAR.bak.$TIMESTAMP"
            echo "[RESUME] 复制 CONTCAR -> POSCAR"
            cp CONTCAR POSCAR
            echo "[OK] 续算准备完成"
        else
            echo "[INFO] 检测到 CONTCAR，如需续算请使用: RESUME=1 run_vasp.sh"
        fi
    fi
else
    echo "[INFO] 无 CONTCAR，非续算"
fi

# ---------------------- 备份输入文件到 snapshots ----------------------
echo ""
echo ">>> 备份输入文件..."
SNAPSHOT_DIR="snapshots/$TIMESTAMP"
mkdir -p "$SNAPSHOT_DIR"

for f in INCAR POSCAR POTCAR KPOINTS; do
    if [[ -f "$f" ]]; then
        cp "$f" "$SNAPSHOT_DIR/"
        echo "    [备份] $f -> $SNAPSHOT_DIR/"
    fi
done
echo "[OK] 输入文件已备份到: $SNAPSHOT_DIR/"

# ---------------------- 移动旧输出文件到 old ----------------------
echo ""
echo ">>> 移动旧输出文件..."
OLD_FILES=(OUTCAR OSZICAR CONTCAR XDATCAR vasprun.xml WAVECAR CHGCAR CHG "$OUT" "$ERR")
moved_count=0

for f in "${OLD_FILES[@]}"; do
    if [[ -f "$f" ]]; then
        fsize=$(stat -c%s "$f" 2>/dev/null || echo 0)
        if [[ $fsize -gt 0 ]]; then
            mkdir -p old
            new_name="${f}.${TIMESTAMP}"
            mv "$f" "old/$new_name"
            echo "    [移动] $f -> old/$new_name"
            ((moved_count++)) || true
        fi
    fi
done

if [[ $moved_count -eq 0 ]]; then
    echo "    [INFO] 没有需要移动的旧输出文件"
else
    echo "[OK] 已移动 $moved_count 个旧输出文件到 old/"
fi

# ---------------------- 记录运行日志 ----------------------
echo ""
echo ">>> 记录运行日志..."
{
    echo "============================================"
    echo "时间: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "目录: $(pwd)"
    echo "NP: $NP"
    echo "EXE: $EXE"
    echo "OUT: $OUT"
    echo "ERR: $ERR"
    echo "RESUME: $RESUME"
    echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
    echo "快照: $SNAPSHOT_DIR"
    echo "============================================"
} >> run.log
echo "[OK] 运行参数已记录到 run.log"

# ---------------------- 运行 VASP ----------------------
echo ""
echo "============================================"
echo "[run_vasp] 开始运行: mpirun -np $NP $VASP_BIN/$EXE"
echo "[run_vasp] stdout: $OUT"
echo "[run_vasp] stderr: $ERR"
echo "[run_vasp] 开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================"
echo ""
echo "[INFO] 监控命令:"
echo "       tail -f $OUT"
echo "       tail -f $ERR"
echo "       aimd_watch.sh"
echo ""

# 记录开始时间
START_TIME=$(date +%s)
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')" >> run.log

# 运行 mpirun，stdout/stderr 分离，同时 tee 到屏幕
# 使用 pipefail 并捕获返回码
set +e
{
    mpirun -np "$NP" "$VASP_BIN/$EXE" 2>&1 \
        | tee "$OUT" \
        | grep -E "(E0|F=|DAV|RMM|Error|error|WARNING|STOP)" || true
} 2> >(tee "$ERR" >&2)

# 获取 mpirun 的返回码（通过 PIPESTATUS）
VASP_RC=${PIPESTATUS[0]}
set -e

END_TIME=$(date +%s)
WALL_TIME=$((END_TIME - START_TIME))

# 转换为 h:m:s
WALL_H=$((WALL_TIME / 3600))
WALL_M=$(((WALL_TIME % 3600) / 60))
WALL_S=$((WALL_TIME % 60))

echo ""
echo "[run_vasp] 结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "[run_vasp] 运行耗时: ${WALL_H}h ${WALL_M}m ${WALL_S}s (${WALL_TIME}s)"

# 记录到日志
{
    echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "运行耗时: ${WALL_H}h ${WALL_M}m ${WALL_S}s (${WALL_TIME}s)"
    echo "mpirun 返回码: $VASP_RC"
} >> run.log

# ---------------------- 错误处理 ----------------------
if [[ $VASP_RC -ne 0 ]]; then
    echo ""
    echo "[ERROR] mpirun 返回非零退出码: $VASP_RC"
    echo "状态: 失败 (rc=$VASP_RC)" >> run.log
    
    # 记录错误日志尾部
    if [[ -f "$ERR" ]]; then
        echo "" >> run.log
        echo "=== $ERR 最后 30 行 ===" >> run.log
        tail -30 "$ERR" >> run.log 2>/dev/null || true
    fi
    
    echo "[INFO] 检查 $ERR 获取详细错误信息"
    echo "[INFO] 检查 run.log 获取完整记录"
fi

# ---------------------- 完成检查 ----------------------
echo ""
if [[ -f OUTCAR ]]; then
    if grep -q "General timing and accounting informations for this job" OUTCAR; then
        echo "[OK] VASP 计算正常完成。"
        echo "状态: 正常完成" >> run.log
        
        # 提取总 CPU 时间
        total_cpu=$(grep "Total CPU time used" OUTCAR | tail -1 || true)
        if [[ -n "$total_cpu" ]]; then
            echo "    $total_cpu"
            echo "$total_cpu" >> run.log
        fi
    else
        echo "[WARN] OUTCAR 存在但未找到正常完成标志，请检查计算是否中断。"
        echo "状态: 未正常完成（可能中断）" >> run.log
    fi
else
    echo "[WARN] OUTCAR 不存在，计算可能失败。"
    echo "状态: 失败（无 OUTCAR）" >> run.log
fi

echo "" >> run.log
echo "[run_vasp] 完成。"
