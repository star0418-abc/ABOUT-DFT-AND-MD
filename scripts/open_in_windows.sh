#!/bin/bash

# Author: Star



# ==============================================================================
# open_in_windows.sh - 在 Windows 资源管理器中打开 WSL 目录
# ==============================================================================
# 用法:
#   ./scripts/open_in_windows.sh outputs/smoke
#   ./scripts/open_in_windows.sh .
#   ./scripts/open_in_windows.sh ~/gel_packmol/outputs
# ==============================================================================

set -euo pipefail

# 颜色
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

if [[ $# -lt 1 ]]; then
    echo -e "${YELLOW}用法:${NC} $0 <path>"
    echo ""
    echo "示例:"
    echo "  $0 outputs/smoke"
    echo "  $0 ."
    echo "  $0 ~/gel_packmol/outputs"
    exit 1
fi

TARGET_PATH="$1"

# 检查路径是否存在
if [[ ! -e "$TARGET_PATH" ]]; then
    echo -e "${RED}[ERROR]${NC} 路径不存在: $TARGET_PATH"
    exit 1
fi

# 转换为绝对路径
ABS_PATH=$(realpath "$TARGET_PATH")

# 检查 explorer.exe 是否可用
if ! command -v explorer.exe &> /dev/null; then
    echo -e "${RED}[ERROR]${NC} explorer.exe 不可用（可能不在 WSL 环境中）"
    echo "请手动打开: $ABS_PATH"
    exit 1
fi

# 转换为 Windows 路径
WIN_PATH=$(wslpath -w "$ABS_PATH" 2>/dev/null)

if [[ -z "$WIN_PATH" ]]; then
    echo -e "${RED}[ERROR]${NC} 无法转换路径: $ABS_PATH"
    exit 1
fi

echo -e "${GREEN}[INFO]${NC} 在 Windows 资源管理器中打开:"
echo "  WSL 路径: $ABS_PATH"
echo "  Win 路径: $WIN_PATH"

# 打开
explorer.exe "$WIN_PATH" 2>/dev/null &

echo -e "${GREEN}[INFO]${NC} 已发送打开命令"

