#!/bin/bash

# Author: Star



# ==============================================================================
# [DEPRECATED] run_full.sh - 已弃用
# ==============================================================================
# 
# ⚠️ 此脚本已弃用！请使用以下替代命令:
#
#   bash scripts/run_workflow.sh --recipe config/recipe.yaml
#
# 更多信息请参见 MIGRATION.md
# ==============================================================================

echo ""
echo "============================================================"
echo "⚠️ [DEPRECATED] run_full.sh 已弃用"
echo "============================================================"
echo ""
echo "请使用新的工作流脚本:"
echo ""
echo "  bash scripts/run_workflow.sh --recipe config/recipe.yaml"
echo ""
echo "============================================================"
echo ""

# 自动调用新脚本
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "${SCRIPT_DIR}/run_workflow.sh" "$@"
