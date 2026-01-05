# 脚本迁移指南

本文档说明 `scripts/` 重构后的新旧脚本映射关系。

## 📋 迁移映射表

| 旧脚本 | 新入口 | 状态 | 等价命令 |
|--------|--------|------|----------|
| `mol2_to_polymer_mol2.py` | `cli/monomers.py` | 保留包装器 | `python -m scripts.cli.monomers --input FILE` |
| `process_monomers_batch.py` | `cli/monomers.py` | 保留包装器 | `python -m scripts.cli.monomers --names EGDA MMA` |
| `recipe_to_counts.py` | 保留原位 | 使用 lib/recipe | 无变化 |
| `make_packmol_from_recipe.py` | 保留原位 | 使用 lib/recipe | 无变化 |
| `build_pegda_network.py` | 保留原位 | 使用 lib/io_mol | 无变化 |
| `gaussian_log_to_pdb.py` | 保留原位 | 独立功能 | 无变化 |
| `generate_topology.py` | 保留原位 | 独立功能 | 无变化 |
| `generate_topology_gaff.py` | 保留原位 | 独立功能 | 无变化 |
| `inject_charges.py` | 保留原位 | 独立功能 | 无变化 |
| `sanitize_pdb_for_ligpargen.py` | 保留原位 | 独立功能 | 无变化 |
| `mol_to_mol2_batch.py` | 保留原位 | 使用 lib/io_mol | 无变化 |
| `run_workflow.sh` | 保留原位 | 统一入口 | 无变化 |
| `run_gmx.sh` | 保留原位 | GROMACS 驱动 | 无变化 |
| `run_htpolynet.sh` | 保留原位 | HTPolyNet 驱动 | 无变化 |
| `run_packmol.sh` | 保留原位 | Packmol 驱动 | 无变化 |
| `run_all_smoke.sh` | **已弃用** | 打印警告后退出 | 使用 `selfcheck.py` |
| `run_full.sh` | **已弃用** | 打印警告后退出 | 使用 `run_workflow.sh` |
| `selfcheck_oligomer.py` | `selfcheck.py` | 合并 | `python scripts/selfcheck.py` |

## 🆕 新增模块

### lib/ - 公共库

| 模块 | 功能 | 用法 |
|------|------|------|
| `lib/recipe.py` | 统一 recipe.yaml 解析 | `from scripts.lib.recipe import load_recipe` |
| `lib/io_mol.py` | 分子文件读写与校验 | `from scripts.lib.io_mol import read_mol2, write_mol2_strict` |
| `lib/naming.py` | 统一命名规则 | `from scripts.lib.naming import get_polymer_output_path` |
| `lib/validate.py` | 文件/物理量校验 | `from scripts.lib.validate import validate_mol2_output` |
| `lib/logging_utils.py` | 统一日志格式 | `from scripts.lib.logging_utils import log_success` |

### cli/ - 命令行入口

| 入口 | 功能 | 用法 |
|------|------|------|
| `cli/monomers.py` | 单体批处理 | `python -m scripts.cli.monomers` |

## 📖 常见用例迁移

### 1. 批量生成寡聚体

**旧命令**:
```bash
python scripts/process_monomers_batch.py --names EGDA MMA TFEMA VA --n_repeat 3
```

**新命令**:
```bash
python -m scripts.cli.monomers --names EGDA MMA TFEMA VA --n_repeat 3
# 或使用包装器（仍然可用）
python scripts/process_monomers_batch.py --names EGDA MMA TFEMA VA --n_repeat 3
```

### 2. 单文件处理

**旧命令**:
```bash
python scripts/mol2_to_polymer_mol2.py --input mol2/EGDA.mol2 --n_repeat 3
```

**新命令**:
```bash
python -m scripts.cli.monomers --input mol2/EGDA.mol2 --n_repeat 3
# 或使用包装器（仍然可用）
python scripts/mol2_to_polymer_mol2.py --input mol2/EGDA.mol2 --n_repeat 3
```

### 3. 自检测试

**旧命令**:
```bash
python scripts/selfcheck_oligomer.py
```

**新命令**:
```bash
python scripts/selfcheck.py
```

## ⚠️ 已弃用脚本

以下脚本已弃用，运行时会打印警告并退出：

- `run_all_smoke.sh` → 使用 `python scripts/selfcheck.py`
- `run_full.sh` → 使用 `bash scripts/run_workflow.sh`

## 🔧 开发者指南

### 添加新脚本

1. 涉及 recipe 读取：使用 `from scripts.lib.recipe import load_recipe`
2. 涉及 mol2 写入：使用 `from scripts.lib.io_mol import write_mol2_strict`
3. 涉及输出验证：使用 `from scripts.lib.validate import validate_mol2_output`

### 命名规则

- 聚合物输出：`P + 原文件名`（如 `EGDA.mol2 → PEGDA.mol2`）
- 输出目录：与输入相同
- 时间戳目录：`outputs/run_YYYYmmdd_HHMMSS`

