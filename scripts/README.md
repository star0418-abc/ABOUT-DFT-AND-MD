# scripts/ 目录结构说明

重构后的脚本目录结构。

## 📁 目录结构

```
scripts/
├── lib/                    # 公共库（不直接作为 CLI 入口）
│   ├── __init__.py        # 模块导出
│   ├── recipe.py          # 统一 recipe.yaml 解析
│   ├── io_mol.py          # 分子文件读写与校验
│   ├── naming.py          # 统一命名规则
│   ├── validate.py        # 物理/文件校验
│   └── logging_utils.py   # 统一日志格式
│
├── cli/                    # 命令行入口（推荐使用）
│   ├── __init__.py
│   └── monomers.py        # 单体批处理统一入口
│
├── selfcheck.py           # 自检脚本
├── recipe_to_counts.py    # 配方计数
├── make_packmol_from_recipe.py  # Packmol 输入生成
├── build_pegda_network.py # PEGDA 交联网络生成
├── gaussian_log_to_pdb.py # Gaussian → PDB 转换
├── generate_topology.py   # 拓扑生成
├── generate_topology_gaff.py  # GAFF 拓扑
├── inject_charges.py      # CM5 电荷注入
├── sanitize_pdb_for_ligpargen.py  # PDB 清洗
├── mol_to_mol2_batch.py   # .mol → .mol2 批转换
│
├── run_workflow.sh        # 统一全流程入口
├── run_packmol.sh         # Packmol 驱动
├── run_htpolynet.sh       # HTPolyNet 驱动
├── run_gmx.sh             # GROMACS 驱动
└── open_in_windows.sh     # Windows 资源管理器打开
```

## 🚀 快速开始

### 批量生成寡聚体

```bash
# 使用新的统一入口（推荐）
python -m scripts.cli.monomers --names EGDA MMA TFEMA VA --n_repeat 3

# 单文件处理
python -m scripts.cli.monomers --input mol2/EGDA.mol2

# 使用 recipe.yaml 配置
python -m scripts.cli.monomers --recipe config/recipe.yaml
```

### 自检测试

```bash
python scripts/selfcheck.py
```

### 完整工作流

```bash
bash scripts/run_workflow.sh --recipe config/recipe.yaml
```

## 📚 公共库使用

### 读取 recipe.yaml

```python
from scripts.lib.recipe import load_recipe, get_oligomer_n

config = load_recipe("config/recipe.yaml")
n = get_oligomer_n(config, cli_override=None, default=3)
```

### 写入 MOL2 文件

```python
from scripts.lib.io_mol import write_mol2_strict, validate_mol2_output

write_mol2_strict(mol, "output.mol2", "PMMA", "PMM")
validate_mol2_output("output.mol2")  # 自动检查存在性、非空、必需段落
```

### 命名规则

```python
from scripts.lib.naming import get_polymer_output_path

# 输入 mol2/EGDA.mol2 → 输出 mol2/PEGDA.mol2
output = get_polymer_output_path("mol2/EGDA.mol2")
```

## ⚠️ 已弃用脚本

- `run_all_smoke.sh` → 使用 `python scripts/selfcheck.py`
- `run_full.sh` → 使用 `bash scripts/run_workflow.sh`

详细迁移指南请参见 [MIGRATION.md](../MIGRATION.md)

