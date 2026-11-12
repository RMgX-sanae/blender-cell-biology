# 先核定来源状态，再设计构象运动

用于同一蛋白在不同伙伴、配体或反应状态中的结构对照。先按 [联网研究](online-binding-research.md) 主动查找条目、assembly、序列映射及原文；然后将实际判断写入项目配置，而不是只凭链名配对。这个入口保留分开的静态端态，不自动插值、不把突变或类似物改成天然底物。

## 翻译终止的实际迁移

配置 [termination-source-states.json](../assets/termination-source-states.json) 来自 [Shao 等的原始研究](https://doi.org/10.1016/j.cell.2016.10.046)，已读取作者公开 PDF 的结果、Figure 7 和终止复合物制备方法，并核对三份实际 mmCIF。其兔核糖体与人释放因子是原实验重构体系，不能改称纯人或纯兔天然循环。

| 来源 | 关键身份和配体 | 可展示的内容 |
|---|---|---|
| [5LZT](https://www.rcsb.org/structure/5LZT) | 人 eRF1 `ii`、人 eRF3a `jj`；`jj:700 GCP` 为 GMPPCP，保留 `jj:701 Mg` | 识别态局部复合物；不是天然 GTP 水解轨迹 |
| [5LZU](https://www.rcsb.org/structure/5LZU) | 人 eRF1 `ii`；没有对应 eRF3/ABCE1 `jj` 聚合物 | 已就位的野生型 eRF1 构象；不能仅凭未建模肽链证明完整释放过程 |
| [5LZV](https://www.rcsb.org/structure/5LZV) | 人 eRF1 `ii` 的作者 183/184 为 Ala；兔 ABCE1 `jj`，保留 `jj:600/601 SF4` | AAQ 阻断构建体中的 ABCE1 结合态；不是已完成释放/回收的天然端点 |

三份 assembly 1 都由恒等操作组成。`jj` 在不同文件中代表不同蛋白；顶层 mutation 字段为空也不能据此判为无突变，必须看实际残基和序列差异。每份 mRNA `hh:49–51` 均核对为 UGA，`hh:100 Mg` 保留。5LZV 未建模 ATP，不能因此宣布生理状态无核苷酸。

显示保留完整的**已建模** eRF1、对应伙伴、P-tRNA 和 mRNA；其他核糖体成分、E-tRNA 与 nascent chain 整体省略。这是局部选择，不是拆散核糖体。mRNA 与埋藏配体可能被完整蛋白遮挡；保存坐标通过不等于界面原子在总览图中可见，需要另做完整单位的观察窗口或独立标记特写。

## 可复用的准备与显示

`prepare_source_states.py` 沿用已有的 `checked_source`、`atom_table`、`select_ligands`、`fit_template` 和 `surface`。新增的是**状态输入核对**和多来源对照：

- 每份来源用 `pdb`、`sha256`、`assembly_id` 绑定；`chains` 指定作者链、精确 description、可选 UniProt/源物种、关键残基/原子、标签和颜色。表达宿主不当作蛋白来源物种。保留 label/entity、sequence-difference 记录和源物种类别。
- `ligand_instances` 精确指定非聚合物身份。复用 helper 的 `owner_role` 只作静态分组；所有配体与该状态的全部链接受同一次刚体变换，不表示已经证明某种生理归属或金属配位。
- 当前仅接受原始 mmCIF 中恒等操作的 assembly。多操作装配须先走专用装配准备；不要删掉操作表来强行通过。
- `alignment` 明确 CA 对应及允许的覆盖/RMSD。旧格式继续支持同一作者编号体系；当前终止配置用 eRF1 N 端 6–135。不同链或编号用下述逐状态映射，不能靠链名或相同数字猜测同源关系。其他对象一起移动；不是逐蛋白居中，没有逐结构域拉伸，也不自动实现柔性运动。
- 共同 `origin_angstrom`、`basis_columns` 只定义显示坐标系。每栏同一尺度，排版仅平移；`build_source_state_comparison.py` 复用既有表面、材质、灯光、字体与合成入口。

```bash
# 使用新项目目录；PYTHON 含 gemmi/numpy/scipy/scikit-image。
for ID in 5LZT 5LZU 5LZV; do
  "$PYTHON" "$SKILL/scripts/fetch_structure.py" "$ID" --output "$PROJECT/sources"
done
"$PYTHON" "$SKILL/scripts/prepare_source_states.py" --sources "$PROJECT/sources" \
  --config "$SKILL/assets/termination-source-states.json" --output "$PROJECT/states"
"$PYTHON" "$SKILL/scripts/test_source_states.py" --sources "$PROJECT/sources" \
  --config "$SKILL/assets/termination-source-states.json" --asset-dir "$PROJECT/states" \
  --output "$PROJECT/source-check.json"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_source_state_comparison.py" -- --asset-dir "$PROJECT/states" \
  --output "$PROJECT/comparison" --font-zh /path/to/Songti.ttf --font-en /path/to/Times-New-Roman.ttf \
  --width 1800 --samples 24 --render
"$BLENDER" -b "$PROJECT/comparison/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_source_states.py" -- --asset-dir "$PROJECT/states" \
  --output "$PROJECT/saved-scene-check.json"
```

`audit_source_states.py` 在普通 Python 中配合 `--sources --config` 比较原始原子身份、坐标、配体和共同拟合；在重开的 Blender 工程中比较表面顶点/拓扑、对象清单、完整刚体变换、配体位置及相机边界。`test_source_states.py` 的来源反例是这个翻译终止配置专用，不是任意新项目自动适用的测试。

### 不同作者链和编号的映射

用 `alignment.reference_chain` 指定参考链，`alignment.maps` 的键必须恰好覆盖全部来源 PDB。每项包含 `chain`、`residue_pairs`、`expected_substitutions`、`min_CA_coverage`、`max_CA_rmsd_angstrom` 和说明对应依据的 `note`；字段意义沿用单体结合的 `fit_template`。`residue_pairs` 每对按“当前来源作者编号 → 参考来源作者编号”排列。普通作者编号必须是整数，带插入码写作 `"184A"`；`184` 和 `"184A"` 是不同残基，`184.0`、`true` 或 `"184"` 不代替整数。参考状态的映射必须保留自己的链与编号。显式格式不能混用旧 `chain/author_residues` 格式。

通过序列、实际坐标和必要的 [SIFTS 对应](online-binding-research.md) 审查映射后写入项目配置；已有 `source_residue_mapping.py` 可提供有来源的对应记录。不能用最佳拟合代替同源性判断。覆盖率须有限且在 0–1 内，RMSD 上限须有限且非负；重复、缺失、共线或含非有限坐标的 CA 不用于拟合。替换清单必须与实际氨基酸差异一致。未出现替换也不证明整个构建体相同，仍需核对完整来源记录。

准备结果保留原 mmCIF 的 `polymer_type` 供蛋白/核酸材质选择，原始来源检查器复核其身份。新的资产入口先检查所声明的状态、链与必需伴随组分清单，再进入 Blender；不能仅根据仍存在的文件哈希宣布组分完整。伴随组分按身份匹配，不依赖 JSON 列表顺序。

2026-11-12 的有限回归复用了既有 2PTN/2PTC：作者链 A → E，保留 184A/188A/221A，各自保留来源钙离子；拟合仍为 0.3351 Å。不同数字编号另用技术坐标验证，不将它冒充新的蛋白研究。终止例新增输入检查共拒绝 25 类反例，88 个既有数组未改变；单体结合例的 22 类反例和 145 个姿态保持通过。原来会被建模入口漏掉的 SF4 组分现已明确拒绝。检查针对来源与交接，不等于验证新的机制或动力学。

### 材质与最终像素

通过 `--profile "$SKILL/assets/style-structure-materials.json"` 选择分角色材质，流程见 [材质与灯光迁移](material-and-lighting.md)。新准备的聚合物按已核定的 `polymer_type` 区分蛋白与核酸，完整结构、原始尺度与各栏共同刚体变换不变。旧准备文件没有该字段时保留通用蛋白表面的兼容外观；要区分核酸材质应从原始来源重新准备，不能按标签猜类型。此入口是清楚的独立来源对照，可用纯色背景；不是细胞环境或机制影片。

本例包含 23,633 个所选聚合物重原子、52 个配体重原子。5LZT N 端拟合 RMSD 为 1.1211 Å，5LZV 为 0.1560 Å；保留这些差异。最初 1.0 Å 静态拟合门槛未通过，改为记录明确的 1.5 Å 对照门槛；没有放宽界面碰撞条件或修改原坐标。

从源码包解压副本和原始坐标缓存重建了 88 个相同数组、11 个完整聚合物显示对象和 52 个配体原子对象。19 个来源/身份/坐标反例被拒绝，保存工程另查变换、配体和相机裁边。仅检查一张最终布局的 Cycles 技术对照图及重建图，未渲染影片；两图存在少量渲染数值差异，不能称逐像素相同。默认深色比较图已检查，其他调色/引擎仍按项目复核。

## 从静态端态继续推进机制

这次迁移证明来源核对和独立端态显示可执行。它没有实现 peptidyl-tRNA 水解、eRF3 GTP 循环、ABCE1 ATP 循环、亚基分离、mRNA/tRNA 释放或回收。此前 9SPI 延伸例也是另一实验体系，不能把它和这些端态直接接成同一份底物的电影。

需要终止/回收动画时，先选定持续构建体和底物，建立肽链/tRNA/水的反应身份、因子占位与物质池，再处理来源支持的局部构象。刚体候选不能自动解释 eRF1 的大幅结构域重排；两个状态看起来可以贴合，也不证明插值路径无穿插。按 [机制深度](mechanistic-depth.md) 和 [反应/工作流迁移](workflow-transfer.md) 验证新增事件后再显示。开发只渲染回答当前问题的关键帧，分发仅包含脚本、配置和来源说明。
