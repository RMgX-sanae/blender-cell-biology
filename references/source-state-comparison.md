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
- `alignment` 明确同一作者编号体系中的 CA 段及允许的覆盖/RMSD。当前配置用 eRF1 N 端 6–135，其他对象一起移动；不是逐蛋白居中，也没有逐结构域拉伸。不同编号/插入码拟合、同源替换和柔性运动需另外建立映射。
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

本例包含 23,633 个所选聚合物重原子、52 个配体重原子。5LZT N 端拟合 RMSD 为 1.1211 Å，5LZV 为 0.1560 Å；保留这些差异。最初 1.0 Å 静态拟合门槛未通过，改为记录明确的 1.5 Å 对照门槛；没有放宽界面碰撞条件或修改原坐标。

从源码包解压副本和原始坐标缓存重建了 88 个相同数组、11 个完整聚合物显示对象和 52 个配体原子对象。19 个来源/身份/坐标反例被拒绝，保存工程另查变换、配体和相机裁边。仅检查一张最终布局的 Cycles 技术对照图及重建图，未渲染影片；两图存在少量渲染数值差异，不能称逐像素相同。默认深色比较图已检查，其他调色/引擎仍按项目复核。

## 从静态端态继续推进机制

这次迁移证明来源核对和独立端态显示可执行。它没有实现 peptidyl-tRNA 水解、eRF3 GTP 循环、ABCE1 ATP 循环、亚基分离、mRNA/tRNA 释放或回收。此前 9SPI 延伸例也是另一实验体系，不能把它和这些端态直接接成同一份底物的电影。

需要终止/回收动画时，先选定持续构建体和底物，建立肽链/tRNA/水的反应身份、因子占位与物质池，再处理来源支持的局部构象。刚体候选不能自动解释 eRF1 的大幅结构域重排；两个状态看起来可以贴合，也不证明插值路径无穿插。按 [机制深度](mechanistic-depth.md) 和 [反应/工作流迁移](workflow-transfer.md) 验证新增事件后再显示。开发只渲染回答当前问题的关键帧，分发仅包含脚本、配置和来源说明。
