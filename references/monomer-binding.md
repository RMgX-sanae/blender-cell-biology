# 从单体结构构建结合候选

用户要求只有单体时也尽力实现可能的结合方式，并默认联网。先执行 [联网结合研究](online-binding-research.md)，主动寻找适用复合物、构象与界面证据，核对物种、构建体、定位、配体和寡聚状态，再选择下面的路线。直接可用的实验复合物优先走 `prepare_interaction.py`；本页用于确实需要组合不同来源单体的情况。不要把两个表面接触的模型称为已证实的结合。

| 输入证据 | 可执行路线 | 输出含义 |
|---|---|---|
| 同源或同构建体的已知复合物 | `prepare_binding_candidates.py --method template`：明确 CA 配对映射、独立刚体拟合、局部姿态搜索 | 模板引导的结合假设；不是该单体组合的实验结构 |
| 界面残基、交联或其他距离约束，缺少可用复合物坐标 | 同脚本 `--method restraints`：由两个界面区域估计朝向，采样 16 个转角，再优化刚体姿态 | 满足所给约束的几何候选；约束可能不唯一，也可能无法同时满足 |
| 两个单体且没有可靠界面信息 | 先检索互作/定位依据，提出多个明确标为探索性的界面假设；分别写入项目约束，再用上述路线检查 | 美术或科学假设的可视化，不是已实现的盲对接预测器 |

模板映射需要显式的源/模板作者残基对、预期替换、覆盖率与 RMSD 上限。不同编号、同源蛋白或结构域只能在审查序列/结构对应后建立这些映射；不要默认编号相同就等价。不自动补环、换侧链、恢复突变或改变主链。保持完整已建模链，模板和候选都不能通过拉伸或切断分子来适配。

当前入口处理两个刚性聚合物伙伴，每个来自一个选定作者链；氢、非聚合物和非 A 替代位置被过滤。需要辅因子、金属、糖链、膜或多链伙伴时，先扩展该项目的实体与障碍物选择，不能把过滤后的两个蛋白当作完整体系。含插入码造成身份歧义的来源会拒绝；多构象和局部柔性还需专门实现。

## 已实际运行的迁移案例

- [1A2P](https://www.rcsb.org/structure/1A2P)：游离野生型 barnase，作者链 A；保留 857 个建模重原子。
- [1A19](https://www.rcsb.org/structure/1A19)：游离 barstar C82A，作者链 A；保留 719 个建模重原子。
- [1BRS](https://www.rcsb.org/structure/1BRS)：作者链 A/D 作模板；barstar 是 C40A/C82A。单体中的 Cys40 保留，未暗改为 Ala。

项目配置为 `assets/binding-monomers.json`。CA 拟合 RMSD 分别约 0.464 Å 和 0.437 Å；直接套入仍有严重重叠，不能作为合格端态。初次尝试要求恢复七个模板残基对的近距离接触，两条刚体搜索均未通过。随后把两类问题分开：用七个 **CA 距离**约束界面位置，原来的侧链接触保留为独立诊断，未恢复的接触仍明确报告失败。CA 距离的 ±0.2 nm 范围是选定的模型容许范围，不是实验误差。

当前配置的模板路线保留 11 个、约束路线保留 4 个聚类后的候选。它们是同一界面附近的不同刚体姿态，**不是 11 或 4 种经过发现/证明的生物学结合模式**。两条路线的约束都来自已知 1BRS 界面，因此是方法迁移对照，不能当作独立预测准确率测试。约束路线运行时不读取模板坐标。

选定候选接入了已有 `interaction_motion.plan`、高斯表面准备、`build_interaction_scene.py`、字体/合成与保存工程检查。各检查 145 个姿态（144 帧及虚拟接缝），整个伙伴只做刚体运动。候选生成不是完整分子对接软件：未包含力场、溶剂、侧链优化、诱导契合、结合能、亲和力或真实时间。

两条路线均从解压源码与经过指纹核对的原始缓存独立重建，各 15 个几何/身份/运动数组完全相同。约束路线的重建目录只有两个单体来源，配置还删除了模板映射。各保存工程已检查；Cycles 关键帧已目视复核，没有为这次迁移渲染完整视频。模板路线的 17 类错误输入、约束路线适用的 14 类，以及 4 个保存场景证据/标签错误均被拒绝。既有 AP2 场景也通过共用检查器的回归检查。这些结果仅覆盖当前示例与参数。

## 从源码重建

沿用 [运行环境](production-workflow.md) 的 Python、Blender 和字体设置。输出全部写到 skill 外的新项目目录。

```bash
python3 "$SKILL/scripts/fetch_structure.py" 1A2P --output "$PROJECT/sources"
python3 "$SKILL/scripts/fetch_structure.py" 1A19 --output "$PROJECT/sources"
python3 "$SKILL/scripts/fetch_structure.py" 1BRS --output "$PROJECT/sources"

"$PYTHON" "$SKILL/scripts/prepare_binding_candidates.py" \
  --config "$SKILL/assets/binding-monomers.json" --sources "$PROJECT/sources" \
  --method template --output "$PROJECT/template" --mode loop --frames 144 --fps 24

# 第二种初始化不读取 1BRS.cif；本例约束仍来自 1BRS，不能称为盲预测。
"$PYTHON" "$SKILL/scripts/prepare_binding_candidates.py" \
  --config "$SKILL/assets/binding-monomers.json" --sources "$PROJECT/sources" \
  --method restraints --output "$PROJECT/restraints" --mode loop --frames 144 --fps 24

"$PYTHON" "$SKILL/scripts/test_binding_candidates.py" \
  --sources "$PROJECT/sources" --prepared "$PROJECT/template" \
  --output "$PROJECT/template/source-and-path-check.json"

"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_interaction_scene.py" -- \
  --asset-dir "$PROJECT/template/interaction" --output "$PROJECT/template/scene" \
  --font-en "$FONT_EN" --style cinematic --engine cycles --width 1280 --samples 24 \
  --partner-colors '#58DCEC' '#FFB660'

"$BLENDER" -b "$PROJECT/template/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_interaction_scene.py" -- \
  --output "$PROJECT/template/saved-scene-check.json"
"$BLENDER" -b "$PROJECT/template/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/template/review" --frames 37,73
```

`--candidate 2` 可在新输出目录选择第二个保留候选；不保证所有输入都有多个候选。`--mode teaching` 切换为分步动作，`loop` 保持循环；时序改变后重新检查。没有通过的候选或路径就保留失败报告，调整证据、构象或算法，不降低穿插门槛凑成片。截图与影片不放进 skill 包。

## 项目输入与检查

复制配置后替换两个来源、SHA256、作者链、准确标签与构建体说明。模板路线还替换模板与显式 CA 映射。约束路线可删除整个 `template` 字段；每条约束给 `anchor`/`moving` 作者残基编号、可选 `anchor_atom`/`moving_atom`、`range_nm` 和来源/假设。普通编号用整数（如 `184`），带字母插入码时用字符串（如 `"184A"`）；映射与约束都区分这两个身份。输出保持整数 `residue_numbers`，有插入码时另存 `insertion_codes`。省略 atom 字段表示两残基所选重原子间的最小距离；CA 约束与侧链最小距离不能混淆。交联上限需按试剂、位点和可达性另定，不能沿用本例数值。

输出 `candidates.json` 保存全部尝试、刚体变换、筛除理由、聚类和映射；`interaction/manifest.json` 保存选中候选、原子来源、端态与路径检查。候选排序的 `geometric_penalty` 只含几何约束罚项及对初始姿态的小偏移惩罚；不同初始化的分数也不是统一的物理量，更不能解释为概率或亲和力。聚类阈值控制显示冗余，不证明新的功能状态。

每个候选端态要求全体原子中心距离通过指定门槛，本例为 0.24 nm；**没有“这是来源接触，所以允许重叠”的豁免**。选中候选的每个运动姿态再以同一门槛检查。该筛查仅排除严重重叠，不等于完整范德华/化学合理性判断。CA 约束通过不表示侧链诊断通过，也不应在画面上额外画成已形成的盐桥或氢键。

测试独立从原始单体核对全部身份与刚体坐标，用稠密全原子距离复核路径，并拒绝拉伸、镜像、未申报突变、错误映射、丢失残基、错误来源、候选冒充测量值和重叠豁免。保存工程另核对几何、父对象、逐帧变换、取景及可读的 `Binding hypothesis` 标识。这个必要的短标签不替代来源文档，也不增加细小的解释字幕。

将候选用于新机制时先保留独立分支：配体占据、核苷酸状态、竞争伙伴、同源/异源寡聚、不同界面和化学计量可能决定不同结合方式。每个分支需要自己的依据、约束与可达路径。柔性变化明显或所有刚体候选失败时，转入有来源的构象集合/专门柔性建模；当前脚本不能通过几何得分判断真实结合是否存在。

## 来源配体与结合态构象的迁移

每个 `anchor`/`moving` 来源可显式添加 `ligand_instances` 与 `ligand_selection_note`，沿用 [配体实例选择](protein-interactions.md#明确保留哪一份配体) 的字段。`owner_role` 必须与这个来源伙伴相同。源配体与蛋白共享拟合、候选及逐帧刚体变换，单独保存原子身份和来源，进入两伙伴之间的全部穿插检查；不会把钙离子并入蛋白表面，也不会用复合物里的另一份离子替代单体的离子。仅支持明确选择的 NonPolymer 重原子；共价修饰、分支糖链、溶剂交换和配体动力学另行建模。

`assets/binding-source-conformer.json` 是一个有意保留来源差异的迁移对照：游离 trypsin [2PTN](https://www.rcsb.org/structure/2PTN) A 链及其 A:480 钙，与网上查得的 [2PTC](https://www.rcsb.org/structure/2PTC) I 链 **结合态 BPTI** 组合，用 E/I 作参考。包含 1,564 + 422 个聚合物重原子和 1 个钙；trypsin 的 `184A`、`188A`、`221A` 保持独立身份。它不是两个游离单体对接成功，也不是测得的混合来源复合物。一般教学任务已有合适 2PTC 时直接沿用其完整结合态；这个对照用于检查需要保留单独来源时的代码迁移。

```bash
python3 "$SKILL/scripts/fetch_structure.py" 2PTN --output "$PROJECT/sources"
python3 "$SKILL/scripts/fetch_structure.py" 2PTC --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_binding_candidates.py" \
  --config "$SKILL/assets/binding-source-conformer.json" --sources "$PROJECT/sources" \
  --method template --candidate 2 --output "$PROJECT/transfer" --frames 144 --fps 24
"$PYTHON" "$SKILL/scripts/test_binding_candidates.py" \
  --sources "$PROJECT/sources" --prepared "$PROJECT/transfer" \
  --output "$PROJECT/transfer/source-check.json"
```

本例要求五对指定原子的接近范围，不能用 CA 位置通过冒充侧链接触通过。两组尝试使用游离 BPTI 4PTI、6PTI 时均未满足既定几何检查，失败报告留在项目；6PTI 的源磷酸也参与了那次尝试。查得结合态后，迁移对照保留了 13 个几何候选；第一名的接近路径仍有穿插，第二名才通过 145 个姿态的严格检查。几何排序不替代路径检查。原子近邻约束通过也不能证明氢键方向、质子状态、催化、亲和力或真实结合路径。记录构象替换，不通过删除失败来源、拉伸模型或改变检查门槛掩盖差异。

该对照已从解压源码和重新在线下载的 2PTN/2PTC 重建：19 个几何/身份/运动数组及候选、运动记录相同；原 Barnase 对照的 15 数组保持相同。源数据检查拒绝 22 类错误配置，另把插入码抹掉或钙移位并重算文件哈希，仍被原始来源检查拒绝；实际保存后把钙挂到另一伙伴的工程也被拒绝。两次保存工程的 145 个姿态通过检查，接近/结合关键帧作了技术复核，重建结合帧最大相差一个 8-bit 通道值。未渲染整片，未把该对照列为自由单体诱导契合或天然相互作用预测。
