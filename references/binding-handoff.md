# 多伙伴竞争占位与完整分子交接

这条路线把已验证的来源核对、刚体配准、运动检查和分子双视图连接起来，处理两个伙伴交替占据同一宿主界面的情况。先按 [联网结合研究](online-binding-research.md) 主动查复合物、物种、构建体、位点与状态；不要因为初始输入只有单体而停止检索。这里的示例有两份实验结合态，但没有一份同时包含两个伙伴的可用装配。

## 来源决定哪些组合可以使用

示例配置 `assets/eif4e-handoff.json` 使用人源 eIF4E、eIF4G1 和 4E-BP1：

| 来源与选择 | 实际保留及处理 |
|---|---|
| [5T46](https://www.rcsb.org/structure/5T46)，assembly 1，作者 A/B | eIF4E A 的已建模 31–206、eIF4G1 B 的 608–642；保留 A301 MGP 帽类似物。排除另一晶体拷贝、甘油和水。不是完整 mRNA 或全长 eIF4G1。 |
| [4UED](https://www.rcsb.org/structure/4UED)，assembly 1，作者 A/B | A 仅用于配准；B 的已建模 49–83 整体迁移，49 Met 是保留并披露的表达标签。该结构没有建模帽类似物。 |
| [5BXV](https://www.rcsb.org/structure/5BXV)，筛查后未使用 | 含小鼠/人混合来源，本例选用两份人源宿主结构；并不表示该条目错误或不可用于其他问题。 |

169 对同作者编号、同残基名的 Cα 显式对应，4UED 宿主到 5T46 宿主的 RMSD 为 1.1201 Å。保留 5T46 宿主与 MGP 固定，只搬运完整 4UED B。因此 BP1/5T46 宿主/帽类似物的组合属于**来源引导的推定组合**，不是某份实验结构已经观测到的三者同态；不得把两份宿主形变合并。

机制依据：[Grüner 等，2016](https://doi.org/10.1016/j.molcel.2016.09.020) 的出版者 highlights/summary、[Peter 等，2015](https://doi.org/10.1016/j.molcel.2015.01.017) 的原始摘要，以及实际沉积坐标、entity、链和 assembly。本次未读取两篇完整正文。前者支持人源 eIF4G 同时利用 dorsal/lateral 界面；不能直接沿用更早果蝇模型中“lateral 面不被 eIF4G 使用”的解释。

对齐后的两个完整结合姿态有 985 对重原子距离小于 2.4 Å，最近 0.0349 Å，故这些特定姿态不能同时保留。这个几何结果不能排除其他部分结合构象，也不能给出竞争速率、亲和力或真实置换路径。

## 准备、检查和显示

先设置本机 `SKILL`、项目目录 `PROJECT`、含 NumPy/SciPy/gemmi/scikit-image 的 `PYTHON`、`BLENDER` 和宋体/Times 字体路径。使用项目新目录；资源与影片留在项目内。

```bash
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 5T46 --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 4UED --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_binding_handoff.py" \
  --sources "$PROJECT/sources" --config "$SKILL/assets/eif4e-handoff.json" \
  --output "$PROJECT/prepared" --mode loop --frames 288 --fps 24
"$PYTHON" "$SKILL/scripts/audit_binding_handoff.py" \
  --sources "$PROJECT/sources" --config "$SKILL/assets/eif4e-handoff.json" \
  --assets "$PROJECT/prepared" --output "$PROJECT/source-motion-audit.json"
"$BLENDER" -b --python-exit-code 1 \
  --python "$SKILL/scripts/build_molecular_sequence_scene.py" -- \
  --asset-dir "$PROJECT/prepared" --display "$SKILL/assets/eif4e-handoff-display.json" \
  --output "$PROJECT/scene" --font-zh "$FONT_ZH" --font-en "$FONT_EN" --render
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_molecular_sequence_scene.py" -- \
  --output "$PROJECT/scene/audit.json"
```

下载后的来源哈希若不匹配，重新审查条目版本；不要仅改期望哈希来绕过审查。配准要求明确的残基表、构建体和替换记录，不能把 Cα RMSD 通过当作侧链或配体已经无冲突。

准备器输出完整重原子身份、每帧坐标、运动路径、来源清单和背景网格；直接使用 [共用分子序列显示](molecular-sequence-display.md)，不用为新案例复制场景生成器。两个视图、材质、字体、渲染器和合成器仍是独立显示参数。背景是完整源蛋白形状的艺术分布，不表示实验测得的浓度或空间组织。

连续模式为第一伙伴离开、第二伙伴结合、第二伙伴离开、第一伙伴返回。`--mode teaching --frames 384` 使用同一几何路线，并在两个结合态停留。为教学模式复制显示 JSON 到项目，将 `captions` 替换为新 `prepared/manifest.json` 的 `captions`，如需检查停留可设 `key_fractions` 为 `[0.05, 0.55]`。字幕区间必须与准备结果一致；不能仅把连续视频改标题就称为步骤教学。

## 包绕界面的运动和共同占位

普通弧线对这两个包绕片段未找到合格路径。`rigid_escape.py` 在局部约束下同时平移、转向完整片段，每个尝试步检查 17 个实际非线性子步，再按平移长度和转角形成路径参数。失败时保留失败记录并改路线或重新核定来源构象，不能放宽原子距离门槛、切断片段或拉伸蛋白。

顺序规划第一伙伴离开宿主与固定配体，再让第二伙伴绕开宿主、配体及第一伙伴的实际等待位置；最终仍须检查整条同步时间线上的**所有组件对**。单独检查“每个伙伴对宿主”会漏掉两个移动伙伴相撞，本例曾实际遇到该问题。

`audit_binding_handoff.py` 不调用准备器、配准函数或运动规划器，从原始 mmCIF 重新读取身份并做独立刚体配准。检查固定宿主/配体、结合端态、片段内部结构、循环闭合、全组件对及帧间四分点。保存工程另由共用检查器核对两个视图的全部坐标、显示变换、完整分子取景、字体及字幕。

2026-11-06 的本例：连续模式 289 个姿态、教学模式 385 个姿态；2,057 个保留重原子，0.24 nm 统一中心距离门槛，无源接触豁免。连续模式采样最近距离 0.242624 nm；整数帧保持源刚体。Blender 对原子坐标做线性帧间插值，半帧刚体拟合最大误差约 0.00173 nm，检查上限 0.005 nm；不要声称任意连续时间都严格刚体或保证无穿插。八项错误输入检查包括同时占位、非刚体形变、帽移动/丢失、标签改号、错端态、循环未闭合和缺帧。

画面检查同时关注**投影遮挡**：三维没有穿插，不代表观众能分清等待者与结合者。调整整体视角和侧视角，让离开的伙伴在二维投影中分开；完整构建体留在两个窗口内，不靠剖切或遮掉部分分子解决。还要检查两位伙伴同时等待的画面，避免投影重叠被误读为二者结合。

默认只渲染两个正式尺寸结合态和两位伙伴都已离开时的状态，以及对新运动有疑问的 1–3 秒短预览。既有实测使用 25 张 854×480 帧、步长 2、12 fps 编码，保留源 24 fps 的原速时间间隔。完整 12/16 秒视频未因此宣称已验收，见 [高效检查和清理](review-and-retention.md)。

## 推广时必须更换的内容

- **相容性与来源**：重新联网核对两个伙伴是否确实共享/遮挡同一界面。可共存的复合物应走共同装配路线；本准备器要求这两个完整已配准结合姿态实际冲突，不适用于把任意两种蛋白强行写成竞争者。
- **共同宿主及伴随物**：选定来源宿主，映射另一来源；当前入口只接受固定宿主伴随配体，并以相同作者链名选择两个来源的宿主；不同宿主链名需先显式扩展映射，不能改写原始身份。移动伙伴自带配体、多操作 assembly、跨膜侧别或第三个移动伙伴需要先扩展资格和全组件检查，不能静默遗漏。
- **运动与画面**：重新选择离开方向、距离与等待位置，检查所有伙伴的同步坐标和投影。局部优化可能失败；通过只是这个构象、门槛和采样条件下的示意路线，不是全局最优对接或分子动力学。
- **机制事件**：此例没有磷酸化、解折叠、RNA 招募、eIF4A 或完整翻译起始。两个源结合片段在游离阶段也保持刚体示意。要加入调控或催化，另建立状态证据和事件/反应账本，不能从占位动画自动推断。

发布只包含本页、脚本和参数；原始下载、NPZ、`.blend`、字体、帧和短片留在项目。解压后的包应能利用已核定原始输入重建本工作流。
