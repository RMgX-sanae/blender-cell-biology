---
name: blender-cell-biology
description: Build editable Blender models and animations of cells, organelles, membrane receptors and cellular processes from microscopy, molecular structures and literature. Analyze reference images and reproduce their visual style through materials, lighting, cameras and selectable compositors. Use for cell biology 3D scenes and mechanism animations.
---

# Blender 细胞生物学建模与动画

把生物学问题转成有来源、可编辑、可复现的 `.blend` 工程及图片/动画。支持清晰教学示意和写实细胞环境；风格、渲染引擎、合成器分别选择。中文交流为默认，标签按交付语言设置。

**默认背景偏好：**图片和视频优先采用有景深虚化的细胞微结构环境，按所处区室选择有身份的结构，并保留前后层次。主反应区、关键分子与文字保持清楚；项目指定纯色、透明底或其他参考风格时按项目选择。实施与可运行示例见 [cellular-environment.md](references/cellular-environment.md)。

**质量优先：**先读 [quality-standard.md](references/quality-standard.md) 和 [production-workflow.md](references/production-workflow.md)。默认使用精细制作入口。旧 `build_scene.py` 仅作机制白模，不得再当作成品示例；“可运行、方向正确”不能替代精细建模、参考风格复现和连续动画的验收。

**分发边界：**用户要求 skill 包不包含生成的视频、渲染帧和项目输出。保留生成脚本、参数、来源与资产获取方法；下载、建模、渲染和检验结果写入 skill 外的项目目录。用 `scripts/package_skill.py` 打包，并从解压副本实际重建代表性产物；仅指向本机旧成品不能证明工作流可复现。

## 接到任务后

1. 明确物种、细胞类型、状态、要解释的过程、观众和交付形式。缺失信息可用明确标注的教学假设开始；若物种/状态会改变机制，先澄清再构建相应部分。不要将哺乳动物、植物、细菌的结构混搭。
2. 与用户确认或说明当前选择的 **表现风格 / 渲染引擎 / 合成方案**。用户给参考图片时，实际看图并按 [style-from-reference.md](references/style-from-reference.md) 识别与复现风格，先出可对比的小样。允许 Blender Compositor、可用的 Fusion/Natron/After Effects 等外部合成器；不要把合成器与渲染引擎混为一谈。见 [rendering-and-compositing.md](references/rendering-and-compositing.md)。
3. **先检索、参考并尽量复用公开资产，再决定哪些部分需要新建。** 按 [public-assets-first.md](references/public-assets-first.md) 检查本地缓存、原始结构/显微数据及可编辑模型；有合适资产时导入、修整、实例化，不从随机团块重新造形。记录物种、状态、装配、尺度、许可与加工过程。用 [sources-and-evidence.md](references/sources-and-evidence.md) 核对机制依据；检索摘要不能替代未读的图例和方法。
4. 写一个简短场景说明：对象及 compartment、尺度、证据级别、事件顺序、视觉简化和检查标准。可复制 [scene-brief.json](assets/scene-brief.json)，删除本项目用不到的字段。
   用户要求“全过程/上游到下游”时，按 [mechanistic-depth.md](references/mechanistic-depth.md) 建立因果事件表、化学状态和物质回收闭环；局部动作的改善不能替代缺失的机制。
   动画形式另行选择：**分步教学**用分段、停顿和同步解释；**连续循环**保持运动、镜头和首尾衔接。它独立于材质风格、渲染引擎和合成器。用户此次选择连续循环时，不自动插入章节卡、停顿或旁白。实现与接缝检查见 [blender-workflow.md](references/blender-workflow.md)。
5. 先做低分辨率静态白模和关键时点分镜，确认生物学拓扑与画面可读性，再加材质、环境拥挤和动画。可在首个模型、动画分镜、最终视频阶段询问满意度；无必要时不要让每一步都等待批准。
6. 实际运行 Blender，抽查关键帧，查看动画文件，再交付工程、生成脚本、来源和简化说明。检查办法见 [validation.md](references/validation.md)。未运行的部分标为“配方/待验证”，不得称为已验证自动化能力。

## 按任务读取

| 任务 | 读取 |
|---|---|
| 精细模型、真实结构资产、连续动画、可复现渲染命令 | [production-workflow.md](references/production-workflow.md)，默认入口 |
| 寻找、选择、缓存和复用网上公开资产 | [public-assets-first.md](references/public-assets-first.md)，所有新对象先走此流程 |
| 参考图、截图、视频关键帧的风格识别与模仿 | [style-from-reference.md](references/style-from-reference.md) |
| 扩展参考风格、几何节点、插画/透明/丝带/密度/显微外观小样 | [style-library.md](references/style-library.md)：18 类可执行起点、来源目录、参数与迁移方法；显微外观是示意，不能冒充真实成像 |
| 组合渲染滤镜、保护文字、主体强调与景深层次 | [render-recipes.md](references/render-recipes.md)：可编辑有序节点链、引擎与 pass 限制、实际画面检查 |
| 细胞、细胞器、膜和分子的背景风格对照 | [cellular-environment.md](references/cellular-environment.md)：按尺度选择环境，浅/深底、独立景深层与高反差分隔线 |
| 将结合/释放动作推广到不同蛋白、肽段与抑制剂 | [protein-interactions.md](references/protein-interactions.md)：明确源构象、保留配体、逐帧路径检查与独立动画模式 |
| 激酶供体进入、磷酸转移、产物释放及下一轮 | [phosphotransfer.md](references/phosphotransfer.md)：PKAc/SP20 专用实例，原子身份守恒、口袋开合与同步双视图；其他酶需重新建立反应与来源 |
| GEF 核苷酸交换与 GAP 促进 GTP 水解 | [ras-switch.md](references/ras-switch.md)：两套独立入口、完整原子身份、同步双视图及产物释放；连续状态循环仍需衔接检查 |
| 独立重建 Arp2/3 与 formin 完整案例 | [portable-actin-films.md](references/portable-actin-films.md)：公开坐标到两条完整视频，不依赖原制作目录；局部接触诊断仍须审核 |
| 分子剖面、遮挡与内部观察窗口 | [section-views.md](references/section-views.md)：用完整单体和前后层次表达 |
| 不同细胞类型、组织极性与基本形态 | [cells-and-organelles.md](references/cells-and-organelles.md) 的细胞表 |
| 细胞器、膜层、骨架与内部连接关系 | 同文件的细胞器表和拓扑检查 |
| 受体、通道、糖链及空间分布 | [receptors-and-distributions.md](references/receptors-and-distributions.md) |
| 分裂、DNA/RNA、翻译、运输等动画 | [process-atlas.md](references/process-atlas.md)，只读相关过程组 |
| 肌动蛋白聚合/解聚、踏车、切断、退火、成束、Arp2/3、封端及微管动态 | [polymer-templates.md](references/polymer-templates.md)：11 个可执行模板 |
| 肌动蛋白上游到下游：Arp2/3 与 formin 两条完整路线 | [actin-full-cycle.md](references/actin-full-cycle.md)：两套各 25 阶段的机制计划、状态闭环、资产选择、配体对照及已渲染的本地完整视频案例；其他情境需按项目制作 |
| formin、cofilin、中间丝、septin、FtsZ、ESCRT-III、clathrin、COPII 装配 | [assembly-storyboards.md](references/assembly-storyboards.md)：可复制的分镜与资产选择模板，按项目实现 |
| 将机制变成 Blender 动画 | [blender-workflow.md](references/blender-workflow.md) |
| V 型 ATPase 的整体与 c 环双视图连续循环 | [vatpase-loop.md](references/vatpase-loop.md)：公开结构准备、可移植生成器、同步检查及渲染命令 |
| PDB/mmCIF、EM/分割网格进入 Blender | [data-import.md](references/data-import.md) |
| 风格、Eevee/Cycles、合成节点和外部合成 | [rendering-and-compositing.md](references/rendering-and-compositing.md) |

## 不可被视觉效果改变的约束

- 每个关键对象/运动标为 `measured`（来自指定实验数据）、`inferred`（来源支持的解释）、`schematic`（教学构造）之一。实验结构模型仍受分辨率、缺失区域和构建体限制，不能称为完整活细胞真值。
- 坐标、单位、物种、膜内外方向、核酸 5′/3′ 和蛋白 N/C 端先确定。全细胞与分子特写通常用不同场景；不要暗中放大蛋白却保留单一精确比例尺。
- 膜是一道拓扑边界。跨膜必须对应孔道、转位、出芽/融合或膜破裂；糖基化结构域通常处于分泌系统腔面/胞外侧，不能随机翻转。
- 分子层面的剖面优先展示前后不同的完整单体、选择性隐藏遮挡物；不直接切开蛋白、核酸、脂质或受体。保持原坐标，省略实体有记录；分解图的移位明确标为示意。见 section-views。
- 受体密度/聚集不是装饰。RNA 表达量、IF 亮度和蛋白结构不直接给出表面拷贝数或纳米簇；无定量依据时用明确标注的假设分布。
- 动画关键帧不是动力学模拟；播放时长不等于生物学时长。布朗运动需说明扩散模型与步长；路径动画不可冒充实验轨迹或 MD。
- 科学数据图沿用朴素科研风格；Blender 场景可有不同艺术表现，但不能用生图伪造显微图、结构证据或数据图。

## 默认精细制作入口

| 内容 | 脚本与范围 |
|---|---|
| PDB/mmCIF 到可追溯分子资产 | `prepare_molecular_asset.py`，Gemmi 解析、nm 坐标、按链表面/原子数据 |
| 公开结构候选检索 | `search_public_assets.py`：RCSB 标题/全文及可选物种检索，返回候选 metadata，不自动认定适用性 |
| 骨架装配与拆解 | `prepare_filament_assets.py` + `build_polymer_scene.py`，复用 8F8R/8F8Q/7AQK/3J6F；模板、命令与边界见 polymer-templates |
| 肌动蛋白化学状态与完整流程计划检查 | `prepare_actin_states.py` + `build_actin_state_comparison.py` 保留 8A2R/S/T 配体并生成静态对照；`audit_mechanism_plan.py` 检查计划覆盖与闭环，不证明全流程成片完成 |
| 膜受体、转录复合物、运输 | `build_molecular_scene.py`：3SN6/5FLM/4HNA；只有 motor 内置过程动画 |
| 核糖体与复制复合体结构 | `build_complex_scene.py`：9SPI/7PFO 完整装配或完整隐藏部分亚基的内部视图，属于静态结构 |
| 翻译延伸的底物循环 | `build_translation_scene.py`：9SPI 衍生 tRNA/mRNA，示意的进入、短肽转移、三核苷酸转位；完整单体分层展示，非完整催化循环 |
| 转录中的一个核苷酸加入 | `prepare_transcription.py` + `build_transcription_scene.py`：5FLM 产物与 CCD CTP，明确标注的底物接近/连接变化/PPi 离开 |
| 放慢转录并继续下一轮 | `prepare_transcription_next.py` + `extend_transcription_scene.py`：23 秒、连续 C→U 两轮，中间停顿和转位；第二轮构象明确标为教学重构 |
| 复制叉和两条新链 | `prepare_replication.py` + `build_replication_scene.py`：完整核苷酸、连续/分段延伸、方向标记；CMG 完整后移为结构参考层 |
| 线粒体 | `build_organelle_scene.py`：内外膜、连续嵴膜、ATP synthase 示意、局部剖切 |
| 动物细胞总览 | `prepare_cell_geometry.py` + `build_cell_scene.py`：核膜—ER 连通腔体、Golgi、细胞器与骨架 |
| 中期后的分离与胞质分裂 | `prepare_division.py` + `build_division_scene.py`：连续收缩/分开、逐帧网格缓存 |
| 重渲染、外部合成交接、视频 | `render_saved.py`、`encode_video.py`；完整命令与限制见 production-workflow |
| 参考图辅助分析与多风格对照 | `analyze_style_image.py`、`prepare_style_asset.py`、`build_style_scene.py`、`style_nodes.py`；可编辑 GN/材质/合成见 style-library |
| 可移植完整肌动蛋白案例 | `build_actin_films.py`，项目路径、字体、运行时、质量档位、分步制作与续跑见 portable-actin-films |
| 可迁移的蛋白局部相互作用 | `prepare_interaction.py` + `build_interaction_scene.py`，来源结合末态与完整配体、自然接近/转向/释放；不自动完成催化或完整通路 |
| PKAc/SP20 磷酸转移与产物释放 | `prepare_phosphotransfer.py`、`prepare_phosphotransfer_motion.py`、`build_phosphotransfer_scene.py`、`audit_phosphotransfer_scene.py`；完整催化实例，不包含 PKA 上游激活 |
| Ras/SOS 核苷酸交换 | `prepare_ras_switch.py`、`prepare_ras_exchange_motion.py`、`build_ras_exchange_scene.py`、`audit_ras_exchange_scene.py`；旧 GDP 释放、新 GTP 结合及 SOS 释放，不包含 GAP 水解 |
| Ras/p120GAP 水解与关闭 | 同一结构准备后运行 `prepare_ras_hydrolysis_motion.py`、`build_ras_hydrolysis_scene.py`、`audit_ras_hydrolysis_scene.py`；水中氧进入 Pi、GDP/Pi 中间态、Pi 释放和构象复位 |
| 可组合滤镜与受保护标签 | `apply_render_recipe.py` + `compositor_stack.py`，节点链保存在新工程，按所需 pass 和引擎检查 |

精细入口的风格档案为 [style-studio.json](assets/style-studio.json)、[style-molecular-reference.json](assets/style-molecular-reference.json)、[style-organelle.json](assets/style-organelle.json)。它们控制几何表现、材质、灯光、镜头和合成分支。受体示例默认采用用户选定的侧视剖面（`--membrane-view side`）；`--environment` 增加有身份和尺度的前/中/后景。保持主结构清楚，环境数量和分布注明是否为假设。

其他细胞与 60 个过程由 agent 按配方、来源和验收标准制作具体工程；不能把目录覆盖范围称为已经渲染完成的 60 个动画。

## 机制白模与旧兼容入口

以下命令的 `SKILL` 与 `BLENDER` 均指向实际安装位置；在项目新目录输出，避免覆盖已有工程。生成器使用标准 `bpy`，不要求联网或付费插件。默认 seed 固定。

```bash
BLENDER='/actual/path/to/blender'
SKILL='/actual/path/to/blender-cell-biology'
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_scene.py" -- \
  --scene receptors --style teaching --engine auto --compositor clean \
  --output /absolute/project/receptors --preview
```

`--scene` 支持 `cell`、`receptors`、`mitosis`、`replication`、`transcription`、`translation`、`vesicle`、`structure`。它们是可编辑起点，**不是所有过程的完整自动模型**。`mitosis` 从中期开始；`vesicle` 只展示沿微管运输，不冒充已实现膜融合。`--style realistic` 只是材质与灯光的另一种起点，真实环境的结构和拥挤程度仍需参考数据。

```bash
# 保存工程并生成全部 PNG 帧；再使用 FFmpeg 编码
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_scene.py" -- \
  --scene vesicle --frames 96 --fps 24 --output /absolute/project/vesicle --animate

python3 "$SKILL/scripts/encode_video.py" /absolute/project/vesicle --fps 24

# 下载实验结构及元数据；不自动假定 asymmetric unit 就是功能装配体
python3 "$SKILL/scripts/fetch_structure.py" 3SN6 --output /absolute/project/structures
```

这些旧入口用于快速核对机制，不用于最终视觉交付。运行细节见 [blender-workflow.md](references/blender-workflow.md)。旧结构示例使用随 skill 附带的公开 1BNA 坐标及 [来源记录](assets/1bna-source.json)。

## 完成标准

交付 `.blend`、可重建脚本/参数、预览图或可播放视频，以及简洁的来源/简化/验证记录。科学正确性、画面可读性和软件运行成功分别报告；漂亮画面不能代替前两项检查。报告实际覆盖的内容与尚未执行的配方，不用“全部支持”掩盖尚需实现的过程。
