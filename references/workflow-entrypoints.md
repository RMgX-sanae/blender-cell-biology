# 按任务选择入口与验证范围

先确定交付是静态结构、局部动作还是完整过程，再选下面的入口。新蛋白或新机制接着读 [workflow-transfer.md](workflow-transfer.md)，只替换必要的数据和科学规则；已有源构象、反应位点或路径不会因更换名称而自动适用。

## 怎样读验证状态

- **独立重建**：指定案例曾从解压后的 skill 和原始输入重建，并检查所列工程、图片或视频；只覆盖该案例、版本和输出条件。
- **局部验证**：所列配置、阶段或显示方法已实际检查，不意味着整个生物过程或全部风格组合通过。
- **既有片段/模型**：已有生成入口和制作样例；该段的范围较窄，此次没有重新完整验收。
- **配方**：说明如何制作，仍需项目实现。不得将其计为已完成动画。

这些状态描述软件与交付证据。来源结构、推定运动和教学简化另按各案例说明标注；独立重建也不等于实验动力学验证。以下整理依据截至 2026-10-29 的制作与复验记录，2026-11-02 核对入口。后续改动影响的案例应重新验证，不能沿用旧状态作保证。

## 分子动作与过程

| 用户需要 | 首选入口与操作 | 当前验证范围 | 需要另做的部分 |
|---|---|---|---|
| 换蛋白/肽段，展示结合、转向、落位和释放 | [局部相互作用](protein-interactions.md)：项目 JSON → `prepare_interaction.py` → `build_interaction_scene.py` | **局部验证**：Barnase/Barstar、MDM2/p53、trypsin/BPTI；外部配置调用 MDM2/p53 的保存工程和关键帧已检查 | 新目标的装配、配体和进入路径；催化不能由结合动画代替 |
| 激酶的磷酸转移与产物离开 | [PKAc/SP20](phosphotransfer.md)：`prepare_phosphotransfer.py` → `prepare_phosphotransfer_motion.py` → `build_phosphotransfer_scene.py` | **独立重建**：24 秒局部催化示例、原子身份及保存工程 | PKA 上游激活；其他激酶的底物、反应端态和运动 |
| 核苷酸交换、水解或两段接续 | [Ras](ras-switch.md)：选择 exchange、hydrolysis 或 cycle 准备入口，使用对应场景与检查器 | **独立重建**：Ras/SOS 与 Ras/p120GAP 的限定示例，包含 48 秒连续状态序列 | 其他 GEF/GAP；整条信号通路；开放物质池不能称为封闭回收循环 |
| 泛素连接及双泛素释放 | [K63 连接](ubiquitin-transfer.md)：`prepare_ubiquitin_transfer.py` → `build_ubiquitin_transfer_scene.py`；释放由准备阶段的 `--release` 选择 | **独立重建**：预装载体系的净键转移和可选释放 | E1 装载、完整泛素回收、蛋白酶体降解；其他链型需另建规则 |
| 从连接流程迁移到共价键水解 | [AMSH-LP](dub-workflow-transfer.md)：`prepare_dub_hydrolysis.py` → `build_dub_hydrolysis_scene.py` | **独立重建**：30 秒已结合底物→净水解→产物离开；水氧身份和受限末端路径 | 底物进入、质子接力、溶剂重新装载和完整酶循环 |
| Arp2/3 或 formin 的上游到回收 | [完整 actin 案例](portable-actin-films.md)：`build_actin_films.py` 按 prepare/build/render/encode/audit 执行 | **独立重建**：两条指定教学路线、21 个独立章节 | 新细胞/调控情境；来源未解析的成核中间态仍是教学构造 |
| 局部丝端增长、拆解、封端或微管变化 | [骨架模板](polymer-templates.md)：`prepare_filament_assets.py` → `build_polymer_scene.py` | **既有片段/模型**：11 个可执行模板；逐项状态以具体项目为准 | 不能将 actin 完整案例的验收外推到全部模板或微管过程 |
| V 型 ATPase 的连续泵送 | [V-ATPase](vatpase-loop.md)：`prepare_vatpase.py` → `build_vatpase_scene.py` → `audit_vatpase.py` | **独立重建**：16 秒双视图循环；完整低分辨率周期及高分辨率关键帧 | 全部催化构象、质子水链、累积 pH 和膜电位 |
| 转录中的核苷酸加入 | [精细制作流程](production-workflow.md) 的转录段：`prepare_transcription.py`、`build_transcription_scene.py`；第二轮用 `prepare_transcription_next.py`、`extend_transcription_scene.py` | **既有片段**：23 秒 C→U 两轮教学示例 | 起始、终止和完整天然中间态；第二轮构象是教学重构 |
| 翻译的 tRNA 进入、肽转移和转位 | [精细制作流程](production-workflow.md) 的翻译段：`build_translation_scene.py` → `audit_translation.py` | **既有片段**：9SPI 衍生的局部底物运动 | 起始、终止、eEF1A/eEF2、GTP 水解和完整催化循环 |
| DNA 复制叉与新链方向 | [精细制作流程](production-workflow.md) 的复制段：`prepare_replication.py` → `build_replication_scene.py` | **既有片段**：连续/分段延伸示意，完整 CMG 位于结构参考层 | 起始、引物处理、连接与终止；不能当作完整复制机制 |
| 中期后的染色体分离与胞质分裂 | [精细制作流程](production-workflow.md) 的分裂段：`prepare_division.py` → `build_division_scene.py` | **既有片段**：四条已复制染色体、连续收缩和整体分开 | 前期、核膜重建、中体、ESCRT 与分子级切断 |
| 囊泡出芽/切断/识别/融合，或其他未列出的全过程 | [过程目录](process-atlas.md)、[装配分镜](assembly-storyboards.md) → [迁移方法](workflow-transfer.md) | **配方**：运输白模只覆盖沿轨道移动 | 需实现改变的膜连通性、区室侧别和物质交接；不要直接复用刚体进出作为膜融合 |

要做连续循环时，先确认本例能否闭合对象状态、物质池和运动接缝。教学片段的末帧跳回首帧并不会自动成为合格循环；完整 actin 案例的分章剪辑也不等于无缝循环。

## 静态结构与空间分布

| 用户需要 | 首选入口 | 当前验证范围与边界 |
|---|---|---|
| 从 PDB/mmCIF 准备可追溯资产 | [公共资产](public-assets-first.md)、[数据导入](data-import.md)：`fetch_structure.py`、`search_public_assets.py`、`prepare_molecular_asset.py` | 多个案例共用并已实际重建。先核对链、配体、单位和功能装配；`generic` 默认过滤不适用于所有反应物 |
| 受体膜环境或静态转录复合物 | [精细制作流程](production-workflow.md)：`build_molecular_scene.py` 的 membrane/transcription | **既有模型**：分别对应 3SN6/5FLM；同脚本 motor 是运输示意，不包含出芽或融合 |
| 曲面上的受体朝向、簇或极性分布 | [受体分布](receptors-and-distributions.md)：`prepare_receptor_distribution.py` → `build_receptor_distribution_scene.py` | **独立重建**：完整 3SN6 复合物、三种静态分布和六种外观检查；密度/聚集是假设，未验证生理丰度或扩散 |
| 完整核糖体/复制复合体及内部观察 | [完整单体观察](section-views.md)、[精细制作流程](production-workflow.md)：`build_complex_scene.py`、`audit_intact_molecules.py` | **既有模型**：9SPI/7PFO；隐藏完整遮挡亚基，不能切开分子，也不能把结构视图称作催化动画 |
| 线粒体或动物细胞总览 | [细胞与细胞器](cells-and-organelles.md)、[精细制作流程](production-workflow.md)：`build_organelle_scene.py`；或 `prepare_cell_geometry.py` → `build_cell_scene.py` | **既有模型**：连续嵴膜、核膜—ER 关系等示意；20 类细胞与 27 类组件表仍是建模配方，不能当作全部已重建 |
| 新细胞形态、细胞器定位或组织极性 | [细胞与细胞器](cells-and-organelles.md) → [迁移方法](workflow-transfer.md) | **配方**：改用目标物种/细胞状态的来源与空间约束；受体表面采样不能直接替代细胞器体积分布 |

## 表现与交付沿用哪些步骤

| 变化 | 沿用入口 | 验证边界 |
|---|---|---|
| 按参考图片换表现 | [看图实施](style-from-reference.md)、[风格库](style-library.md)：`analyze_style_image.py` 辅助统计，agent 观察后调整实际材质、光、镜头和节点 | 已有实例和 18 类可执行起点；统计色板不等于自动恢复材质，不保证任意图片都已验证 |
| 换细胞/细胞器/膜/分子背景 | [分尺度环境](cellular-environment.md)：`apply_scale_environment.py` 及所选资产准备入口 | 细胞与 ER 环境已有独立重建小样；四种景别不是同一细胞的定量连续放大，原生数据与示意三维深度要分开 |
| 换滤镜或保护主体/文字 | [合成配方](render-recipes.md)：`apply_render_recipe.py`、`compositor_stack.py` | Blender 原生可编辑节点已运行；每个新场景仍需检查身份颜色、文字和 pass。Fusion/Natron/AE 目前是交接说明 |
| 换机器或 Python 环境 | [运行检查](production-workflow.md)：`check_runtime.py` | Eevee/Cycles 的资产、双字体、几何、合成、三帧编码/解码已检查；不认证新机制或全部质量档位 |
| 重开、续渲染、编码和检查 | `render_saved.py`、`encode_video.py`、`audit_video.py`；机制检查器见对应案例 | 共用运行步骤已使用；视频可解码不代表机制或视觉合格。外部合成必须实际执行后才能声明完成 |

## 已有验证记录摘要

下列是上述状态的依据摘要，不是要求每次任务都重跑的检查清单。完整制作记录和生成媒体保留在项目目录，skill 只分发生成方法。用户对旧片段的满意反馈不自动替代新项目验收。

- **配置迁移**：MDM2/p53 以项目自己的 ID 和 JSON 调用原生成器；818 个原子、193 个含接缝时点的保存工程检查和关键帧复核。不新增完整影片的验证声明。
- **反应身份共用**：PKA、Ras 水解、K63 接入同一账本后，22 个 NPZ 中的 107 个数组保持不变；六类错误账本被拒绝。AMSH-LP 后续复用该入口并单独定义水氧和断键/成键。
- **完整 actin 案例**：21 章来源资产、保存几何、字体与路径重建；两条影片完整解码并抽查实际画面。适用范围仍为各章声明的教学情境。
- **K63→AMSH-LP**：重新准备的 24 个数组与复核项目一致；721 个准备姿态、两种风格各 720 个保存帧、六个反例和 30 秒完整视频检查；来源几何及远端刚体释放的失败案例没有被略过。
- **静态受体分布**：均匀、成簇、极性三种假设分布，完整复合物和膜叶朝向检查；六个所选外观组合。没有把静态分布升级为扩散验证。
- **共用运行环境与采样**：Eevee/Cycles 实际应用指定采样数；原遗漏采样设置会被检查识别。AMSH-LP 显式 64 次采样重建的四张关键帧与已验收像素一致。

新增结果只更新它实际覆盖的行。例如 DUB 片段通过，应更新 DUB 的范围；不能因此将“泛素—蛋白酶体降解全过程”或其他 60 项目录标成完成。
