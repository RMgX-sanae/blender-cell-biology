# 把已验证工作流用于另一种机制

扩展这个 skill 的首要产物是可执行的迁移方法。新任务先选一个已验证案例作为起点，复用其中适用的资产处理、运动、显示和检查步骤；新增视频只用于回答尚未验证的问题。过程目录用来查找需求，不以条目数或视频数衡量完成度。

先用 [入口与验证范围](workflow-entrypoints.md) 区分可运行案例、已验收范围和仍需项目实现的配方；不要仅按脚本名字选择。

## 先选迁移层级

| 新任务改变了什么 | 优先做法 | 必须重新判断的内容 |
|---|---|---|
| 同一类结合动作，换蛋白、肽段或抑制剂 | 用项目配置调用现有结合生成器 | 功能装配、作者链、配体、源结合状态、进入路径是否能通过 |
| 同一对象，换画面风格或动画节奏 | 保留结构和身份，调整对应风格/时序入口 | 新滤镜的颜色辨认、标签可读性、重新定时后的全部姿态和接缝 |
| 发生另一种化学变化 | 复用资产、运动和显示工具，另写本机制的反应准备与检查 | 底物/产物、原子映射、断键/成键、末端化学、金属与必要的中间态 |
| 将已有阶段接成长过程 | 复用阶段输出，在交接处对齐同一批对象 | 对象身份、坐标、配体状态、区室、位点占有及开放物质池 |
| 新增膜拓扑、晶格或运动约束 | 复用来源处理和渲染流程，单独实现变化的几何规则 | 膜侧别/连通性、轨道或聚合端、单体去向；刚体移动不能代替融合或切断 |

不需要把每个任务都升级为新框架。先使用已有入口；只有当前目标确实需要的新科学内容才新增实现。

用户要求默认联网，在只有单体结构时仍尽力探索可能的结合方式。先按 [联网结合研究](online-binding-research.md) 主动核对可用复合物、位点对应和组装状态；适用实验复合物可以直接复用，需要组合单体时再走 [单体结合候选](monomer-binding.md)。明确保留单体来源、配对映射、构象限制和筛选理由。已有 `prepare_binding_candidates.py` 将模板迁移与约束初始化接入共用运动/场景，并区分实验复合物、模板迁移、约束候选和探索性摆放。几何接触、分数较优或运动能通过都不等于确认结合，更不能据此编造亲和力；当前是刚体候选工作流，不是通用盲对接预测器。

## 现有代码中可以直接沿用的部分

| 已验证工作流中的部分 | 现有入口 | 迁移边界 |
|---|---|---|
| 分子视频的共用运行链 | `check_runtime.py` | 在新环境检查实际资产准备、Blender 读取/渲染及编码/解码；不认证机制、特定滤镜或项目额外依赖，见 production-workflow |
| 来源下载与留痕 | `fetch_structure.py`、`search_public_assets.py` | 下载成功不等于选中了正确构建体；继续核对状态、突变、缺失区和装配 |
| 按选定链准备资产 | `prepare_molecular_asset.py --config` 及 `atom_table`、`xyz`、`surface`，见 [项目配置](molecular-asset-config.md) | 可替换聚合物链/残基和有依据的刚体坐标系；使用者另行选择必要配体，不得因默认过滤而丢失反应物 |
| 从不对称单元恢复源装配 | `prepare_clathrin_assembly.py` 的 `rigid`、`compare_reference`，见 [clathrin](clathrin-assembly.md) | 复用刚体与逐身份/坐标对应检查；源单体分组及实际装配操作另行资格确认。静态重建不证明生长轨迹 |
| 结合/释放弧线 | `interaction_motion.py` 的 `plan`、`matrices`、`contacts` | `plan` 面向给定源结合态的刚体伙伴；`matrices` 只给位姿，不自动保证无穿插或催化正确 |
| 单体起步的结合候选 | `prepare_binding_candidates.py` → 共用结合场景，见 [单体候选](monomer-binding.md) | 明确 CA 映射或界面约束；候选端态及路径不允许来源接触豁免。保留未满足的侧链诊断，几何筛查不是结合预测准确率 |
| 准备数据与工程的对应关系 | `interaction_assets.py` 的 `verify_assets`，由局部结合构建器/检查器调用 | 核对全部实体和运动文件，再将 manifest 指纹绑定到新工程；改过坐标后必须重新准备和检查路径。接口采用局部结合 manifest，其他机制应适配其实际文件清单 |
| 反应端态身份与净键变化 | `reaction_ledger.py` 的 `validate` 或命令行 | PKA、Ras 水解、K63、AMSH-LP 及两轮转录准备共用；各自定义反应，只检查输入原子和键子图，不代替几何检查 |
| 同一批原子的连续反应 | `reaction_ledger.py` 的 `validate_sequence` 或 `--sequence` | 核对阶段交接和独立定义的总反应；固定所选原子集合，不自动处理开放物质池或坐标交接 |
| 原子/网格动画、阶段可见性与视图层 | `build_phosphotransfer_scene.py` 的 `animated_mesh`、`visibility`、`view_layer` | Ras 场景已直接调用这些函数；它们不包含通用反应化学或通用坐标文件格式 |
| 整分子双视图与分层背景 | `build_ubiquitin_transfer_scene.py` 的 `keyed_points`、`envelope`、`background`、`paired_composite` | K63 与 AMSH-LP 共用；背景换入目标蛋白形状，化学和相机取景仍按项目定义 |
| 已检查坐标序列的项目显示 | `build_molecular_sequence_scene.py`、`audit_molecular_sequence_scene.py`，见 [输入约定](molecular-sequence-display.md) | PTP1B 与 AMSH-LP 共享两图原子、标记、背景和保存检查；新机制仍提供经过检查的 nm 坐标、身份、反应账本和显示配置 |
| 字体、材质、环境和合成 | `detail_scene.py`、`style_nodes.py`、`cellular_environment.py`、`compositor_stack.py` | 沿用清楚的主反应区、宋体/Times New Roman 和独立标签层；环境身份与尺度按新区室改选 |
| 对已有画面追加调色 | `apply_render_recipe.py`、`publish_composite_outputs`，见 [合成输出约定](render-recipes.md) | 保留已合成的分屏/背景及独立文字；多视图不能直接借用一个视图的法线/深度，缺少对齐 pass 时拒绝对应配方 |
| 双视图分别处理深度、轮廓与对象强调 | `panel_compositor.py`、`publish_panel_sources`，见 [各视图通道](render-recipes.md) | PTP1B/AMSH-LP 已重建检查；先用各自 Render Layers 的 pass 和原 alpha 处理，再回到原布局并叠加文字。明确每个视图的深度范围与主体选择；不改变反应几何或扩充机制 |
| 保存后渲染与视频检查 | `render_saved.py`、`encode_video.py`、`audit_video.py` | 可以共用执行步骤；场景内的机制检查仍按本例定义，成功解码不能证明科学正确 |

复用这些脚本时把实际 skill 的 `scripts` 目录加入导入路径，不绑定作者机器的绝对路径。先查看函数签名和一个现有调用点，保持 Å/nm、帧数及数组形状一致。不要假定不同机制的 manifest 可以直接互换。

## Agent 实际执行的五步

1. **选择起点并说明相似性。** 根据发生的事件选择，而不是根据蛋白名字：只结合用局部相互作用；更换配体可参考 Ras 交换；化学转移可参考 PKA 或 K63；阶段组合参考 Ras 状态序列；装配与回收参考 actin。
2. **在项目中写一张简短迁移说明。** 列出目标及范围、沿用的脚本/函数、替换的数据、必须新建的科学规则、通过标准。已有来源与结果用文件链接引用，不复制整套旧项目。
3. **先替换最小输入。** 同类结合先试项目 JSON；化学改变先建立可核对的反应端态和身份账本；阶段组合先核对交接。只有这些输入通过后才生成全部运动。
4. **先验证变化最大的部分。** 如新界面的进出路径、另一反应的键变化、另一膜曲率上的朝向。检查每个输出姿态，并加入能击中本次风险的错误输入；不要为了增加测试数量重复检查已确定的内容。
5. **复用显示与交付。** 生成少量关键帧检查主结构、双视图和遮挡；机制通过后再做所需完整视频。从解压的 skill 在项目目录重建，记录实际复用和新增的部分。只在代表性验证通过后把方法写回 skill。

每次迁移记录两个结论：本目标是否通过；哪些步骤可以继续用于下一目标。一个新蛋白成功，不足以证明整类酶都能自动处理。

## 可直接运行的同类迁移

局部结合入口已经把生物选择与路径算法分开。下例从已核对的配置中取 MDM2–p53，放入项目自己的文件，并使用不同项目 ID；生成器及路径代码保持不变。它检验的是从蛋白抑制剂结合流程迁移到结构域–肽段识别的配置接口，不增加一个新机制样例。

```bash
# SKILL、PROJECT、PYTHON 均为当前机器的实际路径；PROJECT 已存在。
"$PYTHON" - "$SKILL" "$PROJECT" <<'PY'
import json, sys
from pathlib import Path
skill, project = map(Path, sys.argv[1:])
entries = json.loads((skill/'assets/interaction-examples.json').read_text())['examples']
entry = next(e for e in entries if e['id'] == 'mdm2-p53').copy()
entry['id'] = 'project-domain-motif'
path = project/'interaction.json'
with path.open('x') as stream:
    json.dump({'version': 1, 'examples': [entry]}, stream, indent=2)
PY

# 先复用来源哈希相符的缓存；没有缓存时再运行下载命令。
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 1YCR --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_interaction.py" \
  --examples "$PROJECT/interaction.json" --example project-domain-motif \
  --sources "$PROJECT/sources" --output "$PROJECT/assets" \
  --mode loop --frames 192 --fps 24
```

接着使用 [局部结合工作流](protein-interactions.md) 中原有的场景构建、保存后检查及渲染命令。新蛋白需要经来源审查后替换 `pdb`、`anchor`、`moving`、配体选择、标签和构建体边界；无需为了新名字复制 `prepare_interaction.py`。`ligands` 中请求但未找到的组分会失败；需要区分具体副本或独立作者链时用 [ligand_instances](protein-interactions.md#明确保留哪一份配体)，写明其刚体伙伴和选择依据。路径失败时查看实际失败原子对，再决定是否需要柔性或其他源构象。

## 可复用的反应身份检查

`prepare_phosphotransfer.py`、`prepare_ras_hydrolysis_motion.py` 和 `prepare_ubiquitin_transfer.py` 在准备时调用同一个检查器，并生成项目文件 `reaction-ledger.json` 和 manifest 内的 `reaction_ledger` 检查结果。PKA 的转移基团、Ras 的水氧归属、K63 的连接位点分别由各自机制定义；通用检查器没有内置这些反应。

```bash
"$PYTHON" "$SKILL/scripts/reaction_ledger.py" "$PROJECT/assets/reaction-ledger.json" \
  --output "$PROJECT/ledger-audit.json"
```

新机制在自己的准备脚本中构造相同格式的字典，调用 `from reaction_ledger import validate`；让检查通过后再生成运动。字段约定如下：

| 字段 | 内容 |
|---|---|
| `version`、`scope`、`bond_scope` | 版本 `1`、本账本包括/省略的对象、检查哪些键；明确重原子或含氢范围 |
| `reactant_atoms`、`product_atoms` | 分别从端态资产提取的 `[唯一原子ID, 元素]` 列表；产物中离开的分子和离子也要列入所选物质范围 |
| `atom_map` | `[反应物ID, 产物ID]` 一一映射；允许产物名称改变，元素必须相同，两个端态的所列原子都必须覆盖 |
| `reactant_bonds`、`product_bonds` | 各自端态 ID 下的无向共价键列表；只列经过判断的子图，写明它不是完整分子图 |
| `break`、`form` | 预期净断键/成键，**均用反应物 ID**；由机制依据独立确定，与端态键差核对，不能直接复制计算出的键差作为预期值 |

这个检查能发现“总数相同但原子重复”、映射遗漏、元素被换成另一种元素、错误反应位点和未声明的额外键变化。它不能发现两个端态同时遗漏了同一批必要原子，也不能判断一个自洽但生物学错误的声明。PKA 与 Ras 仅覆盖反应物/产物重原子及金属；K63 覆盖所选结构的重原子；三者都只提供反应中心键子图。源结构适用性、质子/电荷、键级、全分子价态、中间姿态和实际画面仍需各自检查。

## 将同一批原子交给下一反应阶段

单个阶段通过检查，不能证明下一阶段沿用了它的产物。`validate_sequence` 接受 `{"version": 1, "stages": [第一阶段账本, 第二阶段账本], "net": 总反应账本}`；每份账本均使用上表格式，可有更多阶段。总反应的预期映射和断键/成键须由目标机制独立定义，不能从阶段结果复制出来充当验证。

```bash
"$PYTHON" "$SKILL/scripts/reaction_ledger.py" "$PROJECT/reaction-sequence.json" \
  --sequence --output "$PROJECT/sequence-audit.json"
```

检查器逐一验证阶段，再要求上一阶段产物与下一阶段反应物的原子 ID、元素和所列键完全相同；最后将逐段映射串起来，与独立总反应的端态和映射核对。阶段内部允许有明确映射的名称改变，交接边界必须使用相同 ID。共价中间体可以形成后再消失，净反应因此没有这条键，但中间的原子不能被替换。

这里采用**固定的所选原子集合**。等待反应的水、已离开的产物也要在每段账本中保留身份和相应键；离开画面不等于删除物质。选择只检查反应中心时，明确其他原子未被验证。不同范围的局部账本（如下面的两轮转录）不能直接拼接；需要先定义持续状态及进入/离开的物质池，而不是取交集忽略缺失对象。这个入口不检查坐标连续、质子/电荷、开放池交换，也不证明催化循环已闭合。

最初用 PTP1B 两步净反应的 **23 重原子局部图原型**验证交接；反例中每步都能通过，但改换交接原子、改变中间连接或交换总反应的氧身份时，序列检查会失败。反应依据为 [Brandão 等的原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC2871455/)，组分内连接参考 [RCSB CCD](https://www.wwpdb.org/data/ccd)。随后 [PTP1B 几何迁移](ptp1b-workflow-transfer.md)以持续的 3I7Z 构建体扩展到 2,482 个所选重原子，另行处理局部中间体、Q262 扭转、水位置和碰撞；不能把早先的图检查当作几何已经通过。两版仍不含质子/电荷或完整酶循环。

## 把身份检查用于核苷酸加入

转录的 `prepare_transcription.py` 与 `prepare_transcription_next.py` 复用同一账本检查器，`addition_ledger` 则定义本例的 RNA 链端、NTP 核心和 PPi 规则。方法的迁移发生在原子对应和图检查层；底物选择、断键/成键及来源端态仍是本机制自己的内容。

- [5FLM](https://www.rcsb.org/structure/5FLM) 提供静态牛 Pol II 与合成核酸支架及 C20 产物坐标；[CCD CTP](https://www.rcsb.org/ligand/CTP)、[UTP](https://www.rcsb.org/ligand/UTP) 提供各自理想组分与连接。进入姿态、U21 和转位仍是教学构造。
- RNA 3′ 端进攻与 PPi 释放的机制依据见 [SSRL 对该原始研究的说明](https://www-ssrl.slac.stanford.edu/research/highlights_archive/pol-II.pdf)，原论文为 [Wang 等，2006](https://doi.org/10.1016/j.cell.2006.11.023)。这里只采用净连接规则，不重建 trigger loop、质子转移或金属催化，也不把论文的酵母构象搬到牛 Pol II 模型上。
- 第一轮账本列出 P19 与 CTP 的 52 个重原子，第二轮列出 C20 与 UTP 的 49 个；同一个 `RNA:P:20:*` ID 在第一轮产物和第二轮链端中延续。新底物与 PPi 的 ID 带目标残基位置，以免把两轮游离分子混成一份。
- 两个账本各自只覆盖一个链端和一个底物，早先 RNA、第一轮已离开的 PPi 及细胞物质池不被第二轮重复计入。需要更长的序列或完整回收时，另建这些对象的持续状态和去向。

`prepare_transcription_next.py` 另用 `addition_sequence` 生成 `reaction-sequence.json`，把**现有两轮**扩展为同一套 81 个所选重原子。原来的 52/49 原子局部账本仍保留。扩展不是把局部账本取交集，而是显式补入等待/离开的对象和不参与本轮反应的连接：

| 时点 | 同一套所选原子的状态 | 数量 |
|---|---|---|
| 第一轮开始 | RNA P 链第 19 位的 23 个原子 + CTP 29 个 + 等待的 UTP 29 个 | 81 |
| 两轮交接 | 第 19–20 位的 43 个原子 + 第一份 PPi 9 个 + 同一个 UTP 29 个 | 81 |
| 第二轮结束 | 第 19–21 位的 63 个原子 + 两份具有不同身份的 PPi 各 9 个 | 81 |

第二轮同时保留第 19–20 位已经形成的 RNA 连接，不能因为它落在当前局部反应范围之外就删除。总反应的端态和映射从输入的 RNA/CTP/UTP 数组单独构造，预期两次 α 磷酸断键和两次 3′ 加入；不把阶段输出的键差复制成“独立”预期。第二次净成键的亲核氧在**整个序列初态**仍属于 CTP，因此总反应使用 `CTP_for_P20:O3'`，不能误用尚未形成的 `RNA:P:20:O3'` 初态 ID。

```bash
# prepare_transcription_next.py 完成后，manifest 的 reaction_sequence 记录其文件和哈希。
"$PYTHON" "$SKILL/scripts/reaction_ledger.py" "$PROJECT/two-cycles/reaction-sequence.json" \
  --sequence --output "$PROJECT/two-cycles-sequence-audit.json"
```

这是有限的、明确选定的图范围，不包含第 19 位之前的 RNA、聚合酶、模板、金属、溶剂或氢；不新增进入/终止/回收的机制，也不验证无限开放物质池。原有坐标和动画不因此改变。新目标应从自己的输入建立同样明确的等待/释放状态及独立净反应，而不是复制本例的 81 原子或 C→U 规则。

更换核苷酸或聚合酶时，先核对可反应的末端原子、模板与底物、α 磷酸/离去基团的命名及真实/构造产物。重新得到目标端态后再沿用检查器，不能把“总数仍是 29→20+9”当作原子身份正确的充分证据。

## 跨机制迁移：保留方法，重建改变的科学部分

| 起点 → 目标 | 可以保留 | 需要替换或新增 | 有效的验收问题 |
|---|---|---|---|
| PKA → 另一激酶 | 供体/产物身份账本的方法、可动结构域处理、同步双视图与合成 | 新酶与底物来源、反应原子、必要配体、开合和路径 | 被转移的原子是否来自正确供体，产物是否守恒，蛋白是否被拉伸？ |
| Ras → 另一种小 GTPase | 交换与水解分开的事件设计、旧/新核苷酸身份、阶段交接检查 | 该 GTPase 的 GEF/GAP、端态构象、反应位点与配体映射 | 是否错误地把旧 GDP 当成新 GTP，或直接沿用 SOS 的作用方式？ |
| K63 连接 → 断开共价连接的过程 | 共价图、端态几何和整分子释放的检查方法 | 新反应的水/底物/产物账本、切割位点和必要中间态 | 是保持身份的化学转化，还是倒放旧动画或直接删除对象？ |
| actin → 其他装配过程 | 单体身份、位点占有、接近/释放、回收账本和阶段组织 | 目标的装配界面、晶格、调控因子、化学状态与端行为 | 有无单体凭空出现、穿过邻居、错误占位或未交代回收？ |
| 受体分布 → 其他膜蛋白 | 曲面法向、排除区域、完整复合物间距与分层环境 | 新蛋白的膜锚定、拓扑朝向、尺度和分布依据 | 改变膜曲率后是否翻面或穿插，分布是否仍被正确标为假设？ |

只换膜形态时，先保留已核对的受体资产，用项目 NPZ 替换 `surface`，沿用 [受体分布](receptors-and-distributions.md) 的原准备、场景与保存检查入口。把世界变换应用到膜顶点，受体仍只做刚体变换；封闭膜启用 `require_closed`。用整体反面、局部反面、缺面和镜像输入验证错误能否被发现。几何检查通过仍不表示形态来自显微数据，也不构成新受体机制的验证。

这些是迁移路线，不是目标已实现的声明。去磷酸化需要单独定义水和产物；囊泡融合需要膜拓扑变化；马达需要轨道接触与能量状态。遇到这类变化时新增相应机制步骤，同时继续沿用已有的来源、显示和交付流程。

## 代表性验收

改变教学节奏时可复用 `event_timing.phase_progress(time, intervals, before=...)`：它检查 0–1 时段和调用方声明的先后关系，再给出各事件的进度。Ras 水解和 AMSH-LP 已使用同一入口，但各自决定什么必须先发生；它不推断生物因果关系。`continuous_intervals(intervals, order)` 仅适用于依次发生的事件，按活动时长比例移除观察停顿；不能压缩本来需要重叠耦合的事件，不能把有限反应片段变成闭合循环。AMSH-LP 的字幕由同一时序派生，具体用法见 [两种节奏](dub-workflow-transfer.md)。

两者还共用 `frame_window(start, end, frames)`：按实际 `(frame-1)/frames` 采样把半开时段映射到 Blender 帧号，避免简单四舍五入让字幕比事件早一帧。短于一帧且没有落入采样点的时段返回空帧窗，不能凭空增加一帧。它适用于这两个入口的采样约定；其他输入帧、实测时间或包含重复终点的片段必须先说明自己的时间约定。

用少量互补情况验证方法：一个只换项目配置的同类目标、一个确实改变反应规则的目标、一个会使旧假设失效的反例。已经存在的 PKA→Ras 显示模块复用、Ras 两阶段组合、以及 K63 对 `interaction_motion.matrices` 的调用，可作为跨机制复用的代码实例；它们各自的科学准备和检查仍独立保留。

检查材料应回答：复用了哪些确切入口；新增的代码为何必要；错误输入能否被检出；保存工程与实际画面是否一致；目标的适用边界是什么。尚未实现的过程继续明确标注，但不以重复制作整段视频替代对迁移能力的验证。

## 改变反应规则的实际迁移

[AMSH-LP 水解迁移案例](dub-workflow-transfer.md)从 K63 工作流沿用原子处理、绕键转动、刚体运动、双视图、背景和合成器；新建水解端态、水氧归属与受限末端释放。它还保留了两个有用反例：源晶体的连接几何不能直接当作天然反应物，能通过的近端刚体释放不能自动用于远端通道。这个案例覆盖已结合底物到产物离开；底物进入和溶剂重装载没有实现。
