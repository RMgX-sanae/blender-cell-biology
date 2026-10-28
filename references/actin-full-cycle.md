# 肌动蛋白：Arp2/3 分支网络与 formin 直丝的完整过程

本说明把旧延伸短片扩展为两条独立的制作路线，覆盖上游调控、单体供给、成核、延伸、化学老化、功能输出、终止和回收。两套 JSON 本身是**有来源的分镜与状态模板**，不能传给 `build_polymer_scene.py`。完整案例的独立生成入口现见 [portable-actin-films.md](portable-actin-films.md)，与核苷酸状态对照分别检查，不把计划字段检查当作成片验收。

读取顺序：通用 [mechanistic-depth.md](mechanistic-depth.md) → 本文件 → 选择一套计划。需要两条路线时分别实施，共用素材仍核对构象与配体。

- [Arp2/3 全周期计划](../assets/mechanism-templates/actin-arp23-cycle.json)
- [mDia1/formin 全周期计划](../assets/mechanism-templates/actin-formin-cycle.json)
- [来源索引](../assets/mechanism-templates/actin-sources.json) 与 [资产适用性表](../assets/mechanism-templates/actin-assets.json)

## 共同边界与单体供给

默认是哺乳动物非肌细胞的机制教学综合模型，并非某一种细胞的完整实验重建。具体受体、GEF/GAP、formin、actin 亚型及刺激时程按项目选；来源结构中兔骨骼肌 actin、鼠蛋白、人蛋白或酵母构建体的身份不能被标签抹去。

起始先展示单体池和已有细丝的来源。Tβ4 参与单体封存，profilin 参与单体供给与核苷酸交换；不能把它们都画成永久黏在新生丝上的装饰。Tβ4/actin 结构中存在工程化融合或杂合构建体，不能当作天然全长人复合物。[原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC4217450/)

跟踪一个亚基 `actin_001`：G–ADP → 释放 ADP 并结合池中 ATP → G–ATP → 加入丝 → F–ATP → F–ADP·Pi → F–ADP + 游离 Pi → 经拆解回到 G–ADP → 再交换并进入第二轮。这是选择的一条示范闭环，**不表示 ADP-actin 不能聚合**。ATP 水解、Pi 释放、核苷酸交换是三件事；不能在单体上把 ADP 直接“加磷酸变 ATP”，也不能把 actin ATPase 的产物画成 PPi。[装配与老化](https://www.nature.com/articles/s41586-022-05241-8)、[Pi 释放机制研究](https://www.nature.com/articles/s41594-023-01101-9)

## 路线 A：Rac1–WAVE–Arp2/3

| 镜头组 | 因果内容 | 实施要点 |
|---|---|---|
| 上游与定位 | 项目选定刺激 → Rac 调控 → 膜上 Rac1–GTP → WRC 活化 | GDP 离开/GTP 进入与 GTP 水解分别记录；胞质侧膜锚定和膜外配体不能混淆 |
| 原料与母丝 | 单体准备；已有母丝；WCA、actin 与 Arp2/3 的协同 | 起始母丝是边界条件，要解释来源；从无细丝开始的项目需要另补种子成核 |
| 分支成核 | Arp2/3 活化、母丝结合和子丝起始 | 七亚基复合物完整；固定成熟接头再延长不算重新成核 |
| 延伸与输出 | 子丝倒钩端至少两轮增长，分支网络与膜前缘关系 | 局部分支角用选定结构；网络方向/密度和推动膜的位移另有来源或标示意 |
| 调控与周转 | 封端、稳定/去分支、老化丝的 cofilin 作用、片段拆解 | 同一网络内可并行；切断普通细丝和拆开分支接口不是同一动作 |
| 回收 | CAP/cofilin 相关单体回收、交换、再供料；信号关闭 | 追踪原亚基、核苷酸与释放产物，回到第二轮而非删除后新建 |

WRC 保留五亚基整体；A 位点 Rac 结合与 WCA 解除抑制相关，D 位点结合不能单独当成已激活。不要把 WAVE/Rac 路线替换成 N-WASP/Cdc42 标签后继续播放。[WRC 原研究](https://www.nature.com/articles/s41467-022-33174-3)

7AQK 是已形成分支的模型；7TPT 是另一成熟接头资产。用它们做末态，不把它们冒充 NPF 招募的完整轨迹。研究中的工程化 CapZ–CA 捕获构建体也不能被解释为天然 CapZ 是 NPF。[成熟分支研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC9295785/)、[活化中间态研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC10434305/)

Arp2 与 Arp3 自己的核苷酸状态另记，不能跟着 actin 一起统一变色。GMF 的作用对象是 Arp2/3；cofilin 去分支与切断分别列事件。分支脱离后 Arp2/3 是否留在母丝、是否重新成核取决于条件，不强制所有分子按同一顺序完全解离。[Arp2/3 老化](https://www.nature.com/articles/s41467-024-46179-x)、[去分支实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC9304009/)、[再次成核实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC10816697/)

## 路线 B：RhoA–mDia1/formin

| 镜头组 | 因果内容 | 实施要点 |
|---|---|---|
| 上游与抑制 | RhoA 的定位和 GDP/GTP 开关；mDia1 DID–DAD 自抑制 | 明确选 mDia1；RhoA 结合不自动证明全长蛋白在所有条件下完全活化 |
| 活化与成核 | 膜/辅因子情境，FH2 二聚体参与成核 | 全长调控与结构中 FH2 片段分开；缺失 FH1 是柔性示意，不是已测坐标 |
| 递送与延伸 | FH1 捕获 profilin–actin → 递送至端部 → 亚基加入 → profilin 释放与 formin 重新定位 | 检查结合伙伴是否相容，至少两轮；端部步态选定一个来源模型并说明推断 |
| 功能输出 | 直丝进入交联网络；若讲收缩，另加入 RhoA–ROCK 支路与肌球蛋白/黏附条件 | ROCK 与 mDia1 是 Rho 下游的并行支路；平行成束本身不是收缩性应力纤维 |
| 终止与回收 | 端部竞争/封端、formin 解离、信号失活、老化丝拆解和单体回收 | 端部占有状态、分子身份、核苷酸交换各自连续 |

RhoA 对纯化全长 mDia1 的活化依赖条件，不能把片段实验外推为无条件开关。[全长 mDia1 研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC3416746/)。mDia1 和 ROCK 协同调控细丝组织的证据来自特定实验体系，不能直接当作每种细胞的统一下游。[mDia1/ROCK 原研究](https://www.nature.com/articles/ncb0799_136)

**避免跨亚型拼接成“观测到的全过程”：** 9B27/8RU2 可作为 mDia1 端部状态；9AZP 的 incoming profilin–actin 是 INF2，8RTY 的 ring complex 是裂殖酵母 Cdc12。2024 年两组工作给出相关但不同的动作解释；可并列讲解，也可选一种推断模型，不把混合插值称为已测 mDia1 轨迹。1Z2C 是 RhoC–mDia1/GppNHp，不能标成 RhoA–GTP。[Nature 研究](https://www.nature.com/articles/s41586-024-07637-0)、[Science 研究与条目](https://www.rcsb.org/structure/8RU2)、[Rho–mDia1 调控结构](https://pmc.ncbi.nlm.nih.gov/articles/PMC1356318/)

封端蛋白与 formin 可形成竞争相关的端部复合物，不能规定“任何共占有都不可能”，也不能让封端后不受限制地继续增长。[单分子研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC4660045/)

## 上游接口怎样接到受体

计划给出两个有依据的接口供项目选择，不把它们拼成一个细胞中的实证链：

- PDGF 刺激与 PI3K/PIP3 依赖的 Rac 活化有实验支持；该证据不能直接指定所有细胞都使用 Tiam1。受体、接头与具体 GEF 的连接需按所选细胞补充来源。[PDGF–Rac 原研究](https://www.sciencedirect.com/science/article/pii/S0960982295000807)
- 可选择 thromboxane A2 受体相关的 p115RhoGEF 膜募集，再结合 Gα13 激活 p115RhoGEF 的机制；不同实验层级要分开。不能将 Gα12 与 Gα13 写成所有条件下等价。[膜募集研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC1223324/)、[Gα13–p115RhoGEF 研究](https://pubmed.ncbi.nlm.nih.gov/9641916/)

未选细胞/受体时，可以完成核心机制教学计划，但交付写明“受体与具体 GEF 接口待项目指定”，不要声称已重建某细胞从配体到骨架的全部细节。

## 回收镜头的检查

两套计划还分别展开条件性的 Rac–PAK1–LIMK 或 Rho–ROCK–LIMK 支路，以及 SSH 相关的 cofilin 去磷酸化。Ser3 磷酸化/去磷酸化调节的是 cofilin；不能与 actin 结合核苷酸的水解混为一件事。它们描述选定调控池，不要求每个本来有活性的 cofilin 都先被磷酸化一次。[PAK1–LIMK](https://www.nature.com/articles/ncb0999_253)、[ROCK–LIMK](https://pubmed.ncbi.nlm.nih.gov/10436159/)、[LIMK–SSH 调控](https://pmc.ncbi.nlm.nih.gov/articles/PMC548651/)

切断产生两个完整片段，不直接得到一堆游离球；片段继续解聚才产生单体。cofilin 对丝的装饰、接口切断、与单体结合、释放需区分。CAP 的片段结构可支持局部接口，不代表已有完整 CAP 的所有动作。回收的 monomer ID 进入交换步骤，再由可用单体池进入新的丝端。[CAP 结构与交换](https://www.nature.com/articles/s41467-018-04231-7)、[网络回收实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC10152149/)

## 可运行的配体保留与状态对照

通用 `prepare_molecular_asset.py` 默认只保留聚合物，**会丢掉游离配体、核苷酸、离子**。不能用其输出证明化学过程已建立。下列专用入口读取 8A2R/S/T 的完整模型链 C，并保留同一作者链下的 ADP、Mg 以及 BEF 或 PO4；核验来源哈希、配体身份和重原子数。三个结构均是 F-actin，不能拿一个独立显示的 F 亚基充当 G-actin。

```bash
SKILL='/actual/path/to/blender-cell-biology'
PYTHON='/actual/python/with/gemmi/numpy/scipy/skimage'
BLENDER='/actual/path/to/blender'
python3 "$SKILL/scripts/fetch_structure.py" 8A2R --output /absolute/project/sources
python3 "$SKILL/scripts/fetch_structure.py" 8A2S --output /absolute/project/sources
python3 "$SKILL/scripts/fetch_structure.py" 8A2T --output /absolute/project/sources
"$PYTHON" "$SKILL/scripts/prepare_actin_states.py" \
  --source-dir /absolute/project/sources --output /absolute/project/state-assets
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_actin_state_comparison.py" -- \
  --asset-dir /absolute/project/state-assets \
  --profile "$SKILL/assets/style-studio.json" --engine cycles --compositor clean \
  --font /absolute/path/to/chinese-font.ttf --output /absolute/project/state-comparison --render
python3 "$SKILL/scripts/audit_mechanism_plan.py" \
  "$SKILL/assets/mechanism-templates/actin-arp23-cycle.json"
python3 "$SKILL/scripts/audit_mechanism_plan.py" \
  "$SKILL/assets/mechanism-templates/actin-formin-cycle.json"
```

8A2R 的 ADP–BeF₃ 是 ATP 状态类似物；8A2S 为 ADP·Pi；8A2T 为 ADP。主视图保持完整亚基，下方另列明确放大倍率的配体视图。状态对照没有反应轨迹、催化水或质子转移；Pi 的真正出口与速率还需相应证据。选择风格和合成器不会改变坐标身份。

## 与现有模板衔接

12 秒 `actin-polymerization` 仅验证接近/转向/落位；`arp23-branch-growth` 仅验证成熟子丝延伸。它们可供镜头节奏参考，旧 F-actin 代理、缺配体和旧避碰路径不能直接升级成全过程。换入 G-actin、profilin、NPF 或 formin 后重新验证相邻结构、路径和配体，不沿用旧“无穿插”结论。

项目验收分别写明：两条路线哪些阶段只有计划、哪些已用真实资产构建、哪些有保存的动画、哪些已导出完整视频并检查。完整视频最后应能回答“刺激怎样启动它、亚基怎样加入、能量和产物怎样变化、它做了什么、怎样停止并再次利用”。

## 已制作的本地完整视频案例

2026-10-24 的教学案例保存在原制作环境的项目目录。持有该项目时，先读其中 `观看与工程入口.md`，再按章节打开 `shots/<id>/scene.blend`；大型工程和视频不包含在本仓库或 skill 压缩包内。新环境重建使用 [portable-actin-films.md](portable-actin-films.md)。

- `movies/Arp23_完整循环.mp4`：5 分 36 秒、14 章，从 PDGF/PI3K/Rac 接口，经 WRC、分支成核、连续延伸、状态变化、网络输出、关闭与去分支，到 cofilin/CAP 回收和重新加入。
- `movies/Formin_完整循环.mp4`：5 分 27 秒、13 章，从 GPCR/Gα13/p115RhoGEF/RhoA 接口，经 mDia1 自抑制解除、成核、两轮供料延伸、功能输出与终止，到拆解回收和重新加入。

共有 21 个可编辑章节，1600×900、24 fps，中文合成讲解；中文宋体常规字重，英文和数字 Times New Roman 常规字重。`validation/final-movies.json` 记录实际导出文件的哈希、解码帧数、连续时间戳与音频检查；`validation/视觉检查.md` 记录抽帧检查范围。制作脚本保存在案例的 `scripts/`；其他机器重建前需解决该案例对父目录缓存和已安装 skill 的依赖，不能将本地案例冒称随包附带的通用生成器。

这些视频是跨物种结构资产支持的机制教学综合。FH2 中间步态、G/F 构象过渡和缺失柔性段为标注的教学插值，局部原子中心重叠诊断仍有剩余；它们不构成全原子无碰撞轨迹或分子动力学验证。条件支路部分用文字/关系图表达。保留原计划的 `design_only` 状态；已导出视频的证据以该项目实际文件和检查记录为准。
