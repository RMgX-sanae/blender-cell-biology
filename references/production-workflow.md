# 精细制作入口与实测边界

先完成参考图拆解和生物结构选择，再选脚本。不要把脚本名当作自动完成整个机制的承诺。下列脚本支持 Blender 5.x，在 5.2.3 LTS / Apple M5 Metal 实际运行过；其他版本需要适配节点和小样验证。

## 准备结构资产

先按 [public-assets-first.md](public-assets-first.md) 检索并核对可用坐标、分割或现成模型；已有适用缓存就复用。肌动蛋白与微管的 11 个可执行模板和命令见 [polymer-templates.md](polymer-templates.md)，其他 8 个装配过程的项目分镜见 [assembly-storyboards.md](assembly-storyboards.md)。

在项目自己的虚拟环境安装 `numpy scipy scikit-image gemmi`，不要修改用户其他项目的依赖。Blender 内的生成器只需要自带的 `bpy`、`numpy` 和 skill 脚本；Gemmi 等用于 Blender 外部的资产准备。

```bash
SKILL=/actual/path/to/blender-cell-biology
PROJECT=/absolute/new_project
python3 -m venv "$PROJECT/.venv"
"$PROJECT/.venv/bin/python" -m pip install numpy scipy scikit-image gemmi
python3 "$SKILL/scripts/fetch_structure.py" 3SN6 --format pdb --output "$PROJECT/structures"
"$PROJECT/.venv/bin/python" "$SKILL/scripts/prepare_molecular_asset.py" \
  "$PROJECT/structures/3SN6.pdb" --preset beta2 --output "$PROJECT/beta2"
```

`beta2` 仅用于 3SN6：保留受体 R 中原始编号 <400 的部分及 Gs A/B/G；去掉工程化 T4 融合体、稳定抗体 N 及非聚合物。受体相对膜方向来自 UniProt P07550 TM1–7 端点，属于推定朝向，不是测量的膜平面。缺失 N/C 端和环区不会自动补全。

`transcription` 仅用于 5FLM；`motor` 仅用于 4HNA A/B/K，排除 DARPin。`generic` 用于已核对的其他结构。每个输出保存原始文件哈希、选取规则、单位变换、各链范围和表面参数。聚合物修饰残基按 Gemmi entity 类型保留；氢、游离配体及非 A alternative location 默认排除，若它们对问题必要，应在项目中明确改选。

高斯表面是原子位置的显示包络，不是实验 cryo-EM 密度，也不是精确溶剂排除表面。`surface` 和 `spacefill` 可由风格档案选择；不能把表面平滑或原子半径说成从参考图片测得。

## 分子镜头

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_molecular_scene.py" -- \
  --scene membrane --asset-dir "$PROJECT/beta2" \
  --style-profile "$SKILL/assets/style-molecular-reference.json" \
  --output "$PROJECT/membrane" --width 1920 --samples 64 --render
```

- `membrane`：实验结构衍生的 β2AR–Gs；双叶头部/尾部、蛋白排斥区与侧视剖面。上下叶分别做有排斥距离的面内采样；密度是示意假设。默认 `--membrane-view side`，可选 `arc` 局部弧形剖口。
- `transcription`：5FLM 的静态转录复合物，蛋白与 DNA/RNA 分开。它不是已经执行了核苷酸加入的转录动画。
- `motor`：4HNA 衍生的 α/β-tubulin 与 kinesin 头部，示意的 13 原丝阵列、双头/颈部/茎和囊泡。两头交替保持接触或完成 16 nm 的摆动；整体时间、连接形变、囊泡比例、晶格拼装均为示意。`--animate` 生成 PNG 帧。

### 根据用户反馈增加环境层次

用户明确希望在保持主结构清楚的同时，增加环境层次。可以用 `--environment` 为膜镜头添加前/中/后景；不能简单加满无身份的小球。先准备 1UBQ (`generic`) 和 9SPI (`generic`) 资产，后者用 `--spacing 2.2 --sigma 1.7 --level .6` 制作远景 LOD，再传入 `--crowder-dir /.../ubiquitin --ribosome-dir /.../ribosome_background`。

这一示例包括 3 个受体复合物、2 个外围核糖体和 24 个 ubiquitin 形状。数量和位置属于美术布局，不是特定细胞中测得的组成或分子密度；大分子保持同一 nm 坐标尺度。明确保留主受体轮廓周围的空间，排除明显相交，外围对象允许有意裁切，主对象必须完整。用透视、距离造成的尺寸差、柔和景深和明暗层级表达空间；不要把所有背景对象做成相同尺寸或同样清晰。

用户随后指出方形剖口看着不自然，在两张新对照中选择了第二张侧视剖面，因此将 `side` 设为当前示例默认。弧形剖口的实心补面试验被弃用，不保留为成品或通用规则。其他项目仍按机制可读性和参考图选择视角。

### 核糖体和复制复合体

以 `generic` 准备 9SPI 或 7PFO，使用 `--spacing 1.1 --sigma 1.25 --level .6`，然后运行：

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_complex_scene.py" -- \
  --scene ribosome --asset-dir "$PROJECT/ribosome" \
  --output "$PROJECT/ribosome_view" --reveal --render
```

`ribosome` 的 A/P/E tRNA 与 mRNA 保持原始坐标；`replisome` 使用 7PFO 的链映射。`--reveal` 以完整亚基为单位隐藏前景，形成观察窗口，不移动或剖开分子。具体隐藏对象和解释见 [section-views.md](section-views.md)。完整装配与内部可见版本均为静态实验结构视图，不是翻译或复制催化动画。分子包络在这些大复合体中固定用 surface LOD；不要将 profile 的 spacefill 字段误称为已经对所有链生效。`audit_intact_molecules.py --asset-dir /... --output /.../audit.json` 可在重开工程后逐一比对表面顶点、拓扑与变换，确认完整性。

## 翻译延伸的结构底物循环

`build_translation_scene.py` 使用 9SPI 的 tRNA 与七个已建模 mRNA 残基，制作 A 位点进入、短肽转移、原 E 位点 tRNA 离开、A→P/P→E 和 mRNA 三核苷酸转位。颜色跟踪同一个分子，不随它所在位点重新着色。短肽由三个残基符号增加到四个，新增 C 端随 A-tRNA；N/C 和 mRNA 极性独立标注。

入场使用从左上前景接近的 Bézier 弧线、小幅转向和五次缓动，保留完整 tRNA 与所携 Ile 的连接。路径绕开当前 mRNA；控制点、初始朝向和编排性质写入 manifest。用户对这一修订的反馈是“现在自然多了”，因此作为本例当前的运动参考。它只是这个镜头的可读性方案，换视角或复合体时应重新核对路径，不能直接移植其方向作为生物学规律。

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_translation_scene.py" -- \
  --asset-dir "$PROJECT/ribosome" --output "$PROJECT/translation" \
  --font /absolute/path/to/a/CJK-font.ttf --width 1920 --samples 64 --render

"$BLENDER" -b "$PROJECT/translation/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_translation.py" -- \
  --asset-dir "$PROJECT/ribosome" --output "$PROJECT/translation-audit.json"

"$BLENDER" -b "$PROJECT/translation/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/translation/frames" --frames all
python3 "$SKILL/scripts/encode_video.py" "$PROJECT/translation" --ffmpeg /actual/path/to/ffmpeg
```

提供 CJK 字体时使用中文标注，并将字体打包进本地 `.blend`；未提供时使用英文。按字体自身许可决定工程是否可公开分发。`--style-profile` 仍独立选择风格；浅底自动使用深色标注，合成分支与其他精细入口一致。

这是明确标注的分层观察：40S 全链隐藏，完整 60S 沿观察方向移至活性分子之后，源对象留在工程中；不是核糖体生理解离。tRNA 刚体端点由共享 C4′ 坐标拟合，拟合误差保存在 manifest；中间运动为教学编排。源 mRNA 片段为 CGUAUCC，A/P 位点对应 AUC/CGU；示意短肽不是 9SPI 中观测的肽。连接线表达连通关系，不声称动画每帧都有真实化学键长。

此示例没有 eEF1A/eEF2、GTP 水解、亚基旋转或混合态，不可称为完整翻译催化循环，也不能由此推断 E 位点释放的唯一时序。源结构缺失的 tRNA 区域和 mRNA 两侧不补成假实验数据。完整天然装配的结构视图继续由 `build_complex_scene.py` 提供。

检查除了数值，还必须包括遮挡：从相机向肽链符号发射射线，确认转移与结尾时可见；逐帧检查主要 tRNA 的实际网格边界、刚体变换、肽链连接和 mRNA 的最终三残基位移。核对 A/P/E 分子间最近距离；进场路径也可能穿插，不能只检查静态端点。`audit_translation.py` 另外检查入场 tRNA 与 mRNA 的原子中心距离，用于发现明显交叉；这一检查不是完整碰撞或能量学认证。

## 转录：一次 C 残基加入

用已准备的 5FLM 资产和公开 [CCD CTP.cif](https://files.rcsb.org/ligands/download/CTP.cif)。把后者下载到项目的结构目录，先检查是否已存在，保留来源链接/哈希；不要把理想化学组分坐标称为 Pol II 中观测到的底物。

```bash
"$PROJECT/.venv/bin/python" "$SKILL/scripts/prepare_transcription.py" \
  --asset-dir "$PROJECT/transcription_asset" --ccd "$PROJECT/structures/CTP.cif" \
  --output "$PROJECT/transcription_substrate"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_transcription_scene.py" -- \
  --asset-dir "$PROJECT/transcription_asset" --substrate-dir "$PROJECT/transcription_substrate" \
  --output "$PROJECT/transcription" --font /absolute/path/to/a/CJK-font.ttf --render
"$BLENDER" -b "$PROJECT/transcription/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_transcription.py" -- \
  --asset-dir "$PROJECT/transcription_asset" --substrate-dir "$PROJECT/transcription_substrate" \
  --output "$PROJECT/transcription-audit.json"
```

RNA 的 P:19 G 残基为初始 3′端，增加的 P:20 C 最终采用 5FLM 原始坐标，与模板 T:17 G 配对。CTP 的 29 个重原子分成新增 C 残基的 20 个和 PPi 的 9 个；α 磷酸保留在 RNA 中。预处理以碱基对齐，再拟合四个可转动键的扭转，保留 CCD 起始键长；中间构象插值与接近/离场路径是教学重构。不能把这一过程当作过渡态、金属配位、反应能量或真实动力学的演示。

完整 Pol II 后移，非模板 N 链整条隐藏；该链的 14–25 残基缺口没有被填补或连接。核酸其他已建模原子保留。初末端点、29→20+9 的计量、12 个完整蛋白亚基、画内边界和新增碱基的可见性由 `audit_transcription.py` 检查。仍须查看关键帧和最终视频。

**构建体检查：**2E2H 可提供底物识别参考，但其 R:10 末端缺 O3′，原研究使用 3′去氧阻断延长。不能拿这个末端直接画正常加聚。选用其他结构时，核查修饰、末端原子、非水解类似物及缺失区后再设计事件。

### 放慢并接下一轮 C→U

用户对单轮视频要求“节奏放慢并增加下一个循环”。应延长可读的事件阶段，在对接与连接变化之间、两轮之间留出停顿，再按模板选择下一个底物；不能只重复同一个 CTP。此 5FLM 支架下一位是 T:16 A，因此下一轮用 [CCD UTP](https://www.rcsb.org/ligand/UTP)，新增 U21。它不在实验条目中。

下载 [UTP.cif](https://files.rcsb.org/ligands/download/UTP.cif) 到结构目录，记录来源和哈希。下面的扩展入口针对上面脚本生成、使用中文字体的单轮工程，继承它的材质、灯光、引擎和全部合成分支。每个输出使用新目录：

```bash
"$PROJECT/.venv/bin/python" "$SKILL/scripts/prepare_transcription_next.py" \
  --asset-dir "$PROJECT/transcription_asset" --substrate-dir "$PROJECT/transcription_substrate" \
  --ccd "$PROJECT/structures/UTP.cif" --output "$PROJECT/transcription_two_cycle_substrate"
"$BLENDER" -b "$PROJECT/transcription/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/extend_transcription_scene.py" -- \
  --asset-dir "$PROJECT/transcription_asset" --substrate-dir "$PROJECT/transcription_two_cycle_substrate" \
  --output "$PROJECT/transcription_two_cycles" --render
"$BLENDER" -b "$PROJECT/transcription_two_cycles/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_transcription_cycles.py" -- \
  --asset-dir "$PROJECT/transcription_asset" --substrate-dir "$PROJECT/transcription_two_cycle_substrate" \
  --output "$PROJECT/transcription-two-cycle-audit.json"
```

当前时序为 552 帧、24 fps、23 秒：0–11 秒第一轮，11–13.5 秒转位，13.5–17.5 秒 UTP 进入和停顿，17.5–19 秒连接变化，19–21.5 秒第二个 PPi 离开，末尾保留产物观察时间。这些时间服务于教学可读性。已加入的 C 在转位时改成 RNA 的橙色，金色继续标记当前加入的 U。第二个底物从画外进入；文字保持面向镜头，随核酸端点更新引线；镜头轻微拉远以容纳转向后的完整核酸。

UTP 以 CCD 化学连接为基础，借用源结构 A–T 配对的空间朝向，作几何约束拟合。完整 A16 核苷酸作小幅刚体调整，其两侧 DNA 连接和新增 RNA 连接长度受到约束；整个支架再作一个位置的编排转位。第一轮 C20 保留实验产物坐标，第二轮 U21、A16 加载构象和转位均为教学重构。约束键长并不验证主链连接角、金属催化、能量稳定性或真实转位路径。

`audit_transcription_cycles.py` 检查每轮 29→20+9 重原子计量、保留的完整核酸和 12 条 Pol II 链、552 帧中的 DNA 连接长度、U 的进入/成键附近距离、主体和标签边界以及 U 碱基遮挡。第一轮 PPi 的离场和第二轮 UTP 的入场允许完整对象跨越画框；反应停顿阶段须完整可见。数值通过之后仍需查看实际关键帧，并在编码后检查最终 MP4。

## 复制：连续与分段延伸

```bash
"$PROJECT/.venv/bin/python" "$SKILL/scripts/prepare_replication.py" \
  --output "$PROJECT/replication_prototypes"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_replication_scene.py" -- \
  --asset-dir "$PROJECT/replication_prototypes" --replisome-dir "$PROJECT/replisome" \
  --output "$PROJECT/replication" --font /absolute/path/to/a/CJK-font.ttf --render
"$BLENDER" -b "$PROJECT/replication/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_replication.py" -- \
  --asset-dir "$PROJECT/replication_prototypes" --replisome-dir "$PROJECT/replisome" \
  --output "$PROJECT/replication-audit.json"
```

1BNA 提供四组完整互补核苷酸的局部几何；72 bp 的 CGAT 重复序列、复制叉曲率及时间均为构造。两条亲代链始终保留，每条子双链与一条新链配对。前导新链的索引递增；后随链在每段引物之后按索引递减，二者化学方向均为 5′→3′。每个 12 单位片段含 3 个 RNA 引物符号和 9 个 DNA 单位，是为画面可读性设定的长度。

RNA 引物仅用符号表示，不能拿它讲解核糖化学。前导链的起始位于画外，未显示聚合酶/引物酶、校对、ATP 化学、引物去除和连接。不能称为完整复制全过程。主链连接管表达路径连通，不是每帧精确的共价键长。

7PFO 的 CMG 11 个完整亚基（2–7、C/D/E/F/G）作为一组后移，不切开，也不暗示当前构造的 DNA 路径就是实验 DNA–CMG 界面。其他复合体成员未选入这一 CMG 参考层；完整 7PFO 仍用 `build_complex_scene.py` 查看。`audit_replication.py` 检查亲代身份、互补链共用轨迹、增长次序、引物先决条件、完整亚基、前后层次和标签边界。

两种新示例都先抽查关键帧，再用 `render_saved.py --frames all` 和 `encode_video.py` 完成视频。小内存机器优先顺序渲染；同时运行多个 Cycles 任务若使单帧耗时明显上升，应保留已完成帧，按已保存的同一工程和设置续渲缺失区间。逐段帧必须尺寸/fps/合成设置一致；检查全部 PNG 可读、帧号连续后编码，再完整解码 MP4 并查看关键时点。不要把部分序列当成完成品。

## 细胞器与细胞总览

`build_organelle_scene.py` 使用连续的线粒体内膜、外膜/膜间隙、嵴脊上朝基质的 ATP synthase 二聚体行、环状 mtDNA 和蛋白示意。按图像参考重做了可选材质、灯光、镜头和合成。原始几何改编自用户的线粒体工程，其原件保留。分子大小/数量有夸大，不能在同一图上加精确分子比例尺。

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_organelle_scene.py" -- \
  --style-profile "$SKILL/assets/style-organelle.json" \
  --output "$PROJECT/mitochondrion" --width 2400 --render

"$PROJECT/.venv/bin/python" "$SKILL/scripts/prepare_cell_geometry.py" --output "$PROJECT/cell_asset"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_cell_scene.py" -- \
  --asset-dir "$PROJECT/cell_asset" --output "$PROJECT/cell" --render
```

细胞总览使用共用隐式腔体生成核膜—ER 连通边界，并检查腔体连通分量；Golgi 膜囊保持分开。它是未指定细胞身份的动物细胞示意，不能直接改名为“真实肝细胞”“肿瘤细胞”或“细菌”。核孔、许多细胞器和大部分骨架不在总览里；特写必须按参考另做，不能把缺失藏起来。

## 连续细胞分裂与视频

```bash
"$PROJECT/.venv/bin/python" "$SKILL/scripts/prepare_division.py" \
  --output "$PROJECT/division_asset" --frames 120
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_division_scene.py" -- \
  --asset-dir "$PROJECT/division_asset" --output "$PROJECT/division" --animate
python3 "$SKILL/scripts/encode_video.py" "$PROJECT/division" --ffmpeg /actual/path/to/ffmpeg
```

该动画从中期开始，展示 4 条已复制染色体的 8 条姐妹染色单体分离、分裂沟收缩和整体分开。每个子细胞得到 4 条子染色体；这不是人类核型。连续的 Cassini 曲面序列保持总体积，控制腰部半径趋近零后再打开间隙，避免可见的瞬移和体积跃变。它是可移植的逐帧网格缓存/可见性关键帧，不依赖重新打开 `.blend` 后自动执行 Python handler。核膜重建、染色质去凝聚、中体与 ESCRT 细节尚未包含，不能称为完整有丝分裂全过程。

优先先保存工程，再用 `render_saved.py --frames 1,60,105,106,120` 输出关键帧检查；只有通过后再 `--frames all`。始终查看视频与临界时点，数值体积一致不等于所有拓扑/运动都生物学正确。编码器核对 manifest 的 fps 和完整帧数，不自动重定时。

保存后重新打开运行 `audit_animation.py --kind motor|division --output /absolute/audit.json`，逐帧检查实际评估的姿态/可见性/主体画内范围。马达还需用源坐标和导出的头部位移检查穿插；站立头固定不等于摆动头一定不碰撞。

## 合成器与外部交接

每个新工程保存 `COMP_clean`、`COMP_reference`、`COMP_monochrome` 独立节点组，并设置 fake user，确保未选分支在重开文件后仍保留。构建器或 `render_saved.py` 的 `--compositor` 可选择分支；镜头、灯光和材质由 JSON 独立控制。参考分支先合成背景再调对比，避免透明 premultiplied RGB 被偏移后形成黑边。

```bash
"$BLENDER" -b "$PROJECT/membrane/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/membrane/exr_handoff" --frames 1 --exr
```

该路径关闭成片合成，输出场景线性 32-bit multilayer EXR、Z、diffuse/glossy/emission 和 object Cryptomatte，供 Fusion、Natron 或 After Effects 等可用程序接手；保存 fps、显示变换和 alpha 说明。实际外部软件执行与仅导出交接材料分开报告。EXR 不能再被当作已烘焙 AgX 的 PNG，避免二次显示变换。

## 验收

同时保留：输入与来源、风格观察/假设、运行参数、可编辑 `.blend`、原始渲染、可播放视频、实际发现的问题和修订记录。脚本的运行、几何/时间检查、人的视觉判断分开；当前脚本没实现的 60 类机制配方，仍需要 agent 按 process-atlas.md 制作具体项目，不能宣传为“一键完成 60 段高质量动画”。

## 中文与英文混排字体

用户默认选择中文宋体，英文与数字 Times New Roman。检查实际字体族与字重：macOS 的 Songti.ttc 第一个 face 可能是黑体字重，不能只把 TTC 路径交给 Blender 就假设用了常规宋体。选择宋体常规 face；核验标点、连接号、希腊字母和箭头的字形覆盖，嵌入工程并检查实际渲染。

Blender 单一 Text Curve 需要混排时，可把宋体常规放在 `font`，将 Times New Roman **常规字体文件**放入 `font_bold` 备用槽，再按字符将拉丁文字、数字及其标点的 `body_format.use_bold` 设为 true。这里的布尔值用于选择备用字体，并非要求视觉粗体；要核对备用槽中的实际文件确为常规字重。也可采用独立文本段排版。换字体后重新检查真实像素，不只检查文件名。

完整视频应按固定帧率图像序列编码。完全相同的停留帧可以复用同一张渲染，但序列仍须有全部时间点；不要让静态图片 concat 的默认时间基造成隐式丢帧。导出后解码最终 MP4，核对帧数、连续时间戳、音轨和逐章节抽帧，并分别记录软件检查与视觉检查。
