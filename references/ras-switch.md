# Ras 核苷酸交换与 GAP 水解的分工

这组入口把核苷酸交换单独实现为一个可重建的机制示例：SOS 接近、Ras 开口、旧 GDP 释放、无核苷酸复合物、新 GTP 装载、SOS 离开。两个核苷酸身份不同，原有共价键不变。GEF 不是把结合的 GDP 磷酸化为 GTP。

另一个入口覆盖 **p120GAP 促进的水解**：GAP 接近、反应构象调整、γ-磷酸键断裂、水中氧进入 Pi、蛋白结合 Pi、Pi 释放与 GDP 构象复位。两段可独立重建，也可用下方的组合入口生成同一镜头中的 **GDP→GTP→GDP 状态序列**。它包含核苷酸交换和水解关闭，不是包含膜募集、效应器和所有上游受体的完整 Ras 信号通路。

## 来源与构建体

- [4Q21](https://www.rcsb.org/structure/4Q21)：GDP 状态。为定义一致的 G-domain，取完整的 1–166 残基，明确排除原结构 167–169 残基；这不是全长、带脂质锚定的 Ras。
- [1BKD](https://www.rcsb.org/structure/1BKD)：人 H-Ras 与 SOS1 催化区的无核苷酸复合物。作者链 R 的部分侧链没有可靠坐标；链 S 的已建模区间在 568–1044 内，未解析环不补造。
- [1QRA](https://www.rcsb.org/structure/1QRA)：真正 GTP 的预水解参考，不能和不可水解类似物混称。未解析的 switch-II 侧链位置由残基局部坐标框架推断，逐原子记录。
- [1WQ1](https://www.rcsb.org/structure/1WQ1)：H-Ras 与 p120GAP 催化区，含 GDP/AlF3 过渡态类似物。AlF3 不作为 GTP 或 Pi；新反应中的磷原子必须来自 GTP。该参考不证明动画中的接近或解离轨迹。

SOS 引起的结合口袋变化依据[原始结构研究](https://doi.org/10.1038/28548)。Ras/GAP 的催化解释可参考[复合物结构](https://doi.org/10.1126/science.277.5324.333)及[时间分辨 FTIR 研究](https://doi.org/10.1073/pnas.0712095105)；后者涉及 neurofibromin 和 260 K 实验条件，不能把其数字速率直接搬到 p120GAP 的细胞动画上。

## 重建交换动画

设置实际 `PYTHON`、`BLENDER`、`FFMPEG`、`SKILL`、`PROJECT`、`FONT_ZH`、`FONT_EN` 路径。结构准备需要 Gemmi、NumPy、SciPy、scikit-image；中文使用宋体常规 face，英文/数字使用 Times New Roman Regular。

```bash
for code in 4Q21 1BKD 1QRA 1WQ1; do
  "$PYTHON" "$SKILL/scripts/fetch_structure.py" "$code" --output "$PROJECT/sources"
done
"$PYTHON" "$SKILL/scripts/prepare_ras_switch.py" \
  --sources "$PROJECT/sources" --output "$PROJECT/assets"
"$PYTHON" "$SKILL/scripts/prepare_ras_exchange_motion.py" \
  --asset-dir "$PROJECT/assets" --output "$PROJECT/motion" \
  --mode teaching --frames 576 --fps 24
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_ras_exchange_scene.py" -- \
  --motion-dir "$PROJECT/motion" --output "$PROJECT/scene" \
  --font-zh "$FONT_ZH" --font-en "$FONT_EN" \
  --width 1280 --engine eevee --samples 32 --render
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_ras_exchange_scene.py" -- \
  --output "$PROJECT/scene/scene-audit.json"
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/scene/frames" --frames all
"$PYTHON" "$SKILL/scripts/encode_video.py" "$PROJECT/scene" --ffmpeg "$FFMPEG"
"$PYTHON" "$SKILL/scripts/audit_video.py" "$PROJECT/scene/preview.mp4" \
  --ffmpeg "$FFMPEG" --fps 24 --frames 576 --output "$PROJECT/scene/video-check"
```

`--mode continuous` 去掉交换阶段之间的停留，仍完整保留事件。单独交换从 Ras–GDP 到 Ras–GTP，**不是化学闭环**；不能在最后一帧把 GTP 重置成 GDP。组合入口让交换进入的同一个 GTP 接着水解。首轮离开的旧 GDP 与水解形成的新 GDP 身份不同，状态循环不等于封闭物质循环。

## 重建 GAP 水解动画

复用上面从公开结构准备的 `assets` 和 `sources`，不依赖交换工程或已有渲染结果。两个模式使用相同反应事件：`continuous` 连续推进，`teaching` 增加反应前后停留。默认 576 帧、24 fps 是教学节奏，不是实验速率。

```bash
"$PYTHON" "$SKILL/scripts/prepare_ras_hydrolysis_motion.py" \
  --asset-dir "$PROJECT/assets" --sources "$PROJECT/sources" \
  --output "$PROJECT/hydrolysis-motion" --mode continuous --frames 576 --fps 24
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_ras_hydrolysis_scene.py" -- \
  --motion-dir "$PROJECT/hydrolysis-motion" --output "$PROJECT/hydrolysis-scene" \
  --font-zh "$FONT_ZH" --font-en "$FONT_EN" --width 1280 --engine eevee --samples 32 --render
"$BLENDER" -b "$PROJECT/hydrolysis-scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_ras_hydrolysis_scene.py" -- \
  --output "$PROJECT/hydrolysis-scene/scene-audit.json"
"$BLENDER" -b "$PROJECT/hydrolysis-scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/hydrolysis-scene/frames" --frames all
"$PYTHON" "$SKILL/scripts/encode_video.py" "$PROJECT/hydrolysis-scene" --ffmpeg "$FFMPEG"
"$PYTHON" "$SKILL/scripts/audit_video.py" "$PROJECT/hydrolysis-scene/preview.mp4" \
  --ffmpeg "$FFMPEG" --fps 24 --frames 576 --output "$PROJECT/hydrolysis-scene/video-check"
```

两个场景生成器均可传入 `--environment-dir "$PROJECT/environment"`，加入默认的细胞微结构景深背景；公开资产准备见 [cellular-environment.md](cellular-environment.md)。背景是来源形状的示意摆放，不参与反应区的化学或碰撞判断。

水解入口的专项约束：

- 共 34 个化学重原子：GTP 32 + 水中氧 1 + Mg 1，生成 GDP 28 + Pi 5 + Mg 1。PG、O1G/O2G/O3G 始终来自同一个 GTP，Pi 的第四个氧来自 1WQ1 作者链 R 的水 230。没有 Al/F 原子变成磷酸。省略氢、质子转移、质子化状态、电子结构和其他溶剂。
- 1QRA 的真正 GTP 拟合到 1WQ1 的 GDP 骨架；不为避开 Mg 而任意平移整个 GTP。Pi 末态是理想化四面体，不是已测得的产物晶体。反应键的粗细插值用于表达成键/断键，不是计算得到的键级。
- 1WQ1 支持 p120GAP Arg789 和 Ras Gln61 的催化构象解释；动画没有计算活化能。蛋白结合 Pi 及其释放/构象复位的定性顺序参考上述 neurofibromin FTIR 工作，不能据此宣称 p120GAP 的相同数值速率或逐原子轨迹。
- `ras-gap-route.json` 保存来源哈希、整体接近路径和 Lys949/Arg903 的侧链转角。两个侧链随路径调整，其余 GAP 原子保持刚体几何；这些转角是避免穿插的示意安排，不是实测 rotamer 动力学。全部原子保留，检查侧链与其余 GAP 的内部接触。
- 全时序检查 Ras、GAP、核苷酸骨架、磷酸、水中氧和 Mg 的相关原子对。Mg–磷酸氧的特定短距离配位按源结构处理；不得把全局 0.20 nm 阈值降下来。需要时记录局部平滑的 GAP 刚体绕行和单个 Mg 位移，端点位移与速度均为零，再检查所有伙伴。
- 从保存的工程重新读取形变坐标、双视图变换和反应键；独立核对 CCD 图、34 个身份、反应前连通分量 32/1/1 与产物 28/5/1、Pi 释放时的刚体几何、GDP 终点和逐帧字幕。几何检查通过后仍需查看编码视频，尤其是水中氧遮挡和 Pi 离开方向。

## 重建连续状态序列

`prepare_ras_cycle_motion.py` 使用前两段的准备坐标与验证记录，合成一份连续时序；不拼接已有视频，不输入旧 `.blend`。两段必须来自同一个资产目录、相同 fps 和相同动画模式。下面使用连续模式；把两处 `continuous` 同时改成 `teaching` 可保留教学停顿。默认合计 1152 帧、48 秒。

```bash
"$PYTHON" "$SKILL/scripts/prepare_ras_exchange_motion.py" \
  --asset-dir "$PROJECT/assets" --output "$PROJECT/cycle-exchange-motion" \
  --mode continuous --frames 576 --fps 24
"$PYTHON" "$SKILL/scripts/prepare_ras_hydrolysis_motion.py" \
  --asset-dir "$PROJECT/assets" --sources "$PROJECT/sources" \
  --output "$PROJECT/cycle-hydrolysis-motion" --mode continuous --frames 576 --fps 24
"$PYTHON" "$SKILL/scripts/prepare_ras_cycle_motion.py" \
  --exchange-dir "$PROJECT/cycle-exchange-motion" \
  --hydrolysis-dir "$PROJECT/cycle-hydrolysis-motion" --output "$PROJECT/state-cycle-motion"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_ras_hydrolysis_scene.py" -- \
  --motion-dir "$PROJECT/state-cycle-motion" --output "$PROJECT/state-cycle-scene" \
  --font-zh "$FONT_ZH" --font-en "$FONT_EN" --width 1280 --engine eevee --samples 32 --render
"$BLENDER" -b "$PROJECT/state-cycle-scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_ras_hydrolysis_scene.py" -- \
  --output "$PROJECT/state-cycle-scene/scene-audit.json"
"$BLENDER" -b "$PROJECT/state-cycle-scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/state-cycle-scene/frames" --frames all
"$PYTHON" "$SKILL/scripts/encode_video.py" "$PROJECT/state-cycle-scene" --ffmpeg "$FFMPEG"
"$PYTHON" "$SKILL/scripts/audit_video.py" "$PROJECT/state-cycle-scene/preview.mp4" \
  --ffmpeg "$FFMPEG" --fps 24 --frames 1152 --output "$PROJECT/state-cycle-scene/video-check"
```

组合时保持以下条件：

- 两段交接处 Ras、GTP 骨架、γ-磷酸与 Mg 的坐标必须一致，不能用淡入淡出遮住坐标跳变。统一镜头，两个视图共享同一批形变坐标；γ-磷酸始终来自新进入的 GTP。
- 共 63 个化学重原子：旧 GDP_A 28、新 GTP_B 32、水中氧 1、Mg_A/Mg_B 各 1。GDP_A 保持完整并离开；GTP_B 生成 GDP_B 和 Pi_B，Pi_B 的第四个氧来自水。最终 Ras 回到同一 GDP 构象，但绑定的是另一个 GDP 分子。没有在循环接缝凭空补磷酸或重置分子身份。
- 催化水中氧在 GTP 装载后沿示意曲线进入；它不能提前停在会排斥该位置的 GDP/GEF 构象中。路径须同时通过伙伴接触检查和相机检查，不能从镜头近裁剪面突然出现。远处 GAP 在其到达前、离开后应完整退出画面，不能残留在分隔线旁。
- 合并后重新检查九类对象的全部 36 组伙伴对，并保留源结构允许的特定配位/反应键例外。单独两段通过不等于合并后的旧 GDP、水、SOS、GAP 与新产物之间也通过。
- 最终工程核对 63 个身份、反应前后连通分量、SOS 原始网格与逐帧刚体变换、Ras/GAP 形变、字幕覆盖、无核苷酸阶段的两个视图，以及水中氧的可见进入。成功解码后，仍抽看真实编码帧与交接前后的连续片段。

该入口只完成一次开放底物池中的状态循环。要重复多轮，下一轮必须接收上一轮的 GDP_B，并给新进入的 GTP_C、水和 Mg 分配新身份；不得直接把 48 秒视频首尾重复并宣称分子身份连续。

## 几何、化学与画面验收

1. 状态间使用共同的 1322 个 Ras 原子身份；缺失侧链由 4Q21 的残基局部框架推断，并记录推断清单。共同构建体不包含 1WQ1 的末端 OXT。来源匹配不等于每个过渡原子都被实验观测。
2. SOS 接近、Ras 开口及核苷酸离开必须联动。不得把旧 GDP 停放在已经排斥其磷酸/Mg 环境的无核苷酸末态中。
3. 搜索候选后检查输出的全部帧；严重重叠阈值为 0.20 nm。源结构例外只允许相同原子对最多额外加深 0.01 nm，不删除原子、不降低阈值。几何约束改变口袋中间构象后，可搜索一个端点位移和速度都为零的局部绕行，保持整个核苷酸的共价几何，并记录偏移。整个时序还检查旧、新核苷酸彼此之间的穿插。
4. Ras 构象插值带键长与相邻角度的几何约束，表面从原子重新生成；SOS 的原始表面只做刚体移动。它们不是分子动力学结果，约束也不能证明所有内部非键相互作用都合理。
5. 两个核苷酸和两个 Mg 共 62 个重原子始终保留身份。核对 CCD 重原子集合与连通性；读取 `C1'` 等名称时只去掉 CIF 外层引号，保留糖环原子名称中的撇号。工程检查再次核对边集合和逐帧键长。Mg 跟随对应核苷酸是便于跟踪的显示约定；不据此推断实际水合、解离顺序或扩散复合物。原子球半径用于显示，不表示范德华半径。
6. 整体图与局部图共享时序和化学坐标，局部图省略整条蛋白，不切割单个分子。γ-磷酸的黄色是已有四个原子的显示叠层，不是新增的游离磷酸。
7. 选择能看清核苷酸进出的镜头；只检查三维距离还不够。无核苷酸阶段的两个视图都应没有残留的核苷酸影像；核对编码后视频的进出、标题和整体/局部同步。

该双视图工程使用独立合成节点。滤镜插在视图合成后、Typography 前；不要用单 Beauty 层配方直接覆盖双视图。视频、帧、字体、工程和准备坐标仅放在外部项目目录，分发包保留生成方法与来源。
