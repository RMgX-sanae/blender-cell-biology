# Ras 核苷酸交换与 GAP 水解的分工

这组入口把核苷酸交换单独实现为一个可重建的机制示例：SOS 接近、Ras 开口、旧 GDP 释放、无核苷酸复合物、新 GTP 装载、SOS 离开。两个核苷酸身份不同，原有共价键不变。GEF 不是把结合的 GDP 磷酸化为 GTP。

当前可执行动画对应 **SOS 促进的交换**。GAP 水解还需独立实现 γ-磷酸键断裂、水中氧进入 Pi、蛋白结合 Pi、Pi 释放与构象复位；不能把交换视频称为完整 Ras 循环。资产准备保留 GAP 状态和失败的刚体路径报告，为后续开发提供可检查的起点。

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

`--mode continuous` 去掉阶段之间的停留，仍完整保留交换事件。它从 Ras–GDP 到 Ras–GTP，**不是化学闭环**；不能在最后一帧把 GTP 重置成 GDP。要做循环，先补上 GAP 催化水解与产物释放，再检查循环接缝和核苷酸身份交接。

## 几何、化学与画面验收

1. 状态间使用共同的 1322 个 Ras 原子身份；缺失侧链由 4Q21 的残基局部框架推断，并记录推断清单。共同构建体不包含 1WQ1 的末端 OXT。来源匹配不等于每个过渡原子都被实验观测。
2. SOS 接近、Ras 开口及核苷酸离开必须联动。不得把旧 GDP 停放在已经排斥其磷酸/Mg 环境的无核苷酸末态中。
3. 搜索候选后检查输出的全部帧；严重重叠阈值为 0.20 nm。源结构例外只允许相同原子对最多额外加深 0.01 nm，不删除原子、不降低阈值。几何约束改变口袋中间构象后，可搜索一个端点位移和速度都为零的局部绕行，保持整个核苷酸的共价几何，并记录偏移。整个时序还检查旧、新核苷酸彼此之间的穿插。
4. Ras 构象插值带键长与相邻角度的几何约束，表面从原子重新生成；SOS 的原始表面只做刚体移动。它们不是分子动力学结果，约束也不能证明所有内部非键相互作用都合理。
5. 两个核苷酸和两个 Mg 共 62 个重原子始终保留身份。核对 CCD 重原子集合与连通性；读取 `C1'` 等名称时只去掉 CIF 外层引号，保留糖环原子名称中的撇号。工程检查再次核对边集合和逐帧键长。Mg 跟随对应核苷酸是便于跟踪的显示约定；不据此推断实际水合、解离顺序或扩散复合物。原子球半径用于显示，不表示范德华半径。
6. 整体图与局部图共享时序和化学坐标，局部图省略整条蛋白，不切割单个分子。γ-磷酸的黄色是已有四个原子的显示叠层，不是新增的游离磷酸。
7. 选择能看清核苷酸进出的镜头；只检查三维距离还不够。无核苷酸阶段的两个视图都应没有残留的核苷酸影像；核对编码后视频的进出、标题和整体/局部同步。

该双视图工程使用独立合成节点。滤镜插在视图合成后、Typography 前；不要用单 Beauty 层配方直接覆盖双视图。视频、帧、字体、工程和准备坐标仅放在外部项目目录，分发包保留生成方法与来源。
