# V 型 ATPase：可重建的双视图连续循环

本入口复现膜上整体复合体与 c 环俯视图，包含 ATP→ADP＋Pi 图形、同步转动、两幅图的 H⁺ 运动及左侧向囊泡腔的释放。384 帧、24 fps、16 秒为一个连续周期。支持 Eevee/Cycles 与 clean/reference/monochrome 合成分支；其他场景风格仍按参考图分析后调整材质、灯光和镜头。

包内仅包含源码与说明。视频、PNG/EXR、`.blend`、结构衍生网格和字体保存在独立项目目录；不从旧视频或旧帧恢复动画。

## 输入与来源

- 主体：[RCSB 6VQ6](https://www.rcsb.org/structure/6VQ6)，大鼠脑 V-ATPase，ADP/SidK 结合的静态结构模型。[Abbas 等，Science 2020](https://doi.org/10.1126/science.aaz2924)。
- 省略 SidK 的 Q/R/S 三条完整链，保留来源中已建模的 31 条泵相关聚合物链；不自动补齐未解析的 H 亚基或柔性区。源模型去除 SidK 不等于活性构象实验。
- 九条常规 c 链 Glu139 与 c″ 链 b 的 Glu98 定义载体角位置；[完整牛脑 V-ATPase 研究](https://www.nature.com/articles/s41467-020-17762-9) 亦提供相关位点依据。十个载体与三个 ATP 图形事件表示理想耦联计量。
- 背景复用公开 [1UBQ](https://www.rcsb.org/structure/1UBQ) 的几何，只用于环境层次，不表示测得的胞质组成。
- 中文宋体与英文 Times New Roman 从用户本机合法可用字体提供。使用常规字重 TTF；如原字体为 TTC，应先选取相应 face 并检查字形。字体不随包分发。

## 从空项目开始

需要 Blender 5.x、FFmpeg；Blender 外的 Python 环境需要 `gemmi numpy scipy scikit-image`。以下变量填写实际位置，不能指向 skill 内部输出目录。

```bash
SKILL=/absolute/path/to/blender-cell-biology
PROJECT=/absolute/path/to/new-vatpase-project
PYTHON=/absolute/path/to/venv/bin/python
BLENDER=/absolute/path/to/blender
FFMPEG=/absolute/path/to/ffmpeg
FONT_ZH=/absolute/path/to/Songti-SC-Regular.ttf
FONT_LATIN=/absolute/path/to/Times-New-Roman.ttf

"$PYTHON" "$SKILL/scripts/fetch_structure.py" 6VQ6 --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 1UBQ --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_vatpase.py" \
  --structure "$PROJECT/sources/6VQ6.cif" --output "$PROJECT/assets/pump"
"$PYTHON" "$SKILL/scripts/prepare_molecular_asset.py" "$PROJECT/sources/1UBQ.cif" \
  --preset generic --output "$PROJECT/assets/ubiquitin"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_vatpase_scene.py" -- \
  --asset-dir "$PROJECT/assets/pump" --crowder-dir "$PROJECT/assets/ubiquitin" \
  --font-zh "$FONT_ZH" --font-latin "$FONT_LATIN" \
  --output "$PROJECT/scene" --engine eevee --compositor reference --width 1600 --samples 64 --preview
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_vatpase.py" -- \
  --asset-dir "$PROJECT/assets/pump" --output "$PROJECT/scene/validation"
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- --directory "$PROJECT/scene/frames" --frames all
"$PYTHON" "$SKILL/scripts/encode_video.py" "$PROJECT/scene" --ffmpeg "$FFMPEG"
```

得到 `scene/scene.blend`、`manifest.json`、`frames/` 和 16 秒 `preview.mp4`。浏览器的 video `loop` 属性可持续循环。先检查关键帧再完整渲染；`--width 640 --samples 8` 可做完整周期的快速重建验证，不作为成片质量标准。宽度须至少 480 且为 32 的倍数，以保证 16:9 输出及 H.264 偶数尺寸。

已有公共结构缓存可代替下载，但要保留来源与哈希；不能把旧渲染或预制 `.blend` 当作重建输入。链文件以十六进制链名保存，避免 A/a 在大小写不敏感文件系统上互相覆盖。

## 检查与科学边界

`audit_vatpase.py` 重开工程，检查 31 条主体链、10 条放大环链与源网格一致；逐帧检查定子固定、左右载体相位、整体/标签画内边界和循环接缝。视觉检查需另看实际输出，尤其是入口接近、环上遮挡、出口离场与 16 秒接缝。编码后用 FFmpeg 完整解码并核对 384 帧、24 fps、时间戳连续及关键帧；代码检查不能替代成片检查。

表面为原子坐标衍生的高斯包络；膜面由结构推定。完整单体分层表达，不切开蛋白。右侧仅保留完整 c 环链，整体省略 a 亚基等遮挡物；曲线是半通道功能示意，不是实验孔道表面。左侧结合标记向观察者偏移 2 nm，入口、出口球表达开放池中的物质流；两图结合/释放相位一致，但接近/离场路径按景别分别编排。它们不是单个氢原子的实测轨迹。

转动、ATP 转换和时长均为机制编排，未重建全部催化构象、质子接力水链或原子碰撞能量。循环表达稳态局部泵送，不计算累积 pH、膜电位、缓冲、泄漏或离子补偿；不能靠倒放酸化来复位。

## 轻量打包与复验

```bash
"$PYTHON" "$SKILL/scripts/package_skill.py" --output /absolute/path/to/blender-cell-biology.skill.zip
```

打包器按文本、脚本与结构参考文件的允许类型收集，并排除帧、输出、验证缓存等目录。发布前解压到新位置，用解压副本运行以上命令，至少生成一个完整低分辨率周期，检查来源、帧数、左右同步、解码与关键帧。新生成的验证视频和帧仍在项目目录，不加回 skill。

2026-10-26 已在本机 Blender 5.2.3 LTS 对解压副本实跑：重新下载 6VQ6/1UBQ、重新生成网格与场景，输出三张 1600×900/64 samples 预览，以及 384 帧、24 fps 的 640×360/8 samples 完整循环。完整视频解码、41 个来源网格、双视图相位与循环接缝检查通过，并检查了解码关键帧。验证没有把旧 `.blend`、视频或帧作为输入；低分辨率完整周期用于检验可执行性，不能代替新项目的成片质量验收。验证记录保存在 skill 外的项目目录，不是运行依赖。
