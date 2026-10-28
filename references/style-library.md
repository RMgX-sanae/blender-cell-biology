# 从参考图到可编辑风格库

用于扩展参考图模仿、比较不同表示法，或将选定外观迁移到细胞/细胞器/分子机制工程。先读 `style-from-reference.md` 的看图方法，再按本文件选择实现。风格小样使用真实结构坐标或明确标注的示意细胞；它们不替代过程本身的生物学建模。

入口为 `assets/style-library.json`、`scripts/style_nodes.py`、`scripts/build_style_scene.py`。参考图出处、原图 URL、哈希和再利用边界在 [style-reference-catalog.json](style-reference-catalog.json)。**包内没有参考图、渲染图、动画或 `.blend`**。

## 选择与实现

表中的参数是实现起点，不是从参考图反演得到的原作者设置。参考实际可见的特征与实现假设分别记录；不要因为配色接近就宣布复现完成。

| ID | 参考与观察 | 实施重点 | 重点检查 |
|---|---|---|---|
| `watercolor` | Goodsell：扁平色块、细轮廓、密集完整分子 | 正交、低反射的 emission 色块、物体空间轻微纹理、独立内轮廓与外轮廓合成、完整分子前后叠放 | 不能把表面所有小皱褶画成粗黑线；样例重复 ubiquitin 只作外观测试，不代表胞质组成 |
| `illustrative` | Mol*：原子球面、细边线、较浅的内部 | Geometry Nodes 原子实例，Eevee Shader to RGB 控制明暗，Normal/Alpha 描边 | 原子相接关系与大小不被描边改坏；此材质不能直接用于 Cycles |
| `clay` | ChimeraX 柔光表面：宽软高光与凹处阴影 | 坐标衍生包络、粗糙漫反射、宽 Area 光、低环境填充 | 阴影不能把相邻亚基伪造成同一块或遮住缝隙 |
| `cinematic` | WEHI：暖主结构、冷环境、原子细节和景深 | 主体原子实例、相机空间分布的完整背景分子、主/填/轮廓光、有限 DOF、轻 bloom | 环境要同时覆盖前后层次；关键结合/释放部位保持清楚 |
| `teaching` | Allen：浅背景、单一强调色、柔和表面 | 简洁表面、低饱和辅助颜色、控制高光、Standard 显示变换 | 不复制网站 UI；标签仍用宋体/Times New Roman，不加入细小脚注 |
| `transparent` | ChimeraX calmodulin：淡外壳、清楚的内部丝带 | Transparent/Principled 混合，Facing 控制边缘，内部丝带独立材质 | 保留整个单体；外壳不能盖住内部或出现严重透明排序问题 |
| `ribbon` | ChimeraX 二级结构表示 | 读取原结构 HELIX/SHEET 注释、CA 路径、肽羰基方向；螺旋带/片层箭头/细环 | 注释缺失、链断点、插入码必须核对；平滑带不是原子轨迹 |
| `ball-stick` | 原子与键的局部结构表示 | 小原子实例＋细曲线管 | 样例蛋白键由距离初筛，正式展示化学需 CCD/拓扑核对；不能把示意键当键级证据 |
| `spacefill` | 原子填充表示 | 重原子位置＋常规 VDW 半径、GN 实例、元素或链颜色 | H 默认省略；图示颜色不是分子的真实可见颜色 |
| `density-wire` | ChimeraX fitmap：透明网格＋内部链 | Points to Volume → Volume to Mesh → Decimate → Wireframe，内链独立 | 本示例是原子衍生示意密度，**不是实测 EM 图**；真实图应导入原始体数据 |
| `density-cloud` | 从密度外壳/透明叠加拓展的体积方案 | 原子点生成 live volume、Volume Info 驱动散射与发光 | 这是拓展实现；参考图不提供散射系数。不要用不透明胶状表面假装体积数据 |
| `xray` | ChimeraX 暗场透明外壳 | Facing 薄轮廓、半透明 emission 外壳、亮色内部丝带 | 这里的 xray 是分子可视化名称，不是医学 X 射线成像模拟 |
| `fluorescence` | ChimeraX 展示的 Allen 多通道图 | 独立结构通道、GN 曲线纤维与点实例、少量 blur/glow、固定纹理种子 | 样例通道是示意配色；不从 emission/glow 测表达量或共定位 |
| `sem` | HKU SEM：表面突起、边缘亮度、灰度颗粒 | 真实表面几何、沿表面法线实例化小突起、弯曲细突起、Facing、细 bump、灰度颗粒 | 表面纹理与细胞形态分开核对；示例不是病毒颗粒检测或 SEM 物理模拟 |
| `dic` | MicroscopyU HeLa：亮暗成对，具有方向性 | 高度代理 pass → 两个相反方向的像素位移 → 差分 → 灰背景 | 调方向、位移和增益；代理高度不等于光程，不能据此估计折射率 |
| `phase` | 同一 HeLa 的相差图：边缘亮晕/暗边 | 高度代理＋宽核模糊之差产生 halo，灰度细纹理 | 与 DIC 的方向差分分开；该核未标定光学环板/PSF |
| `histology` | MicroscopyU：粉色胞质、紫色核、细纹理与间隙 | 独立胞质/核材质、多尺度染色纹理、核内示意点、平面投影 | 仅为染色外观起点，不是特定组织或疾病重建，不用于诊断 |
| `tem` | HKU TEM：暗膜、亮腔、密度颗粒 | 完整器官样轮廓与灰度纹理层、膜线独立、固定细颗粒 | 示例是投影式教学外观，不复制真实超薄切片拓扑；不切开单个分子 |

前 12 种使用分子结构；后 6 种使用示意细胞形态，便于比较外观算子。后者仍明显是示意图，**不能宣称达到参考显微照片的形态/光学真实性**。需要更接近具体图片时，优先替换为该项目的分割轮廓、细胞形态与内部结构，再调像差/纹理，不能仅加滤镜。

补充实查参考：[Nikon Small World 的 HaCaT 角蛋白图](https://www.nikonsmallworld.com/galleries/2017-photomicrography-competition/immortalized-human-skin-cells-expressing-fluorescently-tagged-keratin)可见细密、弯曲的黄橙色网络和有纹理的蓝色细胞核。迁移到这种参考时应采用局部共享的纤维方向场、受控的小幅弯曲和更接近参考的景别，避免把随机直线网当作已经匹配。当前示例只检验通道、曲线实例和合成流程。

## 先看图，再用测量辅助判断

用图像查看工具看参考的原图，记录形体、镜头、遮挡和光照。必要时排除 UI/水印区域后测量色板、亮度分位数和边缘比例：

```bash
"$PYTHON" "$SKILL/scripts/analyze_style_image.py" /absolute/project/reference.png \
  --crop 0 0 800 600 --candidate /absolute/project/test.png \
  --output /absolute/project/style-analysis.json
```

不要不经检查照抄示例 crop。统计只描述图像；不能自动得出渲染引擎、roughness、材料类别或三维结构。看图后填写 JSON 内的观察与实现假设，选择最接近的起点，复制风格档案到项目目录修改。`strong_edge_fraction` 会受分辨率、标注和构图影响，不是风格相似度评分。

## 从空项目重建

需要 Blender **5.2**；外部 Python 使用 `gemmi numpy scipy scikit-image pillow`。这些脚本的 API 已在 5.2.3 LTS 实跑，不能据此宣布支持全部 Blender 版本。

```bash
SKILL=/absolute/path/to/blender-cell-biology
PROJECT=/absolute/path/to/new-style-project
PYTHON=/absolute/path/to/python
BLENDER=/absolute/path/to/blender

"$PYTHON" "$SKILL/scripts/fetch_structure.py" 1UBQ --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_style_asset.py" "$PROJECT/sources/1UBQ.cif" \
  --output "$PROJECT/assets/ubiquitin"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_style_scene.py" -- \
  --asset-dir "$PROJECT/assets/ubiquitin" --styles all \
  --output "$PROJECT/comparison" --width 960 --samples 48
```

每种风格得到独立的 `scene.blend`、`frames/frame_0001.png` 和 `build-report.json`。`--styles` 可选一个或逗号分隔的多个 ID；仅生成细胞外观时不需要 `--asset-dir`。`--chain A,B` 或 `--chain all` 选择保持共同坐标系的多条链；是否对应生物装配仍须查来源。准备器默认省略非聚合物，配体/金属/heme 等需要按机制另行导入，不能因此声称模型完整。

`--profiles /absolute/project/style-library.json` 使用项目参数副本。`--engine profile|eevee|cycles` 与 `--compositor profile|clean|ink|filmic|fluorescence|grain|dic|phase|soft` 分开选择。遇到材质的硬性引擎限制应更换材质方案，不静默忽略。原有精细场景的 `clean/reference/monochrome` 与这里的 `COMP_STYLE_*` 是不同节点组。

## 迁移到实际机制工程

1. 先完成该过程的实体、拓扑和状态序列，保留来源坐标。小样中的 ubiquitin 重复阵列和示意细胞不是通用生物场景资产。
2. 从项目脚本导入 `style_nodes`。用 `material()` 创建所需材质，按实体身份赋色；`instance_spheres()` 适用于真实原子坐标或有说明的粒子，`curve_tubes()` 适用于核酸/纤维路径，`density()` 只用于明确标注的示意密度。实际 EM/光学体数据需走 `data-import.md`。
3. GN 实例尽量保留为实例。需要导出网格、逐面编辑或其他网格操作时才 Realize Instances；不要为渲染提前复制数万对象。实例中心、半径和单位属于结构层，材质/轮廓属于外观层。
4. 调光和构图时保留同一实体在各镜头中的语义颜色。背景布局按相机右/上/深度轴组织，不能把所有背景物都放到斜视相机的下方。显示内部采用完整单体分层或淡外壳，不直接切开蛋白。
5. `compositor(scene, profile)` 创建可编辑 `COMP_STYLE_*`。保留 clean 分支，用其定位几何/材质问题。支持原流程的线性多层 EXR、Z、Normal、对象 Cryptomatte；外部 Fusion/Natron/AE 仅在实际可用时执行，不能把导出 EXR 算作运行了外部合成器。
6. 物理模拟、Geometry Nodes Simulation/Repeat Zone、Grease Pencil/Line Art、粒子/毛发、雕刻和外部 Molecular Nodes 等工具应按项目需要选用。当前库没有验证所有这些工具，不要求把无关工具强行加到场景里；新节点/插件先用一个对象试运行并记录版本。

## 动画验证与可复现边界

风格与动画模式独立。机制动画继续使用该过程的完整状态机。这里的 `--animate` 仅是 **4 秒风格稳定性预览**：主体小幅周期转动、荧光示例的点整体作周期偏移；不代表细胞过程、布朗运动或分子动力学。

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_style_scene.py" -- \
  --asset-dir "$PROJECT/assets/ubiquitin" --styles cinematic \
  --output "$PROJECT/motion" --width 720 --samples 32 --frames 96 --animate
"$PYTHON" "$SKILL/scripts/encode_video.py" "$PROJECT/motion/cinematic" \
  --fps 24 --ffmpeg /absolute/path/to/ffmpeg
```

`--render-frames 1,25,49,73,97` 可检查 96 帧周期的首尾等价状态；正式视频只编码 1–96，不能再加入等价的第 97 帧造成重复停顿。核对实际像素、噪点、描边、透明层、DOF 和首尾速度；固定随机种子不能替代观看视频。

生成器用整个周期的结构范围安排分子镜头，避免静帧合适但转动时裁掉末端。重开工程后运行：

```bash
"$BLENDER" -b "$PROJECT/motion/cinematic/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_style_scene.py" -- \
  --asset-dir "$PROJECT/assets/ubiquitin" --output "$PROJECT/scene-audit.json"
```

此检查核对原子坐标、节点连线、颗粒贴图、GN 周期状态与关键姿态画面边界；不替代实际视频解码和视觉对照。

生成器不读取旧帧/旧视频。生成的 grain 图由代码重建并打包在输出 `.blend` 内，既不依赖外部缓存，也不进入 skill。新输出一律使用 skill 外的新目录。发布时运行 `package_skill.py`，从解压副本准备新资产、生成新场景和一个完整短循环，检查解码及节点效果；保留验证报告在项目目录。

## 本轮检查的具体边界

2026-10-28：18 类已在 Blender 5.2.3 实际生成静态图并逐类查看，修正了法线描边过重、暗场环境偏在下半部、透明外壳过厚和颗粒层变黑。分子表面/原子/丝带/网格/体积等表示不同，非仅替换色板。示意细胞外观保留上述局限，不称为实测图重建。保存节点可编辑，软件成功与视觉接近程度分别记录；新项目仍须与其实际参考逐项对照。

Blender 5.2 的具体接口差异：Compositor 数学/渐变/Mix 使用共享 `ShaderNodeMath`、`ShaderNodeValToRGB`、`ShaderNodeMixRGB`；Resample Curve 的模式通过 `inputs['Mode']` 设置。生成图像应先设置 colorspace，再填充 pixels、update、pack，避免颜色空间切换把刚填入的颗粒缓存重置为黑色。不要将这些具体约束推广到未经试跑的旧版本。
