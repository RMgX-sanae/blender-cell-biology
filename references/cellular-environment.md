# 有景深的细胞微结构背景

用户在 2026-10-28 指定图片与视频的默认起点：用多种细胞微结构形成有远近层次的虚化背景。它适用于教学和写实表现；风格参考、纯色底与透明背景仍可按项目覆盖。

先选区室和景别，再决定环境内容。胞质可使用核糖体、适当的骨架和可溶蛋白；膜邻域需要正确朝向的膜、糖链和邻近蛋白；细胞核可使用染色质和核内复合物。不能为了丰富而把核糖体放进核内或把细胞器缩成分子大小。无定量依据的密度、位置及背景漂移记录为示意。

将主反应区、环境与文字分层。近、中、远景应有不同虚化程度、遮挡和对比；保留部分可辨认的轮廓，避免整片随机光斑。背景可以缓慢运动，但不得盖住核苷酸/离子进出、让同屏视图错时或造成文字闪烁。双视图中的局部图可能独立放大，背景不能暗示它和整体图共享精确比例尺。

分隔线采用与局部背景反差明显的颜色：暗底优先浅灰白，亮底优先深灰。最终输出约 2 px 的可读线宽为起点；检查动画全程，不能让线消失在背景结构中。此规则来自用户对双视图中间分隔线的直接反馈。

## 可运行的胞质背景例子

`prepare_cellular_environment.py` 从公开原始坐标准备三种粗略表面：人 ubiquitin [1UBQ](https://www.rcsb.org/structure/1UBQ)、人 80S 核糖体 [9SPI](https://www.rcsb.org/structure/9SPI) 和哺乳动物肌动蛋白形状代理 [8F8R](https://www.rcsb.org/structure/8F8R)。后者来源为兔肌动蛋白，不把它宣称为人细胞中测得的丝；使用完整内部单体及 D→E 变换延伸示意背景丝。实际项目需要不同区室时重新选资产，不能盲目复用这一组合。

```bash
for code in 1UBQ 9SPI 8F8R; do
  "$PYTHON" "$SKILL/scripts/fetch_structure.py" "$code" --output "$PROJECT/sources"
done
"$PYTHON" "$SKILL/scripts/prepare_cellular_environment.py" \
  --ubiquitin "$PROJECT/sources/1UBQ.cif" \
  --ribosome "$PROJECT/sources/9SPI.cif" \
  --actin "$PROJECT/sources/8F8R.cif" --output "$PROJECT/environment"
```

构建 Ras/SOS 场景时默认采用上述环境，在原有命令中加入 `--environment-dir "$PROJECT/environment"`。省略此参数保留简洁纯色版，供指定纯色或对比检查使用。其他生成器由 agent 按同一分层方式接入；不能把此入口称为已对所有场景自动生效。

`cellular_environment.py` 使用完整表面实例，保存来源、变换及示意布局。远景核糖体、中景丝状骨架与分散蛋白分别进入三个 View Layer；`Cellular environment ... defocus` 节点可分别编辑虚化宽度。当前是有明确深度层次的艺术化合成景深，不是标定的相机光圈或显微 PSF。主体、化学对象和 Typography 在背景合成后叠加，不参与背景模糊。

验收静帧时比较有/无背景的同一反应时点；动画检查入场、反应、离场和首尾，确保环境层不挡文字、不掩盖分子路径问题。保存节点、资产准备方法和参数；视频、帧、准备网格仍在 skill 包之外。

## 按尺度制作背景对照

用户要求从细胞、细胞器、膜到分子分别试验环境。`apply_scale_environment.py` 为以下已知生成器的保存场景提供浅色/深色两种可编辑起点；它不是任意 `.blend` 的自动分类器。

| `--scale` | 输入生成器 | 环境内容 |
|---|---|---|
| `cell` | `prepare_cell_geometry.py` → `build_cell_scene.py` | Allen Cell 三个不同间期细胞的完整膜、DNA 与内部标记结构，保持各自来源相对位置 |
| `organelle` | `build_organelle_scene.py` | 显微参考图约束的 ER 分支曲线，带示意管径和三维深度 |
| `membrane` | `build_molecular_scene.py --scene membrane --membrane-view side --environment` | 保留完整远处受体与胞质大分子；在外叶添加示意糖链 |
| `molecule` | `build_ras_exchange_scene.py --environment-dir ...` | 已准备的胞质分子环境；保持双视图与标签独立 |

膜场景需要与源 `.blend` 同目录的 `manifest.json`，糖链锚点从其中已生成的外叶头基位置选择；糖链组成与数量均是示意，不能解读为 3SN6 测得的糖基化。细胞背景从 Allen Cell 公开显示网格准备，保留其原成像及显示方法限制。细胞器的细密 ER 背景用下面的参考轮廓建模流程，不能冒充原子或电镜重构。分子形状沿用各输入场景的来源记录。四种景别分别构图，不能把这组对照称为同一细胞的连续定量放大。

```bash
"$PYTHON" "$SKILL/scripts/prepare_allen_environment.py" \
  --output "$PROJECT/allen-context"
"$BLENDER" -b "$PROJECT/cell/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/apply_scale_environment.py" -- \
  --scale cell --style dark --engine cycles \
  --context-assets "$PROJECT/allen-context" \
  --output "$PROJECT/cell-dark" --width 1280 --samples 32 --render
```

浅色版使用 `--style light`；其他景别改变 `--scale` 和对应输入。Ras 对照可指定 `--frame 352`，具体帧须匹配源场景时序。`--engine` 独立选择 `cycles` / `eevee`；`--samples` 控制 Cycles 采样。每次输出新目录，保存场景、来源哈希、参数与主体网格/变换签名；该签名只检查原对象的基础网格坐标和物体变换，不替代生物学或修改器求值检查。

浅底采用深色文字和分隔线，并加深 Ras/GTP 等颜色；暗底保留浅色文字和分隔线。主体与环境进入独立 View Layer；细胞近景只在画面边缘参与遮挡，主体位于远中景与近景之间。细胞与细胞器相机适度拉远，以留出环境边缘；原主体形状和摆放保持不变。背景专用灯仅属于对应环境层。先检查结构关系、主体轮廓、膜内外及字线反差，再决定是否使用这一风格扩展动画。

用户反馈指出，等距放置的规则细胞、平行波浪线、散落小球及 Voronoi 管网会显得虚假，且不能解释环境中的结构关系。相应构造已从此入口移除；不要继续靠随机椭球、正弦曲线或几何镶嵌图案模拟真实微环境。先检查参考图的轮廓、分支、局部疏密和空间关系，优先使用合适的公开三维模型。细胞背景重视来源形状、不同朝向与遮挡；细胞器背景使用有身份的局部结构。

形态参考包括 [NIGMS/Janelia 超分辨图 5855](https://www.nigms.nih.gov/image-gallery/5855) 和 [Allen Cell 的结构示意与分辨率说明](https://www.allencell.org/visual-guide-tutorial.html)。实际背景网格来自 [Allen Institute for Cell Science — Visual Guide to Human Cells](https://www.allencell.org/visual-guide-to-human-cells.html)，不是从 NIGMS 二维图反推出的真实三维结构。

`assets/allen-environment-sources.json` 记录公开 OBJ 地址、SHA-256、资源页面及使用条件；`prepare_allen_environment.py` 校验源文件并保留顶点与面，输出外置 NPZ。可用 `--source-dir` 复用已下载的命名 OBJ，但仍校验哈希。生成器只做刚体摆放和均匀缩放，配对膜、DNA、内部结构采用同一个变换；场景记录各对象变换。背景在画面中的位置不是前景模型真实测得的邻域，选择 hiPSC 形态代理也不能自动将前景泛化到该细胞类型。

Allen 资源按其[使用条款](https://alleninstitute.org/legal/terms-of-use)用于非商业研究/教学，并保留[来源署名](https://alleninstitute.org/legal/citation-policy)。公开展示这些渲染时，在同页说明或项目来源文档注明提供者、模型页面及显示加工方式；不必向画面塞入小字。源 OBJ、准备网格、参考图片和渲染均不加入 skill 或 Git 仓库。商业项目需另行核对该资源授权或选择适用资产。

## 用显微参考轮廓制作三维 ER

用户在看到来源网格和显微图小样后，选择继续细化三维背景。细胞保留不同来源形状；细胞器避免放大粗分辨率 ER 表面，改用参考图约束的曲线建模。`prepare_reference_er_network.py` 对参考图的局部轮廓进行阈值化、细化和短末端去噪，保留非均匀分支密度与环路；`apply_scale_environment.py` 生成可编辑的 Blender 三维曲线及可变管径，设置体积侧光和环境虚化。

```bash
# 先按下一节下载同一张公开参考图，再准备曲线坐标
"$PYTHON" "$SKILL/scripts/prepare_reference_er_network.py" \
  --image "$PROJECT/er-reference.png" --output "$PROJECT/er-network"
"$BLENDER" -b "$PROJECT/mitochondrion/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/apply_scale_environment.py" -- \
  --scale organelle --style dark --context-assets "$PROJECT/er-network" \
  --output "$PROJECT/mitochondrion-3d-environment" --width 1280 --samples 32 --render
```

曲线 XY 形态参考原图，深度使用固定 seed 的平滑场，管径、光照与布局是教学建模选择；不是从图中恢复的实测 z，也不是定量分割。投影中交叉的分支可能被合并，因此不能据此推断真实 ER 连通拓扑或接触位点。需要真实拓扑或定量尺度时改用适用的三维分割数据。当前细密网络面向背景、小幅镜头视差；大幅环绕应先补足完整体积布局并检查。将其用于别的细胞区室前重新核对结构身份。

## 固定镜头的显微图背景选项

若公开显示网格的分辨率不足，放大后会变成粗厚团块；调色和虚化不能恢复缺失细节。固定镜头可直接把合适的显微参考图作为背景。它保留来源图中的分支、疏密和空隙，背景是二维图层；需要绕行镜头、真实视差或实体相互作用时仍要使用适用的三维资产。

以下示例把 NIGMS 5855 的 ER 图像裁切、单色调色并虚化，在其前方渲染清楚的线粒体教学模型。它不表示该线粒体实际处于原图中的位置，原图的深度颜色、强度和比例尺也不再作为数据使用。相应处理记录在 `environment-reference-er.json`；来源比例尺不转移到示意合成画面。下载文件保存在项目中，原图不随 skill 分发。

```bash
curl -fL 'https://www.nigms.nih.gov/sites/nigms/files/image-and-video-gallery/Figure%203A.png' \
  -o "$PROJECT/er-reference.png"
"$BLENDER" -b "$PROJECT/mitochondrion/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/apply_scale_environment.py" -- \
  --scale organelle --style dark --engine cycles \
  --reference-image "$PROJECT/er-reference.png" \
  --reference-profile "$SKILL/assets/environment-reference-er.json" \
  --output "$PROJECT/mitochondrion-reference-dark" --width 1280 --samples 32 --render
```

此路径替代 `--context-assets`，不能同时传入；目前接入细胞/细胞器生成器，使用源场景第一个 View Layer 作为主体。图像哈希、裁切范围与宽高比会检查；不拉伸参考结构。新工程内嵌图像并保存可编辑裁切、调色和虚化节点。其他参考图须另写项目 profile，核对图像身份、区室、授权和裁切，不能只替换哈希就宣称机制或尺度适用。

参考提供者为 Jennifer Lippincott-Schwartz / HHMI Janelia，来源为 [NIGMS 5855](https://www.nigms.nih.gov/image-gallery/5855)；图库标注 [CC BY-NC-SA 3.0](https://www.nigms.nih.gov/image-gallery)。使用该示例公开展示时，在图片同页或配套来源文件署名，注明裁切、调色、虚化与三维主体合成，并将改编图片按同一许可发布。无需在画面上增加小字。

2026-10-28 已在 Blender 5.2.3 LTS / Cycles 渲染四尺度明暗静帧，并根据用户反馈以来源网格和显微图背景替换规则的程序化环境。该轮验证是背景小样，审美反馈仍待用户确认；不等于四种完整动画已制作或任意参考风格都能自动复现。生成图片仍保存在仓库外。
