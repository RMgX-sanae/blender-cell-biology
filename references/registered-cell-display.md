# 保留同源坐标的细胞与细胞器展示

当目标是展示一个实际来源细胞中某结构的位置，先沿用来源中已经配准的细胞边界、DNA 和标记结构。不要分别居中或归一化这些网格，也不要把不同细胞的细胞器放到同一外壳内后称作实测分布。这将分子工作流中的“资产身份 + 共同变换 + 保存工程核对”迁移到细胞尺度，材料和合成步骤继续复用。

本例使用 [Allen Visual Guide](https://www.allencell.org/visual-guide-to-human-cells.html) 的三个间期单结构模型。提供者[说明](https://www.allencell.org/visual-guide-tutorial.html)指出，这类视图由同一细胞的膜、DNA 和一种 GFP 标记结构组成；显示表面受共聚焦成像分辨率限制。ER、线粒体和高尔基体例子分别来自不同细胞，不是三种结构同时测得的一套完整细胞模型。DNA 区域不等同于核膜。

## 复用入口

```bash
# 新目录；可用 --source-dir 复用已有原始 OBJ，仍核对来源哈希。
"$PYTHON" "$SKILL/scripts/prepare_allen_environment.py" --output "$PROJECT/source-cells"
"$BLENDER" -b --python-exit-code 1 \
  --python "$SKILL/scripts/build_registered_cell_scene.py" -- \
  --asset-dir "$PROJECT/source-cells" --cell mito --style dark \
  --engine cycles --recipe neutral --width 1280 --samples 32 \
  --font-zh "$SONGTI" --font-en "$TIMES" --output "$PROJECT/mito-dark" --render
"$BLENDER" -b "$PROJECT/mito-dark/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_registered_cell_scene.py" -- \
  --asset-dir "$PROJECT/source-cells" --output "$PROJECT/mito-dark/audit.json"
```

`--cell er|mito|golgi` 选择来源三件套；`--style dark|light`、引擎和合成配方分别选择。`--display` 接受项目自己的 `registered-cell-display.json` 副本，可改共同朝向、标题、身份颜色和显示跨度。深底使用亮结构，浅底使用深结构；透明细胞边界提供位置参照，不切掉源网格。显示单位没有另行定标，不加定量比例尺。

准备阶段给每个 NPZ 记录哈希。加载时复核源 URL、原始 OBJ 哈希、准备文件哈希，并将坐标与面和原始 OBJ 解析结果逐项比较；旧准备目录缺少新哈希时重新准备。三件网格保留源顶点和面，使用同一个旋转、平移与均匀缩放。材质、光照、相机、双字体和合成复用现有入口；标题/标签进入独立文字层。

保存检查器检查源身份、完整顶点/面、共同正向均匀变换、画框边界和打包字体。实际图像仍需检查遮挡、透明度、颜色反差与标注；通过数值检查不代表视觉自动合格。

2026-11-03 的独立解包复验保留九个源网格的 18 个数组，实际检查了线粒体明暗、浅色高尔基体及深色 ER 四张 1280×800、Cycles 64 samples 图和保存工程。所用配方为 neutral，以及高尔基体的 crisp-neutral；没有逐一验证全部引擎、配方或动态视角。八个错误输入/保存改动被拒绝，范围包括缓存、原始文件、来源错配、单独移动/拉伸、改顶点和删面。

## 换自己的细胞资料

1. 确认各组件来自同一对象、同一坐标系和相容时点。来源分别对齐过时，需要保存已验证的配准变换及误差，不能用包围盒独立居中代替配准。
2. 重新定义组件身份、单位、分割/显示分辨率和文件哈希。这里的读取器只承认所列 Allen 源文件；新来源需自己的准备适配器，不能只替换哈希或结构名称。
3. 沿用共同变换、原网格核对和相机/材质/合成步骤。保留所有组件共同的形态与位置，再选择完整对象显隐或透明显示。
4. 要生成新的细胞器位置时，另定义胞质可用体积、核与其他区室的排除区、对象尺寸/形状及接触约束，并检查整个对象，不能只检查中心点。膜受体的面积采样和法线朝向不等同于细胞器的体积采样。

本入口只展示已有来源位置，不产生新的空间分布或动力学。不会从粗表面推断膜的闭合/连通性、真实接触、ER 腔、Golgi 膜囊、线粒体嵴、原子结构、浓度或完整细胞器数量。透明度和配色是表现选择，不是测量通道强度。

源文件的边检查已经发现实际缺陷：ER、线粒体和高尔基体例子的细胞边界分别有 14/113/105 条开边和 7/78/67 条非流形边。这里只核对三角面边的邻接数，没有证明其他拓扑性质；显示时忠实保留这些来源面，不擅自补洞。**不能把这些显示表面直接当作封闭胞质体积**，也不能将来源表面的孔洞解释成生理开口。需要体积排除、穿膜或膜融合时换用合格的封闭表面/体积分割，并另行验证。

许可和署名沿用[细胞背景来源说明](cellular-environment.md)：按提供者条件用于非商业研究/教学，公开渲染在同页或配套文档署名 Allen Institute for Cell Science、链接源页面并说明显示加工。原始 OBJ、准备网格、字体、`.blend` 和图片均保留在项目目录，不随 skill 分发。
