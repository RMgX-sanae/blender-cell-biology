# Blender Cell Biology

Blender Cell Biology 指导 AI 智能体把显微数据、分子结构和论文描述的过程做成可编辑的三维模型与动画。你可以用它展示线粒体在整个细胞中的分布，也可以放大两条蛋白链相遇的地方，看清它们怎样结合。

Blender Cell Biology guides an AI agent through building editable 3D models and animations in Blender from microscopy data, molecular structures, and published research. Show how mitochondria are distributed throughout a cell, or zoom in on the site where two protein chains meet.

[中文](#中文) · [English](#english) · [微形互动展厅 / Interactive gallery](https://rmgx-sanae.github.io/blender-cell-biology-gallery/)

[![互动展厅里的鱿鱼血蓝蛋白十聚体 / Squid hemocyanin in the interactive gallery](https://rmgx-sanae.github.io/blender-cell-biology-gallery/assets/exhibition.jpg)](https://rmgx-sanae.github.io/blender-cell-biology-gallery/)

## 中文

### 适合用在什么地方

- **课堂和课件。** 讲 DNA 复制时，把复制叉处正在合成的两条新链放在同一个镜头里，让学生看清前导链连续延伸、后随链分段合成。讲质子泵时，可以让侧视图和转子俯视图同步循环。
- **组会、论文配图和学术报告。** 先用全貌交代复合物，再用同一套坐标放大配体、金属离子或一根化学键。也可以把一张显微 CT 切片与同一标本重建出的三维外形并排展示。
- **科普视频和展览。** 让观众跟着同一批货物从细胞表面进入早期内体。做互动展品时，按亚基分组的模型也方便拆开展示，像展厅里的鱿鱼血蓝蛋白，观众可以亲手把十个亚基分开再合拢。
- **实验室培训。** 根据仪器手册制作装配图和拆解视图，让新人看到移液器按压机构、活塞、弹簧和密封件之间的关系。

### 在 Blender 里能做什么

智能体通过 Blender 的 Python 接口（`bpy`）搭建场景。打开交付的 `.blend` 工程，你可以自己移动相机、隐藏一个亚基、换一套材质，或在时间轴上调整某个动作的快慢。

| 要展示的内容 | 智能体在 Blender 里的做法 |
| --- | --- |
| 细胞、组织和细胞器 | 从同一标本的三维分割数据生成细胞轮廓、核区和细胞器网格，保留它们在标本中的位置和比例。也能参照显微图塑造略有弯曲的细胞，让线粒体内膜折叠成嵴，或按组织结构排列耳蜗毛细胞和支撑细胞。 |
| 膜、受体和空间分布 | 让完整受体沿弯曲的膜排列，胞外端朝外、胞质端朝内。调节密度、每簇大小和聚集位置，就能比较均匀、成簇和单侧富集三种分布。在细胞内部，还可以比较离细胞核远近不同的区域，按实测比例放置细胞器片段。 |
| 蛋白、核酸和结合位点 | 按 PDB 条目给出的装配方式搭建复合物，把不同蛋白链和功能单元分组。总览用分子表面，特写可切换为蛋白带状图、空间填充原子或配体的球棍模型。移开前方遮挡视线的完整亚基，就能露出核心，并留下邻近结构作参照。 |
| 几何与动画 | 用网格和曲线建膜、骨架和器械零件，用几何节点（Geometry Nodes）重复放置脂质等组件。用关键帧安排接近、转向、结合、释放和回位。对有对应原子坐标的多个构象，可用形态键（Shape Keys）做连续过渡。教学片可分步停留，也可做成流畅循环。 |
| 风格、材质和镜头 | 参照图片的配色、表面质感、照明方向和景深，调整材质节点、灯光和相机。教学图可以用正交镜头比较位置，用与背景对比鲜明的线条指出位点。表现细胞环境时，可以让近处的反应位点清晰、远处的微结构自然虚化。 |
| 渲染、合成和交付 | 选择 Eevee 或 Cycles，用 Blender 合成节点调整对比度、灰度、描边或荧光外观。先看关键帧或按原速播放的短预览，再输出高清静帧、图像序列或供后期使用的多层 OpenEXR。FFmpeg 可将图像序列编码为 MP4。 |

遇到新蛋白时，智能体会查阅 RCSB PDB、UniProt 和原始研究，核对物种、复合物、配体和结合位点。手上只有单体结构时，它会寻找同源复合物或实验测得的接触位点，用这些信息建立候选结合构型。项目说明会写明构型的依据、哪些位置来自测量，以及哪些颜色和运动路径是为了讲解而设计的。[蛋白结合与结构来源](references/online-binding-research.md) · [受体在膜上的分布](references/receptors-and-distributions.md) · [细胞内的实测分布](references/volume-distributions.md) · [根据参考图设置场景](references/style-from-reference.md)

### 从哪些作品开始看

- **整体与局部。** 红藻捕光天线可以从数百个亚基放大到两个亚基上的五个色素。酵母剪接体可以从整个复合物放大到 RNA 分支点的 2′–5′ 键。[捕光天线与移液器](examples/three-scale-studies/README.md) · [剪接体与有孔虫 CT](examples/biological-form-studies/README.md) · [血蓝蛋白与小鼠耳蜗](examples/structure-and-function/README.md)
- **讲清一次反应。** PKA 示例追踪一个磷酸基团从 ATP 转移到底物肽。Ras 示例把两件容易混淆的事拆开演示，先看 SOS 促进核苷酸交换，再看 GAP 促进 GTP 水解。[PKA 磷酸转移](references/phosphotransfer.md) · [Ras 开关](references/ras-switch.md)
- **让过程连续播放。** 可运行的项目包括 RNA 聚合酶 II 转录、哺乳动物翻译从起始扫描到核糖体亚基回收、囊泡货物进入早期内体，以及 V 型 ATPase 连续泵送质子。[转录](references/pol2-transcription-cycle.md) · [翻译](references/translation-appearance.md) · [内吞与融合](references/vesicle-topology.md) · [质子泵](references/vatpase-loop.md)
- **准备一节短课。** 可以从一个现成镜头讲起，例如复制叉向前推进、有丝分裂中期以后的染色体分离与胞质分裂，或肌动蛋白丝端的延伸与解聚。过程目录还列出 60 个题材的分镜和制作提示。[复制与细胞分裂](references/production-workflow.md) · [肌动蛋白装配](references/polymer-templates.md) · [过程目录](references/process-atlas.md)

### 使用

把仓库克隆到智能体的 skills 目录，或在本地打开仓库，让智能体从 SKILL.md 开始。使用 Blender 5.x，已有示例在 5.2.3 LTS 和 5.3.0 中运行。处理结构和显微数据时，脚本按需要使用 NumPy、SciPy、scikit-image 和 Gemmi。[打开 SKILL.md](SKILL.md)

```sh
git clone https://github.com/RMgX-sanae/blender-cell-biology.git
```

给膜运输课备课时，可以直接把镜头和交付要求写给智能体。

> 用 Blender 做一段可循环播放的 V 型质子泵课堂动画。左侧从侧面看质子从胞质进入囊泡腔，右侧从上方看 c 环旋转。让膜两侧的质子都清楚可见，背景微结构略微虚化。给我可编辑的 .blend 工程和一段预览。

仓库保存构建脚本、参数和来源说明。下载的数据、可编辑的 `.blend`、静帧和视频保存在你的项目目录。中文标签默认用宋体，英文和数字用 Times New Roman。公开数据与字体按各自的署名和许可条款使用。[选择建模流程](references/workflow-entrypoints.md) · [构建命令](references/production-workflow.md) · [渲染与合成](references/rendering-and-compositing.md) · [打包 skill 源码](scripts/package_skill.py)

## English

### When to use it

- **Teaching and lecture slides.** Show both nascent strands at a replication fork so students can compare continuous leading-strand synthesis with the fragments made on the lagging strand. For a proton pump, synchronize a side view with a top view of the rotor and let the animation loop.
- **Lab meetings, paper illustrations, and research talks.** Introduce a complex in its entirety, then use the same coordinates to zoom in on a ligand, a metal ion, or a single chemical bond. Pair a micro-CT slice with a 3D reconstruction of the same specimen.
- **Science communication and exhibits.** Follow a group of cargo molecules from the cell surface into an early endosome. Models grouped by subunit also lend themselves to interactive exhibits. In the gallery, visitors can pull apart the ten subunits of squid hemocyanin and put them back together.
- **Laboratory training.** Turn an instrument manual into assembled and exploded views. A pipette model shows how the push button, piston, return spring, and seal fit together.

### What the agent can do in Blender

The agent builds scenes through Blender’s Python API (`bpy`). Open the resulting `.blend` file to move a camera, hide a subunit, change a material, or adjust the timing of an action yourself.

| What you want to show | What the agent builds or controls |
| --- | --- |
| Cells, tissues, and organelles | Convert segmentations of the cell surface, nucleus, and organelles from one specimen into meshes that preserve their positions and proportions. For an anatomical illustration, shape gently curved cells from microscopy references, build folds in the mitochondrial inner membrane, or arrange hair cells and supporting cells in the cochlea. |
| Membranes and receptor distributions | Orient each receptor on a curved membrane so that its extracellular and cytosolic ends face the appropriate compartments. Vary density, cluster size, and location to compare uniform, clustered, and polarized distributions. Inside a cell, place organelle fragments in proportions measured at different distances from its nucleus. |
| Proteins, nucleic acids, and binding sites | Build the biological assembly described in a PDB entry and organize Blender Collections by chain or functional unit. Use molecular surfaces for the overview, then switch to protein ribbons, space-filling atoms, or ball-and-stick ligands for a close-up. Hide whole foreground subunits to expose the core while keeping its neighbors in view. |
| Geometry and motion | Build membranes, filaments, and instrument parts with meshes and curves. Use Geometry Nodes to instance repeating components such as lipids. Keyframe an approach, turn, binding event, release, or return. Shape Keys can animate a series of conformations with matched atom coordinates. Stage the events with pauses for teaching or as a continuous loop. |
| Materials, lighting, and cameras | Match a reference image’s palette, surface finish, lighting direction, and depth of field using shader nodes, lights, and cameras. For a teaching diagram, use an orthographic view and high-contrast lines to point out a site. For a cellular environment, keep the reaction in focus and let the more distant structures soften. |
| Rendering, compositing, and delivery | Choose Eevee or Cycles. Use Blender’s compositor for contrast, grayscale, outlines, or a fluorescence-style appearance. Review selected frames or a short clip at playback speed, then render high-resolution stills, image sequences, or multilayer OpenEXR files for post-production. FFmpeg can encode the sequence as an MP4. |

For a new protein, the agent consults RCSB PDB, UniProt, and the original studies to check the species, complex, ligands, and binding sites. When the starting structures are monomers, it looks for a homologous complex or experimentally identified contacts to build a candidate arrangement. The project notes explain the evidence for the arrangement, which positions come from measurements, and which colors or motion paths were chosen for teaching. [Protein structures and binding sites](references/online-binding-research.md) · [Receptor placement](references/receptors-and-distributions.md) · [Measured distributions within cells](references/volume-distributions.md) · [Working from a visual reference](references/style-from-reference.md)

### Examples to explore

- **From a whole structure to a single site.** Zoom from a red algal light-harvesting antenna containing hundreds of subunits to five pigments bound to a pair of protein chains. In the yeast spliceosome, follow the same assembly down to the 2′–5′ bond at the RNA branch point. [Light-harvesting antenna and pipette](examples/three-scale-studies/README.md) · [Spliceosome and foraminiferal micro-CT](examples/biological-form-studies/README.md) · [Hemocyanin and mouse cochlea](examples/structure-and-function/README.md)
- **Follow one reaction.** The PKA example tracks a phosphate group from ATP to a substrate peptide. The Ras example shows two distinct events in sequence, first SOS-mediated nucleotide exchange and then GAP-stimulated GTP hydrolysis. [PKA phosphate transfer](references/phosphotransfer.md) · [The Ras switch](references/ras-switch.md)
- **Mechanisms in motion.** Runnable projects cover RNA polymerase II transcription, mammalian translation from initiation scanning through ribosomal subunit recycling, cargo delivery to an early endosome, and continuous proton pumping by a V-type ATPase. [Transcription](references/pol2-transcription-cycle.md) · [Translation](references/translation-appearance.md) · [Endocytosis and fusion](references/vesicle-topology.md) · [V-type ATPase](references/vatpase-loop.md)
- **Build a short lesson.** Start with a working scene of a replication fork advancing, chromosomes separating after metaphase, the cell dividing its cytoplasm, or an actin filament growing and shortening at its ends. The process atlas gives storyboards and production notes for 60 topics. [Replication and cell division](references/production-workflow.md) · [Actin assembly](references/polymer-templates.md) · [Process storyboards](references/process-atlas.md)

### Quick start

Clone the repository into your agent’s skills directory, or open a local checkout and ask the agent to start with SKILL.md. Use Blender 5.x. Examples have run in Blender 5.2.3 LTS and 5.3.0. Data preparation scripts use NumPy, SciPy, scikit-image, and Gemmi as needed. [Open SKILL.md](SKILL.md)

```sh
git clone https://github.com/RMgX-sanae/blender-cell-biology.git
```

For a membrane transport class, give the agent the camera angles, the movement to show, and the files you need.

> Make a looping classroom animation of a V-type proton pump in Blender. On the left, show protons moving from the cytosol into the vesicle lumen. On the right, show the rotating c-ring from above. Keep the protons visible on both sides of the membrane and soften the background with depth of field. Deliver an editable .blend file and a preview clip.

The repository contains build scripts, parameters, and source notes. Downloads, editable `.blend` files, stills, and videos go in your project folder. Labels default to Songti for Chinese and Times New Roman for English and numerals. The original providers supply the attribution and license terms for data and fonts. [Choose a workflow](references/workflow-entrypoints.md) · [Build commands](references/production-workflow.md) · [Rendering and compositing](references/rendering-and-compositing.md) · [Package the skill source](scripts/package_skill.py)
