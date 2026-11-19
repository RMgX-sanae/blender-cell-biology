---
name: blender-cell-biology
description: Build editable Blender models and animations of cells, organelles, membrane receptors and cellular processes from microscopy, molecular structures and literature. Analyze reference images and reproduce their visual style through materials, lighting, cameras and selectable compositors. Use for cell biology 3D scenes and mechanism animations.
---

# Blender 细胞生物学建模与动画

制作有来源、可编辑、可重建的 Blender 模型、静帧或动画。先写清观众在这个镜头里应能辨认的结构或事件，再选建模入口。形体是否可信、画面是否好看和清楚，都是验收条件。

## 按镜头选择建模方法

| 镜头与问题 | 从哪里起步 | 在画面里核对什么 |
|---|---|---|
| 全细胞，细胞器在哪里、细胞有何极性 | 有同一标本的膜、核和标记时，共用原始坐标与变换。做教学形态时，按参考确定长轴走向、两端厚薄和接触面。较长的细胞或细胞器可轻微弯曲。需要把完整对象放入细胞时使用体素掩膜。[同源细胞](references/registered-cell-display.md)　[细胞形态](references/cells-and-organelles.md)　[体积分布](references/volume-distributions.md) | 先看剪影能否表现细胞形态和极性，再看内部位置及相邻细胞的遮挡与疏密。 |
| 细胞器，内部如何连接 | 从相应区室的显微图、分割或有依据的参数模型重建连续膜、腔与连接处。已有线粒体入口和可调整材质。[形态与拓扑](references/cells-and-organelles.md)　[细胞器材质与灯光](references/material-and-lighting.md#细胞与细胞器的迁移) | 选能看清内外膜和通道的视角，让周围胞质微结构随景深退后。 |
| 膜，受体朝哪边、如何分布 | 核定膜蛋白的生物装配与膜侧，按网格法线放置完整拷贝。记录密度、富集区和聚集尺度的测量依据或示意假设。[受体分布与建场](references/receptors-and-distributions.md)　[结构资产选择](references/molecular-asset-config.md) | 从正面看分布，从侧面看朝向。观察内部时优先布置前后不同的完整单体，保留膜两侧的区室关系。 |
| 分子，哪些伙伴、位点或连接发生变化 | 默认联网核对实验复合物、必要配体及状态。只有单体时查同源模板或可靠界面约束。依任务选择刚体结合、局部拟合或已核定的坐标序列。[结合与组装](references/online-binding-research.md)　[单体候选](references/monomer-binding.md)　[局部端态](references/coordinate-edits.md)　[反应双视图](references/molecular-sequence-display.md) | 总览交代伙伴和方向，特写留给观众需要辨认的位点、底物和连接。保留结构比例与配体身份。 |

跨尺度任务可以用多个镜头。保留同一对象的身份、方向和语义色，为每个镜头单独定单位、取景和背景。放大示意的幅度和省略的结构写入项目说明。[四种尺度的环境](references/cellular-environment.md)

## 从任务到 Blender 工程

1. **确定来源与交付。** 明确对象、物种、区室、状态、观众和交付形式。默认联网查公共资产、显微数据、结构库和原始研究，有合格缓存则复用。记录来源、许可、装配与尺度，将关键几何或动作标为 `measured`、`inferred` 或 `schematic`。若缺少的信息会改变机制，先澄清。正式科学数据图按实际数据绘制。[公开资产与缓存](references/public-assets-first.md)　[证据与简化](references/sources-and-evidence.md)
2. **落实美术方向。** 教学示意和写实表现均可选。有参考图时实际看图，将轮廓或表面表示、主次配色、材质、光源、镜头与合成分别落实到所选生成器。渲染引擎与合成器按项目选择。默认用符合尺度与区室的微结构作景深背景，主体、相邻关键结构和分隔线具有清楚的明度差。用同一构图比较正式尺寸像素。[参考图到场景](references/style-from-reference.md)　[材质灯光](references/material-and-lighting.md)　[渲染与合成](references/rendering-and-compositing.md)
3. **建立可编辑几何。** 按本次镜头选择来源网格、装配实例、参数几何或受约束的候选。保持完整单体、各向同性比例、膜的拓扑和实体身份。固定朝向、单位和坐标后再布光。个体轮廓、微环境和拥挤程度应可解释；分子内部视图优先移开遮挡单体。[精细制作](references/production-workflow.md)　[内部视图](references/section-views.md)
4. **按需要安排时间。** 静态模型直接选视角。分步教学安排可辨的事件并适当停留。连续片段保持接近、转向、落位与释放的衔接。明确要求循环时，检查首尾位置、运动和物质去向。各视图共享事件时点和方向。声称覆盖全过程时，逐项核对上游、执行、能量、产物与回收。演示时间和推定轨迹的依据写入项目记录。[机制深度](references/mechanistic-depth.md)　[运动交接](references/coordinate-motion.md)　[过程与动画](references/process-atlas.md)
5. **重开工程看画面。** 检查正式尺寸的代表帧，看观众能否直接认出主体、连接或运动方向。中文默认宋体，英文与数字默认 Times New Roman。标签用简短自然的文字，画面只保留必要标签，来源脚注和模型限制留在工程说明。遇到遮挡、表面拉伸、同色粘连或过度虚化，回到对应几何、取景或灯光修正。动作问题只取变化区间的原速短预览，需要成片时再渲染全片。[视觉质量标准](references/quality-standard.md)　[抽帧与产物清理](references/review-and-retention.md)

可在项目中复制[场景说明模板](assets/scene-brief.json)，按需填写。静帧可删去时间和事件。在 `entities` 中记录关键对象的 `id`、`compartment`、`evidence` 与 `source_ids`，跨镜头沿用同一身份。每个参考特征写下对应的场景改动，为代表帧写下帧号和要看的现象。`review_frames` 可直接交给保存工程渲染器，`"saved"` 指工程保存时的当前帧。

```bash
"$BLENDER" -b "$PROJECT/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/review" --review-brief "$PROJECT/scene-brief.json"
```

## 查找专门入口与验证范围

| 需要进一步解决的问题 | 读取 |
|---|---|
| 已有哪个生成器，覆盖静态、局部动作还是完整机制 | [入口与实际验证范围](references/workflow-entrypoints.md)　[把既有流程迁移到新对象](references/workflow-transfer.md) |
| 新蛋白的装配拷贝、作者编号、界面交接或化学变化 | [装配身份](references/source-assembly-instances.md)　[结合运动](references/protein-interactions.md)　[多伙伴交接](references/binding-handoff.md)　[反应键级与衔接](references/workflow-transfer.md#反应前后的键级) |
| 分裂、复制、转录、翻译、囊泡或聚合物装配 | [过程目录](references/process-atlas.md)　[囊泡拓扑](references/vesicle-topology.md)　[聚合物模板](references/polymer-templates.md)　[机制分镜](references/assembly-storyboards.md) |
| 改变表现方法、保护标签或接外部合成器 | [风格库](references/style-library.md)　[合成配方](references/render-recipes.md)　[Natron](references/natron-compositing.md) |
| 核查 Blender 构建、复用机制白模或验证视频 | [Blender 工作流](references/blender-workflow.md)　[验证与交付](references/validation.md) |

按用户要求交付工程、图片或视频，保留重建脚本、输入参数、来源与验收记录。发布 skill 时用 `scripts/package_skill.py` 从源码打包，在项目目录复核解包后的入口。字体、准备网格、渲染帧、视频和 `.blend` 留在项目目录。
