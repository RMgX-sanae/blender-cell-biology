# 材质、灯光和结构共同验收

美观是建模 skill 的主要要求。先用有依据的形体建立轮廓、结构关系和完整单位，再按参考选择表面、光照、配色、景深与合成。不能用交杂的棒、粗圆柱束或随机球团代替需要辨认的细胞骨架、膜系统或蛋白复合物；噪声、散射和虚化也不能补救这样的模型。细胞骨架应来自核定的结构装配、分割/断层数据或有文献约束的可解释形态，按景别保留其真实层级，而不是装饰性的管线。

## 可执行的共享入口

[style-structure-materials.json](../assets/style-structure-materials.json) 是可选择的分角色外观起点。它沿用已观察的 WEHI 参考特征：暖色主体、冷色环境、完整分子轮廓与深度层次。它不是原作者材质、光学参数或原始节点的复刻。

`detail_scene.material(..., settings=...)` 接受 `roughness`、`specular`、`subsurface_weight`、`subsurface_radius`、`subsurface_scale`、`bump`、`bump_strength`、`texture_scale` 和 `texture_detail`。旧调用不传设置时保持原外观；未知键、非有限值和越界权重会拒绝。实际解析值记录在 Blender 材质的 `appearance_parameters`，节点保持可编辑。

| 角色 | 起点与用途 |
|---|---|
| `protein` | 柔和实体表面、少量散射和细微凹凸；主要细节仍来自来源坐标 |
| `nucleic_acid` | 较清楚的窄高光、更低散射与凹凸，保留链的轮廓与相邻链色差 |
| `membrane_heads` | 温和反光与散射，帮助双叶边界可读 |
| `membrane_tails` | 较粗糙、低散射，避免把核心画成强反光栅栏 |
| `environment` | 更粗糙、较低反光的冷色背景，降低与主体竞争；保留真实来源形体 |
| `cofactor` | 清楚的小型原子标记，不用发光或过大的散射半径遮住位置 |

这些值是艺术表现参数。当前配方按 `1 Blender unit = 1 nm` 调整；换成细胞/细胞器尺度时重新选择纹理频率、凹凸距离和散射长度。凹凸只接法线，不接网格位移；不能改变源坐标、轮廓或结构比例。空间填充、平滑表面和核酸原子显示各有信息损失，按解释任务选择，而不是给所有对象套同一种几何和塑料材质。

分子场景与来源状态对照入口消费这些参数：

来源状态对照默认采用这个配置；分子场景入口保留原默认外观，使用下方 `--style-profile` 显式选择，以免同时改变旧项目的几何表示。可随时指定项目自己的配方。

```bash
# 膜复合物及来源分子环境。资产准备沿用 production-workflow.md。
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_molecular_scene.py" -- --scene membrane \
  --asset-dir "$PROJECT/assets/beta2" --environment \
  --crowder-dir "$PROJECT/assets/ubiquitin" --ribosome-dir "$PROJECT/assets/ribosome_background" \
  --style-profile "$SKILL/assets/style-structure-materials.json" \
  --engine cycles --compositor clean --width 1600 --samples 32 --frames 1 \
  --output "$PROJECT/material-review" --render

# 独立来源状态；先按 source-state-comparison.md 完成来源检查。
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_source_state_comparison.py" -- --asset-dir "$PROJECT/states" \
  --font-zh "$FONT_ZH" --font-en "$FONT_EN" \
  --profile "$SKILL/assets/style-structure-materials.json" \
  --width 1800 --samples 32 --output "$PROJECT/state-review" --render
```

`build_molecular_scene.py` 按已定义的蛋白、核酸、膜部件与背景角色路由；来源对照入口按 mmCIF 聚合物类型路由。背景实例的完整结构保留原尺寸，位置、数量和景深只是展示安排，不能称为实测细胞组成。

通用反应双视图 `build_molecular_sequence_scene.py` 可通过 `--style-profile` 使用显式组件角色、灯光、语义配色和背景三层合成；配置及保存检查见 [分子序列外观](molecular-sequence-display.md#按参考选择外观)。[PTP1B 外观配置](../assets/style-sequence-ptp1b.json) 是一个具体起点；新反应重新对应组件，取景和结构表示继续由其 display 文件控制。

曲面受体入口 `build_receptor_distribution_scene.py` 可通过 `--style-profile assets/style-receptor-materials.json` 使用同一角色着色器、按链配色、灯光与远中近层虚化。深浅两套配色共享完整源几何与分布；镜头继续由准备配置 `display` 控制。字段、保存工程检查及边界见 [受体外观](receptors-and-distributions.md#参考外观与分层背景)。

## 大型复合体的角色与构图

`build_complex_scene.py` 的 9SPI/7PFO 结构视图也可读取 `surface_roles`。材质角色由来源链的蛋白/核酸身份决定，不能由色板键名猜测：例如复制体中使用 `heads` 颜色的 AND-1 仍是蛋白，核糖体的青色 rRNA 仍用核酸材质。未提供角色配置时保留原材质数值；不把这些表面颜色或散射参数称为实验测量。

`assets/style-complex-replisome.json` 给出 7PFO 的金色主结构、冷色辅助亚基和洋红/米色 DNA 配方。它采用完整源包络与核酸原子显示，没有模仿参考中的原子空间填充或添加无来源的拥挤背景。继承的 `generic` 输入只包含已建模聚合物，排除非聚合物配体、离子和水；这个入口是结构观察，不是完整化学装配或催化模型。

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_complex_scene.py" -- \
  --scene replisome --asset-dir "$PROJECT/replisome" --reveal \
  --style-profile "$SKILL/assets/style-complex-replisome.json" \
  --width 1600 --samples 48 --output "$PROJECT/complex" --render
"$BLENDER" --disable-autoexec -b "$PROJECT/complex/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_complex_scene.py" -- \
  --asset-dir "$PROJECT/replisome" --profile "$SKILL/assets/style-complex-replisome.json" \
  --output "$PROJECT/complex/saved-check.json"
```

该入口限正交镜头。显式 `camera.position/target` 保持所选视角与瞄准点；`scale` 太小会拒绝，而不是自动覆盖或裁掉源结构。共享 `studio` 支持有限数值的 `camera.roll_degrees`，用于让长复合体在横向画面中更清楚；它只滚转镜头，不改变分子坐标。取景保留全部已建模亚基的边界，包括 `--reveal` 隐藏但仍留在工程中的完整亚基。未设置镜头控制时保留原自动取景。

保存检查核对完整组件清单、源表面、核酸原子中心/显示半径、实际分配的角色材质与颜色、指定灯光和取景；距离推断的键、光学真值和最终美观不由此证明。`test_complex_scene.py` 在重开的真实工程中破坏材质分配、源几何、镜头或灯光并确认拒绝，随后恢复；实际画面仍须独立查看。

## 相互作用入口的迁移

相互作用入口 `build_interaction_scene.py` 可用 `--surface-profile assets/style-structure-materials.json` 接入同一组 PBR 角色着色器，适用于 `cinematic` / `teaching`。此参数只读取 `surface_roles`；灯光、相机与合成仍由该入口的风格和取景参数控制。来源 `polymer_type` 决定蛋白/核酸角色，配体用 `cofactor`。`--representation surface|atoms` 选择原有表面或原子实例表示；原子球半径是显示参数，不是实验测量的原子半径。

需要同时实施参考灯光、镜头、身份色或有来源背景时，改用 [相互作用外观交接](interaction-appearance.md) 的 `--style-profile`。新入口有明确字段，复用共用角色着色器、合成配方和景深节点；原 `--surface-profile` 的含义不变。3SN6 实验构建体保持纯底，已有膜侧 AP2 用于检验经过区室核对的原尺寸细胞质背景。

多链伙伴可为各成员设置 `color` 与简洁 `label`，主体与图例同步分色；显式 `--partner-colors` 则统一按两个伙伴着色。逐图检查标签长度与数量，避免复杂装配的图例拥挤。孤立结构可用 `--orbit-degrees` 绕伙伴轴转动镜头与灯光，重新取景并保留全部源坐标；膜场景沿用其已核定的 `membrane.view`。观察被遮挡的界面时先调整视角，保持完整分子。保存工程检查同时核对身份、颜色、共同变换和取景；实际材质参数仍用 `audit_surface_materials.py` 检查。

已有 Arp2/3 与 formin 章节的外观接入见 [肌动蛋白外观迁移](actin-appearance.md)：完整来源模型、示意图元和环境分别使用角色材质，并保护原有字幕与章节切换。该入口沿用正面教学布局，镜头控制及审阅范围有明确边界。

已有 V 型 ATPase 整体/局部循环通过 [双视图参考外观](vatpase-appearance.md) 使用同一角色材质与合成流程；保留完整链、同步质子流和催化高亮，背景限制在原尺寸整体图的胞质侧。

已有翻译底物片段通过 [翻译参考外观](translation-appearance.md) 将核糖体蛋白、核酸和短肽符号分别路由；以完整 60S 作柔化背景，保留 tRNA 身份色、入场运动和动态双语标签。它不补齐尚未实现的翻译阶段。

## 细胞与细胞器的迁移

以下入口也读取 `surface_roles`，但使用针对各自景别调节的配置，不能直接把分子场景中的 nm 长度照搬到细胞场景：

| 入口 | 显式选择的配置 | 角色与表示边界 |
|---|---|---|
| `build_organelle_scene.py` | `assets/style-organelle-materials.json` | 连续膜面用 `membrane_surface`，膜边缘/头基与核心用 `membrane_heads/tails`；蛋白示意、核酸和原子标记分别路由。形体仍是已有的线粒体教学模型，不新增原子结构证据 |
| `build_cell_scene.py` | 同上 | `cell_boundary`、`membrane_surface`、`chromatin`、`protein`，以及复用的线粒体角色；核仁采用粗粒度核内材质，不据此推断成分或光学性质 |
| `build_registered_cell_scene.py` | `assets/style-registered-cell-materials.json` | `cell_boundary`、`chromatin`、`organelle_marker`；保留同源网格与共同变换，粗分辨率标记表面默认关闭 bump，不伪造超微结构 |

旧配置和未传角色设置的调用保持旧材质数值。新配置通过 `--style-profile` 显式选择。角色可以局部覆盖，未配置的角色沿用入口默认；材质的 `surface_role` 记录实际路由。凹凸使用 Object 局部坐标，随原对象/集合实例的缩放改变显示细节大小；散射半径、Scale 和 Bump Distance 是艺术参数，需要在实际尺度、实例大小和灯光下检查，不自动换算为实测组织光学。

同源细胞入口保留 `--style dark|light` 的身份配色，配置可用 `registered_palette` 覆盖 `cell/dna/structure`、用 `background` 覆盖底色，并设置 `lighting/world/exposure/camera`。显式角色材质关闭旧入口的固定自发光，避免发光掩盖体积。相机目前限正交；换视角时按实际源顶点重新取景，显式 `camera.scale` 太小会报错，不静默裁掉来源形体。`--display` 的共同源旋转与相机视角是不同控制；接近相机上轴的角度可能产生明显画面旋转，先查看小样。

浅底透明边界可能在背景合成后暴露采样颗粒，即使 Cycles 已启用渲染降噪也要查看最终像素。新同源细胞配置用 `composited_beauty_denoise: true` 在背景合成后、文字叠加前再做一次可编辑降噪；它适用于此类粗分辨率显示表面，不能作为补回显微细节的方法。旧配置默认关闭；精细原子或纹理特写需重新检查细节损失。

`apply_scale_environment.py` 从已保存场景读取同一配置，将 `environment` 角色与 `environment_palette` 接入来源细胞/参考 ER 背景。`environment_palette` 可覆盖 `cell/dna/structure/er`。可选 `environment_compositor.dark/light` 接受 `background`、三个远/中/近层的 `blur_pixels_1280` 和 `opacity`；模糊随输出宽度缩放，主体仍走清晰层。无这些字段时保留旧合成设置。此控制针对三维环境层，不改变独立显微背板的来源配置，也不重写分子入口的背景模型。

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_organelle_scene.py" -- \
  --style-profile "$SKILL/assets/style-organelle-materials.json" \
  --frames 1 --width 1600 --samples 64 --output "$PROJECT/organelle"
# 按 cellular-environment.md 准备参考 ER 曲线，再接入场景中的环境配置。
"$BLENDER" -b "$PROJECT/organelle/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/apply_scale_environment.py" -- \
  --scale organelle --style dark --context-assets "$PROJECT/er-network" \
  --width 1600 --samples 64 --output "$PROJECT/organelle-context" --render

"$BLENDER" -b "$PROJECT/organelle-context/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_surface_materials.py" -- \
  --profile "$SKILL/assets/style-organelle-materials.json" \
  --require-roles membrane_surface protein nucleic_acid environment \
  --output "$PROJECT/organelle-context/material-check.json"
```

该检查器核对配置与保存节点的权重、纹理、局部坐标连接、无几何位移，以及主着色器是否连接到活动输出；`--require-roles` 用于防止预期角色完全没有接入。它不能证明源结构保真、颜色/灯光正确或画面美观，仍须使用相应几何检查器并实际看图。`test_surface_materials.py` 可在 Blender 内运行，覆盖透明混合的有效路径以及断连、未使用着色器、参数漂移、位移和非法背景合成参数。

不是所有生成器都已接入角色配置；未列出的入口需单独核对。此次扩展不升级原有形体的证据等级：细胞总览和线粒体仍含夸大的教学结构；同源标记表面仍受来源成像/显示分辨率限制。

## 逐项目调节并检查

1. 复制外观配置到项目，观察参考的体积、表面细节尺度、亮暗关系和环境层次。主角/相邻组分/背景分别选色，颜色应跟随实体身份。深底配较亮主体，低亮度冷色用于环境；避免所有对象同样鲜亮。
2. `lighting` 中明确主光、补光、轮廓光的位置、尺寸、颜色和功率。主光塑造体积，补光保留暗处结构，轮廓光分离主体与环境。更改景别时同步检查光源尺寸/距离；不要靠过曝代替对比。世界填充与合成底色是不同设置。
3. 先用中等分辨率查看真实像素。背景不应整齐排列、左右镜像或抢占主体；选与区室和景别相容的完整来源结构。景深不能模糊关键反应位点。静态来源对照可明确采用纯色背景，以免引入无依据的细胞情境。
4. 保存源几何、材质参数、灯光、相机与合成设置后，再提高分辨率和采样。推荐只渲染足以判断当前问题的静帧；只有交付需要时输出高清完整图或影片。高清不会自动修复构图、低质量形体或错误机制。
5. 分别记录：来源与几何是否保持；节点是否使用所选参数；最终像素中哪些轮廓、暗部、相邻实体与背景关系清楚。没有查看的风格、相机和动画时段保持未验证。

中文宋体、英文/数字 Times New Roman；成图不堆来源小字。工程、预览与高清图存项目目录，skill 仅分发生成代码、参数与说明。外部参考只作参考，不把 AI 增强图称为 Blender 渲染或结构验证。
