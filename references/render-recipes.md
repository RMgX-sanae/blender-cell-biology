# 渲染滤镜与合成器组合

先完成结构、摄影和动画，再选择合成；滤镜不能修复错误拓扑、缺失机制或严重穿插。这里的配方使用 Blender 原生可编辑节点；另有已实际执行的 [Natron 调色与独立文字入口](natron-compositing.md)。Fusion/After Effects 等仍按项目确认可用并执行后才能声明支持。

`assets/render-recipes.json` 定义有序滤镜链，`compositor_stack.py` 实际构建节点。材质风格、渲染引擎、动画模式和合成配方相互独立。

项目滤镜链还可显式加入 `{"type":"denoise"}`，对已合成背景的 Beauty 做 Blender Denoise，随后才叠加受保护文字。此操作用于实际观察到的透明边界采样颗粒；当前检查对象是同源细胞的 Cycles 静帧。它不需要数据 pass，也不恢复真实微结构；细密表面及动画的细节损失/闪烁需另行检查，不默认应用到全部场景。

| 配方 | 实际操作与用途 |
|---|---|
| neutral | 中性基准，用来比较滤镜影响 |
| teaching-outline | 法线与轮廓边线、轻微对比，浅色背景 |
| depth-layered | 固定深度区间的雾色与轻微光晕，改善环境层次 |
| soft-cinema | lift/gamma/gain、光晕、轻微暗角和固定颗粒 |
| contact-relief | AO 与法线轮廓，突出表面接触；不把 AO 当作相互作用强度 |
| fluorescent-soft | 亮部光晕与曝光；须配合明确的发光材质，不是显微成像模拟 |
| warm-print | 降饱和、暖色增益与纸面颗粒 |
| cool-duotone | 双色映射，会改变分子身份配色，需重新检查图例 |
| subject-focus | 用 Cycles Object Index 保留指定完整对象颜色，降低其余画面饱和度与亮度 |
| crisp-neutral | 轻微对比、饱和度与曝光调整 |

## 对现有工程应用

```bash
"$BLENDER" -b "$PROJECT/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/apply_render_recipe.py" -- \
  --recipe crisp-neutral --output "$PROJECT/look-crisp" --frame 96 --render
```

输出为新的 `.blend`、PNG 与实际节点/参数记录。原工程保持原样；原目录有 manifest 时保留其中的来源/机制记录，另追加 `display_postprocessing`，记录本次配方及原 `.blend` 哈希，不把旧建模参数冒充新的渲染设置。用 `--recipes /project/recipes.json` 指定项目自己的配方；可以调整滤镜顺序和参数，不必新增脚本。

对已有合成树，工具保留原双视图、背景和取景，将新配方**追加**在其 `Beauty` 后，再叠加原 `Labels`。从同一个基准工程分别生成版本来比较不同配方；对已调色版本再应用会累加效果。原有背景已是图像的一部分，此时配方的 `background` 不替换它；需要换背景/文字明暗时回到场景生成参数。没有合成树的原始场景仍由 `build_stack` 建立基础合成。

已有合成树须提供下面的输出约定；未知树会报错，不能默认取一个 Render Layer 来代替整幅画面。旧工程可用更新后的原生成器重建，或在确知每条支路含义后为项目合成树显式导出这些输出。不要用修改版本标记来跳过检查。

### 跨场景保留布局的输出约定

`publish_composite_outputs(tree, beauty, labels=None, passes=None)` 为当前可编辑节点组增加 `Beauty`（已经完成分屏/背景/取景，尚未叠加受保护文字）、`Labels`（透明底、预乘 alpha 的文字/标注层）和确实可用的 pass；它不改动原 `Image` 输出。Ras 水解/状态序列、K63/AMSH-LP 双视图、局部蛋白结合构建器当前选用的树，以及 `build_stack` 创建的场景提供此约定。项目制作的自定义合成树由 agent 确定这两个实际 socket 后调用，而不是靠节点名字猜测。

- 曝光、对比、饱和度、调色、光晕、暗角、颗粒和双色映射只需合成图像，可接在这些布局之后。
- 描边还需与合成像素对应的 Normal/Alpha；景深雾需 Depth/Alpha；AO 和主体隔离各需对应 pass。原始单视图 `build_stack` 会导出本次配方实际启用的 pass。没有导出时重新准备这些 pass，再应用配方。
- 双视图经缩放、遮罩和混合后，一个原始视图的法线/深度已不对应完整画面。整体输出不宣称这些 pass 已对齐；默认 `--scope layout` 仍会拒绝依赖缺失 pass 的配方。对于已经登记源视图的工程，使用下面的 `--scope panels` 在分屏前分别处理，不能拿左视图的深度给右视图加雾。

主体选择使用实际对象或集合名称，例如 `--recipe subject-focus --engine cycles --subject actin_001`；源合成也必须已有对齐的 Object Index 输出。空集合不算有效选择。不要将膜内外整个群体一起错误选成单个主体。

### 整体图与局部图分别使用自己的通道

`panel_compositor.py` 将指定 Render Layers 的 Image 绑定到原布局的明确输入。`paired_composite` 会登记 `overview` 和 `detail`，因此用当前生成器重建的 PTP1B/AMSH-LP 分子序列、K63/DUB 专用双视图可以提供这一约定。老 `.blend` 或其他自定义树没有登记时要先按其真实支路适配，不能仅填写版本号。当前只接受同一场景的已知 View Layer；其他场景/相机组合需另行适配。

```bash
# 用当前生成器重建后，以中性工程为共同起点。仍保留原完整分子和背景。
"$BLENDER" -b "$PROJECT/base/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/apply_render_recipe.py" -- \
  --scope panels --recipe teaching-outline \
  --output "$PROJECT/outlined" --frame 156 --render

# PTP1B 的示例范围来自本例显示坐标；其他场景改用项目配方。
"$BLENDER" -b "$PROJECT/base/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/apply_render_recipe.py" -- \
  --scope panels --recipes "$SKILL/assets/panel-render-recipes.json" \
  --recipe ptp1b-panel-depth --output "$PROJECT/depth" --frame 156 --render

"$BLENDER" -b "$PROJECT/depth/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/panel_compositor.py" -- --output "$PROJECT/depth/panel-audit.json"
```

配方开头连续出现的 `outline`、`depth_fog`、`ambient_occlusion`、`subject_focus` 在各视图上处理；后续曝光、调色、光晕、颗粒等在原完整布局上处理，然后叠加受保护文字。后半段再次插入依赖 pass 的滤镜会被拒绝，不自动调换顺序。`--panels detail` 仅限制前半段作用的视图；后半段仍作用于完整画面。需要只修改右图时，使用仅含视图滤镜的项目配方，不附加全图滤镜。

深度滤镜必须显式提供 `panel_ranges: {"overview": [near, far], "detail": [near, far]}`。各自使用渲染后的场景坐标；局部图通常放大过，不能把范围解释成共同的分子物理距离。范围是整镜头固定的美术参数，应覆盖关注动作并检查关键时点，不能逐帧自动归一化造成闪烁。光晕接在已合成的画面上，避免被单体 alpha 边界裁掉。

对象强调需 Cycles，且所选的每个视图中都要有实际可见的索引对象。PTP1B 示例可选择 `peptide overview`、`peptide detail`、`phosphoryl overview`、`phosphoryl detail` 四个实体，配合 `ptp1b-substrate-focus`；该例只降低其余主体的亮度、保持饱和度，仍要检查图例对应。不能只在整体图选了对象，就默认局部图存在相同对象。`--subject` 可重复；比较新选择时从同一个中性工程开始。

节点执行“解除预乘 → 在不透明颜色上处理 → 恢复该视图原 alpha → 重新预乘”，再送回原有遮罩和布局。这保留透明边缘的覆盖度，也避免深度雾给空背景铺上不透明颜色。实际原理见 Blender 的 [Alpha Convert](https://docs.blender.org/manual/id/5.1/compositing/types/color/alpha_convert.html)；数据 pass 与模糊/透明度的限制见 [Render Passes](https://docs.blender.org/manual/en/4.2/render/layers/passes.html)。本实现已针对 Blender 5.2.3 检查，并拒绝主相机景深或运动模糊与未模糊数据通道混用；独立背景层的后期虚化可以保留。透明材质和复杂光学效果仍需另做项目视觉检查。

保存检查核对源图与布局输入、每条 pass 的所属视图、所需通道开关，以及 alpha 的转换和恢复。它不替代生物模型检查。分子序列滤镜输出同步保留 `display.json`，可以继续运行 `audit_molecular_sequence_scene.py` 检查完整坐标序列、字体和同步标注。只修改显示时无需重新渲染整段视频。

本次已在 PTP1B 和 AMSH-LP 的已审核分子序列上独立重建并检查双视图；重建复用了已准备的生物学资产与字体缓存，不代表重新下载或重新认证反应机制。中性与原柔和电影配方的对照帧逐像素一致；只修改右视图深度时，左半幅和上下文字区逐像素不变。另检查 PTP1B 的 AO/轮廓、Cycles 对象强调及不同反应时点，AMSH-LP 的轮廓与 1 秒连续局部预览。18 个错误配置被拒绝，包括跨视图深度/alpha、关闭所需通道、错误 alpha 转换、无效范围与缺失主体。未把共享生成器中的其他机制或全部动画组合列为已验证。

### 文字和语义颜色

初建 `build_stack` 时默认将字体对象移入独立的 Typography View Layer，保留其动画和字体，在艺术处理后叠加，文字明暗按背景调整。对已有合成追加配方时直接复用原 `Labels`，不移动字体或重设字号。使用带景深的相机时仍需把字体置于焦平面，或采用独立二维字幕。各项目的 Labels 可能只含字体，也可能含图例色点和分隔线，必须检查实际支路；字体之外的标注不自动得到保护。彩色文字承载变量时，不能把自动中性色字当作科学语义验证。

双色、局部去色、饱和度和调色均可能改变身份颜色。至少对照中性版本，检查 ATP/ADP、亲代/新生链、不同离子和同一跟踪分子是否仍可区分。显微外观只用于教学模型；不得对实验数据图套这些艺术滤镜。

### 软件与视频检查

Blender 5.2 的 Render Layers 插槽名称是 `Ambient Occlusion` 和 `Object Index`，不是旧式缩写。Object Index 主体遮罩目前限定 Cycles；法线、深度和 AO 必须先启用相应 pass。任何新引擎/版本组合都做实际小样。

固定颗粒在生成后打包进本地 `.blend`，不会逐帧随机闪烁；skill 包不包含这些图片。深度范围在构建时固定，运动主体超出该范围会被钳制，摄像机移动较大时应显式给出覆盖全镜头的 `range`。

`render_saved.py --compositor <名字>` 可重新选择保存的 `COMP_STACK_` 分支。重开工程后检查关键时点的透明边缘、光晕吞并相邻粒子、暗角遮挡和文字对比；轮廓闪烁等时间问题用 1–3 秒连续预览检查。完整视频作为交付或待检查问题涉及全片时，再渲染并检查全片。单帧或局部预览成功不代表全部动画已经验证。
