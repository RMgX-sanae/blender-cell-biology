# 渲染滤镜与合成器组合

先完成结构、摄影和动画，再选择合成；滤镜不能修复错误拓扑、缺失机制或严重穿插。当前入口是 Blender 原生可编辑节点，外部 Fusion/Natron/After Effects 仍按项目确认可用并实际执行后才能声明支持。

`assets/render-recipes.json` 定义有序滤镜链，`compositor_stack.py` 实际构建节点。材质风格、渲染引擎、动画模式和合成配方相互独立。

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
- 双视图经缩放、遮罩和混合后，一个原始视图的法线/深度已不对应完整画面。目前双视图不导出这些 pass，所以依赖它们的配方会明确失败；不能拿左视图的深度给右视图加雾。若项目需要，分别处理各视图，或同步合成全部相关 pass 后再声明可用。

主体选择使用实际对象或集合名称，例如 `--recipe subject-focus --engine cycles --subject actin_001`；源合成也必须已有对齐的 Object Index 输出。空集合不算有效选择。不要将膜内外整个群体一起错误选成单个主体。

### 文字和语义颜色

初建 `build_stack` 时默认将字体对象移入独立的 Typography View Layer，保留其动画和字体，在艺术处理后叠加，文字明暗按背景调整。对已有合成追加配方时直接复用原 `Labels`，不移动字体或重设字号。使用带景深的相机时仍需把字体置于焦平面，或采用独立二维字幕。各项目的 Labels 可能只含字体，也可能含图例色点和分隔线，必须检查实际支路；字体之外的标注不自动得到保护。彩色文字承载变量时，不能把自动中性色字当作科学语义验证。

双色、局部去色、饱和度和调色均可能改变身份颜色。至少对照中性版本，检查 ATP/ADP、亲代/新生链、不同离子和同一跟踪分子是否仍可区分。显微外观只用于教学模型；不得对实验数据图套这些艺术滤镜。

### 软件与视频检查

Blender 5.2 的 Render Layers 插槽名称是 `Ambient Occlusion` 和 `Object Index`，不是旧式缩写。Object Index 主体遮罩目前限定 Cycles；法线、深度和 AO 必须先启用相应 pass。任何新引擎/版本组合都做实际小样。

固定颗粒在生成后打包进本地 `.blend`，不会逐帧随机闪烁；skill 包不包含这些图片。深度范围在构建时固定，运动主体超出该范围会被钳制，摄像机移动较大时应显式给出覆盖全镜头的 `range`。

`render_saved.py --compositor <名字>` 可重新选择保存的 `COMP_STACK_` 分支。重开工程后渲染关键时点，并检查完整编码视频的轮廓闪烁、光晕吞并相邻粒子、暗角遮挡和文字对比。单帧滤镜成功不代表所有场景和动画已经验证。
