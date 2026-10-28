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
  --recipe depth-layered --output "$PROJECT/look-depth" --frame 96 --render
```

输出为新的 `.blend`、PNG 与实际节点/参数记录。原工程保持原样。用 `--recipes /project/recipes.json` 指定项目自己的配方；可以调整滤镜顺序和参数，不必新增脚本。主体选择使用实际对象或集合名称，例如 `--recipe subject-focus --engine cycles --subject actin_001`。空集合不算有效选择。不要将膜内外整个群体一起错误选成单个主体。

### 文字和语义颜色

默认将字体对象移入独立的 Typography View Layer，保留其动画和字体，在艺术处理后叠加。文字明暗按背景调整；这样浅底不会保留无法阅读的白字，也不会给文字加颗粒、景深或边缘滤镜。使用带景深的相机时仍需把字体置于焦平面，或在项目中采用独立的二维字幕渲染。彩色文字本身承载变量时，按项目保留独立图例并复查，不能把自动中性色字当作科学语义验证。

双色、局部去色、饱和度和调色均可能改变身份颜色。至少对照中性版本，检查 ATP/ADP、亲代/新生链、不同离子和同一跟踪分子是否仍可区分。显微外观只用于教学模型；不得对实验数据图套这些艺术滤镜。

### 软件与视频检查

Blender 5.2 的 Render Layers 插槽名称是 `Ambient Occlusion` 和 `Object Index`，不是旧式缩写。Object Index 主体遮罩目前限定 Cycles；法线、深度和 AO 必须先启用相应 pass。任何新引擎/版本组合都做实际小样。

固定颗粒在生成后打包进本地 `.blend`，不会逐帧随机闪烁；skill 包不包含这些图片。深度范围在构建时固定，运动主体超出该范围会被钳制，摄像机移动较大时应显式给出覆盖全镜头的 `range`。

`render_saved.py --compositor <名字>` 可重新选择保存的 `COMP_STACK_` 分支。重开工程后渲染关键时点，并检查完整编码视频的轮廓闪烁、光晕吞并相邻粒子、暗角遮挡和文字对比。单帧滤镜成功不代表所有场景和动画已经验证。
