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

两个现有生成入口已经消费这些参数：

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

不是所有生成器都读取 `surface_roles`。新项目先确认所用入口是否实际连接设置，再在保存的 `.blend` 中检查节点和最终渲染；只修改 JSON 不算实施。`build_molecular_scene.py` 按已定义的蛋白、核酸、膜部件与背景角色路由；来源对照入口按 mmCIF 聚合物类型路由。背景实例的完整结构保留原尺寸，位置、数量和景深只是展示安排，不能称为实测细胞组成。

## 逐项目调节并检查

1. 复制外观配置到项目，观察参考的体积、表面细节尺度、亮暗关系和环境层次。主角/相邻组分/背景分别选色，颜色应跟随实体身份。深底配较亮主体，低亮度冷色用于环境；避免所有对象同样鲜亮。
2. `lighting` 中明确主光、补光、轮廓光的位置、尺寸、颜色和功率。主光塑造体积，补光保留暗处结构，轮廓光分离主体与环境。更改景别时同步检查光源尺寸/距离；不要靠过曝代替对比。世界填充与合成底色是不同设置。
3. 先用中等分辨率查看真实像素。背景不应整齐排列、左右镜像或抢占主体；选与区室和景别相容的完整来源结构。景深不能模糊关键反应位点。静态来源对照可明确采用纯色背景，以免引入无依据的细胞情境。
4. 保存源几何、材质参数、灯光、相机与合成设置后，再提高分辨率和采样。推荐只渲染足以判断当前问题的静帧；只有交付需要时输出高清完整图或影片。高清不会自动修复构图、低质量形体或错误机制。
5. 分别记录：来源与几何是否保持；节点是否使用所选参数；最终像素中哪些轮廓、暗部、相邻实体与背景关系清楚。没有查看的风格、相机和动画时段保持未验证。

中文宋体、英文/数字 Times New Roman；成图不堆来源小字。工程、预览与高清图存项目目录，skill 仅分发生成代码、参数与说明。外部参考只作参考，不把 AI 增强图称为 Blender 渲染或结构验证。
