# V 型 ATPase 双视图的参考外观

`build_vatpase_scene.py --style-profile` 将共用角色材质、灯光、独立文字、来源背景和合成配方接入已有连续循环。来源与机制边界沿用 [vatpase-loop.md](vatpase-loop.md)：31 条主体链、10 条完整放大环链、左右同步质子图形及开放底物池。外观调整不改变原网格、平滑修饰器、1.6 倍局部显示倍率或已编排的运动。

[style-vatpase-reference.json](../assets/style-vatpase-reference.json) 是可修改的参考起点：冷色支架、暖色转子、清楚的入口/出口和三层背景。参考用于指导体积、灯光和景深关系；本例保留来源包络，没有把原子球风格当作新结构。按实际参考与用途调整参数，并检查最终像素。

## 输入与作用

- `surface_roles`：`protein`、`membrane_heads`、`membrane_tails`、`chemical_token`、`proton_marker`、`flow_guide`、`environment`，复用 [材质参数](material-and-lighting.md)。ATP/ADP/Pi 与 H⁺ 仍是图形符号。
- `palette`：使用示例中的身份键。同一泵链、转子、质子和入口/出口在整体与局部图保持颜色对应；催化 A 亚基原有微弱高亮节奏保留。
- `lighting/world/background/view_transform/look/exposure`：实际控制光源与显示。面积光的位置、尺寸和功率直接采用章节坐标，不额外缩放。
- `camera.scale_multiplier` 与 `camera.center_offset_nm: [x,z]`：围绕原正面布局缩放/平移正交取景。标签保留既有位置，调整后检查整个循环的主体与文字边界。透视和轨道旋转需要重新设计双视图，入口明确拒绝这些字段。
- `render_recipe`：选择共用 [合成配方](render-recipes.md)，默认 `neutral`。主体和背景经过所选节点，文字独立叠加；配方受引擎与通道支持限制。

不传 `--style-profile` 时，原 `--crowder-dir/--compositor` 调用保持旧外观；新模式使用 `environment/render_recipe`，不混用旧参数。

## 背景的尺度与区室

先用 `prepare_cellular_environment.py` 从核定的 1UBQ、9SPI、8F8R 原始坐标准备原尺寸 nm 资产，方法见 [来源微结构背景](interaction-appearance.md)。传 `--environment-dir`，并在 profile 中声明 `environment.compartment/context_note/context_sources/objects`。

当前版式要求完整背景处于胞质侧，最低点高于膜头基外缘 z=2.30 nm；还须避开整个前景运动范围，并位于整体泵一侧 x≤2 nm。右侧是放大的孤立环，不给它混入同尺度物理背景。检查完整对象，不通过切割或缩放单体使其通过。

环境包含人源核糖体/泛素和兔肌动蛋白代理形状；布局不代表大鼠囊泡周围的实测组成或浓度。换区室、物种或版式时重新选择来源与空间约束。

`environment.compositor` 控制三层不透明度、按输出宽度换算的模糊及背景色。主体和文字保持清晰。无背景项目可删除 `environment` 并省略 `--environment-dir`。

## 生成与核验

沿用已核定泵资产和本机合法可用的宋体/Times New Roman；只改外观无需重复准备结构。输出为 skill 外的新目录。

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_vatpase_scene.py" -- \
  --asset-dir "$PROJECT/pump" --environment-dir "$PROJECT/environment" \
  --style-profile "$SKILL/assets/style-vatpase-reference.json" \
  --font-zh "$FONT_ZH" --font-latin "$FONT_LATIN" \
  --output "$PROJECT/styled" --engine cycles --width 1600 --samples 32

"$BLENDER" -b "$PROJECT/styled/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_vatpase.py" -- \
  --asset-dir "$PROJECT/pump" --output "$PROJECT/styled/validation"

"$BLENDER" -b "$PROJECT/styled/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/styled/keys" --frames 53,145,278
```

工程旁的 `appearance.json` 保存输入快照；它的内容哈希、保存工程与 manifest 必须一致。修改 profile 后在新目录重新构建，不把旧工程或旧图片当作新配置。

保存检查读取实际使用的材质、身份色、包络平滑、光源、相机、背景边界/源几何、合成通路、字体与 385 个时点的催化高亮；原完整网格、定子、左右相位及循环接缝检查继续执行。文字从独立图层读取实际评估坐标，局部图的源网格也纳入画内边界检查。

在带背景的新工程中运行 `test_vatpase_appearance.py --output /absolute/report.json` 可检查无效请求和定向破坏，恢复正例，不保存被修改的工程。这不代表所有光学参数、所有配方或原子级轨迹均已验证。

实际查看少量正式尺寸图，检查 H⁺ 与金色转子的对比、两图身份颜色、膜下运动和文字。循环边界只补必要的原速短预览；全片渲染用于成片交付。源码分发复验保留原始输入重建、完整姿态检查和有针对性的实际像素检查。
