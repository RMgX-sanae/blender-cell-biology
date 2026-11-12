# 相互作用场景的参考外观交接

在已经核定的 `prepare_interaction.py`、结合候选或 AP2 膜侧资产上使用；不改变准备好的完整伙伴、伴随组分、运动矩阵和接触门槛。先按 [参考图片](style-from-reference.md) 分析画面，再选择适用字段。纯底适合孤立实验构建体；有细胞背景时先核对区室、来源物种和表示边界。

旧 `--surface-profile` 仍只读取 `surface_roles`。新的 `--style-profile` 显式接入材质、各 actor 的颜色与标签、灯光、相机、合成配方和可选三维环境；二者不能同时使用。它有自己的输入字段，不能把其他入口的 JSON 原样传入后假定所有字段都生效。未知顶层字段会报错。

## 两个已有来源的起点

```bash
# 完整实验构建体控制：纯底，全部准备姿态留在可编辑工程中。
"$BLENDER" --disable-autoexec -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_interaction_scene.py" -- \
  --asset-dir "$PROJECT/candidates/interaction" --output "$PROJECT/assembly-look" \
  --style cinematic --engine cycles --orbit-degrees 150 \
  --style-profile "$SKILL/assets/style-interaction-assembly.json" \
  --font-en "$FONT_EN" --width 1600 --samples 48

# 已开放 AP2 的局部货物识别：沿用 prepare_ap2_capture.py 的核定膜侧资产。
"$BLENDER" --disable-autoexec -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_interaction_scene.py" -- \
  --asset-dir "$PROJECT/ap2-assets" --output "$PROJECT/ap2-look" \
  --style teaching --engine cycles \
  --style-profile "$SKILL/assets/style-interaction-ap2.json" \
  --environment-dir "$PROJECT/environment" \
  --font-zh "$FONT_ZH" --font-en "$FONT_EN" --width 1600 --samples 48

"$BLENDER" --disable-autoexec -b "$PROJECT/ap2-look/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_interaction_scene.py" -- \
  --output "$PROJECT/ap2-look/saved-check.json"
```

第一组输入来自 [完整伙伴控制](source-assembly-instances.md)：3SN6 保留 T4L、Nb35、P0G，是来源约束候选，不是天然招募。第二组沿用 [膜侧 AP2](protein-interactions.md) 的 6YAF/EMD-10748 准备；它从已开放构象开始，不包含激活、原子级 PIP2、完整内吞或真实动力学。背景资产按 [细胞质环境](cellular-environment.md) 准备，不把它们的存在解释为该实验原位观测。

默认不渲染。需要静态局部检查时用 `--snapshot-frame`，其取景和检查范围只有该姿态；整段工程的保存检查则覆盖全部准备姿态及原有接缝样本。用 Blender `-o /absolute/project/key_#### -f 121` 渲染选定关键帧；`--render` 会输出入口原有五个代表时点，不生成整片视频。

## 外观输入

- `surface_roles` 采用 [已有着色器](material-and-lighting.md) 的参数；聚合物按来源类型选择 `protein/nucleic_acid`，配体选择 `cofactor`。膜为 `membrane_heads/membrane_tails`，未解析的教学连接线另用 `schematic_tether`；背景用 `environment`。这是艺术材质，不是实测光学属性。
- `actor_colors`、`actor_labels` 按实际 actor ID 设置，图例与主体使用相同颜色。`--partner-colors` 与逐 actor 配色不能同时指定。`representation` 为 `surface/atoms`，显式 CLI 表示选项优先；原子球仍是显示半径。
- `membrane_palette` 可选 `heads/tails/tether`。只改变已有完整膜符号及连接线的颜色，不改变厚度、朝向和来源约束。
- `lighting` 是具名面积灯列表，每项含 `name/position/power/color/size`。位置是源场景世界坐标 nm，显式灯光不随 `--orbit-degrees` 暗中旋转。`world/view_transform/look/exposure/background` 控制实际世界光和显示变换。
- `camera` 支持 `type=ORTHO/PERSP`、`position/target`、`scale`、`lens`、`fstop`、`roll_degrees`。位置和目标必须一起提供，不能与非零 `--orbit-degrees` 混用；它们不移动源分子。正交用 `scale`，透视用位置及 `lens`；不接受无视觉作用的投影参数组合。未指定位置时沿用来源视角及完整路径的自动取景。`fstop=0` 关闭主体光学景深，背景仍可单独虚化。
- 取景检查全部可见源对象的包围盒，为标题和图例留出区域；过小比例或错误透视位置明确失败，不通过缩小分子或隐藏组分修补。镜头光学仍为艺术设定。
- `render_recipe` 选择 `render-recipes.json` 中的现有配方，默认 `neutral`；通道和引擎限制沿用共用合成器。字体最后独立叠加。涉及调色的配方必须复查图例，不外推到所有组合。

中文标签要求显式提供宋体 `--font-zh`；中文中的英文、数字及所用希腊字母由 Times New Roman 字体槽显示。纯英文使用 `--font-en`。生成器将选定外观原文件复制到输出的 `style-profile.json`，并记录输入哈希、解析后相机、灯光及参数；后续修改仓库配方不会破坏旧工程的输入快照。

## 有来源的三维背景

`environment` 与 `--environment-dir` 必须同时提供，不能和旧 `--crowder-dir` 混用。当前支持已有 `schematic-cytosolic-environment` 原尺寸 nm 资产；其核糖体、泛素和肌动蛋白形状有明确来源，包含既有兔肌动蛋白代理的边界。其他区室应先建立适用资产，不能把 `compartment` 改名后沿用这一组。

环境中声明 `compartment: cytosol`、`context_note`、`context_sources` 和 `objects`。每个对象提供 `asset`、`depth=far/middle/near`、世界坐标 `position_nm`、`rotation_degrees`；肌动蛋白可用 `repeat_indices` 沿准备文件中核定的重复变换放置完整亚基。所有对象保持源形状与原尺寸，没有投影缩放。

`palette` 按源资产名覆盖颜色；`compositor` 可设置 `background`、三个 `blur_pixels_1280` 和三个 `opacity`。背景底色同时用于文字对比。主体和文字不进入背景虚化；三层合成是艺术景深近似，不是显微 PSF。背景为静态布局，不生成布朗运动。

背景包围盒必须避开该工程全部保留姿态的主体范围；这是保守的空间筛查。存在核定膜面时，背景完整包络必须在胞质侧外叶头部之外。失败应调整有依据的布局，不能拉伸、裁切背景或放宽膜边界。位置、数量和拥挤程度仍是示意，不是实测浓度。

## 验收

`audit_interaction_scene.py` 保留原来的源几何、身份、共同运动、膜约束及全部姿态检查，并增加实际材质分配、身份色和图例、请求灯光与镜头、字体、源背景形状/位置及真实景深连接。背景在其参与渲染的层中求值。保存检查不能替代最终尺寸的画面复核。

`test_interaction_appearance.py --output "$PROJECT/appearance-negative.json"` 在已打开的、带三层背景和中文标签的工程里做实际破坏检查并恢复，禁止保存破坏状态。新目标仍须独立核定输入、检查画面，并从源码包重建；案例通过不证明其他蛋白、任意环境或完整机制。
