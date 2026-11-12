# 肌动蛋白章节的参考外观迁移

`build_actin_films.py --style-profile` 将参考外观接入已有完整案例。它复用角色着色器、合成配方和来源背景，不改变原有结构、机制范围或教学动作。原来不传此参数的调用保持旧外观。准备来源与完整路线见 [portable-actin-films.md](portable-actin-films.md)。

[style-actin-reference.json](../assets/style-actin-reference.json) 是基于已观察参考图的可编辑起点：暖色重点、冷色环境、柔和体积和分层虚化。它保留当前完整链包络，并不复制参考图的原子球表示或宣称取得作者的原始材质。先看项目参考，按实际画面调整，不把这份色板当作所有细胞过程的默认答案。

## 输入与实际作用

- `surface_roles`：`protein`、`cofactor`、`schematic_symbol`、`schematic_linker`、`membrane_symbol`、`environment`。按对象功能路由，同一种身份色可以用于不同表面；着色器参数见 [material-and-lighting.md](material-and-lighting.md)。膜符号和 FH1 供料箭头仍是示意，材质不能升级其结构证据。
- `palette`：覆盖已有语义色键 `teal/blue/gold/pink/green/red/cream/dim/ink`，保留章节内身份对应。图示框内的深色文字单独保护，其他文字沿用指定前景色；最终像素仍须检查对比。
- `lighting`：每盏面积光明确给出 `name/position/power/color/size`，位置和尺寸直接使用当前章节坐标，不暗中再次缩放。`world/background/view_transform/look/exposure` 控制环境与显示。未提供灯光时沿用已有章节光源尺度。
- `camera`：已有正面教学布局支持正交 `scale_multiplier` 和 `center_offset_nm: [x,z]`；标题与步骤随取景重排。透视、轨道和任意镜头旋转需要重新设计文字与框图，此入口明确拒绝，不能接收后忽略。缩放或平移后检查完整运动取景与动态标签。
- `render_recipe`：选择 [render-recipes.json](../assets/render-recipes.json) 中的配方 ID，默认 `neutral`。主体与背景经过所选节点，文字独立叠加；具体配方仍受引擎和所需通道限制。此模式不要同时传旧 `--compositor`。
- `chapters`：仅按章节覆盖 `camera/environment`。示例只给已审阅的 `a03/f04` 配置来源背景，其余章节使用纯底，避免把未核对的环境套到上游框图或放大配体视图。

## 背景、字体与章节交接

使用来源背景时，传 `--environment-dir` 指向 `prepare_cellular_environment.py` 的原尺寸 nm 输出，配置 `environment.compartment/context_note/context_sources/objects`。实例可选择核糖体、肌动蛋白和泛素；旋转、位置及重复编号由项目指定，不能根据装饰需求任意放大。对象必须避开整个主体运动范围。布局、数量和远中近三层虚化是艺术示意，不是实测浓度、动力学或显微点扩散函数。字段和准备命令可参考 [interaction-appearance.md](interaction-appearance.md)。

`environment.compositor` 可设置三层 `blur_pixels_1280/opacity` 及背景色；模糊随构建输出宽度变化。完整来源形状仍保留，背景独立于清晰主体和文字。

中文宋体和英文/数字 Times New Roman 沿用原字体路由并嵌入工程。文字保护保留框内深色字、普通前景字和原定字幕时序。`a03` 的章节淡出同时作用于画面和文字；它是明确的教学镜头切换，不是无缝运动。批量构建时每章重新建立对象和视图层，避免继承上一章被排除的背景或文字。

## 调用与核验

```bash
"$PYTHON" "$SKILL/scripts/build_actin_films.py" \
  --output "$PROJECT" --steps build,audit --shots a03,f04 \
  --quality final --engine cycles \
  --style-profile "$SKILL/assets/style-actin-reference.json" \
  --environment-dir "$PROJECT/environment"

"$BLENDER" -b "$PROJECT/shots/a03/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/review/a03" --frames 116,360,576
```

上述调用要求项目已按完整路线准备来源和动作，且背景已单独准备；第一次运行还需指定实际 Blender、FFmpeg 和字体路径。`final` 使用 1600×900、24 fps、32 samples；关键帧编号应根据实际章节时长选择。

只改外观时复用已核定的准备资产，执行 `build,audit`，无需重复耗时的 `prepare`。来源或动作准备规则改变时再按对应流程重算；发布前的独立原始输入重建与日常调色迭代分开执行。

项目保存带内容哈希的外观快照和背景资产签名。改变外观、背景资产、质量或帧率后需重新构建受影响章节；修改快照或把旧场景/旧帧当作新配置时拒绝。显式提供新配置可进入下一次构建，既有快照不覆盖。删除所有背景配置时也解除旧环境目录的继承。

`audit` 与逐章渲染入口都会重开工程，检查实际角色/身份色、着色器、光源、镜头、背景完整形状、文字层路由、字体和所有时点的标题/步骤；同时保留既有完整网格、姿态和严格接触检查。`audit_appearance.py` 可单独检查外观；`test_appearance.py` 对带背景工程做定向破坏并确认拒绝，之后恢复，不写回工程。配置快照与缓存行为可用普通 Python 运行 `actin_films/test_appearance_config.py`。

检查结果不能替代看图。先查看少量正式尺寸关键帧；只有验证淡出或运动确有需要时再补短区间。迁移验收分别记录来源/动作保持、实际外观、像素可读性及独立重建范围，不能由两个章节外推整套机制或所有风格组合。
