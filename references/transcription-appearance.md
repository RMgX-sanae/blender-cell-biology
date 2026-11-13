# 转录两轮核苷酸加入的参考外观

用于已有 5FLM 单轮 C 加入和 23 秒 C→U 教学片段的美术迁移。先按 [制作流程](production-workflow.md#转录一次-c-残基加入) 准备来源、CTP 和第二轮 UTP；保留原始核酸、完整 Pol II、底物身份、账本及已认可的事件节奏。它不补齐起始、终止、回收或金属催化。

从 [外观配置](../assets/style-transcription-reference.json) 复制项目参数。参考图片提供暖色精细前景、冷色结构背景及克制的景深；配色、灯光和虚化不是实验光学数据。整组 Pol II 后移作为相关结构背景，非模板 DNA 整条隐藏；不切开单体，不添加无身份球团。

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_transcription_scene.py" -- \
  --asset-dir "$PROJECT/transcription_asset" --substrate-dir "$PROJECT/transcription_substrate" \
  --style-profile "$SKILL/assets/style-transcription-reference.json" \
  --font "$FONT_ZH" --font-latin "$FONT_EN" \
  --engine cycles --width 1600 --samples 32 --output "$PROJECT/single"
"$BLENDER" -b "$PROJECT/single/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_transcription.py" -- \
  --asset-dir "$PROJECT/transcription_asset" --substrate-dir "$PROJECT/transcription_substrate" \
  --output "$PROJECT/single/audit.json"
"$BLENDER" -b "$PROJECT/single/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/extend_transcription_scene.py" -- \
  --asset-dir "$PROJECT/transcription_asset" --substrate-dir "$PROJECT/transcription_two_cycle_substrate" \
  --output "$PROJECT/two-cycles"
"$BLENDER" -b "$PROJECT/two-cycles/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_transcription_cycles.py" -- \
  --asset-dir "$PROJECT/transcription_asset" --substrate-dir "$PROJECT/transcription_two_cycle_substrate" \
  --output "$PROJECT/two-cycles/audit.json"
"$BLENDER" -b "$PROJECT/two-cycles/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/two-cycles/keys" --frames 97,181,325,457,517
```

`FONT_ZH` 为宋体常规，`FONT_EN` 为 Times New Roman 常规。省略 `--font` 会使用英文标签；新外观必须提供 `--font-latin`。第二轮继承单轮工程已经打包的字体、标签语言、材质、灯光与合成，不需要重新指定字体。修改风格时从单轮重新生成，避免把新请求写在未改变的工程旁边。字体只随项目工程保存，不分发进 skill。

| 配置 | 实际对象及作用 |
|---|---|
| `surface_roles.protein` / `palette.heads` | 12 个完整已建模 Pol II 蛋白链 |
| `surface_roles.nucleic_acid` | DNA、RNA、进入的核苷酸核心及已有的连接表示 |
| `surface_roles.phosphate` / `palette.accent` | 两轮 β/γ 磷酸与离开的 PPi、示意断键连接 |
| `palette.dna2/dna1` | 模板/非模板 DNA；非模板链保持整条隐藏 |
| `palette.rna/protein` | RNA 与当前进入的核苷酸；`protein` 是沿用的色板键名，CTP/UTP 仍采用核酸材质 |
| `palette.ink` | 独立叠加的标签和极性引线 |
| `lighting/world/exposure/view_transform/look` | 实际灯光、世界填充和显示变换，不追加旧的固定补光 |
| `camera.scale_multiplier/center_offset_nm` | 初始正交尺度及相机 right/up 平移；第二轮保留原有轻微拉远 |
| `context_blur_pixels_1280` | 完整 Pol II 背景层的轻度虚化，按输出宽度缩放；0 为清晰 |
| `render_recipe` | 共用合成配方，背景、主体和保护文字分别处理 |

第一轮 C20 在转位时改成 RNA 颜色，进入色继续指示 UTP；两轮 PPi 保持自己的身份。不要让材质变更覆盖这段颜色动画。默认没有来源脚注、解释性小字或文字描边；证据边界保留在 manifest 和配套文档，画面只留过程标签、极性与物质关系。

新外观不接收旧 `compositor`、`representation` 或任意相机方向，命令行不同时使用 `--compositor`。放大、平移后仍要检查完整主体和文字边界，不能自动删去原子以适应镜头。旧的无角色配置保留既有调用方式；旧中文工程的两轮扩展继续可用。

翻译与转录共用 `layered_molecular_appearance.py` 的材质、字体与合成操作，各自保留分子角色和机制规则。`audit_transcription.py` 和 `audit_transcription_cycles.py` 对新工程额外调用外观检查：实际角色/着色器、语义颜色、C20 颜色衔接、镜头、灯光、背景路由、字体及每帧字幕事件和边界。检查分层场景时，使用对应渲染层的依赖图求值对象；被当前层排除的原对象可能仍带着旧帧状态，不能据此判断实际位置或可见性。

RNA 的 3′ 标签置于离开的 PPi 路径另一侧，并用短引线指向末端。`audit_transcription_annotations.py` 把引线起点逐帧对照来源 O3′ 坐标和既有的连接、转位事件；单轮转为两轮时也从文字层读取实际求值位置，不能捕获预览帧的陈旧位置。

`test_transcription_appearance.py` 在内存中破坏第二轮颜色动画、字幕路由等实际状态，并要求检查器拒绝；不会保存被修改的工程。角色、路径和计量检查分别进行，不能用外观通过代替原子和反应账本核对。

实际像素另查前景与背景的主次、核苷酸/PPi 的可辨认性、标签与结构的分离，以及两轮之间的颜色衔接。少量关键帧配合已保存运动的比较即可验证纯外观迁移；改变运动时再检查相应短预览。合成配方的运行支持不等于每个配方都已完成审美验收。
