# 让不同反应共用分子序列显示

`build_molecular_sequence_scene.py` 将已经检查的原子坐标序列和项目显示 JSON 转成双视图工程。它复用原有原子动画、整分子包络、源结构背景、字体和合成器；不推断反应、生成运动或修补碰撞。PTP1B 与已有 AMSH-LP 数据用于检查这个接口，原有机制专用生成器继续可用。

## 准备数据的约定

在 skill 外的项目资产目录提供以下文件。已有准备格式不相同时，在项目中显式适配，不能只改字段名就宣称单位和身份正确。

| 输入 | 内容 |
|---|---|
| `manifest.json` | `passed` 和 `checks.passed` 为真；`frames`、`fps`、`units: "nm"`、`scope`；每个组件的 `atom_ids` 和文件名到 SHA-256 的 `data_sha256` |
| `motion.npz` | 每个组件名对应有限坐标数组 `(frames+1, atom_count, 3)`，单位 nm；全部所选原子 ID 跨组件唯一 |
| 可选键文件 | NPZ 的 `bonds` 数组，使用对应组件的零起点原子索引；由项目说明是化学键还是显示导引线 |
| `background_<name>.npz` | 与显示配置相符的完整源结构形状，`vertices` 为 nm、`faces` 为三角面；来源和美术分布另行说明 |
| 反应身份账本 | 有化学变化时提供 `reaction-ledger.json` 或 `reaction-sequence.json`，列入哈希清单；保存检查器将原子集合与显示数据核对 |

渲染帧 1…F 对应 `time=(frame-1)/F`；第 F+1 个坐标样本用于端点检查，不渲染。端点存在不表示无缝循环。旧准备数据没有单位字段时，显示配置须有 `units: "nm"` 及说明转换依据的 `unit_note`；AMSH-LP 的依据是其准备脚本将内部 Å 坐标乘以 0.1，不能据此猜其他旧缓存的单位。

准备数据可提供 `captions`（`[start,end,text]`，0–1 半开时段）。有该字段时，显示字幕的时间边界必须相同，文字可以翻译；没有该字段的旧数据，项目需自行核对字幕和事件。PTP1B 从事件时间生成这些边界，修改节奏后不能继续用旧显示时序。

## 显示配置改什么

从 [PTP1B 显示配置](../assets/ptp1b-display.json) 复制项目文件，只改变有理由改变的部分。

- 两个 `views` 各有正交右手 `basis_columns`、`center_nm`、水平 `span_nm` 和 `actors`。左图须保留所有组件；两图使用同一批原子序列，只做刚体视图变换和均匀缩放。`keep_visible` 指定必须完整落在取景框内的主体；离开的产物仍保留在数据中。
- 组件可用 `surface`、`sticks` 或 `spheres`。线棒来自 `bond_file` 或显式 `bonds`；`radius_nm`、`bond_radius_nm` 为显示尺寸。PTP1B 的 `enzyme_context.npz` 是完整 Cα 骨架导引线加选定残基的化学键，**不是完整共价图**。原始原子点全部保留，取景窗口裁的是画面，不对分子作布尔切割。
- `colors.dark/light` 保持实体语义身份；`color_key` 可选用同一风格色板内的强调色。默认让主体与背景有较高明度对比，背景形状压暗并虚化。改变风格、灯光或合成器后重新检查实际像素，不能仅以十六进制颜色判断。
- `highlights` 和 `bond_annotations` 追踪指定原子坐标，在前景层加粗显示；每个视图只在自己的窗口和声明时段内显示标记。它们帮助看清反应位置，不表示新的原子或实测键级。`legend`、主标签和字幕用宋体/Times New Roman，不加来源小字。
- `background_components` 指向项目的背景形状；数量、景深和摆位是美术设计，不推断细胞丰度或真实结合伙伴。

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_molecular_sequence_scene.py" -- \
  --asset-dir "$PROJECT/assets" --display "$PROJECT/display.json" \
  --output "$PROJECT/scene" --font-zh "$FONT_ZH" --font-en "$FONT_EN" \
  --style dark --recipe neutral --engine eevee --width 1280 --samples 24 --render
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_molecular_sequence_scene.py" -- \
  --output "$PROJECT/scene/audit.json"
```

检查器核对准备文件/显示配置指纹、原子身份、全部保存坐标、视图变换、显示边、标记位置与可见性、字体和字幕。它在内存中关闭生成几何的 viewport modifiers 来检查原始动画点；不会保存这一审计状态。表面、线宽、遮挡和颜色必须另外看真实渲染。深/浅色、Eevee/Cycles 和合成配方独立选择；一个已验收组合不能替其他组合通过。

从保存工程继续使用 `render_saved.py`、`encode_video.py` 和 `audit_video.py`。所有 NPZ、工程、帧和视频保留在项目目录；skill 只保留生成器、参数及来源说明。
