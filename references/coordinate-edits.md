# 分子局部补建、全装配核对与同视角对照

已有实验结构需要补一个缺失原子或调整少数坐标时，先确认缺失原因、构建体化学和适用模板。保留原装配作为独立对照，拟合前为每个可移动或可补建的位点写出依据和权限。相邻肽链、配体、金属和水也进入原子清单。预期有整段缺失、共价反应或装配伙伴整体运动时，另行建立相应的构建或运动模型。[位点与来源核对](site-link-checks.md)

## 给新对象建立检查

工作单位为 Å。每个模型原子的身份依次为 `[label_asym_id, operator, auth_asym_id, auth_seq_id, insertion_code, component, atom]`。使用字符串，空插入码用 `""`。装配操作使用原 mmCIF 的名称，组合操作采用 `(X0)(2)` 这类写法。旧六列原子清单仅适用于全部链都只受编号 1 的恒等操作作用的装配。[装配实例与作者编号](source-assembly-instances.md)

1. 从指定版本的原始 mmCIF 和 CCD 核对缺失位置。固定 `assembly_id`、原文件 `source_sha256` 和所用模型。当前读取一个沉积模型内占有率大于零、altloc 为空或 A 的重原子。多个沉积模型、特殊替代构象和同位素氢需单独处理。
2. 写项目自己的 `policy.json`。除来源字段外，`moves` 每项含 `id` 和 `max_displacement_A`，`additions` 每项含 `id` 和 `element`。两个列表均可为空。添加位置必须属于来源装配里已有的残基实例。按文献和局部结构决定允许的位移，只移动被点名的拷贝。
3. 对需要保留方向的化学关系，可在 `preserve_planes` 加入 `{"id":"ester","plane":[…],"site":[…],"max_height_change_A":…}`。用三个不同原子的完整身份定义有向平面，顺序依次为原点、第一方向、第二方向，`site` 为第四个来源原子。比较该原子相对于两个状态各自平面的有符号高度。原结构本身略偏离平面时，保留该偏差。结合具体构型另外检查成键、糖环扭转、手性、溶剂和全装配接触。跨残基连线需依据 CCD 或 `_struct_conn`。
4. 分发前对已写入磁盘的坐标运行独立检查。准备一个无 pickle 的 NPZ，`atom_ids` 是模型的字符串身份表，`elements` 逐行对应，`model_A` 是所有模型原子的 N×3 坐标，`source_A` 按相同顺序跳过新增位点。源与模型均在原装配坐标系，支持 float32、float64 和整行重排。固定原子核对时计入存储精度。

```bash
"$PYTHON" "$SKILL/scripts/audit_coordinate_edits.py" \
  --source "$PROJECT/sources/target.cif" --candidate "$PROJECT/candidate/coordinates.npz" \
  --policy "$PROJECT/policy.json" --output "$PROJECT/candidate/audit.json"
```

核对程序从原 mmCIF 的装配生成表和 atom_site 重建源位置。报告记录来源、装配、坐标文件指纹、每个实际移动或新增的原子、固定原子数和有向平面测量。缺失伴随组分、身份重复、额外移动、行错配和越界位移均会报错。对新蛋白沿用身份与坐标检查，重新审查自己的化学模型和容差。网格简化、视图中省略原子应发生在此步骤之后。

## 在 Blender 比较两个状态

`view.json` 选择实际需要观看的局部。根据位点在三维中的关系选视角，使关键原子在投影里也能分开。两侧使用同一相机、比例和材质，面板只能整体平移。新加原子只出现在模型侧。选择结构片段时在项目说明中记录两端的截取位置。

| 项目字段 | 内容 |
|---|---|
| `source_sha256`、`assembly_id` | 与刚核过的原始装配一致 |
| `sites` | 名称到 `{"id":[七个身份字段],"role":"ribose"}` 的映射。每种原子身份只能入选一次 |
| `roles` | 每个角色的 `color` 用六位 sRGB 色码，`radius_A` 为球棍表示的绘制半径。`surface_role` 可指定 profile 中的 `protein`、`nucleic_acid` 等材质，细调粗糙度、纹理与次表面散射。未指定时按角色同名查找，再用 profile 的基础材质。两侧相同位点共用球体网格和材质 |
| `bonds` | 每条键给两个已选名称。使用审查过的连接图，两端都存在才绘制。水、金属配位和空间接近应单独解释 |
| `toward_camera`、`up` | 两个不平行的三维向量，指向观察者和画面上方。显示入口将其正交化，保留同一个右手坐标系 |
| `titles`、`legend` | 两个简短的状态名称。可选图例每项含角色 `role` 与可读的 `text`，只标出需要辨认的结构或差异 |

```bash
"$BLENDER" -b --factory-startup --disable-autoexec --python-exit-code 1 \
  --python "$SKILL/scripts/build_coordinate_edit_view.py" -- \
  --candidate "$PROJECT/candidate/coordinates.npz" --audit "$PROJECT/candidate/audit.json" \
  --view "$PROJECT/view.json" --output "$PROJECT/paired-view" \
  --font-zh "$FONT_ZH" --font-en "$FONT_EN" --engine cycles --compositor clean --render
```

默认保留可编辑的原子、键、两个面板的刚体变换与已选位点身份。可换用 `--engine eevee` 或其他有依据的灯光材质 `--profile`。左右面板共用项目指定的正交视角和画面范围。所选 profile 的材质、灯光和色彩设置仍生效。合成分支可选 `clean`、`reference` 或 `monochrome`。使用单色分支时，另行核对两个状态是否仍能靠位置和形状辨认。验收时重新打开保存工程，查看关键原子、邻近水、图例及中心分隔线的投影和对比。表格和来源记录留在项目内。

从磁盘打开工程后，可用相同输入批量核对两侧的原子身份、位置、显示半径、键端点、面板变换、角色着色器和嵌入字体。自定义外观时，两条命令传同一个 `--profile`。

```bash
"$BLENDER" -b "$PROJECT/paired-view/scene.blend" --disable-autoexec --python-exit-code 1 \
  --python "$SKILL/scripts/audit_coordinate_edit_view.py" -- \
  --candidate "$PROJECT/candidate/coordinates.npz" --audit "$PROJECT/candidate/audit.json" \
  --view "$PROJECT/view.json" --font-zh "$FONT_ZH" --font-en "$FONT_EN" \
  --output "$PROJECT/paired-view-check.json"
```

## 已核过的一种局部化学重构

下列入口固定用于细菌翻译终止的 A76 底物。9NO7 提供沉积的 2′-脱氧底物、3′ 酯键和邻近水。9MTP 含 3′ 酰胺类似物，仅提供局部核糖的构象参照。模型补入一个 2′ 氧，保留酯键和周围原子身份。首次运行准备两份 mmCIF 及其来源记录，所需 CCD 通过现有下载器核定并缓存。[9NO7](https://www.rcsb.org/structure/9NO7)、[9MTP](https://www.rcsb.org/structure/9MTP)、[原始研究](https://doi.org/10.1126/science.ads9030)

```bash
python3 "$SKILL/scripts/fetch_structure.py" 9NO7 --output "$PROJECT/sources"
python3 "$SKILL/scripts/fetch_structure.py" 9MTP --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_termination_substrate.py" \
  --sources "$PROJECT/sources" --out "$PROJECT/candidate"
"$PYTHON" "$SKILL/scripts/audit_termination_substrate.py" \
  --sources "$PROJECT/sources" --candidate "$PROJECT/candidate" --output "$PROJECT/independent.json"
```

生成器先尝试固定水氧，必要时仅允许指定水氧在预定范围内调整。输出前独立检查四个手性中心、五个糖环扭转角、跨残基酯键及其有向平面、原 Mg–OP2 连接和周围非键接触。实测来源中的酯氧有符号高度为 +0.351 Å。早期候选虽保住三条距离，高度却多出 0.222 Å。当前静态候选的高度差为 0.008 Å，148,436 个来源重原子中有 148,429 个保持原位。报告保留新增氧和被移动水氧的建模身份。催化水归属、质子转移及肽链释放路径仍需各自的证据和验证。

将这个静态模型或其他核定结构接入多阶段动画时，使用[坐标运动交接](coordinate-motion.md)检查第一帧、后续状态继承和两帧之间的实际距离。
