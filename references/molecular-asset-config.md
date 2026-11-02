# 用项目配置更换分子资产

更换膜蛋白时，先替换经过核对的结构选择和坐标系，再沿用包络生成、曲面分布、材质、合成与保存工程检查。`prepare_molecular_asset.py --config` 接受项目自己的 JSON，无需为每个蛋白添加一个硬编码 preset。它处理**已选聚合物的静态显示资产**；配体、水、离子、反应端态和运动仍需相应的机制准备入口。

## 配置约定

以 [aqp1-asset.json](../assets/aqp1-asset.json) 为可运行示例，复制到项目后再修改。使用 `--config` 时不能同时指定非 generic 的 `--preset`。

| 字段 | 要求 |
|---|---|
| `id` | 项目选择的名字，字母/数字/下划线/连字符；会写入 manifest 的 `preset`，供分布配置匹配 |
| `source_id`、`source_sha256` | RCSB 条目 ID 和**实际输入坐标文件**的 SHA256；assembly 文件名与条目 ID 分开记录 |
| `chains` | 有序的 `{"chain":"A"}` 列表；名称以 Gemmi 读出的作者链为准。assembly 可生成 `A-2` 等新名称，不能假设只有 A/B/C/D |
| `residue_ranges` | 每条链可选的 `[[first,last], ...]`，作者残基整数编号、两端包含。不设则保留该链全部已解析聚合物原子；重叠范围不重复原子 |
| `selection_note` | 为什么选这些链/残基，哪些辅助组分被省略，以及构建体限制 |
| `orientation.origin_angstrom` | 源坐标中的三维原点，单位 Å；用于膜蛋白时说明如何确定膜中面 |
| `orientation.basis_columns` | 3×3 矩阵，列依次是源坐标里的局部 X、Y、Z 单位轴。必须正交且行列式为 +1；不允许拉伸、剪切或镜像 |
| `orientation.note`、`evidence_urls` | 朝向的推断方法、膜侧别及来源；数学上正确的旋转矩阵并不能证明生物朝向正确 |
| `limitations` | 项目特有限制，追加到输出 manifest；配置全文及其哈希也会保留 |

局部坐标计算为 `local_nm = (source_angstrom - origin_angstrom) @ basis_columns * 0.1`。对所有链使用同一个变换，保留装配关系。不要逐链居中后把它们误称为原复合物。

该入口读取模型 0，按 entity 选 polymer，保留非氢、正占据及空白/A altloc 原子；残基整数范围包含该编号的所有插入码变体，不提供插入码单独筛选。需要其他模型、占据策略、非聚合物配体或更细身份映射时，先按 [数据导入](data-import.md) 核对并使用专门准备步骤，不能把默认过滤结果直接当反应物。

## AQP1 的迁移示例

[RCSB 1J4N](https://www.rcsb.org/structure/1J4N) 的 assembly 1 为牛 AQP1 四聚体；不对称单元只有一条聚合物链，另一个软件装配是八聚体。这里明确选择 assembly 1。下载后的四条链为 A、A-2、A-3、A-4，保留各 1,852 个已解析聚合物重原子。

配置里的膜坐标系根据 [UniProt P47865](https://rest.uniprot.org/uniprotkb/P47865.json) 六个完整跨膜区段的首尾各四个 Cα 推定；先核对作者编号与序列一致，再根据各端所在区室取平均。+Z 为胞外方向。具体区段、源哈希和方法保存在 JSON，不能把这组原点和矩阵搬到另一份结构文件。

```bash
# PYTHON、SKILL、PROJECT、BLENDER 为本机路径；使用新的输出目录。
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 1J4N --assembly 1 --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_molecular_asset.py" \
  "$PROJECT/sources/1J4N-assembly1.cif" --config "$SKILL/assets/aqp1-asset.json" \
  --output "$PROJECT/assets/aqp1"
"$PYTHON" "$SKILL/scripts/prepare_receptor_distribution.py" \
  --asset-dir "$PROJECT/assets/aqp1" --config "$SKILL/assets/aqp1-distribution.json" \
  --output "$PROJECT/aqp1-distribution"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_receptor_distribution_scene.py" -- \
  --distribution-dir "$PROJECT/aqp1-distribution" --output "$PROJECT/aqp1-scene" \
  --font-zh /path/to/Songti-Regular.ttf --font-en /path/to/Times-New-Roman.ttf \
  --style dark --recipe neutral --width 1280 --samples 32 --render
"$BLENDER" -b "$PROJECT/aqp1-scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_receptor_distribution_scene.py" -- --project "$PROJECT/aqp1-scene"
```

本配置只定义 `uniform`；需要成簇或极性时，按 [分布说明](receptors-and-distributions.md) 在项目配置中定义对应 `fields`，再检查。24 套四聚体、椭球尺寸、密度和颜色均是显示选择，不是测得的细胞形态或丰度。四种颜色只区分对称拷贝。

膜侧和背景要一起审查。本例向外显示胞外面，因此上述命令没有使用现有的胞质背景；不能因为默认背景方便，就把核糖体摆到胞外。若需要环境，改用与该视角/区室一致的来源资产。按项目可切换浅色材质与 `teaching-outline` 等合成器，几何和结构身份不随滤镜改变。

此包络省略水和去污剂，狭窄孔道还可能被高斯表面封住；不能据此展示水分子路径、孔径、质子排斥或门控。这些新机制需要重新建立通路和相应验证。

## 更换为下一个目标时

1. 核对物种、构建体、状态与 assembly；列出实际链名/已解析残基，与序列编号对齐。
2. 写项目选择及膜坐标依据。来源文件哈希变化时重新核对，不能只更新哈希绕过错误。
3. 准备资产，比较原子数、身份、坐标和链间距离。检查空选择、错链、镜像/拉伸矩阵确实失败。
4. 在分布配置同步 `source_id`、`source_preset=id`、链、颜色、图例、+Z 侧别及环境。蛋白保持刚体，膜形态与蛋白轴分别处理。
5. 用原有构建/检查器重开保存工程，查看至少侧向和整体画面；新的膜曲率或朝向需要重新检查嵌入位置。静态模型通过后再决定是否实现结合、门控或运输。

这样可迁移的是明确输入的资产处理和显示链路；选对生物装配、膜中面及功能状态仍是 agent 必须完成的来源判断。
