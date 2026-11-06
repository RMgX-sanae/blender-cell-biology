# 从联网证据到可建模的位点连接

当画面要显示具体位点之间的接触或共价连接时，先完成 [联网研究](online-binding-research.md)，再检查实际原子。`audit_site_links.py` 复用来源/组装资格检查和 CCD 下载缓存，读取所选原子身份、距离及 `_struct_conn` 注释。它不从距离自动推断氢键、反应或亲和力，也不取代原研究。

## 用法与输入

```bash
python3 "$SKILL/scripts/fetch_structure.py" 9MTP --output "$PROJECT/sources"
python3 "$SKILL/scripts/audit_site_links.py" \
  --sources "$PROJECT/sources" --config "$SKILL/assets/site-link-examples.json" \
  --case rf1-analogue --output "$PROJECT/site-links.json"
```

首次缺少 CCD 文件时，脚本主动从 RCSB 下载；已有文件使用 SHA256 核定缓存。结构或 CCD 更新导致哈希不同，应重新审查后更新项目配置，不能跳过核对。配置示例还提供 `rf1-deacylated`（下载 9MTT）和 `trypsin-bpti`（下载 2PTC）。依赖沿用普通准备脚本的 gemmi/numpy/scipy 环境，无需 Blender。

新项目复制配置中的一个 `case`，只保留需要核对的内容：

- `source` 沿用 [来源状态配置](source-state-comparison.md) 的结构哈希、assembly、作者链及身份。当前只接受单模型、恒等操作装配中的聚合物位点；其他装配先显式准备，金属/游离配体使用现有配体选择与机制检查。
- `sites` 用语义 ID 指定作者链、作者残基号、component、原子名、元素。无插入码的残基号用整数，有插入码用 `184A` 这样的字符串；label 链和作者链不混用。所选 component 对应的 CCD 哈希放入 `ccd_sha256`。
- `relations` 为两个 site ID 指定经文献和化学身份审查的 `covalent` 或 `noncovalent`、Å 距离区间，以及 `annotation: required` 或 `inspect`。区间是本项目筛查范围，不是通用化学判据；不能为了通过而扩大区间。
- `interpretation` 和 `evidence` 记录这条连接的意义、来源链接及适用范围。`noncovalent` 只表示本配置未将它当作共价键，不能自动标成氢键或盐桥。

同一残基内部的共价注释取 CCD 键表；跨残基使用原始 `_struct_conn` 中匹配的作者身份、组分、原子和保留构象。其他对称操作或替代构象的连接另列，不当成当前坐标的连接。没有连接注释并不单独否定共价键；使用 `inspect` 时仍需人工审查。源注释与预期非共价关系冲突、原子身份不符或距离超范围时返回失败；保留冲突和测量结果，不自动修复坐标。

## 实际案例与边界

| 来源与原子对 | 实测距离 | 采用的解释 |
|---|---:|---|
| 9MTP：1z:Cys7:C — 1x:8AN76:N3′ | 1.333 Å | 源注释支持合成酰胺连接；不能标成天然 C—O 酯键 |
| 9MTP：1x:8AN76:O2′ — 1z:Cys7:C | 3.594 Å | 指定原子对的空间邻近，不是已经形成的键 |
| 2PTC：E:Ser195:OG — I:Lys15:C | 2.675 Å | 抑制剂界面的近接触，不画成共价酶—底物键 |
| 2PTC：I:Cys5:SG — I:Cys55:SG | 1.973 Å | 保留原有二硫键 |

[8AN](https://www.rcsb.org/ligand/8AN) 的 3′位置为氮，[MEQ](https://www.rcsb.org/ligand/MEQ) 是甲基化谷氨酰胺。[2025 年原研究](https://doi.org/10.1126/science.ads9030) 用类似物捕获状态，并提出经 3′→2′酰基迁移后水解的机制；具体路径不能当作已逐步观测的轨迹。9MTT 的 A76 是另一来源中的天然腺苷，其去酰化状态不提供丢失肽链的坐标。N→O 的类似物重构必须发生在明示的模型准备阶段，不能通过反应账本把氮变成氧。

另在 9MTP/9MTT 的 1A:5MU1915:O2 与 1w:Arg286:NH1 之间读到约 1.30 Å 距离及 `covale` 记录。这里将其作为未解决的来源冲突，不把注释当作正常 RNA—蛋白接触的证明，也不据此断言整项结构无效。9MTP 局部检查通过仅覆盖上表及甲基化位点，**不表示整个核糖体的几何或天然终止机制已通过**。

本入口已有两种体系的源文件回归、错误身份/元素/距离/连接注释反例，以及重新保存源文件后的变异检查。它产出 JSON 检查记录；不生成运动、原子补建或 Blender 场景。化学端态成立后，再接入 [反应身份检查](workflow-transfer.md#可复用的反应身份检查) 和对应机制的几何准备，不直接拼接这几个来源。
