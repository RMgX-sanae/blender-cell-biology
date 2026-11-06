# 联网核对结合位点、构象和组装

默认联网。这是 agent 接到新相互作用任务后的主动研究步骤，不要求用户先提供结构编号。目的在于找到当前物种、构建体和反应状态下有依据的组合方式。实际坐标、方法、图例和组装说明决定如何建模；关键词命中或数据库中的关联关系不自动给出正确三维接触。

## 从对象身份到可用构型

1. **确认身份。** 用蛋白名、别名、基因名和物种在 [UniProt](https://www.uniprot.org/) 核对 accession、isoform、结构域、成熟链边界及关键修饰。把目标位点的编号体系记下来；缺少足以区分物种或状态的信息时，先研究共用部分，再向用户澄清会改变结论的选择。
2. **主动查复合物与状态。** 在 [RCSB PDB](https://www.rcsb.org/)、[PDBe](https://www.ebi.ac.uk/pdbe/) 搜索两个伙伴及其别名/序列，结合 bound、apo、active、inactive、substrate、inhibitor、assembly 等目标词。名称检索不足时用序列或结构相似搜索扩展同源模板；记录同源、工程突变和截短差异。必要时在 [EMDB](https://www.ebi.ac.uk/emdb/) 找相关密度及配套模型，并读条目对应原始研究。已有同状态缓存可复用；新状态不能只换标签。
3. **核对组合方式。** 阅读 biological assembly、链和 entity、亚基数目、源操作矩阵及原研究。非对称单元可以只有部分装配或包含多份装配，见 [RCSB 组装说明](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies)。同一条目列有多个 assembly 时逐个判断适用性；不把全部晶体邻居都并入复合物。[Complex Portal](https://www.ebi.ac.uk/complexportal/) 可帮助核对参与者、组成和功能，同时区分其人工整理记录与预测记录；同属一个复合物不证明每对成员直接接触。
4. **映射位点。** 用 [PDBe–SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/overview.html) 的序列—结构对应，再与实际 mmCIF 的作者链/label 链、作者残基号、插入码、残基名和原子名核对。`184` 与 `184A` 是不同身份。不要用统一加减一个编号来替代逐位点映射。缺失原子、零占有率、altloc、突变与类似物单独记录；不能给未解析的侧链标出精确实验接触。
5. **确认界面与状态条件。** 在同一已核定坐标系中检查谁接触谁、哪些区域被遮挡、哪些位点可能互斥。核对底物/抑制剂、核苷酸、金属、辅因子、桥接分子、磷酸化、膜侧别和锚定关系。明确每一份非聚合物属于哪个刚体伙伴；对接时也要参与穿插检查。单个距离阈值不是氢键、盐桥、亲和力或催化能力的证明。结合态与未结合态差异明显时保留两个来源状态，而不是拉伸游离单体去贴合界面。
6. **选择实现路线。** 优先复用适用实验复合物及其 assembly。只有同源/部分复合物时，审查序列和结构映射后做明确标注的模板迁移；只有位点/交联等证据时构建受约束候选；仍缺少依据时保留多个探索分支。每种不同界面、化学计量或配体状态单独建来源记录和检查，不把同一界面附近的小幅姿态变化算作新结合模式。

## 使用现有检索与生成器

```bash
# 候选检索，不是自动资格确认；分类 ID 应先核实。
python3 "$SKILL/scripts/search_public_assets.py" 'trypsin inhibitor' --title \
  --taxon 9913 --limit 5 --output "$PROJECT/research/pdb-candidates.json"

# 读过候选条目、原研究和 assembly 后才选坐标。
python3 "$SKILL/scripts/fetch_structure.py" 2PTC --output "$PROJECT/sources"
```

这个小工具只检索 RCSB 实验条目，返回数量有限，排序不保证首项最适用；agent 继续用官方网页/序列搜索及原文解决剩余问题。[RCSB 搜索示例](https://www.rcsb.org/docs/search-and-browse/advanced-search/search-examples) 也包含按组装化学计量筛选的方法。检索失败时保留错误并改用官方入口，不把网络错误写成“没有这种结构”。下载 helper 保存检索时间、URL、结构元数据与 SHA256；确认后的缓存供后续建模和重建使用。

- 合适的实验结合态：沿用 [局部相互作用](protein-interactions.md) 的项目选择表、`prepare_interaction.py` 和保存检查。
- 不同来源单体或同源替代：沿用 [单体候选](monomer-binding.md)，明确被替换的构象/构建体及其推断性质。
- 涉及化学变化或全过程：转到对应机制工作流；一个结合结构不会自动给出反应顺序、所有过渡态或完整循环。
- 比较同一蛋白的不同来源状态：使用 [来源状态对照](source-state-comparison.md) 的项目配置，核对链的实际蛋白身份、关键残基、类似物和辅因子后再生成分开的静态模型。翻译终止实例中，同名 `jj` 链可分别代表 eRF3 或 ABCE1，顶层 mutation 字段也没有直接列出 AAQ；不能依赖这两个字段自动拼接机制。
- 要画具体位点的接触/共价键：用 [位点连接检查](site-link-checks.md) 把审查后的作者原子身份、CCD 元素/键表、实际距离和连接注释核对起来。类似物、突变和冲突注释不能直接进入天然反应动画。

## 留下可复核的决定

在项目 `binding-source-review.json` 或已有来源表中记录实际需要的内容：检索词、日期和入口；目标身份与状态；考虑过的 accession/assembly/链；canonical 位点到作者残基/插入码的对应；参与者与化学计量；配体归属；最终采用/排除原因及原文链接；从来源坐标到模型的变换；实验结构支持的部分与另行构造的部分。成功与失败的来源状态都保留，避免只留下最终能通过的模型。

例如 trypsin–BPTI：公开的 [2PTN](https://www.rcsb.org/structure/2PTN) 是 trypsin 单体且含钙；[4PTI](https://www.rcsb.org/structure/4PTI) 与 [6PTI](https://www.rcsb.org/structure/6PTI) 提供不同晶型的 BPTI；[2PTC](https://www.rcsb.org/structure/2PTC) 提供实际复合物。应先核对后者的 E/I 伙伴与装配，再决定是直接展示结合态，还是研究特定游离来源的候选。6PTI 的磷酸组分不能默认为天然抑制复合物的必要成分。游离构象的搜索失败也不否定已知互作，只说明所用构象、刚体假设、约束和搜索未给出合格结果。

本例实际读取了 2PTC 的 assembly 记录：assembly 1 使用恒等操作、二聚体，由作者及软件定义；assembly 2 使用四个操作、八聚体，由软件定义。本次局部 E/I 对应采用前者，不把后者的全部晶体相关拷贝自动加入。另下载并核对 [2PTC 的 SIFTS 残基映射](https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/2ptc.xml.gz)：

| 2PTC 作者身份 | UniProt 身份 |
|---|---|
| E:184A Gly | P00760:Gly187 |
| E:184 Tyr | P00760:Tyr188 |
| E:195 Ser | P00760:Ser200 |
| E:216 Gly | P00760:Gly217 |
| I:15 Lys | P00974:Lys50 |

这些是本条目的核对实例，不是供别的蛋白套用的编号偏移。脚本的普通关键词检索实测返回了其他相关条目，定向网页检索和条目读取才确定 2PTC；因此不得把返回的前几项或最高分直接当作最终选择。

联网不意味着所有相互作用都有唯一已知构型。若结构/文献冲突，按物种、状态和实验条件分支；若只能访问摘要，注明访问层级。不得把模型预测、数据库预测组合或人工摆放升级为实验结构。用户明确离线或实际网络不可用时，说明本次证据缺口并用已核定缓存继续可做的部分。
