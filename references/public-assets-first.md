# 先找公开资产，再决定建模范围

用户明确要求优先参考、复用网上公开资产来提高效率。每个新对象先查看项目已有资产和来源记录；已有经核对的适用缓存时直接复用，出现缺口再检索公开数据或成品模型。已有可用结构时，优先导入、修整、实例化；只为缺失部分或特殊表达补建几何。不要把“从零写球、棒、随机团块”当作默认起点。

## 一次可执行的资产选择

1. 用对象名、物种、装配/结合状态和景别检索。准备少量候选即可，不为凑数量下载完整数据库。分清想找的是原子坐标、密度、分割表面、可编辑模型、轨迹，还是仅供观察的图片。
2. 实际打开候选条目、缩略图或模型。记录来源、作者、accession、物种/构建体、状态、尺度和许可。搜索标题相近不等于合适；例如泛搜“septin tetramer”可能返回其他蛋白的四聚体，应当排除。
3. 选择“可直接复用 / 需轻度修改 / 仅作参考 / 不适用”。同时考虑科学适配、可编辑性和处理成本。全套原始显微数据未必比一个已分割 ROI 高效；大型原子复合体也未必比按链简化的表面适合远景。
4. 保存原始数据及哈希，再生成派生资产。把筛选链、去除配体、修补缺失区、平滑、降面、变换矩阵和单位转换写进 manifest。相同来源与处理参数直接复用缓存，使用共享 mesh/Collection Instance；颜色和镜头变化不应重新提取原子表面。
   链标识可能区分大小写，而 macOS 的文件系统可能不区分。遇到 6VQ6 这类同时含 `A` 与 `a` 的结构时，采用十六进制等可逆、安全的文件名映射，并在 manifest 保留原链名；检查每条链的原子数、网格和路径唯一性，避免资产静默覆盖。
5. 在本项目中检查尺度、朝向、拓扑及关键身份后再装配。可看、可下载、可修改、可再分发是不同状态；不确定许可时可以保留出处作参考，转用许可明确的资产制作。

找不到合适资产时，简短记录检索词、已查入口及缺口，然后制作明确标注的教学模型。无需为重复确认已经授权的建模而停工；涉及新账号、付费下载或对外发布时按实际授权范围处理。

## 按对象选择入口

| 对象/用途 | 优先入口 | 复用方式与检查 |
|---|---|---|
| 蛋白、核酸、配体、骨架、膜复合体 | [RCSB PDB](https://www.rcsb.org/)、[PDBe](https://www.ebi.ac.uk/pdbe/) | 选 biological assembly、链、核苷酸/药物状态；缺失区、融合蛋白和辅助因子另记。原子模型可转共享表面，不必重新雕刻 |
| 大复合体密度、细胞内原位结构 | [EMDB](https://www.ebi.ac.uk/emdb/)、[EMPIAR](https://www.ebi.ac.uk/empiar/) | 先看预览、分辨率、体素、阈值和文件大小；已有拟合 PDB 时优先结合使用。密度阈值表面不是完整原子表面 |
| 细胞形态、细胞器分布、分裂时序 | [Allen Cell 数据](https://www.allencell.org/data-downloading.html)、[BioImage Archive](https://www.ebi.ac.uk/bioimage-archive/) | 优先单细胞/ROI、已有分割和元数据，避免先下载整套大数据。区分原图、分割和预测模型；逐数据集查许可 |
| 亚细胞定位和组织表达 | [Human Protein Atlas](https://www.proteinatlas.org/humanproteome/subcellular) | 可作定位参考；IF 信号不能直接决定受体拷贝数、分子表面密度或纳米簇 |
| 无实验坐标的蛋白形状 | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | 标为预测，读取 pLDDT/PAE；不要把低置信区、未知装配和配体状态当作确定结构 |
| 细胞环境配方、分子 packing | [cellPACK 官方仓库](https://github.com/mesoscope/cellpack) | 查 recipe 的 ingredients、来源结构、尺度、数量和分布假设；代码许可不能替代每个外部资产许可。本轮未执行 cellPACK |
| 已做好的 Blender/glTF/OBJ 模型 | 作者实验室/论文数据仓库；[Sketchfab 下载入口](https://sketchfab.com/developers/download-api/guidelines) | 优先有科学来源、组件分离且许可适用的文件；检查单项作者、许可、下载限制、单位和面数。只有视觉外观时不可用作结构证据 |
| 材质、光和镜头参考 | 已授权的用户参考、原作者公开图片/视频 | 实际看图后走 [style-from-reference.md](style-from-reference.md)；不因模仿风格就重分发原图或照搬其结构结论 |

PDB/EMDB/EMPIAR 数据的开放依据见 [官方公开数据声明](https://www.ebi.ac.uk/pdbe/about/public-data-access-statement)。论文插图、网页图片、商标与第三方模型分别查自己的条款。其他入口的证据边界见 [sources-and-evidence.md](sources-and-evidence.md)。

## 先检索，再下载

`search_public_assets.py` 只检索 RCSB 实验结构，不声称覆盖所有资产库。默认全文检索；加 `--title` 要求标题包含给出的每个词，降低噪声。分类 ID 确知时可用 `--taxon` 检索该分类谱系（包括后代）；混合物种构建体仍需读每个 entity 的 metadata。返回候选、分辨率、组分、物种、assembly IDs 和引文，不自动决定资产适用性。

```bash
python3 "$SKILL/scripts/search_public_assets.py" 'free barbed end' --title \
  --limit 5 --output /absolute/project/research/actin-candidates.json

# 查过条目、构建体和 assembly 后下载；沿用已有来源与缓存。
python3 "$SKILL/scripts/fetch_structure.py" 8F8R \
  --output /absolute/project/sources
```

普通全文词匹配仍可能返回不相关条目；这不是“已找到正确模型”。接口短暂失败时先看错误、有限重试或用官方网页检索，不把错误响应缓存成坐标。输出已存在时脚本拒绝覆盖，供 agent 检查后复用。

## 项目中留下什么

可复制 [asset-selection.json](../assets/asset-selection.json)，保留实际需要的字段。记录候选为何入选或排除，区分 `coordinates`、`density`、`segmentation`、`mesh`、`reference_image`；引用来源并保留加工过程。脚本不执行下载包中的代码，不要求为普通模型导入安装未核实的插件。

建议缓存键包含 accession/version、assembly、chain selection、坐标哈希、表面参数与单位。镜头近时保留精细资产，远时用共享简化实例；不要把近景全部原子复制成几十万个独立 Blender 对象。向其他机器交付时使用相对引用/打包依赖，并验证缺少本机缓存时能按来源重建。

## 本轮可复用资产

| 条目 | 实际用途 | 易误用处 |
|---|---|---|
| [8F8R](https://www.rcsb.org/structure/8F8R) | 肌动蛋白重复形状、末端方向参考；内部 D→E 配准形成延伸变换 | ADP F-actin；入场单体形状代理不等于真实 G-actin 或 ATP 水解过程 |
| [8F8Q](https://www.rcsb.org/structure/8F8Q) | 完整 CapZ 两亚基及结合终态 | 原研究的兔 actin / 人 CapZ 混合构建体；不默认为同物种天然整体 |
| [7AQK](https://www.rcsb.org/structure/7AQK) | 母丝、七亚基 Arp2/3 接头、子丝的相对位置 | 9 Å 原位模型；不能用粗分辨率骨架解释精确侧链化学 |
| [3J6F](https://www.rcsb.org/structure/3J6F) | 三个原丝局部窗口，完整 αβ 二聚体 | 14 原丝参考的 GDP 晶格；不是现成的通用 13 原丝管或 GTP 帽状态 |

以上坐标与来源 metadata 均已实际获取，并用于模板。其他扩展过程的候选入口在 [assembly-storyboards.md](assembly-storyboards.md)，其状态单独标明，不能混称都已制成动画。
