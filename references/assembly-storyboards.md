# 常用装配过程的分镜模板

这里的 8 个 JSON 是可复制的制作说明，包含参与者、事件顺序、公开资产候选、建模/镜头计划、检查项和限制。它们需要 agent 按项目实现，不能传给 `build_polymer_scene.py`，也不代表已完成 8 段动画。已能直接构建的骨架过程见 [polymer-templates.md](polymer-templates.md)。

Formin/cofilin 的完整上下游需要 [actin-full-cycle.md](actin-full-cycle.md) 的两条机制计划；本页 formin 延伸与 cofilin 周转只是局部模块。尤其不要把 INF2 或 Cdc12 的 profilin 结构改名当作 mDia1 的直接实验证据。

## 使用方式

复制选定 JSON 到项目，先解决 `species_selection`，按 [public-assets-first.md](public-assets-first.md) 核对候选或换用更匹配的资产。为每个 actor 指定实际来源、单位、装配和 evidence；给 `stages` 加上已决定的事件时点、对象 ID 与拓扑前后状态，再开始生成 Blender 工程。阶段时长由讲解需要决定，不是实验反应速率。

`asset_candidates` 是已经查看条目/metadata 的检索起点；本轮没有把这些候选全部下载、清理并渲染。来源条件不合适时应更换，不要为了复用一个 accession 而改变用户的生物学问题。显微图、密度和作者补充视频可补充构象/时序证据。

| 模板 | 资产起点与需要补足的部分 |
|---|---|
| [Formin 延伸](../assets/assembly-storyboards/formin-elongation.json) | [8RU2](https://www.rcsb.org/structure/8RU2) 的端部结合结构；[原研究](https://doi.org/10.1126/science.adn9560)。FH1/profilin、头部步态和核苷酸过程分别选择来源 |
| [Cofilin 周转](../assets/assembly-storyboards/cofilin-turnover.json) | [6VAO](https://www.rcsb.org/structure/6VAO)、[5YU8](https://www.rcsb.org/structure/5YU8)；[结构研究](https://doi.org/10.1038/s41467-018-04290-w)。装饰丝终态不等于整个切断动态 |
| [中间丝装配](../assets/assembly-storyboards/intermediate-filament-assembly.json) | [8RVE](https://www.rcsb.org/structure/8RVE) 及相关四聚体；[原研究](https://doi.org/10.1038/s41594-024-01261-2)。二聚体→四聚体→ULF 的具体早期中间体需另外核对 |
| [Septin 装配](../assets/assembly-storyboards/septin-assembly.json) | [12ZN](https://www.rcsb.org/structure/12ZN) 是线虫候选；如做酵母，可从[四聚体研究](https://www.nature.com/articles/s42003-023-05734-w)追溯。该文明确区分实验四聚体与含预测接口的八聚体；不能把晶体中的非生理接口当正常丝 |
| [FtsZ 踏车](../assets/assembly-storyboards/ftsz-treadmilling.json) | [3VO8](https://www.rcsb.org/structure/3VO8) 为 S. aureus GDP 结构；[构象切换研究](https://doi.org/10.1128/mBio.00254-17)。5MN4 是截短突变构建体；结构不提供完整 Z 环或膜锚定体系 |
| [ESCRT-III 膜重塑](../assets/assembly-storyboards/escrt-iii-remodeling.json) | [6TZ9](https://www.rcsb.org/structure/6TZ9)、[6TZ4](https://www.rcsb.org/structure/6TZ4)；[原研究](https://doi.org/10.1038/s41594-020-0404-x)。CHMP1B/IST1 膜管不能替代所有内体膜颈拓扑与 Vps4 阶段 |
| [Clathrin 被覆](../assets/assembly-storyboards/clathrin-coat-assembly.json) | [6WCJ](https://www.rcsb.org/structure/6WCJ) 的局部 cage 顶点；[原研究](https://doi.org/10.1126/sciadv.aba8397)。完整 triskelion、适配蛋白和颈切断需补资产与分镜 |
| [COPII 被覆](../assets/assembly-storyboards/copii-coat-assembly.json) | [6ZGA](https://www.rcsb.org/structure/6ZGA)、[8BSH](https://www.rcsb.org/structure/8BSH) 是内层候选；[膜上结构研究](https://doi.org/10.1038/s41467-021-22110-6)。另找 Sec13/31 外层，不能把内层复制成笼 |

## 共用建模能力

- **重复装配：** 从相邻结构单元拟合刚体变换，记录间距、旋转与误差；不同构象不能只换颜色。使用共享 mesh，完整组分父子 rig 保存相对关系。
- **交换与踏车：** 固定亚基 ID 与来源；保持观察区内老亚基的位置，明确新增/释放池。整丝平移不代表踏车。
- **切断与退火：** 在亚基接口处分离或连接；保持单体网格完整。检查两侧片段的极性、端部配准和身份守恒。
- **曲面上的网络：** 按弧长、局部法线与真实接口摆放完整单元；检测自交，区分膜形变与蛋白构象变化。
- **状态转换：** 从两个已核对状态建立分镜；插值与人工补形标为教学重构。缺失中间态不意味着允许瞬移、缩放蛋白或穿过整片膜。
- **可复用风格：** 所有模板继承当前项目的参考图分析、材质、灯光、镜头和合成选择，不绑定单一艺术风格。先保住结构可读性，再增加前后景。

这些模板细化已有过程，不应与 60 个过程配方简单相加作为独立成片数量。项目完成时分别报告“配方已提供”“工程已运行”“实际画面已检查”。
