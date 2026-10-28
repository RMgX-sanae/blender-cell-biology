# 从公开数据到 Blender 资产

## PDB/mmCIF

优先保存原始 mmCIF、元数据和 biological assembly 文件；legacy PDB 适合小结构的轻量起点，有链名/原子数等格式局限。`fetch_structure.py` 支持传统四字符 ID、可选 `--assembly N`，保留坐标哈希、下载地址、时间、原始 API 元数据。它不会替 agent 选择正确 assembly、清理构建体或验证配体。

操作顺序：

1. 在 RCSB 条目核对物种、突变、截短/融合蛋白、配体、实验方法与分辨率。记下 biological assembly 编号及其中辅助组分。PDB 测量约束下的坐标也是结构模型。
2. 选链/组分并记录删改。缺失柔性区或糖链不要凭空补成“实验结构”；若为连续外形而补模型，单独标 inferred。
3. 默认可用 `prepare_molecular_asset.py`：在项目 Python 环境用 Gemmi 读取 PDB/mmCIF，按 entity 类型选聚合物，转为各链原子数组/高斯显示包络，再由 Blender 加载。具体筛选、构建体和依赖见 [production-workflow.md](production-workflow.md)。也可使用实际安装的 Molecular Nodes 或 ChimeraX/PyMOL；检查本机 API，不编造插件命令。
4. 设定尺度。原子坐标通常为 Å；记录导入后的统一缩放。新旧资产统一坐标轴，蛋白本身不要非均匀缩放以“塞进膜”。
5. 近景选择 atoms/bonds/cartoon/surface；远景降细节但保留来源。空间填充球不等于分子溶剂排除表面，Cα 折线不等于所有原子。
6. 膜蛋白查 UniProt 拓扑和 OPM/PDBTM/相关实验；将疏水区与膜对齐，确认胞外/胞内域；不能用包围盒长轴代替膜法线。
7. 多份结构实例化共享网格，用对象自定义属性记录来源/chain/role；不要为每个受体复制几十万独立原子 mesh。

旧 `structure` 白模读取附带的 1BNA ATOM 行；该轻量解析器不能当作通用清理器。精细资产准备按 Gemmi entity 处理，保留被标作 polymer 的修饰残基；默认仍不处理多构象占据加权、配体/水和缺失区重建，需按任务选择。来源：[Gemmi molecular models](https://gemmi.readthedocs.io/en/latest/mol.html)。

3SN6 assembly 1 是 β2AR–Gs 的实验构建体复合物，含工程/辅助组分。若只展示受体，应清楚说明选取哪些链及为什么；不要把全部链统一摆到膜内。

## 显微图、分割与表面重建

- 2D 图片通常不能唯一确定 3D；没有断层/Z-stack 时，只能做基于图的近似模型并标明隐含假设。
- OME-TIFF/OME-Zarr：读取 channel、time、axes 顺序和 voxel spacing。`ZYX` 存储顺序与 Blender `XYZ` 坐标需显式转换。Z 常与 XY 不等距。
- 分割标签是测量派生结果：记录人工/算法、版本、阈值和修订；预测标签不是原始像素。阈值改变会改变膜接触/孔/连通性。
- marching cubes 应使用真实 voxel spacing，处理世界坐标/原点；平滑/Decimate 前保存原网格。修复非流形时先定位成像/分割原因，不能盲目封掉真实孔洞。
- 薄膜若低于体素分辨率，图像并未提供两叶结构。远景渲染可加可视厚度，但记下该厚度是表现参数。
- 检查导入前后的包围盒、体积/面积（仅在有意义时）、朝向和至少两个正交视图。与原始切片叠加核对；单张漂亮透视图不够。

数据处理例：在已有 Python 环境中用 `tifffile` 读取 OME 元数据，`skimage.measure.marching_cubes(mask, spacing=(dz,dy,dx))` 生成网格，再显式把 ZYX 顶点重排为 XYZ。对标签逐对象处理并保留 ID。先核对依赖版本，不在 Blender 内临时 pip 一大套分析依赖。

## EM 密度 / cryo-ET

先区分 density map、particle coordinates、segmentation 与 fitted atomic model。核对 map 的 axis order、origin、voxel spacing、threshold 和 handedness。阈值改变的表面不是同一个“唯一真实分子边界”。地图与拟合模型要同坐标，不应分别随意居中后再声称精确吻合。

Molecular Nodes 支持哪些 `.map/.star` 路径以当前文档/安装为准。大型体数据先做 ROI 和低分辨率预览；不要在未评估内存时加载完整电镜数据集。公开库的下载大小与许可分别核对。

## 轨迹 / 模拟

检查 topology 与 trajectory 的原子索引一致、时间间隔、周期边界和单位。先 unwrap/align（并保留处理记录），再导入 Blender；相邻帧的跳跃可能是 PBC，不是分子瞬移。渲染插值不能当作真实 MD 新采样点。未知物理模型的关键帧动画只能称 schematic。
