# 从表面分布推广到体积内放置

这个入口在明确的三维体素区域内放置完整对象，并把刚体矩阵交给 Blender。它沿用受体分布工作流的来源记录、固定 seed、拒绝采样、完整对象排斥和独立检查；曲面法向改为体积包含约束。数值布局与保存工程检查分开进行，不代表细胞器动力学或实测分布。

## 输入和范围

`scripts/volume_distribution.py` 提供 `sample`、`validate`、`ball_inside` 和 `enclosing_sphere`。输入包括：

- `mask`：ZYX 顺序的三维布尔数组，True 是明确允许的区域，数组外部禁止进入。
- `spacing_xyz`、`origin_xyz`：XYZ 体素间距与第一个体素中心的坐标。坐标、几何、间距和所有 gap 必须使用同一个长度单位。
- `points`：完整对象的实际几何顶点。只给原子中心、骨架端点或省略了几何节点生成部分的网格，不能证明最终外形的范围；应导出完整评估几何或另行建立包含这些外形的保守包围体。
- `count`、`seed`、`boundary_gap`、`pair_gap`：明确的示意数量、随机性与间隙，不由基因表达量或画面观感推算真实丰度。

输出 `transforms` 是可直接作用于原始顶点的 4×4 刚体矩阵；没有镜像或非均匀拉伸。`centers_xyz` 是包围球中心，**不一定等于矩阵的平移项**：原始模型的局部包围中心可能不在原点。将布局放入 Blender 后，细胞、对象与源掩膜使用同一显示变换，并检查保存工程中的实际几何；不能分别重新居中破坏配准。

## 为什么不能只检查中心或顶点

[SciPy 的 EDT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.distance_transform_edt.html) 计算到背景采样点的距离，可指定各轴间距；它不直接返回到体素实体边界的精确距离。本模块在可用区域外补一层 False，以处理裁剪边缘，再作以下保守推导：

`可用距离下界 = 到最近禁止体素中心的距离 − 体素半对角线 − 候选点相对该采样中心的偏移长度`。

只有这个下界大于对象包围半径加边界间隙时才接受。包围球包含顶点及它们构成的直线边、平面三角面和实体；因此也覆盖细长或非凸模型可能跨过的障碍，而不是只检查若干顶点。对象间用两个完整包围球加 `pair_gap` 排斥。这个方法偏保守，可能拒绝实际上能放下的姿态；装填失败时明确改变已声明的数量、模型或策略，不悄悄缩小对象或间隙。

`validate` 不复用 EDT 判据：它核对刚体矩阵，并直接计算包围球到周围禁止体素盒的距离及数组边界距离，检查对象之间的间隙。仍然只证明**给定体素掩膜**内的包含关系，不能验证源分割准确性或未知细胞组分。静态端点通过不意味着连接它们的运动路径通过；多个时间点也不能冒充同时存在的独立对象送入布局排斥检查。

## 可复现的真实分割案例

[Allen hiPSC 数据集](https://open.quiltdata.com/b/allencell/packages/aics/hipsc_single_cell_image_dataset) 的 [README](https://allencell.s3.amazonaws.com/aics/hipsc_single_cell_image_dataset/README.md) 定义了单细胞分割、通道名称与体素间距元数据。配置 [volume-cell-230745.json](../assets/volume-cell-230745.json) 固定了一个 CellId 230745、元数据状态 M0、TOMM20 标记案例及两个下载文件的哈希；这不是代表总体的随机细胞样本。提供方使用条件见配置中的链接，原始数据不随 skill 分发。

```bash
"$PYTHON" "$SKILL/scripts/prepare_volume_case.py" \
  --inputs "$PROJECT/inputs" --output "$PROJECT/prepared" --fetch
"$PYTHON" "$SKILL/scripts/test_volume_distribution.py"
```

准备需要 NumPy、SciPy、tifffile。已有匹配缓存可省略 `--fetch`。该选定文件约 437 kB，元数据只取 HTTP Range 指定的前 128 KiB；服务器不接受 Range 或文件指纹不同就停止，不下载完整大型元数据或直接改哈希绕过检查。

这个适配器只接收明确 ZCYX、0/255 的分割及**各向同性**元数据间距。当前间距为 0.108333 µm；这不是光学分辨率。别的各向异性来源先核实轴顺序，再调用通用 API，不能从此例猜测元数据轴序。坐标原点是局部裁剪的首体素中心，没有声称恢复全视野配准。

可放置区域明确为 `cell & ~dna & ~structure`，不采用为展示制作的 roof 通道。原始三个掩膜完整保留；本例结构与 DNA 掩膜有 3,142 个重叠体素，记录为分割冲突，不解释成线粒体进入细胞核，也不修改源掩膜掩盖冲突。排除这些已知掩膜后仍有未覆盖的细胞组分，剩余区域不等于生理上全部空闲胞质。

示例的六顶点形体和 64 个位置仅用于检查范围，不表示真实细胞器形状或数量。`source-masks.npz` 保存原掩膜、派生可用掩膜及尺度；`layout.npz` 保存探针顶点、矩阵和数值下界；`manifest.json` 保存来源、配置、冲突、哈希及检查范围。这些文件均留在项目外部目录。

实际模型的项目代码复用同一接口，替换完整几何及有依据的目标区域：

```python
from volume_distribution import sample, validate
layout = sample(available_mask_zyx, spacing_xyz, origin_xyz, mesh_vertices_xyz,
                count, seed=seed, boundary_gap=boundary_gap, pair_gap=pair_gap)
check = validate(available_mask_zyx, spacing_xyz, origin_xyz, mesh_vertices_xyz,
                 layout['transforms'], boundary_gap=boundary_gap, pair_gap=pair_gap)
```

旋转、数量和排斥形成的是示意布局；拒绝采样也会改变实际空间分布，不能称为均匀生理分布、热平衡或扩散模拟。优先保留源细胞器的真实位置；只有问题需要生成假设分布时才使用放置模块。

## 数值准备的验证

2026-11-06，从独立解包的 skill、空输入目录和上述两份公开文件重建；下载合计 568,499 字节。与已有合格缓存的重建比较，两个 NPZ 的 13 个数组及 manifest 一致。64 个完整包围球通过独立体素盒检查；三张正交数据切片上的球截面仅用于直观核对。六项行为测试包含各向异性、非零坐标原点、局部中心不在原点的模型、越界、跨障碍、拉伸/镜像、重叠和不可放置输入；另外确认错误文件哈希和错误细胞状态不能生成准备输出。

上述六顶点探针阶段没有生成视频或保存 Blender 工程；完整模型接入另按下面的流程检查，不能用数值探针替代实际几何。

## 把完整模型接入 Blender

`blender_volume.py` 负责导出、放置和保存后检查；`prepare_volume_layout.py` 用导出的真实外形替代六顶点探针。当前一次布局放置同一模型的多个静态副本，材料、标签、镜头和细胞背景沿用项目自己的显示流程。

1. 准备并保存一个完整的单网格模型。可以先通过细分、加厚等修改器建模；导出使用 [Blender 的评估几何](https://docs.blender.org/api/5.3/bpy.types.Depsgraph.html)，并冻结为三角面快照。集合、子对象、粒子、未实现的几何节点和材质位移需要先落实为完整网格。视口/渲染的修改器开关、细分层级和螺旋步数必须一致；任意其他渲染专用效果仍由项目先烘焙、核对。
2. 显式说明原模型每个单位对应多少 µm；例如 1 unit = 1 nm 时传入 `0.001`。对象世界矩阵只能含正的统一缩放及刚体变换，不能靠非均匀缩放拉长蛋白。模型自身局部中心不必在原点。材质不随几何快照搬运，之后按项目风格另行设置。
3. 导出器另启一个无自动执行脚本的 Blender 进程，重开源文件、在声明帧比较完整几何与世界矩阵，再记录源文件哈希。`bpy.data.is_dirty` 无法可靠检测 Python 修改，不能用它代替实际比较。源文件与当前模型不同就先保存预期模型，不把未保存形状归给旧文件。
4. 先计算布局，再放置到与源区域共用的坐标框架，重开保存工程检查。

```bash
"$BLENDER" -b "$PROJECT/body-source.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/blender_volume.py" -- export \
  --object "Complete body" --micrometers-per-unit 1 --output "$PROJECT/body-export"

"$PYTHON" "$SKILL/scripts/prepare_volume_layout.py" \
  --region "$PROJECT/prepared" --body "$PROJECT/body-export" \
  --output "$PROJECT/layout" --count 40 --seed 230745 \
  --boundary-gap-um 0.05 --pair-gap-um 0.1

"$BLENDER" -b "$PROJECT/context.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/blender_volume.py" -- place \
  --layout "$PROJECT/layout" --anchor "Source micrometer frame" \
  --output "$PROJECT/placed.blend"

"$BLENDER" -b "$PROJECT/placed.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/blender_volume.py" -- audit \
  --layout "$PROJECT/layout" --output "$PROJECT/saved-volume-audit.json"
```

`--region` 是前面 `prepare_volume_case.py` 的合格输出，不能将不封闭的显示 OBJ 直接充当三维掩膜。`context.blend` 中的 `Source micrometer frame` 是项目创建的空物体，其世界矩阵将源区域的 **XYZ µm 坐标**统一映射到显示坐标。细胞、排除区域与新放置物必须沿用它，不能各自居中、缩放。若源背景已经在同一 µm 坐标下，空物体可用单位矩阵；需要旋转或整体缩放时，由一个共同矩阵处理。来源显示网格与掩膜是否对应，仍需该数据适配器自己的检查。

`place` 保留全部冻结顶点和面，以共享网格和刚体矩阵创建副本，不自动加模型外的装饰、路径或动画。项目 Python 可直接调用 `place(layout_dir, anchor)`，给返回对象赋材质后保存；同一场景不重复追加第二套布局。修改器已烘焙进导出几何，保存后的新位移/变形、非均匀拉伸、脱离共同坐标系、隐藏/增删副本或缓存变化会使检查失败。

保存检查读取**实际评估网格与世界矩阵**，用实际浮点坐标计算完整包围球，再进行体素盒包含和副本间隔检查；不只相信准备阶段的矩阵。姿态比较允许 Blender 浮点存储误差，实际包含仍用保存后的几何范围。原始掩膜、导出体、布局和工程之间保留文件指纹，所有生成资产留在项目目录。它不验证背景网格的生物学来源、相机可读性、材质位移以外的渲染效果或运动路径，仍需对应检查及必要的一张实图。

### 完整模型接入的验证

2026-11-06，在 Blender 5.2.3 LTS 中以一个人为建模的完整膜样物体检查接口：原网格 482 个顶点，经细分、加厚后导出 3,972 个顶点和 7,936 个三角面；包括内外两面，不切开模型。模型形状、大小、厚度和 40 个副本的数量都是测试参数，不是该细胞测得的囊泡。

从独立解包的 skill 重开源工程、重新导出和放置后，14 个几何/布局/掩膜数组与开发结果相同。重开工程验证全部 40 个实际外形；三份背景掩膜等值面的顶点、面和共同变换另查，全部主体在镜头范围内。初次建景的 1280×800、Cycles 24 samples 静态图已查看，背景的阶梯来自源分割分辨率，不用平滑改变检查掩膜；没有渲染动画，也未把重建后未重复渲染的工程称为像素复验。

四个导出反例覆盖未保存的几何、视口/渲染细分不一致、未实现的几何节点和对象级材质位移；七个保存场景反例覆盖拉伸、隐藏、顶点改动、追加加厚、仅渲染生效的修改器、脱离父对象和共同框架移动；三项文件反例覆盖区域/模型缓存变化及只更新局部哈希却保留旧上游记录。`export`、`place`、`audit` 命令均实际运行。此次只完成静态完整几何接入，没有新增真实细胞器丰度、跨区室运动或囊泡运输机制。
