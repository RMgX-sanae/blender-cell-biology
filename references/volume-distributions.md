# 从表面分布推广到体积内放置

这是**数值准备入口**：在明确的三维体素区域内放置完整对象，输出刚体矩阵供 Blender 使用。它沿用受体分布工作流的来源记录、固定 seed、拒绝采样、完整对象排斥和独立检查；曲面法向改为体积包含约束。当前验证的是静态布局，不是细胞器动力学、实测分布或新完成的 Blender 场景。

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

## 已验证的范围

2026-11-06，从独立解包的 skill、空输入目录和上述两份公开文件重建；下载合计 568,499 字节。与已有合格缓存的重建比较，两个 NPZ 的 13 个数组及 manifest 一致。64 个完整包围球通过独立体素盒检查；三张正交数据切片上的球截面仅用于直观核对。六项行为测试包含各向异性、非零坐标原点、局部中心不在原点的模型、越界、跨障碍、拉伸/镜像、重叠和不可放置输入；另外确认错误文件哈希和错误细胞状态不能生成准备输出。

这次没有生成新视频或保存 Blender 工程。下一次使用真实模型时，仍需核对完整评估几何、显示变换与保存工程，不能把上述六顶点数值探针升级为已验收的细胞器场景。
