# 风格、渲染与合成分开选择

用户明确偏好：按项目选择教学或写实表现；渲染风格可以选择，可用不同合成器完成；也可以根据图片确定风格并模仿实施。不要固定唯一的材质、色板、背景或 compositor。

## 三个独立决定

| 决定 | 可选方案 | 什么时候用 |
|---|---|---|
| 表现风格 | 清晰教学、写实分子/细胞环境、插画/描边、显微外观示意、参考图模仿 | 由观众、图示任务和用户偏好决定 |
| 渲染引擎 | Eevee、Cycles；项目已有的可用引擎 | Eevee 快速预览/教学；Cycles 用于需要的光传输、透射、细致材质。不是“Cycles 自动科学准确” |
| 合成器/流程 | Blender Compositor、Fusion、Natron、After Effects 或项目已有流程 | 看可用性、节点/图层需求、后期协作；不是要求全部安装 |

仅在当前环境实际存在且可用时使用外部软件。不为一个简单静态图自动安装整套后期工具。FFmpeg 编码是交付步骤，不等于完整的节点式合成器。

## 风格配方

### 清晰教学

背景明暗按项目选择；默认加入低对比、有景深虚化的细胞微结构环境，见 [cellular-environment.md](cellular-environment.md)。使用低反射材质、适量柔和阴影、有限语义色；正交或弱透视镜头。连续膜外壳可用教学开口或局部透明；分子特写通过完整单体的前后层次与整件隐藏露出内部，见 [section-views.md](section-views.md)。保持 membrane/cytosol/lumen 可区分。必要标签独立，缩到最终投影尺寸仍能读。避免在主体上堆放发光、颗粒、景深或高饱和背景。

### 写实细胞环境

先用数据决定结构和拥挤，再加材质。可采用曲面不规则性、实例化蛋白/聚合物、柔和散射和局部照明；噪声纹理的幅度应低于要表达的真实结构尺度。没有数据时将环境密度标作示意，不把随机小球当真实细胞质。按渲染成本使用实例、LOD 和局部镜头。

PBR “真实”指光照模型，不保证分子表面、数量或动力学正确。背景默认可用景深；关键反应区保持清楚，运动模糊须在确认机制可读后再使用。严重透明遮挡优先改剖面，而非把所有膜变成玻璃。

### 显微外观

可以用通道、有限分辨率、PSF 卷积、采样和噪声生成教学成像模拟；说明哪些参数有校准、哪些仅为视觉效果。单纯 emission+glow 不是荧光显微镜的物理模型；不能从这种渲染测表达量/共定位。EM 外观要区分真实密度/断层重建和灰度美术表现。

### 图像风格复现

按 [style-from-reference.md](style-from-reference.md) 从实看图片拆出形体、材质、光、镜头、色调与后期。保持生物场景不变，比较相同构图的小样。不要只套 LUT 就声称完成模仿。

## Blender 合成器起点

扩展库另外提供 `COMP_STYLE_clean/ink/filmic/fluorescence/grain/dic/phase/soft`，与 live Geometry Nodes 和材质一起使用。构建、独立选择 engine/compositor、迁移和适用边界见 [style-library.md](style-library.md)。DIC/相差节点处理的是示意高度代理，不能替代经校准的成像模型。

精细生成器的 `--compositor` 有 `clean`、`reference`、`monochrome`，分别保存为 `COMP_clean`、`COMP_reference`、`COMP_monochrome`。先将 Render Layers 的前景 Alpha Over 到所选背景，再做参考对比或灰度处理，最后输出。对比偏移在背景合成前作用于 premultiplied RGB 会制造负值/黑边，应避免。材质、灯光与镜头由独立 JSON 控制，详见 [production-workflow.md](production-workflow.md)。

旧白模 `build_scene.py` 使用另一套 `none/clean/contrast/glow` 参数；不要把它们传给精细生成器。glow 是美术效果，不是荧光成像模型。新建分支应保留现有工程的节点组。

Blender 5.x 可使用 `scene.compositing_node_group` 和 `NodeGroupOutput`；较早版本使用 `scene.use_nodes/scene.node_tree/Composite`。Alpha Over 等节点的输入名称也可能变化。运行时查询 RNA/输入名称，不照抄旧教程后声称兼容。本次主验证环境为 5.2.3；其他版本需小样验证。

## 外部合成器交接

1. 提前确认 color management、fps、frame range、分辨率和 alpha 约定。输出场景线性多层 OpenEXR：Combined/RGBA、需要的 diffuse/specular/emission、Z、Cryptomatte 或对象 mask；只开所需 pass。
2. 主体、膜、关键 cargo、环境和标注可分 View Layer；所有层共享同一相机和时间。相互遮挡/投影需要 holdout 或 shadow catcher，不能各自渲染后直接相加。
3. Fusion/Natron 等节点式流程：读取 EXR→线性色彩管理→按 mask 处理→合并→显示变换→标注→输出。AE 项目则明确 32-bit/线性工作流程与导入 alpha。具体控件按已安装版本核对。
4. 避免把 AgX/显示 LUT 烘焙两次。预览 PNG/MP4 通常是显示变换后的结果，不能等同线性 EXR。
5. 标明 straight/premultiplied alpha；边缘黑线/白边优先检查 alpha/颜色空间，不直接扩张 mask 改变膜轮廓。
6. 保留原始 EXR 与最终 MP4/PNG，输出可编辑 `.blend` 加可用外部工程。若外部软件不可用，提供图层与准确的交接说明，不能称外部合成已运行。

`render_saved.py --exr` 提供实际导出路径。Blender 5.x 先将 `image_settings.media_type` 设置为 `MULTI_LAYER_IMAGE`，再设置 `OPEN_EXR_MULTILAYER`；它与旧版直接设 file_format 不同。为交接启用 interleaved 存储，仍须在目标软件核对通道和色彩空间。接口依据：[Blender 5.0 API 变更](https://developer.blender.org/docs/release_notes/5.0/python_api/)、[ImageFormatSettings](https://docs.blender.org/api/5.0/bpy.types.ImageFormatSettings.html)。

## 性能与测试

先在 512–960 px、少样本下检查；Eevee/Cycles 的最终路径都实际渲染一个代表帧。Cycles GPU 先探测设备并局部配置，不修改全局偏好。记录机器、版本、渲染耗时和设置；不能根据单帧时间保证所有帧。遇内存压力先降低实例/局部细节与分辨率，不牺牲必要生物结构。

图片展示在聊天中，视频可用本地绝对路径内嵌。询问满意度应指出可调项，如“更哑光/更通透”“更简洁/更多环境细节”“参考图色调是否接近”，而非重复泛问“可以吗”。
