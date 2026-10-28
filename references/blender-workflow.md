# Blender 执行与动画工程

精细模型的实际命令和依赖在 [production-workflow.md](production-workflow.md)。本文件解释工程与动画约定；旧 `build_scene.py` 仅用于机制白模。

## 环境探测

找用户现有 Blender，而非直接下载安装。可依次检查 PATH、常用安装目录和已有项目配置；macOS 也可能在用户 Applications 或 Steam 中。记录 `--version`。实际执行 `--factory-startup` 仅重置当前后台进程，不写用户偏好。

```bash
"$BLENDER" --version
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_organelle_scene.py" -- \
  --output /absolute/project/mitochondrion --width 960 --render
```

脚本参数放在 `--` 后；`--python-exit-code 1` 让 Python 错误成为失败退出码。打开已有工程时不使用构建器清空场景：先复制工程，按对象/集合做局部修改。

## 建模工程约定

按结构职责命名：`Membranes`, `Organelles`, `Receptors`, `Cargo`, `Guides`, `Annotations` 等，简单示例无需强制每组都存在。保存实体 ID 和 compartment，动画时同一货物不要变成另一个对象却沿用错误身份。导入的原始资产放单独 collection；派生网格和实例可分开。

只在确有需要时采用 Geometry Nodes。重复脂质/受体/核糖体适合实例；连续核膜外壳的教学开口可用简单 mesh/Boolean，分子单体则保持完整，按 [section-views.md](section-views.md) 组织前后层次与遮挡。节点有稳定随机种子与命名输入；大场景不提前 Realize Instances。必要时分级：overview proxy → organelle detail → molecular structure，镜头切换时说明尺度。

程序化脚本优先使用 `bpy.data`；依赖上下文的 `bpy.ops` 明确设置 active object/mode。后台运行必须实际测试，不能把 UI 中成功的步骤直接当作 headless 可行。物体存在、保存成功、实际成像是三个不同检查。

## 动画实现选择

先选展示形式，再选技术实现。`stepwise`（分步教学）可拆分步骤、设置观察停顿、配合讲解；`continuous-loop`（连续循环）以同一镜头持续呈现过程，按用户要求选择是否配音，不自动添加阶段卡或解释性脚注。此选择不替代单独的视觉风格和合成器选择。

连续循环需保存一个完整周期，并检查首帧与下一周期首帧的实际姿态、材质和光照相同；最后输出帧与首帧之间的位移应相当于普通相邻帧，不要重复终点造成停顿。粒子回收在画外或不可见区完成，不能在画内突然换位置。反应物消耗、产物释放和转子周期要一致。增长、分裂、酸化程度等不可逆状态不能为接缝而倒放或跳回：可选择稳态局部过程，另行说明未模拟累积量。

可移植的连续循环入口见 [vatpase-loop.md](vatpase-loop.md)：从公开 6VQ6 与 1UBQ 重新准备几何，构建整体泵及 c 环视图，检查双视图 H⁺ 运动并渲染。生成脚本随 skill 分发，成片、帧和 `.blend` 不打包。原制作环境保留了完成案例，但它不是运行可移植入口的依赖。

| 需求 | 起点 | 特别检查 |
|---|---|---|
| 刚体位移/转动 | 位置/四元数 keyframes、约束或曲线 | 不穿膜；父子空间变换正确；方向从轨道切线定义 |
| 同拓扑形变 | Shape keys、Armature/Lattice | 每个 key 顶点数与顺序相同；不要用顶点错序冒充形变 |
| DNA/RNA/肽增长 | 按方向排列曲线点、bevel_factor 或显式逐段增长 | 端部连接和 5′/3′、N/C 方向；整个曲线均匀缩放不是聚合 |
| 大量重复组件 | Geometry Nodes instances、稳定 IDs | seed、每帧数量和身份；不要每帧重随机产生闪烁 |
| 出芽/分裂/融合 | 分阶段网格、形变+切换、缓存 | 连通性改变是拓扑事件；检查切换前后及中间帧 |
| 显微时序 | 保留测量帧间隔的网格/体积序列 | 不平滑出来源没有的事件；明确时间重采样 |
| 分子构象 | 已对齐的同原子序列结构/真实轨迹 | 对应原子一致；线性 morph 是示意，可能产生碰撞/不合理键长 |
| 物理/生化动力学 | 专门模拟器输出→缓存→Blender | 模型参数、边界、单位和验证先于美术动画 |

示意快慢可选 LINEAR/BEZIER，但分裂/绑定/释放事件不能被插值提前触发。对象可见性开关用 CONSTANT。Blender 新版本 Action slots/layers 与旧版 fcurves API 不同；修改 keyframes 后检查实际评估帧。

定向接近时分别设计位置路径、朝向和时间曲线。根据入口与障碍物选择直线或 Bézier 路径；需要对准的刚体可在接近阶段小幅转向，到位前平缓减速。五次 smoothstep 可使起止速度与加速度为零；携带的货物/氨基酸应用同一刚体变换。逐帧查看中间路径、接触前减速和停靠后的短暂停留，检查附近所有关键对象，不能只比较起点和终点。路径与快慢属于教学编排时明确标注；实测轨迹和扩散模型按其原有时间与运动约束处理。

用户要求放慢时，分别调整接近、对接停顿、连接变化、释放和产物观察时间。重新评估保存工程中的实际曲线，不能只改文档中的秒数或降低输出 fps。动作时段、字幕和高亮要同步；对接前尚未连接，释放后离场对象保持完整，最好真正离开镜头后再隐藏。只有确实不变的停顿帧才可复用图像，运动阶段保留原定帧率。

“增加下一个循环”需要推进实体身份和状态：核酸按下一模板碱基选择底物，先前产物成为下一步的起点，必要时展示转位和活性位点重新空出。不要把上一轮原样循环播放，却在字幕里称为继续聚合。新状态若不在源结构中，明确标注为重构；检查新的末端连接、配对、原子/货物去向和可见性。当前 C→U 双循环示例见 production-workflow。

标注文字保持面向相机。分子转向后，先把真实锚点变到世界坐标，再添加屏幕方向的偏移并更新引线；不要把用于前景显示的深度偏移也随分子旋转，否则文字会倾斜或漂离原位置。转位、相机移动后重新检查完整主体和标签的画内边界。用颜色表示反应状态时，统一并记录规则，例如“金色表示当前加入的核苷酸，整合后采用 RNA 的橙色”；不能让颜色变化暗示实体身份改变。

## 输出与编码

先生成关键帧预览；满意后按完整帧范围输出 PNG/EXR 序列，便于断点恢复、重合成。构建器默认避免覆盖非空输出目录；修订用新目录，已有场景改动按项目约定进行。

```bash
# 从已保存的工程渲染；脚本创建新的 frames 目录并配置本进程 GPU
"$BLENDER" -b /absolute/project/scene.blend --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory /absolute/project/frames --frames all

python3 "$SKILL/scripts/encode_video.py" /absolute/project --fps 24 \
  --ffmpeg /actual/path/to/ffmpeg
```

`-a`/`-f` 放在渲染配置参数之后。生成器的 `--animate` 会建 frames 目录。不要直接在尚未验证的长任务中只输出一个 MP4；中断后不易恢复。编码器检查从 `0001` 开始的连续 PNG 序列，不覆盖已有视频。

不要混淆播放 fps 与生物时间。字幕可从事件表生成 SRT，在后期独立合成；中文字体从本机可用字体选择并确认许可。基础模板用英文短标签以避免不确定字体依赖；项目交付按用户语言更换并实际检查。

## 当前脚本的范围

默认精细制作链在 [production-workflow.md](production-workflow.md)；它包括结构资产准备、实际 β2AR–Gs/Pol II/tubulin-kinesin、连续内膜线粒体、连通核膜—ER、连续胞质分裂、图像风格控制和 EXR 导出。过程配方不等同于成片，查看 [validation.md](validation.md) 区分已执行与待实施。

以下是旧白模的历史范围，不能据此评价新的成片质量：

- 标准库+bpy 构建 8 种起点；不安装插件、不调用付费服务。
- 随附 `1BNA.pdb` 是实际下载的公开结构，完整来源在 assets/1bna-source.json。
- `fetch_structure.py` 的 mmCIF/assembly 下载、结构元数据保留单独验证；下载并不表示该结构已导入 Blender。
- `--style-profile` 由 agent 在看图后填写，生成器应用色板/背景/roughness/合成选项；复杂镜头与光照模仿仍需项目脚本编辑。
- 完整测试记录见 [validation.md](validation.md)。仅为教学原型；其他细胞类型和 60 个过程根据配方实施后逐项验证。
