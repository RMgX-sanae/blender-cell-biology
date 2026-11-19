# 囊泡出芽与融合：把区室交接做成可检查的几何

要拍清囊泡从哪里来、内容物又能进到哪里，可在 Blender 中选两段局部动画。第一段从膜内陷到膜颈收紧、载体切断并离开。第二段从载体靠近受体区室到融合孔张开，货物开始进入受体腔。三份液相货物、一个教学示意跨膜受体和配体各自保持身份。另有静态快照可检查整条教学路径上的特定时点。

## 先限定生物学情境

[vesicle-uptake-scope.json](../assets/vesicle-uptake-scope.json) 保存事件、原始研究、证据适用范围和未实现阶段。本例是哺乳动物内吞到早期内体的教学综合，不能称作同一细胞中观测到的一条完整轨迹。

- [AP2/网格蛋白重构研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC4333214/)支持膜与货物参与的激活及被覆出芽；不说明细胞中其他因子都可省略。选择一个成功摄取事件来讲解，不把“先捕获货物”升级成所有小窝起始的必需顺序。
- [Dynamin2 招募](https://pmc.ncbi.nlm.nih.gov/articles/PMC4230619/)和 [Auxilin/GAK–Hsc70 去被覆](https://pmc.ncbi.nlm.nih.gov/articles/PMC7054993/)分别限定切断与释放后的机器步骤；本原型只实现前后几何，没有装配这些蛋白。
- [2026 年 Rab5 时序研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC13388671/)在所测人细胞中把 Rab5 招募放在去被覆之后，涉及 GEF 与脂质变化。不能给刚形成的小窝直接加上 EEA1 阳性标签；该结果不提供本例所需的全部融合机器身份。
- [EEA1 实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC5142606/)支持不对称系留及柔性变化。实现时保留蛋白轮廓长度，用弯曲改变端到端距离。[分离内体实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC2592678/)区分对接与融合，不能把相互靠近直接算成融合。
- [同型早期内体的 SNARE 研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC1413832/)不能直接确定 SUM159 新摄取载体与内体之间的全部 SNARE 身份；[NSF/α-SNAP 研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691687/)也不等于完整回收及 ATP 计量。这里保持未解决状态。旧 P10 是酵母液泡体系，不能把 HOPS 移入本例冒充 EEA1 路线。

事件表的 `after` 是本例的讲解顺序，允许生物学招募和脂质变化重叠。`uncoat/rab/recycle` 等时间区间仍是待接入的安排，不代表对应蛋白、核苷酸或回收过程已经实现。

## 复用方法

1. 复用 `event_timing.py` 定义讲解进度，但另写区室关系：切断前胞外连通载体腔；切断后分离；对接时两腔仍分离；开孔后载体与受体腔连通，始终不通胞质。
2. 用同一标量场生成两侧膜边界和水相/膜/胞质掩膜。当前用平面、球和胶囊体的并集，是无量纲教学几何。双边界不代表已解析两个分子级单层；无脂质面积守恒、真实动力学或半融合中间态。
3. 同时检查**水相和膜的连通分量**。水相分离时，外侧膜仍可能错误相切；只检查腔内容物会漏掉这个问题。分别用 6 邻域与 26 邻域检查水相和膜，并核对具名区室种子；胞质另查 6 邻域连通。
4. 固定内容物身份，并要求它们处于载体所属的连通腔。仅判断“还在水相”会放过货物提前跳进另一封闭内体。球形探针用 `volume_distribution.ball_inside` 检查整体体积，另用解析场的 Lipschitz 界核对；受体锚点、胞质尾及配体侧别分别检查。
5. 独立检查导出网格的分量、Euler 特征、边面关联、相邻绕序、整体朝向及场等值偏差。有限视野中的质膜有裁切边界；只允许其开边位于声明的视野边界。这不是任意网格的完整流形/自交证明。
6. 将选中的完整网格交给 Blender，重开保存工程比较实际评估顶点、三角面、矩阵和探针几何。已有材质、灯光与合成入口可沿用；不把验证中的透明膜外观当作参考风格已经验收。

## 制作可播放的局部动画

先确定镜头要解释哪一对区室在何时连通。切断时镜头近乎平行于供体膜，使新出现的间隙在投影里也可见。开孔时让两个腔的连线横过画面。若前后两份货物互相遮挡，可绕这条连线抬高机位，保留孔道的投影长度，让观众能跟随每份货物通过。半透明膜用于观察内部，同时保留货物、配体和受体的完整轮廓。

将 `SKILL`、`PYTHON` 和 `BLENDER` 设为本机绝对路径，将 `PROJECT` 设为 skill 外的项目目录。准备程序使用 NumPy、SciPy 和 scikit-image。几何节点动画已在 Blender 5.3 上运行。命令用 `scission` 制作切断段，改成 `fusion` 即制作开孔段。

```bash
EVENT=scission
"$PYTHON" "$SKILL/scripts/prepare_vesicle_topology.py" --motion "$EVENT" \
  --frames 97 --fps 24 --output "$PROJECT/$EVENT"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/blender_vesicle_topology.py" -- animate \
  --input "$PROJECT/$EVENT" --output "$PROJECT/$EVENT.blend"
"$BLENDER" --disable-autoexec -b "$PROJECT/$EVENT.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/blender_vesicle_topology.py" -- audit \
  --input "$PROJECT/$EVENT" --output "$PROJECT/$EVENT-check.json"
```

打开保存的工程，播放时间轴，或在三个事件标记间跳转。选中 `Membrane` 可在几何节点修改器中编辑载体位置、尺寸、膜颈半径、孔径及连通开关。节点从同一距离场生成膜的两面，允许网格顶点数随连接关系变化。货物、配体和受体是独立的关键帧对象。修改位置或时序时，在准备脚本里同步修改相应几何与角色，再生成新的准备目录和工程。打开已保存的工程即可播放，无需帧缓存或启动 Python 回调。

默认的 97 帧按 24 fps 播放。切断在第 65 帧，开孔在第 41 帧。切换前后各看一帧实际画面，再用覆盖变化的一小段原速预览判断可读性。重开检查从实际 Blender 网格核对边界连通分量和两腔之间的通路，同时检查货物的位置、完整半径和受体轴向。这两个局部区间使用教学时间。开关代表两种已定义的连通状态，脂质重排和分子中间态尚未解析。

## 提取整条路径的静态快照

外部 Python 需 NumPy、SciPy、scikit-image；Blender 端只加载准备好的数据和自身 NumPy。输出都放到 skill 外的新目录。`time` 是 0–1 的编排进度，不是实验秒数。

```bash
SKILL='/actual/path/to/blender-cell-biology'
PROJECT='/absolute/project/vesicle-topology'
PYTHON='/actual/path/to/python'
BLENDER='/actual/path/to/blender'

"$PYTHON" "$SKILL/scripts/prepare_vesicle_topology.py" --output "$PROJECT/keys"
"$PYTHON" "$SKILL/scripts/test_vesicle_topology.py" --input "$PROJECT/keys"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/blender_vesicle_topology.py" -- build \
  --input "$PROJECT/keys" --time 0.8 --output "$PROJECT/pore.blend"
"$BLENDER" --disable-autoexec -b "$PROJECT/pore.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/blender_vesicle_topology.py" -- audit \
  --input "$PROJECT/keys" --output "$PROJECT/saved-audit.json"
"$BLENDER" --disable-autoexec -b "$PROJECT/pore.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/review" --frames 1 --width 1280 --samples 24
```

静态准备入口检查 48 个粗网格状态及 7 个较细网格状态，写出 11 个关键时点、22 个膜边界网格和哈希清单。网格键为 0、0.24、0.34、0.35、0.43、0.62、0.77、0.8、0.88、0.97、1。`build` 每次制作一个选定时点的可编辑快照。改变路径、场、厚度或探针尺寸后，选受影响的时段重查。

独立解包重建中，11 个文件的 99 个数组完全一致；16 个拓扑/身份/网格错误样例被拒绝。保存工程另外检查实际几何与受体的局部轴方向，避免均匀短杆翻转后外形相同而漏检；10 个保存场景/输入错误样例被拒绝。开孔状态检查了 Cycles 关键帧，透明边界用于调试连通关系；探针在投影中可能遮挡，此图不是分子身份教学视图或最终风格验收。没有渲染完整视频。

## 迁移边界与后续制作

可以迁移的是“区室图 → 完整货物身份 → 同一场的膜/腔网格 → 独立检查 → 少量视觉检查”。COPII、分泌或其他融合路线需重写供受体、膜侧别、参与者和证据，不只是换标签。

要将局部镜头接成完整内吞路线，还须交代被覆及切断机器、去被覆及能量消耗、目标路线的融合参与者、系留柔性和回收。迁移到其他膜事件时，先选两侧水相、应该保持身份的货物、跨膜朝向以及可分辨通道的机位。若来源支持另一种形状，可替换几何节点中的距离场，沿用持久货物、连通开关和切换两侧的检查。只有需要交付成片时才完整导出。[高效审阅](review-and-retention.md)

[AP2–货物膜侧结合](protein-interactions.md)已开始补分子演员，但仍是独立的平膜局部动作，尚未与本例的出芽、切断及融合进度合并。它使用重构体系的核心/短肽，不能直接替代本例理想化受体的全部结构或声明其已完成整条路径。
