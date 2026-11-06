# 囊泡出芽与融合：把区室交接做成可检查的几何

此入口是**几何原型**：从质膜内陷、切断、运输、对接、开孔到三个固定身份的内容探针进入受体腔。它补上了刚体运输模板没有的膜连通性变化；尚不是完整网格蛋白内吞机制或可交付的连续循环。一个理想化跨膜受体及其配体用于检查膜侧别，不是实验受体结构。

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

## 运行

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

准备入口检查 48 个粗网格状态及 7 个较细网格状态，写出 11 个关键时点、22 个膜边界网格和哈希清单。网格键为 0、0.24、0.34、0.35、0.43、0.62、0.77、0.8、0.88、0.97、1。Blender 入口每次只制作一个静态可编辑快照，避免为几何验证保存整套大场景或渲染完整影片。少量关键时点不证明中间每个状态或连续运动都已通过；改变路径、场、厚度或探针尺寸后须重查受影响时段。

独立解包重建中，11 个文件的 99 个数组完全一致；16 个拓扑/身份/网格错误样例被拒绝。保存工程另外检查实际几何与受体的局部轴方向，避免均匀短杆翻转后外形相同而漏检；10 个保存场景/输入错误样例被拒绝。开孔状态检查了 Cycles 关键帧，透明边界用于调试连通关系；探针在投影中可能遮挡，此图不是分子身份教学视图或最终风格验收。没有渲染完整视频。

## 迁移边界与后续制作

可以迁移的是“区室图 → 完整货物身份 → 同一场的膜/腔网格 → 独立检查 → 少量视觉检查”。COPII、分泌或其他融合路线需重写供受体、膜侧别、参与者和证据，不只是换标签。

本例仍缺被覆及切断机器、去被覆与能量账本、目标路线的融合机器身份、系留柔性、完整回收和连续 Blender 事件同步。连续形变与拓扑切换须另外生成缓存或分阶段网格，并检查切换两侧的实际画面；不能对不同顶点数直接做 shape-key 插值，也不能靠隐藏蛋白掩盖未实现的步骤。需要成片交付时再按 [高效审阅](review-and-retention.md) 的规则补齐连续预览和完整导出。

[AP2–货物膜侧结合](protein-interactions.md)已开始补分子演员，但仍是独立的平膜局部动作，尚未与本例的出芽、切断及融合进度合并。它使用重构体系的核心/短肽，不能直接替代本例理想化受体的全部结构或声明其已完成整条路径。
