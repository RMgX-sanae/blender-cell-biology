# 从已核定的分子构型接入动画

当下一阶段只移动装配中的少量原子，并且这段运动保持原有共价连接时，可以直接从坐标审查通过的模型开始。适用于口袋中的局部重排、受约束的接近，以及多个片段之间的坐标交接。每个新对象重新核对允许变化的键、构象和周围障碍物。[局部补建与坐标核对](coordinate-edits.md)

## 准备与检查

1. 使用上一步的 `coordinates.npz` 和 `audit.json`。新阶段始终在该装配的 Å 坐标系中工作。完整七列原子身份包含装配操作、作者编号及插入码。来源记录要能追到具体复合物、缺失区域、伙伴和位点。
2. 给每段运动写一个无 pickle 的 NPZ。`atom_ids` 是该段移动列的身份，`positions_A` 是形状为“采样数 × 原子数 × 3”的浮点数组，`frames` 为严格递增的展示帧号。相邻阶段共用一个边界帧。省略帧号时，从上一阶段结束处起每个样本占一帧。多个阶段可以移动不同的原子，暂时退出移动清单的原子停留在刚到达的位置。
3. 在 `stages.json` 指定 `candidate_sha256`、`source_sha256`、`assembly_id`、`bonds` 和 `stages`。每个 stage 需要唯一的 `id`、相对路径 `file`、完整身份列表 `moves` 和 `max_step_A`。可用 `sha256` 锁定已查过的运动文件。连续呈现时可另设 `max_entry_velocity_change_A_per_frame`，它比较接缝两侧的离散速度，也包含新开始移动或刚停下的原子。
4. `bonds` 用成对的完整原子身份列出所选运动涉及的共价图。根据已核定的 CCD、聚合物连续性及跨残基连接补齐移动原子的邻域，包括固定的连接伙伴及成键隔一个原子的邻居。移动一个有共价系链的肽段时，保留那条跨组分的连接。化学拓扑改变的时刻须另作事件核对。[反应原子身份与净键变化](workflow-transfer.md#可复用的反应身份检查)
5. 按本项目已有的化学依据填写 `minimum_nonbonded_A` 和 `max_bond_change_A`。可用 `component_floors_A` 为特定组分与对方元素指定更高的接触下限。两端的规则都适用，取较高值。全装配中已核定的水、离子和其他固定原子参与搜索。

```bash
"$PYTHON" "$SKILL/scripts/audit_coordinate_motion.py" \
  --candidate "$PROJECT/candidate/coordinates.npz" --audit "$PROJECT/candidate/audit.json" \
  --plan "$PROJECT/motion/stages.json" --output "$PROJECT/motion/motion-audit.json"
```

检查先按身份恢复每个阶段的起点，再检查每段插值中移动与固定原子、两个同时移动原子的最近距离。两个样本之间即使端点都安全，中途也可能穿过障碍。共价键还检查插值段内最短和最长键长。这些结论适用于保存的笛卡尔线性插值。改用曲线缓动、骨骼形变或改变样本时间后应重新检查。角度、手性和具体反应的化学约束随项目单独核对。

用 PTP1B 的 Asp48 核对过这个问题。沉积的 A、B 两种构象中，CG–OD2 的键长分别为 1.251 Å 和 1.250 Å。按原子名直接从 A 匀速移到 B，这条键在中点会缩到 0.358 Å。检查器会拦下这条路径。替代构象的占有率、两个羧酸氧的对应关系，以及从一种构象到另一种构象的实际运动，需要分别判断。[3I80 实验结构](https://www.rcsb.org/structure/3I80)

## 在现有同视角工程中播放

构建和重开检查沿用静态入口，都加上同一对参数 `--motion-plan "$PROJECT/motion/stages.json"` 和 `--motion-audit "$PROJECT/motion/motion-audit.json"`。左侧保持来源位置，右侧使用经过核定的模型和运动。原子与半键端点保留原生 Blender 关键帧，两侧使用覆盖整段动作的统一正交取景。构建时可用 `--snapshot-frame 25 --render` 只查看一帧，存盘工程仍包含整个时段。重开检查会取全部关键时刻及每段中点，也会检查键连线与固定参照物。实际查看少量画面，判断尺度、取景和相邻颜色能否让人辨认动作。[同视角构建与重开命令](coordinate-edits.md#在-blender-比较两个状态)

## 翻译终止中已经接通的第一段

从 9NO7 的 2′-脱氧实验底物补建 A76 后，源结构的酯键、肽链和周围水都保留了身份。下面的命令令尚未成键的 2′ 氧接近酯羰基。驱动角为 C3′–O3′ 单键扭转，肽链通过邻近主链单键的补偿转动避开核糖体内部的原子。配方核对完整的 30 个已建模肽原子，并在整个装配中筛查运动。附近水氧在这个阶段保持原位，其催化作用尚未确定。[9NO7](https://www.rcsb.org/structure/9NO7)、[原研究 Figure 6 的机制提议](https://doi.org/10.1126/science.ads9030)

```bash
"$PYTHON" "$SKILL/scripts/prepare_termination_approach.py" \
  --candidate "$PROJECT/candidate" --sources "$PROJECT/sources" \
  --mode continuous --frames 49 --out "$PROJECT/motion"
```

`continuous` 从起点平滑进入再减速，`teaching` 为起点与终点留出观察时间。两种呈现经过同一组扭转构型，展示时间可以独立修改。默认路径中 2′ 氧到羰基碳由 3.035 Å 减至 2.753 Å，2′ 氧–羰基碳–羰基氧夹角由 122.9° 变为 110.7°。运动终点仍保留原酯键。这段仅展示反应前的一种几何可能性。形成四面体中间体、后续酰基迁移、水解、肽链释放与回收各需新的状态证据与验证。
