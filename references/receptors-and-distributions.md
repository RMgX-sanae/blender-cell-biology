# 受体外形、朝向与空间分布

先查生物身份，再生成位置。来源入口 S05–S09；3SN6 只是一个 β2AR–Gs 实例，不是任意 GPCR 的通用几何。表达、定位、表面密度和纳米簇是四个不同问题。

## 可运行的曲面分布实例

`prepare_receptor_distribution.py` 在 Blender 外使用 NumPy/SciPy；`build_receptor_distribution_scene.py` 在 Blender 内读取准备结果。默认实例把 3SN6 的 R/A/B/G 蛋白部分放在示意内体膜上，输出固定 48 套完整复合物的均匀、成簇或极性分布。它补充了原有点位采样工具，并不自动推断任意受体的膜方向或生物密度。

```bash
python3 "$SKILL/scripts/fetch_structure.py" 3SN6 --output "$PROJECT/sources"
python3 "$SKILL/scripts/prepare_molecular_asset.py" "$PROJECT/sources/3SN6.cif" \
  --preset beta2 --output "$PROJECT/assets/beta2"
python3 "$SKILL/scripts/prepare_receptor_distribution.py" \
  --asset-dir "$PROJECT/assets/beta2" --config "$SKILL/assets/receptor-distribution.json" \
  --distribution clustered --output "$PROJECT/distribution-clustered"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_receptor_distribution_scene.py" -- \
  --distribution-dir "$PROJECT/distribution-clustered" --output "$PROJECT/scene-clustered" \
  --font-zh /path/to/Songti-Regular.ttf --font-en /path/to/Times-New-Roman.ttf \
  --environment-dir "$PROJECT/environment" --style dark --recipe neutral --width 1280 --render
"$BLENDER" -b "$PROJECT/scene-clustered/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_receptor_distribution_scene.py" -- --project "$PROJECT/scene-clustered"
```

每次输出到新目录。`--distribution` 还支持 `uniform`、`polarized`。`--style light` 改变材质/灯光/底色；`--recipe teaching-outline` 等选择已有合成节点链，文字独立保护。`subject-focus` 要求 `--engine cycles`。有色图例也必须实际复查，双色合成不自动保留链身份配色。默认是静态空间分布；没有生成或验证受体扩散、激活及内吞动画。

背景目录由 [cellular-environment.md](cellular-environment.md) 的 `prepare_cellular_environment.py` 从 1UBQ、9SPI 和 8F8R 生成。该实例按原始 nm 大小放置完整核糖体、泛素和肌动蛋白形状，再分别虚化远、中、近层；数量和位置是场景构造。`--environment-dir` 可省略以输出纯底；默认制作优先提供有身份的细胞质背景。

### 修改配置时保留的约束

- `assets/receptor-distribution.json` 保存物种/工程构建边界、来源、膜朝向、数量、种子、三种强度场、排除区、镜头及图例。先复制到项目再修改。3SN6 是工程改造、混合物种、核苷酸空缺的 Gs 复合物；省略 T4 融合部分、稳定纳米抗体及非聚合物配体不等于天然完整人源复合物或受体失活。此入口保留所选解析蛋白部分，不补缺失环、糖链或脂锚。
- 默认 `local_positive_z=inward`：准备资产的 +Z 胞外方向在内体上指向腔内，Gs 在膜外的胞质。用于质膜时，重新定义所处区室、法线和镜头；不能只把标题改为“细胞”。[3SN6](https://www.rcsb.org/structure/3SN6) 提供结构，[Irannejad 等 2013](https://www.nature.com/articles/nature12000) 提供内体 β2AR/Gs 信号背景；本配置的形状、数量和分布参数不来自该研究。
- `surface.kind=mesh` 可读取配置文件旁的 NPZ `vertices`（N×3 世界坐标 nm）和整数 `faces`（M×3）。必须先应用非均匀缩放并核对三角形绕序；生成器按实际三角形面积采样，并使完整复合物局部 Z 对齐指定方向的面法线。`ellipsoid` 是明确标记的教学形态，不是显微重建。
- `fields` 是固定总数下的**相对强度**：均匀提案、带背景项的高斯热点，或平面一侧富集。热点距离使用世界空间欧氏距离，并非曲面测地距离。硬核拒绝会改变最终强度；不能称为严格泊松过程或直接拟合表达量。`excluded_spheres` 约束锚点；完整复合物还要通过全重原子间距检查。
- 计数装不下时明确失败，不自动减少数量。0.20 nm 是严重几何穿插筛查阈值，不是范德华半径、结合能或生物学验证。默认三个例子实际最短距离均更大；任何新输入重新运行数值检查并看图。
- 两层完整脂质由随机面积采样及最小间距构造，分别朝向膜芯。连续膜芯网格补足展示间隙，保守省略靠近蛋白的整片小三角形；不切蛋白或脂质。没有脂质组成、力场平衡或膜动力学含义。
- 换受体时重新限定结构、链、来源构象、膜中心和 +Z 方向，并为所有链提供 `chain_colors`、对应图例与区室标签。普通 `generic` 导入不会识别跨膜方向；应先完成来源有记录的刚体定向，不能直接套用 β2AR 的选择规则。

验证分三层：`test_receptor_distribution.py` 检查面积统计、固定种子、朝向、不可装填和完整分子排斥；准备记录保存全复合物重原子最短距离及原文件哈希；重开 `.blend` 的审计核对未变形的源网格/拓扑、实例矩阵、脂质朝向、字体和构图。最后查看实际浅/深色图及选定合成结果，并从解压包和新下载坐标独立重建。此工作流不证明细胞形态、受体丰度或聚集机制。

Blender 5.2 的 `FLOAT_VECTOR` 属性批量读取使用 `np.float32` 缓冲区；实测 `float64` 缓冲会误读字节，造成错误的朝向报警。检查器记录单位方向向量的最大误差，容许旋转转换的浮点误差，但实际翻转必须失败。

## 组件选择

| 类别 | 外形/拓扑 | 不应默认的内容 |
|---|---|---|
| 单跨膜受体（如特定 RTK） | 胞外域、跨膜螺旋、胞内域；方向由序列/结构核对 | 所有 RTK 都由配体诱导两个孤立单体相遇；有些存在预组装 |
| GPCR | 通常 7 跨膜束；N/C 端与胞内耦联面明确 | 任意七根棒就是目标受体；所有 GPCR 激活都具有同样幅度/寡聚状态 |
| 配体门控/电压门控通道 | 装配亚基、贯穿膜的通路和门 | 孔道等于持续开放；电流不等于每帧一个离子 |
| 转运蛋白/泵 | 朝向交替或循环构象 | 双侧永久打通的孔；ATP/离子计量必须查对应家族 |
| 整合素 | α/β 亚基、胞外结合端与胞内骨架连接 | 所有局部粘附都均匀分布；活化构象仅靠拉长整个蛋白 |
| 钙黏蛋白 | 跨细胞接触面的对向结合、胞内连接 | 分布到所有自由膜面，或让其连接到 ECM |
| TCR/BCR/Fc/细胞因子受体 | 各自多亚基装配、辅助分子和配体方向 | 免疫突触中的簇是静态同心环且适用于所有阶段 |
| GPI 锚定蛋白/外周膜蛋白 | 锚定在特定叶/表面，不一定跨膜 | 给每个膜相关蛋白补一根跨膜螺旋 |
| CD44 等糖基化表面蛋白 | 胞外域、跨膜区、胞内尾及明确糖基位点 | RNA 量代表表面 CD44 数量；受体结合代表完整内吞/降解 |

## 分布模型选择

- **测量坐标**：保留单位、定位精度、检测效率、去重/聚类方法；不凭渲染把定位点变成确切分子数。
- **均匀独立分布**：作为零假设；曲面按面积采样，不能均匀采样 UV/纬度后声称面积均匀。
- **排斥分布**：给足迹/最小距离，通过拒绝采样或 Poisson disk 施加；所需密度装不下时报告不可能，不暗降数量。
- **成簇分布**：明确 cluster centers、簇宽、各簇拷贝数/强度、背景分量、随机种子；没有数据则标 hypothesis。
- **极性/区域分布**：给 apical/basolateral、接触面或权重场 mask；密度为单位真实面积的强度，先归一化。
- **梯度/热点**：保存 λ(x) 的定义、来源与边界；不要让相机视角改变生物密度。
- **动态分布**：区分扩散、结合/解离、围栏、主动转运和内吞移除。关键帧插值不会自动保总数或守恒。

`distributions.py` 可复用：`sphere_points` 是球面面积均匀；`patch_points` 提供有边界的均匀/成簇/极性假设和排斥；`triangle_points` 对世界坐标三角形面积加权。`build_scene.py` 的受体模板是平面单跨膜代理，可换成核对过的结构资产。

## 曲面采样与法线

1. 获取评估后的网格，应用世界变换后计算三角形面积；非均匀缩放下不能直接使用原网格面积。
2. 按面积选三角形 a,b,c。令 s=√U, t=V，点 p=(1−s)a+s(1−t)b+stc，U,V 独立均匀。
3. 用正确外向的法线 n；顶点法线需经 inverse-transpose 变换再归一化，或直接用变换后三角面的叉积。校验法线是否指向胞外/指定腔面。
4. 先把受体局部 +Z 定为胞外方向，再将 +Z 旋到 n。随机旋转只绕 n；任意 Euler 随机会把胞外域翻入胞质。
5. 位置为膜中面 p，亲疏水跨膜段与膜厚对应。胞外域大于跨膜段，排斥足迹应覆盖实际外形；不能只检查锚点。
6. 非凸曲面需要局部邻域/测地距离或 BVH 碰撞检查。两张相邻但不同膜的 Euclidean 距离不能误当同膜近邻。

物理量 `N ≈ ∫_A λ(x)dA`；观测 N 与期望 N 需区分。已知计数时固定 N 抽样，不再额外做 Poisson 波动。未知密度的模板输出显式说明 `N` 是演示参数。球变椭球会改变面积元，不能只缩放球上点后声称椭球均匀。

## 扩散和聚集动画

平面二维自由扩散可用 Δx,Δy~Normal(0,2DΔt)，此处第二参数为方差，传给 NumPy `normal` 的标准差应为 √(2DΔt)，因此二维 MSD=4DΔt。如果 D 是 µm²/s 而场景单位是 nm，先完整换算。边界选择反射/周期/吸收并说明。曲面随机步需要局部切平面、投影与偏差检验；简单投影通常不是严格曲面布朗运动，不能作为定量模拟。

结合/解离概率可用 `1-exp(-k*Δt)`，但双分子结合率涉及浓度/体积或面积，不能直接把三维二阶速率当二维每帧概率。教学动画可用手工时序，标为 schematic。需要实验动力学时使用经验证的模拟器输出，再在 Blender 播放。

验证：数量与排除区域、分区实际面积、法线方向、最小间距/相交、seed 重现；可视化抽样俯视/侧视。若要验证聚集性，可用面积正确的 null model、pair correlation/Ripley K 并说明边界修正；单张“看起来成簇”的图片不是统计证据。
