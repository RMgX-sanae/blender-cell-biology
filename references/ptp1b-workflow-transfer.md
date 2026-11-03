# 从净键变化扩展到共价中间体：PTP1B

本例验证同一磷酰基依次交给底物、酶中间体和 Pi 时，能否保持身份、几何和双视图一致。范围为**预先结合的底物 → 磷酰基转移至 Cys215 → 肽离开 → Gln262 转向及水进入 → 中间体水解 → Pi 留在结合位置**。没有 Pi 释放、底物进入、WPD 开合、质子接力、ER 锚定或完整催化循环。时间、路径和中间插值是教学构造，不是测得的动力学。

## 沿用的方法与必须重建的部分

| 沿用入口 | PTP1B 必须另做的判断 |
|---|---|
| `atom_table`、刚体叠合和共价键提取 | 选择一份持续存在的酶/肽构建体；不同突变体和末端不能整段互换 |
| `reaction_ledger.validate_sequence` | 两步、交接和独立净反应；等待的水与离开的肽始终保留身份 |
| `chemical_samples`、`rotate` | 磷酰基几何、Cys 局部参考和 Gln262 内部扭转，不拉伸整个蛋白 |
| `interaction_motion.matrices`、`phase_progress` | 肽释放及水进入的可通过路径与阶段次序，逐姿态检查 |
| [共用分子序列显示](molecular-sequence-display.md) | 项目镜头、局部残基选择、标记和有主次的高对比色板；科学适配器不包含 Blender 画面布局 |

`ptp_geometry.py` 和 `prepare_ptp1b.py` 是这一组选定来源的适配器，不能靠替换 PDB ID 自动推广为其他磷酸酶。更换机制时保留接口和检查方法，重新建立目标的反应物、来源资格及坐标规则。

## 来源如何使用

- [3I7Z](https://www.rcsb.org/structure/3I7Z) 是第一步过渡态**类似物**结构。持续模型选其解析的作者链 A 2–298（2,426 个重原子）和 B 988–993 的 DADEYL（51 个重原子，含 OXT）。保留原坐标，不把钒改名成磷。底物磷酸基根据类似物方向另建四面体几何，属于解释性重构。
- [1A5Y](https://www.rcsb.org/structure/1A5Y) 是捕获磷酸半胱氨酸中间体的 Q262A 构建体。仅把 CSP215 的磷酰基放入持续模型的局部 Cys 坐标系；不导入整个突变酶。其 mmCIF 还列出序列差异，不能把它当作完整野生型端态。
- [3I80](https://www.rcsb.org/structure/3I80) 是第二步过渡态类似物。核心叠合后，用其 Gln262 参考限定两个侧链扭转；水氧位置由类似物的轴向 O4 相对 Gln 酰胺的关系构造。这个水氧**不是被识别出的实验溶剂水**；最终 Pi 也不是将钨酸盐改名所得。
- 两步机制和 Q-loop/WPD 的区别依据 [Brandão et al., 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2871455/)；中间体来源另见 [Pannifer et al., 1998](https://doi.org/10.1074/jbc.273.17.10454)。已读 2010 摘要、相关正文、Figure 2/5 图例和状态比较；没有声称完成全部方法、补充材料或全部 wwPDB validation reports 审阅。
- [1PTU](https://www.rcsb.org/structure/1PTU) 的 C215S 和酰胺化肽末端不直接替换上述构建体。筛查时不同来源的 WPD 偏差可达数 Å，不能用全酶线性 morph 抹平。
- 背景使用持续 PTP1B 及人 [1UBQ](https://www.rcsb.org/structure/1UBQ)、[1YCR](https://www.rcsb.org/structure/1YCR) 作者链 A 的形状；它们是有来源的拥挤环境示意，不表达 PTP1B 复合物组成或测得的细胞丰度。

原始 mmCIF 和下载元数据有 SHA-256 约束。结构文件版本变化时重新资格审查，不为绕过报错直接更新哈希。主结构与背景源文件均在项目外部缓存，不随 skill 分发。

## 具体检查与失败边界

全程保留 2,482 个重原子：酶 2,426 + 肽 51 + 磷酰基 4 + 水氧 1。两阶段及独立总反应核对同一库存，只声明反应中心键子图，不包含氢、质子化、电荷和键级。

右侧只强调 Cys215、Gln262、底物 Tyr992、磷酰基和水氧；颜色对应下方图例。整酶 Cα 折线曾在放大图里形成难以辨认的网络，用户明确要求修改，因此当前默认去掉这些导引线。所选残基的化学连接仍来自原始键表；全酶、六残基底物肽及运动坐标没有删改。局部显示的残基仍属于左侧完整分子，不表示游离氨基酸，也不表示这几处残基构成全部催化机制。

准备检查覆盖全部输出姿态的分子间接触、Q262 移动侧链自接触、聚合物键长与相邻键角保持、固定酶原子和磷酰基几何。2.38 Å 是本例重原子中心筛查阈值，明确排除反应键及 Q262 的 1–2/1–3 共价邻居；不是能量或生物物理真实性判据。624 帧设置对应 625 个含末端的准备样本。

局部拟合给出的 Q262 扭转约为 χ1 +21.16°、χ2 −108.00°，相对参考的四个移动原子 RMSD 约 0.327 Å；其余 2,422 个酶原子固定。水位点相对全局叠合的类似物氧约移动 0.739 Å。这些是本次重构的量，不应称为文献报告的运动幅度。

失败尝试也约束范围：直接移入早期底物磷酸几何会产生短接触；选择同终点的另一侧链旋转方向可能在中间穿插；尝试的 Pi 刚体释放弧线未通过接触检查。后者不证明所有刚体路径都不可能，但足以拒绝把这段动画标为完整循环。未来释放需要自己的来源、构象和路径验证，不用隐藏碰撞或放宽阈值来补齐章节。

## 重建与改变节奏

Python 需要 NumPy、SciPy、Gemmi、scikit-image；使用 Blender 5.2、FFmpeg 和本机宋体/Times New Roman。变量均替换为实际绝对路径，PROJECT 为 skill 外的新项目目录。已有哈希相符的公开缓存可以复用。

```bash
for pdb in 3I7Z 1A5Y 3I80 1UBQ 1YCR; do
  "$PYTHON" "$SKILL/scripts/fetch_structure.py" "$pdb" --output "$PROJECT/sources"
done
"$PYTHON" "$SKILL/scripts/prepare_ptp1b.py" \
  --sources "$PROJECT/sources" --config "$SKILL/assets/ptp1b-two-step.json" \
  --output "$PROJECT/assets" --frames 624 --fps 24
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_molecular_sequence_scene.py" -- \
  --asset-dir "$PROJECT/assets" --display "$SKILL/assets/ptp1b-display.json" \
  --output "$PROJECT/scene" --font-zh "$FONT_ZH" --font-en "$FONT_EN" \
  --style dark --recipe neutral --engine eevee --width 1280 --samples 24 --render
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_molecular_sequence_scene.py" -- --output "$PROJECT/scene/audit.json"
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- --directory "$PROJECT/scene/frames"
"$PYTHON" "$SKILL/scripts/encode_video.py" "$PROJECT/scene" --ffmpeg "$FFMPEG"
"$PYTHON" "$SKILL/scripts/audit_video.py" "$PROJECT/scene/preview.mp4" \
  --fps 24 --frames 624 --ffmpeg "$FFMPEG" --output "$PROJECT/scene/video-check"
```

调整节奏时复制两个 JSON 到项目中，修改反应配置 `timing`，重新准备并将输出 manifest 的 `captions` 更新到项目显示 JSON；保留需要翻译的主标签。五个事件必须按配置的先后发生。字幕边界由事件推导，旧字幕边界会被显示入口拒绝。去掉观察停顿仍只是连续有限片段，不会自动成为可循环过程。

## 实际验证范围

来源与运动的既有重建：2026-11-03 从 source-only 包的解压副本和原始结构缓存重建，23 个既有 NPZ 数组相同，旧显示的五张深色关键帧逐像素相同；625 个姿态、两图与标记、身份、字体和字幕通过。旧影片的整酶骨架特写后来被用户否定，不能将这些技术检查当作该构图获认可的依据。

当前局部显示另从解压副本重建，复用哈希相符的准备数据，不重算或改动反应运动。保存工程再次核对 2,482 个重原子和每图 625 个姿态，原坐标最大误差约 0.0000059 nm，字体、字幕、选定化学边与标记通过。五张深色关键帧中四张逐像素一致，另一张仅一个像素存在最大 1/255 的通道差异；不声称五张全部像素相同。1280×720、24 fps、Eevee 24 samples、neutral 的新版 26 秒影片完整解码 624 帧，抽取 15 帧与原 PNG 比较，并实际查看初态、中间体、底物离开、Gln262 转向、水进入和终态的六张编码帧；浅色初态和终态另看图。用户认可的是新关键帧的信息量，完整新版影片尚待用户反馈。没有执行原生实时播放，不将解码与抽帧检查写成实时逐秒观看。

同一显示代码另读取既有 AMSH-LP 数据，检查 721 个姿态、2,569 个原子和一个关键帧，未重做 DUB 影片。八项错误输入/运动反例被拒绝：原始来源变化、准备缓存变化、错误单位、旧字幕时间、事件逆序、同终点的错误旋转方向、非均匀视图缩放及显示边被改动。上述事实只支持所列机制范围与输出组合，不认证任意磷酸酶、所有合成器或生理动力学。
