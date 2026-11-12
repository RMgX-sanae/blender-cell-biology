# 有来源端态的 Cα 构象路径

当已有受约束的 Cα 路径，需要把它可靠地迁移到可编辑 Blender 场景时使用此入口。它**消费经过审查的路径，不从两张静态结构自动预测运动**。按项目选择并另行生成候选坐标，再用本入口检查来源、连接与连续线性路径。需要整体刚体或分段铰链转动时沿用 [clathrin 工作流](clathrin-assembly.md)；旋转曲线与其直线弦不是同一路径，不能用本检查器替代旋转路径检查。

## 已有代表输入的边界

[6YAE](https://www.rcsb.org/structure/6YAE) 提供缓冲液中的 AP2 核心，[6YAF](https://www.rcsb.org/structure/6YAF) 提供膜和货物环境中的另一构象。[原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC7375805/) 使用重组构建体；不要将跨物种构建体写成某一种天然细胞的完整过程。

`assets/ap2-ca-conformation.json` 对接已有的 508 点构造路径，固定保留 6YAE 所选 A/B/M/S 链的 1,717 个已建模 Cα。由原始 mmCIF 与 SIFTS 重新确认 1,711 个端态对应；M136–141 的开放端位置是推断。6YAF 另有 51 个未匹配的已建模核心残基，以及单独 P 链的六残基货物肽，均在配置中逐项说明未纳入此固定身份路径的原因。它不是完整 6YAF 复合体，也不展示货物结合、PIP2 化学、AP2 完整激活或内吞。

这份源配置绑定已有候选路径的哈希，**skill 包不包含该轨迹，也没有宣称已将原项目的约束求解器通用化**。复用已有工程时提供其受审查的 `motion.npz`；新对象必须提供自己的路径和来源配置，不能只改蛋白名称，不能用旧视频或 `.blend` 反向充当生成输入。

## 来源与路径交接

准备命令在 Blender 外运行，需要 NumPy、SciPy、Gemmi。已有哈希相符的来源缓存可直接复用；新目标按 [联网结合与组装](online-binding-research.md) 核定构建体、端态和必要伙伴。

```bash
"$PYTHON" "$SKILL/scripts/prepare_ca_conformation.py" \
  --sources "$PROJECT/sources" --path "$PROJECT/authored-motion.npz" \
  --config "$PROJECT/conformation.json" --output "$PROJECT/prepared"
```

候选 NPZ 中 `keys` 为 `(N,2)` 的作者链名与完整残基编号字符串，`points` 为 `(K,N,3)` 的 nm 坐标；顺序与第一来源所选链的已观察 Cα 顺序相同。插入码保留为字符串。输入中的其他端态数组不能替代从原始来源重新读取的坐标。

以 AP2 配置为可执行格式参考，但重新填写目标的来源哈希、链/UniProt 映射、对齐链、插入残基依据、两侧未匹配范围、遗漏伙伴及几何限制。第一来源整体刚体对齐到第二来源的指定对应链，路径首尾分别与原始端态核对。遗漏必要配体或伙伴会改变要讲述的机制时，应使用保留这些组分的全原子/复合体工作流，不能把 Cα 核心视图当成完整替代。

`source_edges` 从 mmCIF 的 polymer label sequence 建立连续已建模连接，不用作者编号相减猜邻接。缺失区域不自动连成肽链；只有明确声明的缺失段端点距离上限参与检查。

`ca_path_constraints.inspect` 对每段线性坐标计算非局部点对的最小距离和相邻点距离的极值，并检查声明的半空间、缺失段端点跨度与最大步长。局部接触排除沿已建模连接图确定。源接触例外仅属于相同点对，且只允许真实观察端态参与；推断端态不能授权例外。AP2 示例用 0.32 nm Cα 接近阈值和 0.005 nm 源对余量，不能解释为全原子范德华半径，也不能直接作为其他精度来源的默认生物学标准。

这些检查不验证原子侧链、手性、能量、完整包络碰撞或真实时间尺度。线性路径通过也不证明真实分子选择此路径。原来的候选生成方法和约束依据仍需作为项目输入保留。

## 风格、保存场景与高效预览

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_ca_conformation_scene.py" -- \
  --asset-dir "$PROJECT/prepared" \
  --display "$SKILL/assets/display-ap2-ca.json" \
  --style-profile "$SKILL/assets/style-ca-conformation.json" \
  --font-zh "$FONT_ZH" --font-en "$FONT_EN" \
  --output "$PROJECT/scene" --width 1600 --samples 32 --engine cycles
"$BLENDER" --disable-autoexec -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_ca_conformation_scene.py" -- \
  --asset-dir "$PROJECT/prepared" --display "$SKILL/assets/display-ap2-ca.json" \
  --style-profile "$SKILL/assets/style-ca-conformation.json" \
  --output "$PROJECT/saved-audit.json"
"$BLENDER" --disable-autoexec -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/keys" --frames 1,121,241
"$BLENDER" --disable-autoexec -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/review" --frames 97:145:2 --width 800 --samples 12
"$PYTHON" "$SKILL/scripts/encode_review.py" "$PROJECT/review" --ffmpeg "$FFMPEG"
```

上述帧号属于 AP2 的 240 区间有限路径。所有候选采样都保留，按累计位移和少量均匀时间分量排序，再采用平滑起止的展示时序；不删掉几乎静止的求解结点来规避 Blender 的时间键合并。它是有限的连续展示，不能直接标为循环。

`display` 控制观察方向、景深、Cα 显示包络与可选膜环境，`style-profile` 为每个链及环境实体指定角色材质、语义颜色、灯光和显示变换。参考图的暖色主体、冷色环境与空间层次进入实际节点；Cα 包络仍是粗粒度形体，不冒充参考图的原子空间填充细节。两个膜叶独立布置完整脂质符号，不切开蛋白或脂质。当前膜是与来源坐标方向一致的平面片段，脂质位置、尺寸和化学结构都是示意；避让只筛查路径结点，不能宣称连续包络—膜碰撞已全部通过。

保存后检查全部形态键、源身份和时间，评估所有结点、结点中点及整数帧；核对真实生效的密度链路、材质、配色、灯光、字体及少数完整表面取景。数值检查之外仍需查看正式尺寸关键帧与短预览。包中只保留脚本、参数、来源说明和测试，轨迹、工程、字体与渲染保存在项目中。
