# clathrin 源被覆与局部单体接纳

此入口生成可编辑的静态源装配，以及完整已建模单体在三邻居晶格片段中的接纳/释放示意。[3IYV](https://www.rcsb.org/structure/3IYV) 是牛来源的 7.9 Å Cα 模型；[原研究](https://www.nature.com/articles/nature03079)使用体外装配被覆。重链已建模 1–1630，轻链仅有 95–164；保持这些缺失范围，不自行补出完整天然序列。显示包络来自 Cα 高斯场，不是全原子表面或实验密度。局部运动不代表完整被覆生长、去被覆或内吞。

## 来源修正：作者链名不保证完整单体

此前 1XI4 的静态重建证明了源装配坐标对应，但错误地把同名重链全部片段当成一个连续三脚单体。实际检查发现其 838/839 与 1279/1280 相邻编号之间存在约 14–15 nm 跳跃。旧单体视图和以这些单体进行的路径搜索不能继续作为完整单体的证据。

默认改用同一研究配套的 3IYV；其条目明确按完整已建模三脚单体组织。保留源 Cα 坐标、缺失范围和长步长，不人工拉直连接。新检查拒绝空间不连续的链，并记录大于 4.5 Å 的粗略 Cα 步长；3IYV 最大约 9.46 Å，10 Å 的门槛只用于排除跨区域片段，不证明正常肽键几何或全原子精度。旧准备清单缺少这一检查时需要重新准备。

## 沿用与新增的部分

沿用 `fetch_structure.py`、`prepare_molecular_asset.py --config`、共享网格实例、材质/相机/合成器及 `render_saved.py`。新建的内容是源单体分组、装配操作、完整装配对应关系和保存后的身份检查。三种源单体共 18 份链网格，经 12 个源操作生成 36 只已建模三脚单体、216 条链；没有按球面随机摆放蛋白，也没有用正六边形代替源被覆。

`prepare_clathrin_assembly.py` 的 `rigid` 和 `compare_reference` 可用于其他有明确链身份的源装配：先检查变换，再与独立下载的完整装配逐坐标、逐身份、一一对应，不能只比较质心或总数。平移单位由 Å 转为 nm。当前主入口的 Cα 范围和 `check_groups` 中重链末端、轻链归属规则专用于 3IYV；换蛋白需重新定义分组和来源，不能只改名称。

## 实际构建

已有哈希相符的来源缓存可复用；缺少文件时再下载。所有输出放在 skill 外，`PYTHON` 需有 Gemmi、NumPy、SciPy 和 scikit-image。

```bash
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 3IYV --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 3IYV --assembly 1 --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_clathrin_assembly.py" \
  --sources "$PROJECT/sources" --output "$PROJECT/assets"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/blender_clathrin_assembly.py" -- build \
  --assets "$PROJECT/assets" --output "$PROJECT/cage.blend"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/blender_clathrin_assembly.py" -- build \
  --assets "$PROJECT/assets" --output "$PROJECT/unit.blend" --show 1-ABC
"$BLENDER" --disable-autoexec -b "$PROJECT/cage.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/blender_clathrin_assembly.py" -- audit --output "$PROJECT/cage-audit.json"
"$BLENDER" --disable-autoexec -b "$PROJECT/unit.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/blender_clathrin_assembly.py" -- audit --output "$PROJECT/unit-audit.json"
"$BLENDER" --disable-autoexec -b "$PROJECT/cage.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- --directory "$PROJECT/cage-review" --frames 1
```

`--show` 只选择完整单体；其他单体仍保留在工程中。`--highlight` 指定强调的单体，金色为其已建模重链，粉色为其已建模轻链片段，其他单体为蓝色。不会剪断长腿。当前模型无文字，说明写在交付说明中。可用既有 `--compositor clean/reference/monochrome` 分支；单色分支不再承担颜色身份区分。

`--style-profile "$PROJECT/style.json"` 可替换项目配色、背景、灯光、粗糙度和表面纹理；从 `assets/style-molecular-reference.json` 的完整结构复制后修改。此入口用 `palette.protein/dna1/alpha` 分别表示强调重链、强调轻链、其他单体，并以正交相机检查完整形状。先按参考图分析这些参数，再实看少量静帧，不改变源坐标。实际使用的配置保存在 `.blend` 的 `style_profile` 中。

保存检查核对全部 216 条链、所属单体、源操作矩阵、共享网格、显示选择和完整可见链的取景。生成资产和 manifest 均绑定哈希。源分组、文件变更或模型改动后重新准备；保存检查不证明生理曲率、原子接触或装配动力学。

此前 1XI4 的重建记录仅保留为历史数值复现证据，完整单体资格已由上述检查撤回。3IYV 的新来源检查增加空间连续性和实际旧来源反例。文件坐标的细小对应误差只说明算术重现，不是该低分辨率实验的精度。

```bash
"$PYTHON" "$SKILL/scripts/test_clathrin_assembly.py" \
  --sources "$PROJECT/sources" --assets "$PROJECT/assets" --output "$PROJECT/source-negative-tests.json"
```

## 从刚体路径迁移到有来源依据的局部构象

修正单体身份后，原有刚体弧线仍未找到合格路径；不能直接径向穿入、逐个显示或切短长腿来掩盖这一问题。联网阅读 [Morris 等的原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC7100586/) Results 的 cage structural modelling 与 Figure 5 后，采用其分段拟合范围附近的 840、1130、1280 作为粗粒度转动支点；12°、17°、18° 为本例受该研究启发的构造上限。该研究使用猪来源被覆，不能把这些上限称为牛蛋白的实测运动限制，更不能把静态构象差异当成观测到的装配轨迹。研究还修订了旧 Cα 模型的部分螺旋取向，因此本例不用于宣称精确适配蛋白接触位点。

`articulated_ca_motion.py` 保留每段刚性、相邻 Cα 源距离及完整身份，用整单体平移/转动和有限的分段转动搜索释放路径，再反向用于接纳。`prepare_clathrin_growth.py` 是 3IYV 专用适配器；运动单体 `1-ABC`，邻居 `1-DEF/2-ABC/6-ABC` 均完整保留已建模的三条重链和三条轻链片段。轻链片段随近端主体移动，未给它们增加虚构伸长。

接触筛查覆盖运动单体与固定邻居，以及自身距离可变的不同刚性片段；同一刚性段内部和同链相隔不超过三个残基的局部邻接另行处理。0.4 nm 是 Cα 接近诊断阈值，不能解释为全原子排斥半径。仅源末态已经存在的相同 Cα 对可沿用源距离，并允许至多 0.01 nm 的额外接近；换位点不能继承例外。数值通过不证明能量可行、结合亲和力或生理路径唯一。

```bash
"$PYTHON" "$SKILL/scripts/prepare_clathrin_growth.py" \
  --assets "$PROJECT/assets" --config "$SKILL/assets/clathrin-growth.json" \
  --output "$PROJECT/growth" --frames 240 --fps 24 --mode teaching
"$PYTHON" "$SKILL/scripts/audit_clathrin_growth.py" \
  --sources "$PROJECT/sources" --config "$SKILL/assets/clathrin-growth.json" \
  --motion "$PROJECT/growth" --output "$PROJECT/growth-raw-audit.json"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/blender_clathrin_growth.py" -- build \
  --motion "$PROJECT/growth" --output "$PROJECT/growth.blend" \
  --font-zh "$FONT_ZH" --font-en "$FONT_EN" --engine eevee
"$BLENDER" --disable-autoexec -b "$PROJECT/growth.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/blender_clathrin_growth.py" -- export --output "$PROJECT/saved-samples"
"$PYTHON" "$SKILL/scripts/audit_clathrin_growth.py" \
  --sources "$PROJECT/sources" --config "$SKILL/assets/clathrin-growth.json" \
  --motion "$PROJECT/growth" --saved-samples "$PROJECT/saved-samples" \
  --output "$PROJECT/growth-scene-audit.json"
"$PYTHON" "$SKILL/scripts/test_clathrin_growth.py" \
  --sources "$PROJECT/sources" --config "$SKILL/assets/clathrin-growth.json" \
  --motion "$PROJECT/growth" --output "$PROJECT/growth-negative-tests.json"
"$BLENDER" --disable-autoexec -b "$PROJECT/growth.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- --directory "$PROJECT/keyframes" --frames 1,65,121
"$BLENDER" --disable-autoexec -b "$PROJECT/growth.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- --directory "$PROJECT/short-review" \
  --frames 41:89:2 --width 854 --samples 8
"$PYTHON" "$SKILL/scripts/encode_review.py" "$PROJECT/short-review" --ffmpeg "$FFMPEG"
```

`--mode loop` 改为持续往复时序；在新输出目录准备、构建和检查。它是局部接纳/释放的教学循环，不是包含 Hsc70/ATP 的去被覆机制。默认独立结构视图便于查看这个体外被覆片段；新项目需要细胞环境时再按区室加入有身份的景深层，不能把无膜示意冒充细胞中的完整被覆。

Blender 保留完整源表面，通过共享的几何节点给表面与 Cα 检查点施加相同转动。顶点按最近源 Cα 归属片段，接缝处只是粗粒度显示插值；相邻 Cα 距离不变也不代表肽链扭转角与全原子几何正确。保存后的实际四分之一帧采样另查源身份、刚性变换、分段转动、接触、首尾、固定邻居和完整表面取景。采样不是连续时间证明。临时采样数组可在检查并记录哈希后清理，当前来源和可重建工程保留。

当前教学实例的保存场景检查覆盖 961 个姿态、5,100 个运动 Cα 和 15,300 个固定 Cα；最大转动约 9.32°，没有全局拉伸。实际检查了三个正式尺寸关键帧与 2.083 秒原速抽帧短预览，未渲染完整影片。AP2/PIP2、膜形变、全晶格装配顺序、Hsc70、ATP、切断和货物运输仍需按目标机制重新接入，不能把该牛来源体外模型直接写成既定 SUM159 细胞路线。

迁移到别的柔性装配时，可沿用约束搜索、源对例外、时序与保存后采样方法；必须重新联网确定单体身份、状态、柔性段、支点/自由度、伴随配体和适用的接触规则。当前通用运动函数不自动预测这些生物学条件，也不支持直接套用本例的牛源分组或角度值。
