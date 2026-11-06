# 从源装配变换生成 clathrin 被覆

此入口完成可编辑的静态源装配，不宣称完成被覆生长或内吞。[1XI4](https://www.rcsb.org/structure/1XI4) 是牛来源的 7.9 Å Cα 模型；[原研究](https://www.nature.com/articles/nature03079)使用体外装配被覆。重链已建模 1–1630，轻链仅有 95–164；保持这些缺失范围，不自行补出完整天然序列。显示包络来自 Cα 高斯场，不是全原子表面或实验密度。

## 沿用与新增的部分

沿用 `fetch_structure.py`、`prepare_molecular_asset.py --config`、共享网格实例、材质/相机/合成器及 `render_saved.py`。新建的内容是源单体分组、装配操作、完整装配对应关系和保存后的身份检查。三种源单体共 18 份链网格，经 12 个源操作生成 36 只已建模三脚单体、216 条链；没有按球面随机摆放蛋白，也没有用正六边形代替源被覆。

`prepare_clathrin_assembly.py` 的 `rigid` 和 `compare_reference` 可用于其他有明确链身份的源装配：先检查变换，再与独立下载的完整装配逐坐标、逐身份、一一对应，不能只比较质心或总数。平移单位由 Å 转为 nm。当前主入口的 Cα 范围和 `check_groups` 中重链末端、轻链归属规则专用于 1XI4；换蛋白需重新定义分组和来源，不能只改名称。

## 实际构建

已有哈希相符的来源缓存可复用；缺少文件时再下载。所有输出放在 skill 外，`PYTHON` 需有 Gemmi、NumPy、SciPy 和 scikit-image。

```bash
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 1XI4 --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 1XI4 --assembly 1 --output "$PROJECT/sources"
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

2026-11-06 从解压脚本和已核对的原始来源缓存独立重建，144 个数组及装配 manifest 相同；总览和单体保存工程通过。检查了深色总览/单体及浅色单体静帧，没有完整视频。单体重渲染只有 13 个像素存在至多一个 8-bit 通道值差异，未声称位级相同。10 个来源/分组反例和 7 个工程破坏反例被拒绝。文件坐标的细小对应误差只说明算术重现，不是该低分辨率实验的精度。

```bash
"$PYTHON" "$SKILL/scripts/test_clathrin_assembly.py" \
  --sources "$PROJECT/sources" --assets "$PROJECT/assets" --output "$PROJECT/source-negative-tests.json"
```

## 刚体结合动作的迁移边界

本次实测使用 `interaction_motion.plan` 搜索 20 nm 接近路径、120 帧教学时序及 0.4 nm Cα 接近诊断阈值。目标 `1-ABC` 在三邻居 `1-DEF/2-ABC/6-ABC` 和两邻居 `1-DEF/2-ABC` 两种固定环境中，搜索所得最好路径仍分别新增或加深 143、52 个位点的重叠。源末态本身的相同原子对例外仍保留，未放宽阈值让路径通过。这只否定本次测试路径作为可用动画，不证明所有刚体路径不存在，也不等于全原子碰撞结论。

因此当前可用产物是完整源被覆和单体模型；不要用逐个显示单体、径向穿入或剪辑跳变冒充已验证装配。下一步需选择有依据的开口晶格顺序或构象/腿部变化，重新检查邻接和交接。AP2 核心、未解析铰链、膜曲率及去被覆仍需单独衔接；不能把该牛来源体外模型直接写成既定 SUM159 细胞路线。
