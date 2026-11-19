# 翻译局部底物场景的参考外观

这一入口给已有 9SPI 教学场景接入共用角色材质、参考灯光、合成与宋体/Times 混排，不增加机制阶段。先按 [制作流程](production-workflow.md#翻译延伸的结构底物循环) 准备合格输入；完整来源、隐藏 40S、后移 60S、短肽符号和运动的证据边界继续有效。

9SPI 含大小写不同的链名。使用修复后的准备器生成 83 份独立链文件，并按 [链身份检查](data-import.md#链身份与文件系统大小写) 直接对照原始 mmCIF；不能沿用旧的 60 文件缓存。材质与动画重放通过不替代来源核对。

### 衔接起始、终止和回收镜头

用同一套分子坐标摆放来自不同标本的镜头，让观众在切换时认出小亚基一侧。以已有 9SPI 延伸工程的 `manifest.json` 为参照，先对齐完整的 18S 序列，再用双方已解析的同种碱基 C4′ 位置配准。密码子使用每份来源自己的作者编号。6YAL 的兔源起始复合物带有人 β-珠蛋白 mRNA。链 3 的 AUG（作者编号 52–54）和链 1 的起始 tRNA 可以同时建模。准备器逐碱基核对正占有率的核糖和碱基原子。5LZT 的 UGA 识别态中可见新生肽 54–68。另一标本 5LZU 提供野生型 eRF1 就位状态。5LZV 的 ABCE1 来自 eRF1(AAQ) 停滞构建体，亚基分离应依据另有功能证据的教学阶段。[β-珠蛋白起始复合物](https://www.rcsb.org/structure/6YAL) · [终止结构研究](https://pubmed.ncbi.nlm.nih.gov/27863242/) · [ABCE1 回收实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC2951834/)

```bash
# 在同一 sources 目录保留以下五份原始 mmCIF 和相邻的 .source.json。
# 对缺少的条目使用 fetch_structure.py（9SPI、6YAL、5LZT、5LZU、5LZV）。
"$PYTHON" "$SKILL/scripts/register_translation_states.py" --sources "$PROJECT/sources" \
  --elongation-manifest "$PROJECT/translation/manifest.json" --output "$PROJECT/registered"
"$PYTHON" "$SKILL/scripts/prepare_molecular_asset.py" "$PROJECT/sources/6YAL.cif" \
  --config "$PROJECT/registered/6YAL.json" --spacing 2.2 --sigma 2.6 --level 0.8 \
  --output "$PROJECT/48s-initiation"
```

每个条目得到可直接传给聚合物准备器的链选择和刚体坐标系。`registration.json` 还保留相同坐标系下的伴随组分原子，单位为 nm。组装镜头时一并放入对应的核苷酸类似物、Mg 或 ABCE1 铁硫簇。起始或终止标本中的缺失肽段、动力学和水解构象继续以各自证据为准。此入口产出构建输入，镜头运动和产物去向另行编排。[条目、角色和适用边界](../assets/translation-route-sources.json)

### 延伸特写的外观

从 [翻译外观配置](../assets/style-translation-reference.json) 复制项目配置。它参考已经查看的 WEHI 图片中的暖色精细前景、冷色分子背景和深度层次；这些艺术参数不是原作者节点或测量的光学性质。保持 tRNA 的分子颜色身份，不能按当前所在 A/P/E 位点重新染色。

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_translation_scene.py" -- \
  --asset-dir "$PROJECT/ribosome" \
  --style-profile "$SKILL/assets/style-translation-reference.json" \
  --font "$FONT_ZH" --font-latin "$FONT_EN" \
  --engine cycles --width 1600 --samples 32 --output "$PROJECT/translation"
"$BLENDER" -b "$PROJECT/translation/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_translation.py" -- \
  --asset-dir "$PROJECT/ribosome" --output "$PROJECT/translation/audit.json"
"$BLENDER" -b "$PROJECT/translation/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/translation/keys" --frames 1,80,150,185
```

`--font` 是宋体常规；`--font-latin` 是 Times New Roman 常规，角色/合成配方调用必须提供后者。省略 `--font` 时使用英文标签及 Times New Roman。字体打包进项目工程，不随 skill 分发。新配置不生成底部说明小字；来源和分层简化保留在 manifest 与文档中。

| 字段 | 实际作用 |
|---|---|
| `surface_roles.protein` | 已建模核糖体蛋白链的材质 |
| `surface_roles.nucleic_acid` | rRNA、三份 tRNA、mRNA 和其示意连接的材质 |
| `surface_roles.peptide_symbol` | 三至四残基短肽符号及示意连接；不是原子级肽转移化学 |
| `palette.alpha/beta/heads` | 大/小亚基 rRNA 与核糖体蛋白的语义颜色；40S 仍整条隐藏 |
| `palette.protein/rna/accent` | 名称沿用旧色板，分别跟踪 At/Pt/Et 三个 tRNA 身份；三者仍用核酸材质 |
| `palette.dna1/dna2` | 名称沿用旧色板，区分 mRNA 的 A 位点密码子和其余已建模残基，不表示 DNA |
| `palette.peptide/amino_acid` | 已有短肽与进入的 Ile 符号 |
| `palette.ink/secondary_ink` | 主标签与极性引线；不提供时根据底色选择 |
| `lighting/world/exposure/view_transform/look` | 真正生成的光源、世界填充与显示变换；不额外插入旧的硬编码补光 |
| `camera.scale_multiplier/center_offset_nm` | 保持原视角，调整正交取景尺度与相机 right/up 方向平移 |
| `context_blur_pixels_1280` | 完整 60S 背景层的显示虚化，按输出宽度缩放；0 为清晰 |
| `render_recipe` | 共用合成配方，文字最后独立叠加，背景单独柔化 |

完整 60S 作为相关结构背景，固定的原几何、nm 尺度和后移量保持不变，不添加无身份球团。背景、底物和字体分层只是合成安排，不能作为新的空间位置、解离状态或实验景深证据。

新配置不能包含旧 `compositor`、`representation` 或任意相机方向。这个镜头已对固定视角的路径与遮挡作过审核；换观察角度需重做路径和标签检查，不能悄悄改变模型或套用原通过记录。旧无角色/无新配方的调用保留旧行为；新调用用 `render_recipe`，不与命令行 `--compositor` 混用。

`audit_translation.py` 在新工程中同时运行 `audit_translation_appearance.py`：检查实际角色分配、材质节点、身份色、灯光/镜头、背景与文字路由、打包字体和所有帧的可见标签边界。原有结构、mRNA 位移、完整 tRNA 与连接检查继续执行。`test_translation_appearance.py` 会在内存中修改真实工程，确认错误角色、颜色、着色器、相机、字体或合成路由不能通过；它不保存破坏性修改。

实际渲染另行检查前景、背景和小型短肽符号的主次，以及动态文字和入场路径。固定源码与源输入重建后比较完整网格、身份、动画连接和姿态。该验证不证明真实翻译动力学、完整原子碰撞、起始、终止、GTP 化学或回收，也不代表所有合成配方已经逐一审美验收。
