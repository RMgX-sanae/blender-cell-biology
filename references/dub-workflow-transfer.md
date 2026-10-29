# 从 K63 连接工作流迁移到 AMSH-LP 水解

这个案例检验“反应规则改变后，原有流程还能复用多少”。它从已经结合的实验双泛素构建体开始，展示水氧并入远端泛素羧基、K63 异肽键断开、两个产物离开。**没有实现底物进入、逐步质子转移、溶剂重新装载或完整泛素回收通路；不能称为完整催化循环。**

## 迁移时保留和改变的内容

| 原工作流的部分 | 本例实际沿用 | 新增内容及原因 |
|---|---|---|
| 结构与原子资产 | `atom_table`、`xyz`、`surface`、`peptide_bonds` | AMSH-LP 构建体选择、Glu292 恢复、两个 Zn 的区分与水的位置 |
| 共价几何 | `ubiquitin_geometry.index/rotate`；绕键转动和整分子刚体运动 | 水解端态与 24 个局部扭转/刚体参数；不能倒放连接动画 |
| 反应身份 | 原有 `reaction_ledger.validate` | 水中的氧映射为远端 Gly76 的 OXT；断 C–NZ、成 C–O，元素和原子身份保持 |
| 释放运动 | `interaction_motion.matrices` 和原有五次缓动 | 近端可用刚体路径；远端末端必须先改变扭转，再离开催化通道 |
| 分子显示 | K63 的 `keyed_points`、`envelope`、`annotation_line`；PKA 的 `animated_mesh/visibility` | 右图显示两条完整泛素的原子/键，整体和细节共同追踪同一个水氧 |
| 环境和合成 | 同一个 `background`、`paired_composite`、`compositor_stack.build_stack` | 背景输入换成本例完整蛋白形状；继续使用受保护字体层和同步双视图 |
| 渲染交付 | 原有 `render_saved.py`、`encode_video.py`、`audit_video.py` | 保存后逐帧核对本例的身份、坐标、产物分离和水氧标记 |

新增的 `dub_geometry.py` 与准备脚本是**这个结构组合的科学适配器**，不是通用水解酶求解器。`assets/dub-hydrolysis.json` 保存经过检查的控制姿态、时序和字幕；它不是渲染缓存。换蛋白时不能只改名字和坐标文件，必须重新判断构建体、反应及可通过的路径。

控制参数使用 Python 的零起点切片：`0:6` 为远端 74–76 的 φ/ψ，`6:10` 为近端 Lys63 的四个 χ，`10:13`/`16:19` 分别为远端/近端的整分子旋转向量，`13:16`/`19:22` 为对应平移，`22:24` 为远端 Arg74 的两个 χ。角度与旋转向量用弧度，内部坐标和平移用 Å；准备输出才换算为 nm。对这些参数的调整必须重新通过接触和共价几何检查。

## 来源和必要修正

- [2ZNV](https://www.rcsb.org/structure/2ZNV)，1.60 Å：选择作者链 A 的 AMSH-LP 265–436、B 远端泛素 1–76、C 近端泛素 1–77。A 是 E292A 失活构建体；B 的 K63R 和 C 的 D77 延伸保留并明确记录。它们不是完整野生型蛋白组合。
- [2ZNR](https://www.rcsb.org/structure/2ZNR)，1.20 Å：提供天然 Glu292 局部侧链、催化 Zn 及水氧 A:517 的参考。仅恢复 Glu292 侧链，使用结合态酶骨架。催化金属按 His347、His349、Asp360 的配位原子定位；与 2ZNV 中已有的结构 Zn 分开。
- [Sato et al., 2008](https://www.nature.com/articles/nature07254) 的 [Figure 4](https://www.nature.com/articles/nature07254/figures/4) 支持金属结合水攻击和异肽键水解的机制解释；[勘误](https://www.nature.com/articles/nature07515) 将 Figure 3c 的 Glu324 更正为 Asp324。[2014 Sst2 研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC4033627/)提供产物/局部开合的比较依据，但酵母 Sst2 产物结构不能冒充人 AMSH-LP 的实测端态。

直接叠合不能作为可用模型：所选 2ZNV 的 Gly76 C—Lys63 NZ 距离为 1.736 Å；整个 apo 酶叠入结合复合物会出现约 0.715 Å 的严重接触；将参考水直接局部叠合，水氧距羰基碳只有约 1.303 Å。准备脚本因此保持结合态识别界面，用共价扭转调整底物端态，并显式重建配位水。所得反应物 C–N 约 1.342 Å，产物 C–水氧约 1.250 Å，属于几何约束下的解释性重构。

当前已读范围是 2008 摘要、Figure 4 图和图例、补充 Table S2 与 S4 图/图例、勘误，以及 2014 机制/产物/局部开合相关正文与图例；没有阅读付费墙后的 2008 全文、全部补充图或完整 wwPDB validation reports。需要原子级论证或新构建体时继续补读，不把这些缺口改写为全面审查。

## 可复用的失败处理方法

先测试原有刚体路径。近端泛素通过；远端泛素的末端与催化通道相撞，且完整底物的进入路径也失败。应先定位失败原子和受限部位，再判断哪些共价扭转可以动。本例只为远端产物补充末端扭转路径，保留每条内部键及相邻键角，水氧随末端局部坐标系一起移动；不拉伸表面、不隐藏穿插、不反复增加无关全局随机路径。

通过控制点并不等于通过插值：初版柔性路径在控制点间仍穿插。增加必要的局部控制后，用密集采样检查中间段，再按最大原子移动距离重定时，最后重新检查实际输出的每个姿态。解离次序、路径及秒数均是教学设计，不能解释为测量的动力学。

底物进入路径的失败仍然保留为范围缺口。不能因为产物能离开，就声称底物能以反向路径进入；两者的共价状态和需要通过通道的对象不同。

## 重建命令

使用 Python 环境中的 NumPy、SciPy、Gemmi、scikit-image，以及 Blender 5.2、FFmpeg。`SKILL`、`PROJECT`、`PYTHON`、`BLENDER`、`FFMPEG` 和两个字体路径替换为本机实际路径；输出放在新项目目录。字体来自用户本机，skill 不分发字体。

```bash
# 已有与准备脚本哈希相符的缓存时直接复用；否则下载公开原始文件。
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 2ZNV --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/fetch_structure.py" 2ZNR --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_dub_hydrolysis.py" \
  --sources "$PROJECT/sources" --output "$PROJECT/assets" --frames 720 --fps 24

"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_dub_hydrolysis_scene.py" -- \
  --asset-dir "$PROJECT/assets" --output "$PROJECT/dark" \
  --font-zh "$FONT_ZH" --font-en "$FONT_EN" \
  --style dark --recipe neutral --width 960 --samples 64 --render

"$BLENDER" -b "$PROJECT/dark/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_dub_hydrolysis_scene.py" -- \
  --output "$PROJECT/dark/scene-audit.json"

# 查看关键帧后，从保存工程渲染和检查真实编码视频。
"$BLENDER" -b "$PROJECT/dark/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- --directory "$PROJECT/dark/frames"
"$PYTHON" "$SKILL/scripts/encode_video.py" "$PROJECT/dark" --ffmpeg "$FFMPEG"
"$PYTHON" "$SKILL/scripts/audit_video.py" "$PROJECT/dark/preview.mp4" \
  --fps 24 --frames 720 --ffmpeg "$FFMPEG" --output "$PROJECT/dark/video-check"
```

使用 `--style light` 另建浅色工程；合成配方通过 `--recipe` 独立选择。上面的 64 次 Eevee 采样对应已验收示例的实际保存设置；降低采样后需重新检查细线和运动噪声。不要把一个配方在 K63 中通过，当作它在新结构上已经通过。默认双视图有不同放大率，水氧及共价键的画面标记有所加粗，背景是源蛋白形状的艺术性分布，不表达该细胞中的实际丰度或相同景深尺度。

准备检查包括 2569 个所选重原子的身份映射、两种 Zn、末端键长度/平面/角度、移动部分的自接触、分子间接触、内部键和相邻键角保持。2.38 Å 是本例的重原子中心接触筛查阈值；只对明确的配位原子和反应键邻域作例外。它不代替能量、键级、电荷或质子路径检查。几何通过后，还要检查保存工程和实际像素中的水氧可见性、整分子释放及文字。
