# RNF4 / Ubc13 / Ube2V2：K63 双泛素的净键转移与释放

这个入口默认制作一个 **24 秒、双视图、在结合产物处结束**的教学片段：已装载泛素的 Ubc13 接近 RNF4 RING 二聚体，Ube2V2 定位受体泛素，供体 Gly76 与受体 Lys63 接近，硫酯连接转为异肽键。准备时加 `--release`，默认扩展为 **36 秒**，继续展示 E2–Ube2V2 成对离开以及完整双泛素产物释放。它不包含 E1 装载、泛素回收或蛋白酶体降解，不是完整泛素化循环。首尾物质状态不同，不能把文件直接反复播放就称为无缝反应循环。

适用于说明一种共价蛋白连接，补充结合/解离、磷酸转移和 GTP 水解示例。这里的公开结构、原子映射和构象参数只适用于所列构建体；不能换个蛋白名称就推广为另一种 E2/E3 或另一种泛素链。

## 来源与模型边界

| 来源 | 本例用途与限制 |
|---|---|
| [5AIT](https://www.rcsb.org/structure/5AIT)、[Branigan 等，2015](https://doi.org/10.1038/nsmb.3052)、[开放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4529489/) | 3.4 Å 交联构建体。author A 为连接起来的 RNF4 RING 结构域二聚体，B 为 Ubc13，C 为供体 Ub，D 为 Ube2V2。受体取 author F 的晶体对称邻居：`-x+y,-x,z+1/3`，晶格平移 `[1,1,0]`。不能把整个 asymmetric unit 当作一个天然复合物。 |
| [1J7D](https://www.rcsb.org/structure/1J7D) | author B 的未突变人 Ubc13，恢复 Cys87 / Lys92 的侧链身份与局部几何。保留 5AIT 主链，按 N/CA/C/CB 局部拟合；这是有来源的重构，不能把原结构 Lys87 的 NZ 改名为硫原子。 |
| [2JF5](https://www.rcsb.org/structure/2JF5) | K63 双泛素产物的连接身份参考：供体 B:Gly76 C 至受体 A:Lys63 NZ。自由产物不能直接刚性塞入 5AIT：整体构象不同，会穿入酶。其受体末端 Gly76 未建模，不能据此丢失本例的原子身份。晶体中的 Cd/Co/Mg/Cl 也不作为本反应的生理辅因子。 |
| [1UBQ](https://www.rcsb.org/structure/1UBQ) | 游离泛素末端具有 OXT 的身份参考。本例选择游离受体 Ub；在原交联构建体 G76 C/O/CA 平面补出 OXT。OXT 的位置为几何推断，未宣称测到 5AIT 中的游离末端构象。 |

5AIT 混合大鼠 RNF4、人 Ubc13/Ube2V2 和牛泛素；所画 RNF4 不是全长蛋白。保留其 **4 个结构锌离子**。脚本固定审核过的坐标文件 SHA-256，公开文件版本改变时需重新审核来源，不自动静默替换。

“完整展示”指保留所选结构全部已解析原子、不对模型做剖切，不代表补全了天然全长序列。这里实际保留 Ub 2–76、Ubc13 4–150、Ube2V2 6–145、连接 RING 构建体 author 131–259；未解析的 Ub Met1 等末端没有另行补造。用于末端特写、分子总质量或全长序列问题时，应先补齐并审核这些缺失，不能沿用这里的显示简化。

源结构中的交联并非天然硫酯键；晶体受体 Lys63 NZ 与供体 G76 C 相距约 12.49 Å。文章也指出仍需构象调整。因此反应构象明确为 `inferred`：蛋白主体刚性移动，局部用真实键轴上的侧链/主链扭转，不直接拉伸蛋白表面。Ube2V2 Ser32 与受体 Gly47 的接触用于约束定位；构象不是实验轨迹或 MD。

只画净连接变化，不指定一个已经证实的 Lys63 脱质子途径。相关机制仍有讨论；例如 [2025 年计算研究](https://doi.org/10.1021/acs.jpcb.5c01486)讨论了不同于简单 Asp119 一般碱图式的解释。这里不画确定的 Asp119→Lys63 质子接力。

释放段同样为 `inferred`。保留 E2–Ube2V2 的伙伴关系，允许 E2 Asp119、Ube2V2 Ser32 和受体泛素 Gln62 的局部侧链绕实际键轴转动；它们在离开后恢复初态。两个泛素随后作为保持异肽键的整体离开 RNF4。[Ubc13 构象动态研究的摘要](https://pubmed.ncbi.nlm.nih.gov/25423605/)支持讨论活性位点附近的柔性，但没有给出这里的释放路径或侧链角度。本例不把这些人为编排的角度、先后顺序或播放时间归为该研究的测量结果，也不声称重建了其主链门控运动。

## 化学身份与检查

- 供体始终为 `donor`，受体始终为 `acceptor`。断开 `e2:87:SG — donor:76:C`，形成 `donor:76:C — acceptor:63:NZ`。不能把连接接到受体 Gly76、K48 或另一个 Ub 身上。
- 初态硫酯 C—S 约 1.82 Å，末态异肽 C—N 约 1.34 Å。供体无游离末端 OXT；受体保留自由末端 OXT。4476 个选定重原子（含四个结构 Zn）全程保留。
- 这个局部转移不消耗 ATP，不放出 Pi、PPi 或水。氢与溶剂/缓冲液未显式建模；不要凭本片推断完整质子计量。
- 所有蛋白伙伴以及非 RING 蛋白与 Zn 的接触逐帧检查。成键前后按实际连接图排除 1–2 / 1–3 共价邻居；这是化学拓扑的例外，不允许把周围任意原子列入免检名单。过渡片段允许两种端态的反应键邻域，不能把它解释为确定的基元反应中间态。
- 局部移动原子与自身其他非键邻居也检查；原有键长及相邻两键端点距离保持不变，避免“没有穿插但蛋白被拉长”。静态 RING/Zn 配位保持原位。
- 端态还检查硫酯羰基平面、反应键角与产物酰胺的平面/反式几何。只满足 C—N 距离并不等于正确成键；把受体绕新 C—N 轴转 90° 的负对照必须失败。
- 释放段继续检查新 C—N 键与酰胺几何、新增可动侧链的自身接触，以及 E2–Ube2V2 伙伴界面 Cα 对的距离守恒。把两个泛素分开移动，或只移动 Ube2V2 而丢失伙伴关系，必须被负对照检出。最终产物与其他结构的重原子间距需大于 8 Å，这只是画面分离标准，不是生物学解离阈值。
- 2.38 Å 是本例重原子中心的几何筛查阈值，不是完整 van der Waals、溶剂、自由能或动力学判据。
- K63 链有情境依赖的功能；不能将其一律画成自动送入蛋白酶体。

## 重建

在 skill 外使用新项目目录。准备步骤使用 NumPy、SciPy、Gemmi、scikit-image；场景生成器使用 Blender 5.x（本例实际检查 5.2.3）。中文传入宋体常规字体，英文传入 Times New Roman。字体只打包进项目 `.blend`，不分发进 skill。

```bash
PYTHON='/actual/python/with/numpy-scipy-gemmi-skimage'
BLENDER='/actual/path/to/blender'
SKILL='/actual/path/to/blender-cell-biology'
PROJECT='/absolute/new/k63-project'
SONGTI='/actual/path/to/Songti-Regular.ttf'
TIMES='/actual/path/to/Times New Roman.ttf'

for pdb in 5AIT 1J7D 1UBQ 2JF5; do
  "$PYTHON" "$SKILL/scripts/fetch_structure.py" "$pdb" --output "$PROJECT/sources"
done

"$PYTHON" "$SKILL/scripts/prepare_ubiquitin_transfer.py" \
  --sources "$PROJECT/sources" --output "$PROJECT/assets" --release --fps 24

"$PYTHON" "$SKILL/scripts/test_ubiquitin_geometry.py" \
  --sources "$PROJECT/sources" --release --output "$PROJECT/negative-controls.json"

"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_ubiquitin_transfer_scene.py" -- \
  --asset-dir "$PROJECT/assets" --output "$PROJECT/scene" \
  --font-zh "$SONGTI" --font-en "$TIMES" \
  --style dark --recipe neutral --engine eevee --width 1280 --samples 24 --render

"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_ubiquitin_transfer_scene.py" -- \
  --output "$PROJECT/scene-audit.json"

"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- --directory "$PROJECT/scene/frames"

"$PYTHON" "$SKILL/scripts/encode_video.py" "$PROJECT/scene" \
  --fps 24 --ffmpeg /actual/path/to/ffmpeg
```

`assets/ubiquitin-transfer.json` 保存可重建的扭转、刚体位姿、接近弧线和播放节奏参数，不含生成网格或图像。PCHIP 插值后重新计算全蛋白原子，再从每帧原子生成显示表面。修改参数或帧数后重新运行几何检查；失败则停止生成场景，不靠取消检查来继续。

不传 `--release` 时保持原有 576 帧的结合产物示例；传入后默认 864 帧。`--frames` 可明确指定其他总帧数，仍需检查全部输出帧。两种选择都可以另行调整渲染引擎与合成器，不等于无缝循环模式。

## 画面与验收

两侧共享泛素坐标但各自取景、放大。左侧保留完整选定复合物，右侧完整省略 E2、UEV 和 RING，展示两个完整泛素；不切蛋白或只显示断开的残基片段。镜头平滑收拢，保留进入画面的整分子。释放路径在离开结合区域后平滑转向侧方；最终同时检查三维原子间距与正面投影的分离，避免“空间中分开、画面中仍重叠”。这条弧线为教学安排。Gly76—Lys63 使用高反差连接标记，线宽为教学可读性放大，不代表原子半径。标题、图例、关键标签清楚，来源细节留在项目文档，不塞进视频小字。

默认背景使用本例公开结构的完整 Ub/E2/UEV 形状，分层虚化；位置、投影大小和分子密度都是艺术安排，不代表某种细胞或生理丰度。它们不会进入前景的反应账本。

可选 `--style dark/light`，渲染引擎 `--engine eevee/cycles`。本入口支持图像级合成器 `neutral`、`soft-cinema`、`fluorescent-soft`、`warm-print`、`cool-duotone`、`crisp-neutral`，文字绕过滤镜。依赖深度、法线或对象索引的其他配方需先为双视图合并相应 pass，不能直接声称全部适用。改色后重新核对分子与图例的颜色对应。

`cool-duotone` 的实图会明显压缩五种分子的色差，适合低饱和外观探索；依靠颜色辨认参与者的教学片段默认选 `neutral` 或 `crisp-neutral`。如需保留双色调，应增加逐个主体的清晰名称或其他身份提示，再检查可辨性，不能仅凭渲染成功就通过教学验收。

交付前实际检查至少进入、定位、转移中、结合产物以及所选释放阶段，并检查编码视频。检查器重新打开 `.blend`，核对所有输出帧的点坐标、原子身份、双视图共享坐标、等比例取景、完整主体边界、字体和阶段文字。产物连接的引线逐帧跟随实际投影，不能在分子移动后仍指向旧位置。负对照包含错误化学状态、错误 K48 位点、分子穿插、拉长键和释放时丢失连接；漂亮图片不替代这些检查。

**仍待扩展：**E1 装载、完整泛素化循环、其他链型、底物连接、DUB 回收和降解过程。这些需要各自的新来源与新验收，不能用本例的局部释放替代。当前释放路径的支持范围是所列构建体与参数下的连续几何检查；未建立生理速率、能量或唯一解离次序。
