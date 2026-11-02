# 蛋白质相互作用：从实际结合构象推广动作

先判断要表达的作用，再选择资产与事件。蛋白结合、催化转化、共价修饰、核苷酸交换、聚合和跨膜运输不是同一种动作；不能只更换两个蛋白的名字。

推广到其他目标时先用 [工作流迁移方法](workflow-transfer.md) 区分配置替换与机制改变；其中给出了 `--examples` 项目配置的实际用法，以及可以从 PKA、Ras、K63 和装配流程复用的入口。

## 已提供的局部结合入口

`assets/interaction-examples.json` 是可扩展的明确选择表：PDB、作者链、固定/移动伙伴、保留配体、标签和构建体边界。`prepare_interaction.py` 保留完整已建模链，单独保留指定配体，生成来源锚定的动作；`build_interaction_scene.py` 保存可编辑模型、全部帧关键帧与完整分子环境。

```bash
python3 "$SKILL/scripts/fetch_structure.py" 1BRS --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_interaction.py" \
  --example barnase-barstar --sources "$PROJECT/sources" \
  --output "$PROJECT/assets" --mode loop --frames 192 --fps 24
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_interaction_scene.py" -- \
  --asset-dir "$PROJECT/assets" --output "$PROJECT/scene" \
  --font-en /installed/fonts/Times-New-Roman.ttf --style cinematic --render
```

完整视频及实际工程检查：

```bash
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_interaction_scene.py" -- \
  --output "$PROJECT/scene/scene-audit.json"
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/scene/frames" --frames all
python3 "$SKILL/scripts/encode_video.py" "$PROJECT/scene" --ffmpeg "$FFMPEG"
python3 "$SKILL/scripts/audit_video.py" "$PROJECT/scene/preview.mp4" \
  --ffmpeg "$FFMPEG" --fps 24 --frames 192 --output "$PROJECT/scene/video-check"
```

更改帧率或帧数时同步修改审计参数。中断后的 PNG 渲染可在相同命令上加 `--resume`；它核对工程和渲染设置，并重新渲染损坏帧。不同工程或设置使用新目录。

准备阶段在 manifest 的 `data_sha256` 中记录每个实体 NPZ 和 `motion.npz`；构建器先核对完整文件清单及哈希，再使用原来的路径检查结果。手工改变坐标、替换运动文件或删掉一条哈希记录都会失败，应从已审查输入重新准备。新 `.blend` 还保存准备 manifest 的指纹，重开检查时核对它，防止把旧工程与另一套准备记录混在一起。哈希只证明这些文件之间的一致性，不证明来源或机制正确。

旧缓存没有这组哈希时，先重新运行 `prepare_interaction.py` 再生成新工程。对已有旧 `.blend`，检查器仍可比较几何和姿态，但会将 `prepared_file_integrity.verified` 明确标为 false；不能把这种旧记录称为文件指纹验证通过。

字体、渲染帧、准备后的网格及视频都留在项目外部目录，不包含在 skill 分发包中。结构准备使用装有 Gemmi、NumPy、SciPy 和 scikit-image 的 Python 环境。

当前示例：Barnase–Barstar（1BRS）、MDM2–p53 片段（1YCR）、Trypsin–BPTI（2PTC），以及用于检验动作限制的 Ras–RasGAP（1WQ1）。可在项目 JSON 中添加新蛋白，而不修改路径生成器；仍必须核对正确装配、作者链、构建体、状态、缺失区域和配体。

### 明确保留哪一份配体

`ligands` 是必需保留的 CCD 组分列表：在已选 `anchor`/`moving` 作者链中保留所有匹配的非聚合物，并随该链的刚体伙伴运动。任何请求的组分未找到、或匹配残基没有可保留重原子时，准备阶段在输出前失败；不能静默省略后继续。空列表表示本配置没有选取非聚合物，不证明来源不含这些组分。

需要指定某一份配体、插入码，或配体使用独立作者链时，改用 `ligand_instances`。例如在 2PTC 的项目配置中保留其他字段，替换配体选择为：

```json
{
  "ligands": [],
  "ligand_instances": [
    {"chain": "E", "residue_number": 462, "insertion_code": "", "component": "CA", "owner_role": "anchor"}
  ],
  "ligand_selection_note": "2PTC 作者链 E:462 的钙离子保留在 trypsin 固定伙伴中；依据条目配体及原子记录。"
}
```

`chain`、整数 `residue_number`、可选单字符 `insertion_code` 和 `component` 合起来确定源残基；`owner_role` 只能是 `anchor` 或 `moving`。所属伙伴由项目依据结构关系确定，不能由链名或距离自动断言。不得同时使用非空 `ligands` 与 `ligand_instances`，同一源残基不得重复归属。输出 manifest 记录作者残基身份与保留的原子名称；原坐标仍与蛋白共用同一原点。保存工程检查逐帧核对每个实体的所属伙伴及整体矩阵，防止配体脱离父对象或被另行平移。

这个入口仅选择模型 0 的 `NonPolymer` 实体，并保留非氢、正占据、空白/A altloc 原子；水、分支糖链、其他 altloc、缺失原子及共价连接需专门处理。保留配体坐标不等于确认其生理作用，也不自动生成配体成键、解离或催化。`prepare_molecular_asset.py` 的静态聚合物入口仍不包含非聚合物；不要混用两个入口的配置。

- [1BRS](https://www.rcsb.org/structure/1BRS) 原研究描述 C40A/C82A barstar 构建体；记录其身份，不把瞬间教学往返当作真实结合/解离速率。
- [1YCR](https://www.rcsb.org/structure/1YCR) 是 MDM2 结构域与 p53 肽段，不能凭这一结合镜头声称展示了完整泛素化和降解。
- [2PTC](https://www.rcsb.org/structure/2PTC) 中 BPTI 是抑制剂，不能作为普通蛋白水解底物剪断。
- [1WQ1](https://www.rcsb.org/structure/1WQ1) 含 GDP/铝氟过渡态类似物。该构象的刚体释放路径可能无法通过界面检查；保留失败报告并改用有来源的构象动作，不能删除配体或放宽阈值掩盖穿插。

## 动作和检查

`interaction_motion.py` 在给定实验结合末态周围搜索接近方向，先完成主要转向，再接触；释放使用另一侧弧线。它不预测结合位点或亲和力。`loop` 连续往返；`teaching` 分为初态、接近、结合停留、释放和末态。模式改变后重新检查，不能沿用另一套时序的通过记录。

粗采样只用于排序候选，最后检查每一个输出姿态。默认 0.20 nm 是原子中心严重重叠诊断阈值，不是范德华接触模型。源末态例外只允许相同原子对最多额外加深 0.01 nm；其他位置出现同样数量的重叠仍然失败。构建器拒绝未经修正的失败路径。深度包埋界面可能需要柔性环/侧链/结构域变化或新的结构状态，不能保证所有蛋白都能用刚体弧线完成。

同一移动伙伴及其配体共用刚体变换，身份保持不变。保存第 N+1 帧的闭环姿态，导出 1…N 帧，检查虚拟末态与第 1 帧相同。源表面是坐标衍生显示包络，不是实验密度；未解析区域不自行补全。

工程审计检查保存后的完整源几何、身份、每帧姿态和相机边界。相机按整个运动的模型包围盒取景，不能只按结合态取景。完整视频解码和时间戳检查之后，仍需查看编码后的接近、接触、释放和首尾帧；低于严重重叠阈值的检查不证明分子动力学合理。

相互作用特写提供 `cinematic`、`teaching`、`clay`、`illustrative`、`fluorescence` 外观；这些外观与 `loop`/`teaching` 动画模式分别选择。荧光外观是发光标记与形状提示，不是光学成像模拟。原子描边外观要求 Eevee。两伙伴使用可区分的颜色及对应图例；文字独立合成，避免描边、颗粒和景深掩盖文字。更换颜色滤镜后重新核对图例与分子颜色。

构建器选用的合成树还导出 `Beauty`/`Labels`，可按 [合成配方](render-recipes.md#对现有工程应用) 用 `apply_render_recipe.py` 追加调色。文字经过独立支路，图例色点随分子一起处理；没有宣称导出 Normal/Depth 等 pass，所以需要这些 pass 的配方仍会拒绝。2PTC 的独立重建检查中，原五张关键帧和追加 neutral 的画面逐像素一致；soft-cinema 改变画面，独立文字缓冲逐像素一致，两种保存版本各 193 个姿态通过检查。该结果不代表全部外观组合或完整调色视频已验收。

## 继续扩展其他作用时

磷酸转移已有独立的 [PKAc/SP20 实例](phosphotransfer.md)，包含供体、产物及原子身份变化。它与局部结合入口分别实现，不把结合视频重新命名为催化视频。

[Ras/SOS 核苷酸交换及 Ras/p120GAP 水解](ras-switch.md) 提供独立与组合入口：旧 GDP 离开，另一个 GTP 进入，再由同一 GTP 的 γ-磷酸与水中氧形成 Pi、释放产物并回到 GDP 构象。组合入口核对交接坐标，重新检查所有伙伴关系，并保留旧/新核苷酸的不同身份；这是开放底物池中的状态循环。

[RNF4/Ubc13/Ube2V2 的 K63 连接](ubiquitin-transfer.md) 提供共价蛋白连接的专用入口：从已装载泛素的 E2 开始，构建净键转移；可选追加伙伴解离与完整双泛素释放。新异肽键在释放期间保留，最终同时检查原子间距与画面分离；这个实例不包含 E1 装载、泛素回收或蛋白酶体降解。

| 作用 | 需要额外实现的机制内容 |
|---|---|
| 激酶/磷酸酶 | 底物状态、供体与产物、磷酸转移/水解、释放；不能把抑制肽当正常底物 |
| GEF/GAP | GDP 释放/GTP 结合与 GTP 水解分开，保留核苷酸和金属/类似物身份 |
| 蛋白水解 | 正确底物与切割位点，完整片段及去向，不能使用抑制剂复合物冒充反应轨迹 |
| E1/E2/E3 | 泛素活化、转移、连接类型、底物识别；结合结构不自动支持整条催化链 |
| 马达与通道 | 接触/闸门状态、能量与方向、膜侧别、货物/离子连续性 |
| 聚合/拆解 | 成核与已有端延伸区分、位点占有、完整单体身份、化学老化与回收 |

这些额外作用需要专门的事件与状态实现；当前局部结合生成器没有自动完成它们。新机制完成与否以实际工程和完整视频验证记录为准。
