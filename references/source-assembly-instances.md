# 显式装配实例与多链伙伴

用于重复链、非恒等装配、多个蛋白共同构成一个运动伙伴，或配体必须随指定实例移动的项目。先按 [联网核对](online-binding-research.md) 选定适用装配、构建体和必要组分。`source_assembly.py` 在普通 Python 中准备；Blender 读取核定后的坐标和身份。

## 共同输入

保留原有 `pdb`、`sha256`、`selection_note`，增加：

- `assembly_id`：原始 mmCIF 的装配编号字符串。
- `assembly_reference.sha256`：独立下载的 `PDB-assemblyN.cif` 哈希。其 `.cif.source.json` 必须对应同一条目和装配。`complete: true` 要求覆盖参考中的全部实例；明确选择局部组分时用 `false` 并说明省略范围。
- `chains`：每个完整已建模聚合物实例分别声明 `instance_id`、作者 `chain`、mmCIF `label_chain` 和 `operator`。实例名使用字母、数字、下划线或连字符，忽略大小写后仍唯一。`description` 如提供，必须与原始 entity 描述一致；显示名用 `label`，不能冒充来源描述。
- `ligand_instances`：除上述实例字段外，指定整数 `residue_number`、可选 `insertion_code`、`component` 和 `owner_instance`。必须与其已选聚合物 owner 使用同一个源操作；不同拷贝的配体有不同实例名。保留空列表也应明确说明为何无需伴随组分。

作者链名和 label 链名常常不同。例如 3SN6 的作者 G 对应 label C。不要改写作者残基编号来区分不同拷贝；用实例身份区分它们。只保留非氢、占有率大于零且 altloc 为空/A 的原子，插入码单独保存。

`operator` 的单操作形式为源 ID，如 `"2"`；组合形式为 `"(X0)(2)"`。可从 `SourceAssembly(path, assembly_id).operations` 查看该装配实际允许的 `(label_chain, operator)`，不能给所有链套用同一张全局操作列表。操作列表、整数范围与括号组的笛卡尔积按 [wwPDB 定义](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_pdbx_struct_assembly_gen.oper_expression.html) 展开，最右侧先作用。Gemmi 0.7.5 的装配读取会忽略组合项，因此这里直接读取原始操作表。

源操作应用一次，随后所有组分接受共同的拟合与显示变换。保留源矩阵精度，不做归一化、拉伸或逐链重新居中。独立参考按全部原子身份对应，允许原子行重排；最大坐标分量差仍采用已有装配检查的 0.001 Å 门槛。该门槛核对坐标重建，不证明几何或生物学正确性。粗 CA 模型仍须标明分辨率和缺失区域。

当前支持完整聚合物，以及明确选择的 NonPolymer/Water 残基。分支糖链、共价连接重建、溶剂交换和未建模区域需要另外处理。必需组分不受支持时应保留缺口，不能删掉它使流程通过。

## 接入状态对照

在 [状态对照](source-state-comparison.md) 的每个 state 中采用上述字段。`chains` 仍提供该入口原有的描述、颜色、来源物种和位点等信息。

逐状态 CA 映射保留作者 `chain`，额外提供 `instance_id`；参考用 `alignment.reference_chain` 加 `reference_instance_id`。全部所选聚合物和配体先接受各自源装配操作，再接受该状态的同一次拟合。旧恒等装配输入仍可使用原格式；多操作输入必须显式选择实例并提供独立参考。

## 接入结合候选

在 [结合候选](monomer-binding.md) 的 `anchor` / `moving` 中分别声明多链伙伴。每个伙伴保留一个共同坐标系、一个候选刚体变换和一条运动轨迹；文件按成员分开保存，避免把不同分子连成一张人工表面。

每条距离约束增加 `anchor_instance` / `moving_instance`，其余作者残基、原子名、距离范围与依据保持原语义。重复编号不能跨链自动合并。

模板方式在 `template.anchor` / `template.moving` 中给出同样的来源及实例选择，并提供 `maps` 列表。每项声明 `instance_id`、`template_instance_id` 及原有 `residue_pairs`、替换清单、覆盖率、RMSD 门槛。全部聚合物必须一一对应；先核对各成员映射，再联合拟合，并重新检查共同变换下每个成员的 RMSD。当前模板拟合使用蛋白 CA，其他类型的配准需要适用方法。

所有聚合物与配体参与跨伙伴端态和逐姿态检查；候选没有来源接触豁免。端态通过后仍须检查接近路径，候选排名不保证可达性。失败时保留报告，检查其他有依据的候选，或重新审查构象和允许自由度；不能用删除组分或放宽门槛替代。

运行 `audit_binding_candidates.py --sources "$PROJECT/sources" --prepared "$PROJECT/candidates" --output "$PROJECT/source-check.json"` 可独立核对成员身份、共同变换、约束和分块稠密全原子距离。保存场景再用原有 `audit_interaction_scene.py` 检查完整几何、共同父对象、实例元数据和逐姿态矩阵。最终画面仍须另行检查形态、配色与可读性。

技术夹具 `test_source_assembly.py` 验证组合顺序、生成器作用范围、重复实例、插入码、伴随组分和两种入口的交接。它不是真实蛋白或新的结合证据。来源核对、候选几何、路径检查和生物机制依据应分别报告。

## 真实结构的有限迁移控制

[binding-assembly-control.json](../assets/binding-assembly-control.json) 复用 [3SN6](https://www.rcsb.org/structure/3SN6)：完整的已建模 R 链（含 T4L 融合部分）和 P0G 激动剂作为一个伙伴，A/B/G/N 四条链作为另一个共同移动的伙伴。它保留五条聚合物及全部 27 个激动剂重原子，共 10,274 个重原子。原始来源是已结合、含稳定抗体的实验构建体，不能称为两个游离伙伴的独立预测、纯人天然体系或生理招募过程。

CA 配对全部保留作者身份，三条几何约束的范围是原结构 CA 距离周围人为设置的 ±0.02 nm 控制区间，不是实测的接触容差或氢键证据。第一个候选端态通过，但接近路径仍有穿插；第二个候选通过 97 个离散姿态的严格检查。第二个候选相对模板的整体位移 RMS 约 0.199 nm，仍须保持 `inferred` 标识。有限方向搜索失败不证明不存在可达路径，离散姿态通过也不是连续时间碰撞或分子动力学证明。

```bash
python3 "$SKILL/scripts/fetch_structure.py" 3SN6 --output "$PROJECT/sources"
python3 "$SKILL/scripts/fetch_structure.py" 3SN6 --assembly 1 --output "$PROJECT/sources"
"$PYTHON" "$SKILL/scripts/prepare_binding_candidates.py" \
  --sources "$PROJECT/sources" --config "$SKILL/assets/binding-assembly-control.json" \
  --method template --candidate 2 --frames 96 --fps 24 --output "$PROJECT/candidates"
"$PYTHON" "$SKILL/scripts/audit_binding_candidates.py" \
  --sources "$PROJECT/sources" --prepared "$PROJECT/candidates" --output "$PROJECT/source-check.json"
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_interaction_scene.py" -- \
  --asset-dir "$PROJECT/candidates/interaction" --output "$PROJECT/bound-view" \
  --font-en "$FONT_EN" --style cinematic --engine cycles --width 1600 --samples 48 \
  --surface-profile "$SKILL/assets/style-structure-materials.json" \
  --representation atoms --orbit-degrees 150 --snapshot-frame 49 --render
```

省略 `--snapshot-frame` 可保存全部准备姿态的动画工程，开发时不加 `--render` 即可避免渲染全部默认关键帧。来源原子、完整伙伴、姿态和相机分别用相应检查器核对，图像只承担其实际展示范围的视觉检查。
