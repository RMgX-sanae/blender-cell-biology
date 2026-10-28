# 肌动蛋白与微管：可执行模板

入口是 `build_polymer_scene.py` 与 [polymer-templates](../assets/polymer-templates/actin-polymerization.json) 中的 JSON。使用公开结构提取的完整分子表面、共享网格和可编辑关键帧。默认 24 fps，普通示例 12 秒，微管缩短示例 14 秒、恢复示例 24 秒；这都是讲解时长，不是生物学速率。Blender 5.x 原生 `bpy`，无需付费模型库或专用插件。

这些是局部动作模块。用户要求肌动蛋白从上游到下游的全过程时，转到 [actin-full-cycle.md](actin-full-cycle.md)：分别实现 Rac1–WAVE–Arp2/3 和 RhoA–mDia1 路线，加上化学状态、调控及回收。旧 12 秒延伸视频不能作为全过程验收证据。

## 模板选择

| JSON 名称（位于 `assets/polymer-templates/`） | 可直接构建的行为 | 边界 |
|---|---|---|
| `actin-polymerization.json` | 倒钩端连续加入 3 个完整亚基 | 从预成丝开始，不含成核、ATP 水解或 G→F 构象转换 |
| `actin-depolymerization.json` | 尖端依次释放 3 个完整亚基 | 是选定的净解聚示例，不表示只有尖端能解聚 |
| `actin-treadmilling.json` | 倒钩端加入、尖端释放，核心亚基保持原位 | 不依靠整根丝平移；不计算临界浓度或实际交换速率 |
| `actin-severing.json` | 亚基接口处分离两个完整片段 | 不切开单体；没有自动补入 cofilin/gelsolin 机制 |
| `actin-annealing.json` | 已配准极性的完整片段退火 | 有序相遇示意，不表示任意碰撞都能结合 |
| `actin-bundling.json` | 三条完整丝转向、靠拢、平行成束 | 保持丝长；未添加 fascin 等交联因子的结合过程 |
| `arp23-branch-growth.json` | 完整母丝/Arp2/3 接头不动，已有子丝继续延伸 | 7 个 Arp2/3 亚基与原始关系保留；不包含 NPF 激活和初始成核 |
| `actin-capping.json` | 完整 CapZ 异二聚体到达来源结合终态 | 兔 actin / 人 CapZ；绑定构象的接近示意，不含全部构象变化 |
| `microtubule-growth.json` | 三原丝局部窗口内逐次加入 αβ 二聚体 | 其他原丝整体省略；不是孤立的三原丝“微管” |
| `microtubule-catastrophe.json` | 末端弯离并释放完整 αβ 二聚体，观察区缩短 | 灾变后行为示意，不计算灾变概率、GTP 水解或裂解能量 |
| `microtubule-rescue.json` | 先缩短，再由具有新 ID 的二聚体恢复生长 | 释放和加入不是同一批对象瞬移；恢复时机人为编排 |

拓展到 formin、cofilin、FtsZ、septin、中间丝和被覆装配，使用 [assembly-storyboards.md](assembly-storyboards.md)。不能把一种聚合物的模板仅换名称就用于另一种体系。

## 准备一次资产，复用于多个场景

先执行 [公开资产选择](public-assets-first.md)。这批模板是以下四个已核对条目的示例，不是任意 PDB 的自动解释器。`fetch_structure.py` 保留原始坐标与 metadata；重复运行时先复用已有文件，不覆盖。

```bash
SKILL='/actual/path/to/blender-cell-biology'
BLENDER='/actual/path/to/blender'
PYTHON='/actual/python/with/numpy/scipy/scikit-image/gemmi'

python3 "$SKILL/scripts/search_public_assets.py" 'free barbed end' --title \
  --output /absolute/project/research/actin-candidates.json

python3 "$SKILL/scripts/fetch_structure.py" 8F8R --output /absolute/project/sources
python3 "$SKILL/scripts/fetch_structure.py" 8F8Q --output /absolute/project/sources
python3 "$SKILL/scripts/fetch_structure.py" 7AQK --output /absolute/project/sources
python3 "$SKILL/scripts/fetch_structure.py" 3J6F --output /absolute/project/sources

"$PYTHON" "$SKILL/scripts/prepare_filament_assets.py" \
  --source-dir /absolute/project/sources --output /absolute/project/polymer-assets
```

输出 `actin/`、`cap/`、`branch/`、`tubulin/` 的原子/表面文件与来源 manifest。1 Blender unit = 1 nm。单条肌动蛋白表面只创建一次，所有重复单体共享 mesh；材料槽按 object 区分身份。螺旋变换由 8F8R 相邻内部链 D→E 配准得到，本轮升高约 2.781 nm、Cα 配准 RMSD 约 0.201 Å；并非凭视觉猜测间距。

## 构建与风格选择

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_polymer_scene.py" -- \
  --template "$SKILL/assets/polymer-templates/actin-polymerization.json" \
  --asset-dir /absolute/project/polymer-assets \
  --profile "$SKILL/assets/style-molecular-reference.json" \
  --engine cycles --compositor reference \
  --font /absolute/path/to/chinese-font.ttf \
  --output /absolute/project/actin --render
```

中文标签需要可用中文字库；若环境没有，指定已有合适字体或把项目模板标签改成目标语言，不下载来源不明的字体。字体仅打包进本项目，不随 skill 再分发。

复制 JSON 后可改字幕、背景开关、fps、时长、1–3 轮交换及肌动蛋白初始 6–18 个亚基。更多亚基/更多循环需同步调整镜头、对象池、时间线并检查，不能仅增加数字。微管生长、缩短、恢复的时长下限分别为 `2 + 4.8 × cycles`、`3.6 + 4.8 × cycles`、`3.6 + 9.6 × cycles` 秒。局部原丝依次加入或离开，以保持接触清楚；这个先后顺序为教学编排，不是实测顺序。

改用 `style-studio.json`、`--engine eevee --compositor clean` 可制作清晰教学版本。`clean/reference/monochrome` 三套合成分支随工程保存，材质/光照与引擎独立选择。参考图的具体模仿仍走 [style-from-reference.md](style-from-reference.md)；选一个预设不等于自动复现任意参考。外部合成可沿用 `render_saved.py --exr` 的场景线性交接。

## 重渲染与视频

```bash
"$BLENDER" -b /absolute/project/actin/scene.blend --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory /absolute/project/actin/frames --frames all

python3 "$SKILL/scripts/encode_video.py" /absolute/project/actin \
  --ffmpeg /actual/path/to/ffmpeg
```

`scene.blend`、`manifest.json` 与重建参数一同保留。模板运动已经烘焙为对象/rig 关键帧，没有外部 driver 或逐帧 handler；用户可以在 Blender 中继续编辑。参考构象、生成的重复单元、路径和画面时长的证据级别分别记录。

## 变更后检查实际保存工程

```bash
"$BLENDER" -b /absolute/project/actin/scene.blend --python-exit-code 1 \
  --python "$SKILL/scripts/audit_polymer_scene.py" -- \
  --output /absolute/project/checks/scene.json \
  --poses /absolute/project/checks/poses.npz

"$PYTHON" "$SKILL/scripts/audit_polymer_contacts.py" \
  --asset-dir /absolute/project/polymer-assets \
  --poses /absolute/project/checks/poses.npz \
  --output /absolute/project/checks/contacts.json
```

第一个检查保存后的所有帧、亚基 ID、刚体尺度/手性、事件端点、字幕边界和合成分支。第二个查重原子严重穿插，默认阈值 0.12 nm；不是范德华能量检查，不覆盖氢、配体、被省略组分或完整物理机制。通过这些检查后仍要查看实际画面，尤其近接触、分支、释放与前后两轮衔接。

`assets/polymer-motion.json` 是针对这些来源坐标编排的接近/离开路径。离体阶段先转向，接触阶段采用局部路径；它是几何排障结果，不是 MD 轨迹或真实反应通路。更换结构、单体构象、装配方向或相邻对象时必须重新检查和调整，不能直接沿用原路径宣称仍无穿插。

## 生物学与表示限制

- 8F8R 提供 ADP F-actin 形状，入场单体使用同一完整形状作代理。未展示的 G-actin 构象与核苷酸状态不能靠颜色补成“真实机制”。
- 7AQK 是 9 Å 细胞内分支结构模型；其表面是粗分辨率骨架派生包络。母丝与 Arp2/3 接头不被重新随意摆放，新增子丝部分才是重复与运动示意。
- 3J6F 为 GDP 晶格，来源对应 14 原丝参考。模板保留三条完整原丝的局部窗口，整体省略其余原丝；不会把蛋白单体切开或把窗口误标为完整 13 原丝微管。
- 生长/缩短、踏车、切断、退火和成束是不同事件。网格数量、内部亚基位置、端点与片段归属分别检查，不能用统一的“整体缩放”动画替代。
