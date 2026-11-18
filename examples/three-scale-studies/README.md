# 三种尺度的生物建模

三个独立的 Blender 工程研究三种不同的结构关系。红藻藻胆体使用同一套分子坐标观察天线装配与色素。中心硅藻用可穿透的孔、壁和筛膜表现硅质壳的层级。空气置换移液器保留同一套零件，用装配、拆解和活塞局部视图说明机械关系。

| 工程 | 在 Blender 的 Scene 菜单中选择 | 输入与表示 |
|---|---|---|
| `01_phycobilisome` | 核心显露、完整天线、色素阵列、局部蛋白对 | 红藻 *Porphyridium purpureum*，PDB 6KGX 生物学装配 1。每个亚基及每个色素可单独选择。 |
| `02_diatom` | 双壳与环带、近端壳面、七个孔室的分层视图 | 参照 *Coscinodiscus granii* 的 SEM 和 AFM 形态进行参数化建模。双壳场景第 1 帧为对合，第 60 帧为拉开。 |
| `03_pipette` | 零件拆解、装配外观、活塞和密封组件 | 参照 Gilson PIPETMAN Classic P200 的部件图。每个视图使用同一套零件数据，标注、材料和灯光可编辑。 |

每个成品目录中的 `model.blend` 可以独立打开，PNG 是从相应的保存工程重开后渲染的。预设为 1920 × 1280。选择相机视图或按 F12 可重新出图。分子模型的原始身份记录在对象属性及 `prepared/inventory.json`。详细来源与建模判断见 [SOURCES.md](SOURCES.md)。

## 从源码重建

用 Blender 5.3.0、Python 3.13.9、NumPy 2.3.5、SciPy 1.16.3、scikit-image 0.25.2 与 Gemmi 0.7.5 进行过制作。先把 `blender-cell-biology` 放在本机或克隆到一个目录，再运行本项目的 `rebuild.py`。项目使用该 skill 的角色材质、几何和分子包络函数。已核对的基础版本是 `55112b6e9e71dafb2fb5fa95bd1f243ba127d1ee`。移液器图使用用户提供的 `Songti-SC-Regular.ttf` 和 `Times-New-Roman.ttf`。提供一个单独的字体目录，文件名可使用相应字体的本地链接。移液器工程会嵌入它实际使用的两种字体。

```sh
python3 -m venv .venv
.venv/bin/python -m pip install numpy scipy scikit-image gemmi
.venv/bin/python rebuild.py all \
  --blender /path/to/blender \
  --skill /path/to/blender-cell-biology \
  --fonts /path/to/font-directory \
  --out /path/to/output
```

将 `all` 换成 `antenna`、`diatom` 或 `pipette` 可以只做一件。藻胆体首次重建会下载并校验 6KGX，也可以用 `--source-cache /path/to/structures` 指向已有的 `6KGX.cif.gz`。源码压缩包只含项目脚本与说明。建模输入和成品保存在项目侧。只做低成本检查时可选 `--width 900 --samples 10`。当前机器支持 Cycles Metal；`SHOWCASE_CPU=1` 可选择 CPU。
