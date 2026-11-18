# 剪接体、有孔虫与超薄切片机

这里有三个独立的 Blender 工程。每个工程从整体转向一个可核对的局部。输出目录中的 `model.blend` 可直接打开，在顶栏选择 Scene 切换视图。保存的工程可以单独渲染。

| 工程 | 可以看什么 | 依据与取舍 |
|---|---|---|
| 酿酒酵母剪接体 | 整个复合体、五条已解析 RNA 与三个分支定位因子、A70 与 G1 的连接 | 使用带 Prp18/Slu7 的 7B9V C 复合体。从 mmCIF 保留全部正占有率主构象重原子及身份。蛋白表面是显示包络。局部保留有序底物片段和两个镁离子。缺失的 RNA 跨段保留空缺。 |
| 有孔虫 *Globigerinella adamsi* | 完整壳体、外壁孔隙、来自同一标本的中央 CT 层 | GE_GAD_2 原始扫描的体素边长 1.6827 µm。完整壳与局部共享网格，来源灰度切片按体素中心配准。象牙色和侧光用于读出壳壁。分割阈值影响细孔和薄壁。 |
| 超薄切片机工作头 | 同一个样品夹头、树脂块、钻石、刀槽和水面切片带 | 按公开实物照片制作可编辑的机械教学模型。夹头、刀口、水面之间的方向有照片依据，尺寸与薄膜反光属表现选择。镜头覆盖工作头。 |

## 从源码重建

需要 Blender 5.x、Python 3.11 或更新版本、NumPy、SciPy、scikit-image、Gemmi、Pillow，以及 Times New Roman 常规字体。只有分支键标签使用字体。可直接使用仓库中的本目录，或使用独立源码包。仓库内复用现有的 `studio.py` 和 `views.py`，独立包包含这两个小文件。

```bash
python3 rebuild.py --out /absolute/new-project \
  --skill /absolute/blender-cell-biology \
  --blender /absolute/path/to/blender \
  --font-en /absolute/Times-New-Roman.ttf
```

第一次运行下载 7B9V 和约 434 MB 的公开 CT 体积。PANGAEA 如需从磁带调出文件，会先给出等待响应。已有原始数据可加 `--source-cache /absolute/sources`，目录中使用 `7B9V.cif.gz` 和 `GE_GAD_2.raw` 两个文件名。准备程序检查来源哈希。程序将角度、照明、材质与各自的相机及合成器保存在三个工程中，禁用自动脚本后重开并渲染。默认 1920 × 1280，快速重建可加 `--width 600 --samples 6`。

`source-check.json` 逐原子比对 7B9V、核对分支端点，并逐像素检查原始 CT 层。三个 `saved-check.json` 核对重开工程中的网格、字体或影像以及切片与水面的相对位置。个别静态视图只描述对应状态。深色分子背景取自工作室布光，CT 数据视图单独使用灰度。

[公开资料与建模方法](SOURCES.md)
