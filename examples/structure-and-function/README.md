# 血蓝蛋白、耳蜗与共聚焦显微镜

三个独立的 Blender 工程都可以切换 Scene 查看整体和局部。构建器从公开结构、人工标注或仪器光路出发。坐标来源、配色含义和教学简化见 [来源说明](SOURCES.md)。

| 项目目录 | 在画面中观察 |
|---|---|
| `01_hemocyanin` | 太平洋褶柔鱼血蓝蛋白的十个亚基构成桶壁。金色 FU-g 和珊瑚色 FU-d* 交错占据腔内。第二个镜头在原位保留六条完整的已建模亚基。 |
| `02_corti` | 小鼠耳蜗中段的毛细胞、柱细胞隧道和 Deiters 细胞为教学重建。第二个镜头放大其中同一外毛细胞的三阶静纤毛和顶端连接。第三个镜头单独显示 P21 小鼠耳蜗顶端一个真实标注视野中的 40 个毛束。 |
| `03_confocal` | 微透镜盘与针孔盘同轴。青色激发光通过透镜与针孔，样品发出的暖色荧光经针孔上行，由二向色镜折向右侧的中继光学系统和相机。 |

## 重建

使用 Blender 5.3 和具有 `numpy scipy scikit-image gemmi pillow` 的 Python。将下例中的路径换成本机路径。默认通过 HTTPS 取回原始输入并核对摘要。耳蜗项目只读取公开档案中所用标注及校准信息的两个字节范围，总量约 0.4 MB。

```bash
python rebuild.py \
  --blender /absolute/path/to/blender \
  --skill /absolute/path/to/blender-cell-biology \
  --out /absolute/path/to/new-projects
```

默认输出 1920 × 1280 静帧。快速检查可加 `--width 600 --samples 8`。已有核定来源可用 `--source-cache /absolute/path/to/sources`，其中可放 `6R83.cif.gz`、`p21-labels.bin`、`p21-calibration.xml`。程序先生成工程并关闭 Blender，再禁用自动执行脚本重新打开，检查可编辑网格并渲染所有镜头。工程默认打开整体场景，各 Scene 有自己的相机和合成节点。

| 检查文件 | 核对内容 |
|---|---|
| `source-check.json` | 按 mmCIF 的序列编号重算功能单元原子数和质心。逐层解码 PNG 标注，按原 CZI 再读一次体素尺寸。 |
| 各工程的 `saved-check.json` | 显露内部时保留整条亚基，全部静纤毛的根部落在所属细胞顶面，顶端连接的两端接触相邻高度的纤毛，显微原始标注与网格一一对应。光学工程检验全部透镜与真实孔口的配准以及二向色镜的反射方向。 |

材料为教学配色。工程是固定结构视图，时间线上未制作氧合、听觉转导或扫描成像动画。
