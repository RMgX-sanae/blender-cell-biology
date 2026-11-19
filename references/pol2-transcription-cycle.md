# 从启动子跟到两段 RNA 的去向

适用于哺乳动物中采用切割和 poly(A) 加尾的蛋白编码基因。观众在总览里跟随同一条珊瑚色 RNA，在特写里看它的 3′ 链端连续加入两个核苷酸。剪切后保留的上游产物继续加尾。金色下游片段仍连着 Pol II，新露出的 5′ 端成为 XRN2 的入口。

| 主时间轴（含首尾，24 fps） | 镜头与来源 |
|---|---|
| 1–144 | 人源闭合及初始转录复合物按 RPB1 核心配准，先展示因子和完整复合物，再隐藏挡住转录泡的完整亚基。5IY8 的坐标里只有两条 DNA。可见的短 RNA 主链及六个碱基位置借用配准后的 5FLM 链 P，只代表教学位置。[5IY6](https://www.rcsb.org/structure/5IY6) · [5IY8](https://www.rcsb.org/structure/5IY8) · [人源起始复合物研究](https://doi.org/10.1038/nature17970) · [TFIIH 与打开启动子](https://pmc.ncbi.nlm.nih.gov/articles/PMC33249/) |
| 145–216 | 切到同一教学身份的基因总览。启动子标记留在原处，聚合酶带着转录泡前进，RNA 延长。蛋白表面复用下一镜头的 12 条源链。 |
| 217–768 | 原有 23 秒两轮延伸原速播放。第一轮加入 C，第二轮加入 U，两份 PPi 分别离开。牛源 Pol II 与合成核酸支架提供局部坐标。[5FLM](https://www.rcsb.org/structure/5FLM) · [两轮延伸构建和证据](production-workflow.md) |
| 769–1032 | 跳到基因的 3′ 加工区。CPSF73 催化段标示剪切，poly(A) 从上游 RNA 的新 3′ 端长出。XRN2 的薄荷色环沿下游 RNA 从 5′ 向 3′ 前进。随后 Pol II 与 DNA 分开，CTD 磷酸标志离开，交代 FCP1 和聚合酶再利用。[CPSF73 催化段](https://www.rcsb.org/structure/2I7V) · [切割与加尾](https://pubmed.ncbi.nlm.nih.gov/35177537/) · [XRN2 与终止](https://pmc.ncbi.nlm.nih.gov/articles/PMC5830926/) · [终止中的协同调控](https://pmc.ncbi.nlm.nih.gov/articles/PMC6938672/) · [FCP1 与再利用](https://pubmed.ncbi.nlm.nih.gov/10385623/) |

基因总览压缩了距离和时间。末段的剪切位点属于教学基因，和 5FLM 支架上的 U21 分处两个尺度。XRN2 的环表示活性位置，CPSF73 的靠近与聚合酶释放是教学编排。此路线覆盖一种基础转录途径；剪接、加帽、输出以及其他终止途径由项目另行取舍。结构校验值、能量来源与产物去向见[来源配方](../assets/pol2-transcription-cycle.json)。

## 重建可编辑工程

令 `SKILL` 指向解包后的 skill，`PROJECT` 指向项目目录，`BLENDER` 指向 Blender 可执行文件。`PY` 使用安装了 NumPy、SciPy、scikit-image、Gemmi 的项目 Python。先沿用[两轮核苷酸加入的工程](production-workflow.md)，保留 `transcription_two_cycles/scene.blend` 和同目录的 `manifest.json`。构建命令可省略 `--render`。

```bash
for id in 5IY6 5IY8 5FLM 2I7V; do
  "$PY" "$SKILL/scripts/fetch_structure.py" "$id" --format cif --output "$PROJECT/pol2-sources"
done
"$PY" "$SKILL/scripts/prepare_pol2_cycle.py" \
  --sources "$PROJECT/pol2-sources" --output "$PROJECT/pol2-prepared"
"$BLENDER" --disable-autoexec -b "$PROJECT/transcription_two_cycles/scene.blend" \
  --python-exit-code 1 --python "$SKILL/scripts/build_pol2_cycle.py" -- \
  --prepared "$PROJECT/pol2-prepared" --output "$PROJECT/pol2-cycle.blend"
```

主场景的序列编辑器可以直接播放和编辑四段镜头。切换到各源场景可改几何、灯光、材质、相机和关键帧。构建器沿用延伸工程的标签语言和字体，旧工程可传入本机的 `--font-zh` 与 `--font-latin`。查看起始 RNA 和剪切产物时，可用 `render_saved.py --frames 127,899` 并指定输出目录。

要把已有特写接进另一条机制，先定一个能跟踪的实体，记录进入、离开特写时的位置、方向和携带物。在各源场景保留该实体的标识与颜色。把有证据的两端交给 `scene_chapters.compose(name, chapters, reference)`，每段提供标题、场景和包含首尾的源帧区间。各镜头使用一致的帧率和画幅。若跨标本或跳过生物学阶段，在镜头交界说明实际对应关系。检查首次改变归属的物质即可判断该衔接是否成立。
