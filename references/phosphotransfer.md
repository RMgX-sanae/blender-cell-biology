# 激酶磷酸转移：PKAc 与 SP20 示例

用于把“结合动作”扩展为真正发生化学状态变化的教学工程。当前入口专门对应 PKAc/SP20：供体进入、底物结合、磷酸转移、产物侧链转向、磷酸化肽释放、金属离子与 ADP 释放，再进入下一轮。它不包含 cAMP 生成或 PKA 全酶激活；不能声称已经展示完整信号通路。

## 来源与身份

- [4IAC](https://www.rcsb.org/structure/4IAC) 给出含 **AMP-PCP** 的复合物。它不是 ATP 的实验坐标。脚本以该类似物为位置参考，用 ATP Chemical Component Dictionary 的原子身份和键连接重建供体，并记录几何调整。
- [4IAD](https://www.rcsb.org/structure/4IAD) 给出 ADP、磷酸化 SP20 与镁离子的产物构象。PKAc 为小鼠来源，SP20 是由人 PKI 序列改造的 20 残基肽；不是天然全长靶蛋白。
- [1J3H](https://www.rcsb.org/structure/1J3H) 的开放构象用于指导口袋开合。构建体、修饰与缺失区域不同，匹配不到的原子位置由局部位移场推断；不补造未解析原子。
- [2015 年原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC4505467/) 提出了产物侧链转向、磷酸化肽释放、M1 离开、ADP–M2 释放的机制解释。本例用这一顺序组织教学时间，不把它当作实测单分子轨迹或普适速率。

反应中的磷原子始终保留同一 ID：ATP 的 γ-磷酰基连接到底物 Ser21 的氧。该重原子示意不包含氢、质子化、溶剂或电荷平衡模拟；不能由此解释质子转移的具体路径，也不能增加一个游离 Pi 作为磷酸化的常规产物。

## 从公开坐标重建

先设置实际可用的 `PYTHON`、`BLENDER`、`FFMPEG`、`SKILL`、`PROJECT`、`FONT_ZH` 与 `FONT_EN`。项目必须在 skill 目录外；准备和场景输出使用新目录。Python 需要 NumPy、SciPy、Gemmi、scikit-image，Blender 自带 Python 用于场景构建。当前实际运行版本为 Blender 5.2.3。

`FONT_ZH` 指向宋体常规 face，`FONT_EN` 指向 Times New Roman Regular。若只有 TTC，先用 fontTools 选择正确的常规 face 并导出到项目；不要把集合的第一个字重当作默认。某些字体没有上标加号，本例使用可读的 ASCII `Mg2+`，避免缺字方框。

```bash
for code in 4IAC 4IAD 1J3H; do
  "$PYTHON" "$SKILL/scripts/fetch_structure.py" "$code" --output "$PROJECT/sources"
done
"$PYTHON" "$SKILL/scripts/prepare_phosphotransfer.py" \
  --sources "$PROJECT/sources" --output "$PROJECT/assets"
"$PYTHON" "$SKILL/scripts/prepare_phosphotransfer_motion.py" \
  --asset-dir "$PROJECT/assets" --output "$PROJECT/motion" \
  --cycles 2 --frames-per-cycle 288 --fps 24
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_phosphotransfer_scene.py" -- \
  --motion-dir "$PROJECT/motion" --output "$PROJECT/scene" \
  --font-zh "$FONT_ZH" --font-en "$FONT_EN" \
  --engine eevee --width 1280 --samples 32 --render
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/audit_phosphotransfer_scene.py" -- \
  --output "$PROJECT/scene/scene-audit.json"
"$BLENDER" -b "$PROJECT/scene/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/render_saved.py" -- \
  --directory "$PROJECT/scene/frames" --frames all
"$PYTHON" "$SKILL/scripts/encode_video.py" "$PROJECT/scene" --ffmpeg "$FFMPEG"
"$PYTHON" "$SKILL/scripts/audit_video.py" "$PROJECT/scene/preview.mp4" \
  --ffmpeg "$FFMPEG" --fps 24 --frames 576 --output "$PROJECT/scene/video-check"
```

增加 `frames-per-cycle` 可放慢完整流程；同步修改视频帧数检查。`--engine cycles` 可以选择另一渲染引擎，仍需实际渲染小样。当前工程使用专门的双视图合成节点；不要直接用单个 Beauty 层的合成配方覆盖它，否则会丢失视图组合。需要新滤镜时把它插在视图合成之后、Typography 叠加之前，并重新检查颜色与文字。

## 动作与画面检查

1. 资产准备检查六条接近/释放路径；完整时序再检查各组分与激酶、各组分彼此之间的接触。默认 0.20 nm 为原子中心严重重叠诊断阈值。源结构例外只适用于相同原子对，并限制进一步加深，不能靠保持相同重叠总数通过。
2. 肽链过渡用键长与相邻键角的几何约束，避免笛卡尔插值压塌末端。非桥接磷酸氧的对应关系是推断的等价标记映射，记录在 motion manifest；不声称测到了它们的运动。
3. 激酶整体用于构象示意：开放参考缺失片段用连接与相邻角度约束进行几何修整，修整范围和偏移量写入来源记录。表面逐帧从全部原子位置重建，避免用局部位移拖拽旧网格而拉长三角面。整体构象插值不是可供定量解释的原子动力学轨迹；右侧反应组分才按原子与键直接显示。
4. 188 个重原子构成一组底物、ATP 和金属离子；反应后身份守恒。第二轮使用新的独立底物包，上一轮产物退出视野，不能把画面里的 ADP 无解释地变回 ATP。
5. 整体图保留完整激酶，低对比度外形上叠加底物以便观察；局部图省略整条激酶，展示完整 SP20 与反应组分。两图共享坐标和时间，但放大倍数不同。没有对单个分子做剖切。
6. 工程检查读取重新打开的 Blender 文件，检查全部帧的坐标、可见性、身份、标题和双视图同步。绝对形态键时间存在浮点误差，脚本逐帧绑定实际 key time，并单独报告坐标误差；这不改变接触诊断阈值。
7. 完整解码后检查真正编码的视频：供体是否能跟踪、磷酸是否落到正确残基、两个金属位点是否可分、进出镜头是否遮挡标签、下一轮是否连续。数值检查通过不代替这一步。

## 推广到其他蛋白时

可以复用来源校验、原子映射、供体/受体账本、路径搜索、双视图与字幕框架。换成另一种激酶、磷酸酶或 ATPase 时，必须重新建立来源构象、反应图、原子对应、金属/辅因子状态和进入/释放路径；不直接套用 SP20 的残基编号、产物顺序或口袋开合。该入口尚不是任意酶的自动催化模拟器。

视频、帧、字体、`.blend`、准备后的网格与坐标缓存仅保留在项目中。skill 分发保留这些生成器、参数与来源获取方法。
