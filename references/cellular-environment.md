# 有景深的细胞微结构背景

用户在 2026-10-28 指定图片与视频的默认起点：用多种细胞微结构形成有远近层次的虚化背景。它适用于教学和写实表现；风格参考、纯色底与透明背景仍可按项目覆盖。

先选区室和景别，再决定环境内容。胞质可使用核糖体、适当的骨架和可溶蛋白；膜邻域需要正确朝向的膜、糖链和邻近蛋白；细胞核可使用染色质和核内复合物。不能为了丰富而把核糖体放进核内或把细胞器缩成分子大小。无定量依据的密度、位置及背景漂移记录为示意。

将主反应区、环境与文字分层。近、中、远景应有不同虚化程度、遮挡和对比；保留部分可辨认的轮廓，避免整片随机光斑。背景可以缓慢运动，但不得盖住核苷酸/离子进出、让同屏视图错时或造成文字闪烁。双视图中的局部图可能独立放大，背景不能暗示它和整体图共享精确比例尺。

分隔线采用与局部背景反差明显的颜色：暗底优先浅灰白，亮底优先深灰。最终输出约 2 px 的可读线宽为起点；检查动画全程，不能让线消失在背景结构中。此规则来自用户对双视图中间分隔线的直接反馈。

## 可运行的胞质背景例子

`prepare_cellular_environment.py` 从公开原始坐标准备三种粗略表面：人 ubiquitin [1UBQ](https://www.rcsb.org/structure/1UBQ)、人 80S 核糖体 [9SPI](https://www.rcsb.org/structure/9SPI) 和哺乳动物肌动蛋白形状代理 [8F8R](https://www.rcsb.org/structure/8F8R)。后者来源为兔肌动蛋白，不把它宣称为人细胞中测得的丝；使用完整内部单体及 D→E 变换延伸示意背景丝。实际项目需要不同区室时重新选资产，不能盲目复用这一组合。

```bash
for code in 1UBQ 9SPI 8F8R; do
  "$PYTHON" "$SKILL/scripts/fetch_structure.py" "$code" --output "$PROJECT/sources"
done
"$PYTHON" "$SKILL/scripts/prepare_cellular_environment.py" \
  --ubiquitin "$PROJECT/sources/1UBQ.cif" \
  --ribosome "$PROJECT/sources/9SPI.cif" \
  --actin "$PROJECT/sources/8F8R.cif" --output "$PROJECT/environment"
```

构建 Ras/SOS 场景时默认采用上述环境，在原有命令中加入 `--environment-dir "$PROJECT/environment"`。省略此参数保留简洁纯色版，供指定纯色或对比检查使用。其他生成器由 agent 按同一分层方式接入；不能把此入口称为已对所有场景自动生效。

`cellular_environment.py` 使用完整表面实例，保存来源、变换及示意布局。远景核糖体、中景丝状骨架与分散蛋白分别进入三个 View Layer；`Cellular environment ... defocus` 节点可分别编辑虚化宽度。当前是有明确深度层次的艺术化合成景深，不是标定的相机光圈或显微 PSF。主体、化学对象和 Typography 在背景合成后叠加，不参与背景模糊。

验收静帧时比较有/无背景的同一反应时点；动画检查入场、反应、离场和首尾，确保环境层不挡文字、不掩盖分子路径问题。保存节点、资产准备方法和参数；视频、帧、准备网格仍在 skill 包之外。
