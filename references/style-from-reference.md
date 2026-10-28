# 根据图片识别风格并在 Blender 复现

用于用户提供图片、截图、视频关键帧或公开作品链接并说“按这个感觉/风格做”。这是 skill 的主要工作模式之一，不需要先给用户一份技术参数问卷。先实际查看图片，不能只凭文件名、网页标题或 image-search 自动描述来推断。缺图时请求上传或有效链接；有本地文件时用图像查看工具。

## 1. 拆出能观察到的视觉特征

按下面维度做简短记录。每项写“观察”“实现假设”“待确认”，不要从单张图声称知道原始节点图、渲染引擎或准确 roughness。

| 维度 | 观察什么 | 可实施的 Blender 控制 |
|---|---|---|
| 形体语言 | 写实表面/简化块体/原子球/线框、细节频率、曲面起伏 | 表现层级、结构数据、Subdivision/Bevel、微细 bump；不改生物实体尺度 |
| 配色 | 背景、主体、强调色、饱和度、语义颜色 | Principled Base Color、色板映射；保留同一实体在所有镜头的一致颜色 |
| 材质 | 哑光/塑料/金属/蜡感、透射、薄膜感、局部光泽 | roughness、IOR、transmission、SSS、coat；逐项试渲染，不把未知属性说成测量值 |
| 光照 | 主光方向、阴影软硬、边缘光、填充、环境亮度 | Area lights 尺寸/位置/能量、world、分层照明 |
| 镜头 | 正交/透视线索、视角、构图占比、景深、景别 | orthographic scale / focal length、相机距离和 DOF；先锁定构图再调材质 |
| 边缘/体积 | 黑轮廓、接触阴影、雾、边缘透明度 | Line Art/Freestyle 或合成边缘、AO/真正阴影、体积；不把过度 AO 变成假膜缝 |
| 后期 | 对比、曝光、色调、glow、颗粒、背景融合 | Compositor 或外部合成器；使用与来源分离的美术效果层 |
| 标注 | 标签位置、字号、线条、是否有剖面/图例 | 独立标注层、文字面向相机；参考 UI 导航和水印通常不属于模型风格 |

语义颜色默认保持一致。若颜色特意表示反应状态变化，应在场景说明中规定映射并在所有循环中一致使用；这不同于无意的色闪或身份混淆。

如果参考是荧光/EM 图，先辨认成像方式及其局限，再决定“模拟该外观”还是“从该数据重建”。艺术化显微风格成片标为 simulated/schematic，不能伪装成新增实验图。

## 2. 转为可实施的风格档案

需要插画描边、透明外壳、丝带、球棍、体积、SEM/DIC/相差等扩展时，读取 [style-library.md](style-library.md)。其中 18 类起点有实际生成器和 live 节点；[来源目录](style-reference-catalog.json) 保留原图链接，图片本身不分发。`analyze_style_image.py` 可测色板/亮度/边缘，仍由 agent 实际看图并确定语义特征，不能把统计结果当原始材质参数或自动识别结论。

保存 `reference URL/path + 允许使用的范围 + 观察 + 假设 + 参数版本`。精细制作优先复制 [style-molecular-reference.json](../assets/style-molecular-reference.json)、[style-studio.json](../assets/style-studio.json) 或 [style-organelle.json](../assets/style-organelle.json)：

```bash
"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/build_molecular_scene.py" -- --scene membrane \
  --asset-dir /absolute/project/beta2 --style-profile /absolute/project/style-profile.json \
  --output /absolute/project/style-test --render
```

由 agent 实际看图后填写档案；脚本本身不自动识别图片。精细档案控制 `representation`（surface/spacefill）、语义色板、roughness/bump、灯位/尺寸/强度/颜色、camera、view_transform/exposure 和 compositor。物种、结构、膜拓扑和动画仍由场景脚本决定。完整单体分层、连续外壳的教学开口、透明材质或轮廓线需在项目脚本实施，不能只改色板就宣布匹配。

新档案颜色输入为 **sRGB 十六进制**，`detail_scene.py` 在设置材质/灯光前转成线性值。旧 `reference-style.json` 仅给旧白模 `build_scene.py` 使用，它采用线性 RGB 数组，不能混用。色调映射也影响成片；保存 view transform。用户参考图中的高光色不能直接当材质基础色。

镜头景深以场景单位工作。分子场景的美术 f-stop 仅用于视觉匹配，不能声称恢复了原参考的物理镜头。先确认主结构全部在画内，再设置景深。

## 3. 迭代，而不是猜一次就结束

1. 锁定同一构图、分辨率和几何，先做 512–960 px 的静态测试。
2. 比较主体尺寸、背景亮度、配色、轮廓软硬和阴影。每次优先修正两三项差异，避免同时更改全部参数。
3. 若有两种合理解释，给两个小样，例如“较哑光”和“更通透”。询问用户偏好可与其他独立建模并行。
4. 核对膜/腔/受体尾部是否仍可读。若参考透明效果遮蔽了关键结构，提出完整单体分层或双版本，不切开分子；连续细胞/细胞器外壳的教学开口另行注明。
5. 动画再检查色闪、噪点漂移、景深遮挡、发光溢出与文字闪烁；单张图匹配不代表整段动画匹配。
6. 保存最终风格档案、节点、灯光和相机，下一镜头复用；用户不满意时记录具体差异并修改。

## 实查案例：Allen Visual Guide

参考：[官方教程](https://www.allencell.org/visual-guide-tutorial.html)的 [公开截图](https://www.allencell.org/uploads/8/1/9/9/81996008/step8-1_orig.png)。此次实际查看该图：中央使用浅背景、红色线粒体、淡紫色核、细胞背景弱化，结构有柔和高光和阴影。截图中央存在播放器遮挡，不能据此补出被遮住的细节。

由此提炼“浅背景 + 单一红色强调 + 低饱和辅助结构 + 柔和表面”的风格；不复制其网站 UI、圆形导航或水印。`reference-style.json` 是这一视觉提炼的参数起点，**不是 Allen 原始材质参数，也不是 hiPSC 的结构重建**。Allen 原图仅作链接参考，不随 skill 分发；具体再利用遵守其条款。

## 实查案例：WEHI 分子动画静帧

实际查看 [WEHI p53 静帧](https://www.wehi.edu.au/wp-content/uploads/2023/11/p53-animation-still.jpg) 和 [DNA replisome 静帧](https://www.wehi.edu.au/wp-content/uploads/2023/06/dna-replisome-bjork.jpeg)，来源页面为 [WEHI-TV](https://www.wehi.edu.au/wehi-tv/)。前者可观察到暗蓝环境、金色主蛋白、洋红/米色核酸、青色远景和细密分子轮廓；后者的环境更拥挤、暖色更强。不能从静帧推断原节点图或精确参数。

`style-molecular-reference.json` 保存前者的观察与实施假设。用户反馈“结构更清楚了，但需要更丰富的环境层次”后，膜示例加入有真实来源的大分子形状、距离/透视/景深层次和暗部填充，并将整齐脂质阵列改为有最小距离的面内采样。它模仿视觉特征，不把 p53 的生物身份移植给受体；参考原图不随 skill 分发。

第二轮反馈认为原剖口不自然，实际对照后选择了第二张侧视剖面。当前膜示例默认 `--membrane-view side`；这个选择适用于本例，不要求所有细胞器都侧视。实心补面试验显得像板材，已从生成器移除。保留失败小样作为本项目记录，不将其推荐为风格模板。
