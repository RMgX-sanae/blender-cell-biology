# Blender Cell Biology

用于指导 agent 制作可编辑的 Blender 细胞生物学模型与机制动画。以公开结构和文献为依据，支持教学示意、细胞环境、参考图片风格复现及可选合成器。

[**打开微形互动展厅**](https://rmgx-sanae.github.io/blender-cell-biology-gallery/)，旋转和拆解血蓝蛋白的十个亚基，或移开超净工作台的前窗，查看内部部件。

[![微形互动展厅与血蓝蛋白模型](https://rmgx-sanae.github.io/blender-cell-biology-gallery/assets/exhibition.jpg)](https://rmgx-sanae.github.io/blender-cell-biology-gallery/)

从 [SKILL.md](SKILL.md) 开始，按交付类型在[入口与验证范围](references/workflow-entrypoints.md)选用已核对的流程。具体构建命令见 [精细制作工作流](references/production-workflow.md)。源码当前以 Blender 5.x 为主要运行目标，已在 Blender 5.2.3 LTS 上运行代表性案例。

## 使用

将本仓库克隆到 agent 可发现的 skill 目录，或把仓库路径交给 agent。建模、下载和渲染输出另放到项目目录。资产准备按对应流程建立独立 Python 环境，主要依赖为 NumPy、SciPy、scikit-image 和 Gemmi；编码视频需要 FFmpeg。字体由运行环境提供：中文默认宋体，英文默认 Times New Roman。

## 能力与边界

- 细胞、线粒体、膜受体和分子结构场景，按尺度及区室选择环境。
- V 型 ATPase、PKA 磷酸转移、Ras/SOS 核苷酸交换、Ras/p120GAP 水解及若干核酸、骨架与运输案例。
- 分步教学和连续播放、参考风格分析、可编辑材质、灯光及合成节点。
- 来源记录、对象身份与几何检查、保存场景检查、渲染和解码检查。

每个入口的实现范围和未完成部分以对应说明为准。过程目录和分镜不等于全部机制已经实现；Ras 交换与 GAP 水解可独立生成或合并为开放底物池中的状态循环；部分完整肌动蛋白案例仍有接触诊断需要审核。教学路径、插值和示意密度不代表分子动力学或实验轨迹。

[三个独立建模工程的重建脚本](https://github.com/RMgX-sanae/blender-cell-biology/tree/main/examples/three-scale-studies)展示了红藻藻胆体装配、硅藻壳的分级孔隙和空气置换移液器拆解。每个工程包含多个可切换视图。仓库仅提供源码和来源说明，成品由 Blender 在项目目录生成。

## 仓库内容

本仓库保存脚本、参数、来源记录及说明，不保存生成的视频、渲染帧、`.blend`、准备网格或字体。公开坐标等第三方数据的来源与使用条件按各资产记录；不因本仓库公开而改变其原有条件。生成可分发 skill 包使用 `scripts/package_skill.py`，输出目录必须位于仓库外。
