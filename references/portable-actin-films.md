# 独立重建 Arp2/3 与 formin 完整教学案例

`scripts/build_actin_films.py` 从公开结构重新准备资产、生成 21 个可编辑章节、渲染、编码并检查两条路线。输入是原始坐标和生成脚本，不需要旧 `.blend`、网格、视频、帧或原制作目录。两条路线仍采用 [actin-full-cycle.md](actin-full-cycle.md) 的特定教学综合情境；修改细胞、受体、亚型或机制时重新核对来源。

## 运行

准备 Python 依赖：`numpy scipy scikit-image gemmi fonttools Pillow`。Blender 使用自己的 Python。需要 Blender 5.x（当前实际测试版本为 5.2.3）和 FFmpeg。项目目录必须在 skill 外。

```bash
"$PYTHON" "$SKILL/scripts/build_actin_films.py" \
  --output "$PROJECT" --blender "$BLENDER" --ffmpeg "$FFMPEG" \
  --font-zh /installed/fonts/Songti-Regular.ttf \
  --font-en /installed/fonts/Times-New-Roman.ttf \
  --quality draft --steps prepare,build,render,encode,audit
```

中文默认宋体常规字重，英文和数字默认 Times New Roman 常规字重。macOS 可自动发现这两种系统字体，并将宋体集合中的指定字形提取到项目目录；其他平台传入实际安装路径。字体只打包进本地项目工程，不放进 skill 压缩包。

- `draft`：640×360、12 fps、8 samples，用于完整时序与技术检查。
- `review`：960×540、24 fps、16 samples，为默认档位。
- `final`：1600×900、24 fps、32 samples，仍需按实际画面增加采样。
- `--narration none` 为默认，保留可读主标签与独立 SRT，不插入自动旁白。`macos` 使用本机 `say` 与可选 `--voice`；`provided` 读取项目 `audio/<shot-id>.wav`。音频必须实际听检。无需旁白不代表删除机制阶段。
- `--engine eevee|cycles` 与 `--compositor clean|reference|monochrome` 独立选择。新滤镜叠加见 [render-recipes.md](render-recipes.md)。

可按 `prepare`、`build`、`render`、`encode`、`audit` 分步运行，后续调用从项目配置读取路径。`--shots a03,f04` 可定位章节；只有所有章节完成后才能合成整条路线。`--preview` 在构建时渲染三个关键时点。

## 可复现性与续跑

来源坐标逐份保存官方元数据和 SHA256；资产准备缓存检查输入及输出哈希。准备步骤中途失败的生成目录会改名保留，完整缓存才复用。来源或预处理变化时使用新项目目录。章节渲染按已保存工程的哈希分开缓存，完成的重复状态只渲染一次；重新执行 `render` 可补齐未生成的帧。编码用图像序列保持固定 fps，不用插帧填补缺失动作。

本入口重新计算链配准、G/F 示意形变、Rac/WRC 与 DID/DAD 接近方向和 FH2 动作参数。不能把某个旧缓存 JSON 当作无法解释的固定真值。详细命令、返回码、日志和软件信息保存在项目 `validation/`。

## 检查边界

`audit` 检查保存工程的完整基底网格、字体、全部时点变换、部分关键章节的原子中心严重重叠，以及实际 MP4 的帧数、时间基准、连续时间戳和图像尺寸。自动生成逐章视频抽帧图；agent 必须实际查看。接触诊断有警告时，应定位修正或明确记录，不能把“脚本退出成功”称为无穿插。

目前原案例中 FH2 中间动作、部分回收与端部装配仍存在局部接触诊断警告；模板保留完整链和来源构象，没有删掉相撞原子来使检查通过。它是带明确简化的机制教学工程，不是分子动力学，也不宣称原子级连续催化轨迹。两条全周期计划的 `design_only` 描述计划本身；成片是否完成，以当前项目实际视频、检查与视觉记录为准。

交付留在项目目录。skill 只分发生成代码、参数、来源规则和说明。
