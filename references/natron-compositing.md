# Blender 到 Natron 的实际合成

将同一套已核定画面换到外部合成器，不重建分子或改变机制。入口复用 [合成输出约定](render-recipes.md)：`Beauty` 是包含分屏、背景和构图的底图，`Labels` 是独立文字层，原 `Image` 用于比较。当前实现在线性 Rec.709 中调色底图，再用预乘 alpha 合并文字；最后交回原版本 Blender，执行一次原 AgX/显示变换。

**范围：**实际验证了 Blender 5.2.3 与 Natron 2.5.0 macOS x86_64（Apple Silicon 上使用已有 Rosetta）。这是 Grade + Over 的可编辑 `.ntp` 工程，不是所有 Natron 滤镜、外部 Z/Normal/Cryptomatte 或其他软件的兼容认证。Natron 使用自己的 OCIO 配置做线性到线性的恒等转换，不加载 Blender 5.2 的新版 AgX 配置；Natron Viewer 的外观不能代替最终 PNG。其他版本先做一个代表帧。

外部软件已有可用安装时优先复用。需要为项目准备运行环境时，从 [Natron 官方发布](https://github.com/NatronGitHub/Natron/releases/) 获取适合系统的版本，记录 URL、版本及哈希；不把软件包放入 skill，不修改系统色彩配置。软件不可用时保留 EXR 交接，明确哪些步骤没有执行。

## 1. 导出少量正式尺寸帧

`BLENDER`、`NATRON` 是实际可执行文件，`SKILL`、`PROJECT` 是绝对路径；输出目录必须不存在。下面的时点属于 PTP1B 示例，其他过程按事件和构图变化选帧。

```bash
"$BLENDER" -b "$PROJECT/source/scene.blend" --python-exit-code 1 \
  --python "$SKILL/scripts/export_composite_layers.py" -- \
  --output "$PROJECT/layers" --frames 1,194,613 --width 1280 --samples 24
```

生成三套 RGBA 32-bit ZIP EXR、原合成 PNG 和 `handoff.json`，保存实际帧号、尺寸、fps、预乘 alpha、源工程哈希及 Blender 色彩资源哈希。`Beauty`、`Labels`、`Image` 必须对应原树的最后一次 Labels-over-Beauty；未知接法、自定义 OCIO 或曲线显示变换先单独核定。导出只修改内存副本，不重存源 `.blend`，并暂停副本中已有的 File Output 节点。

## 2. 准备配方与可编辑工程

令 `NATRON_OCIO` 指向已安装 Natron 的线性 Rec.709 配置。本次 macOS 包的位置是 `Natron.app/Contents/Resources/OpenColorIO-Configs/blender/config.ocio`；不能将不同配置里同名的 `scene_linear` 当作同一色域。准备任务文件：

```bash
python3 - "$PROJECT/layers/handoff.json" "$PROJECT/natron" \
  "$NATRON_OCIO" "$PROJECT/natron-job.json" <<'PY'
import hashlib, json, sys
from pathlib import Path
handoff, output, config, destination = map(Path, sys.argv[1:])
job = {
    'handoff': str(handoff.resolve()), 'output': str(output.resolve()),
    'project_format': 'HD_720',
    'grade': {'multiply': [1.12, 1.03, 0.93], 'gamma': [0.96, 0.96, 0.96]},
    'ocio_config': str(config.resolve()),
    'ocio_sha256': hashlib.sha256(config.read_bytes()).hexdigest(),
}
destination.write_text(json.dumps(job, indent=2)+'\n')
PY

CELLBIO_COMPOSITE_JOB="$PROJECT/natron-job.json" \
  "$NATRON" -t "$SKILL/scripts/natron_composite.py" < /dev/null
```

先用 `multiply`、`gamma` 都为 `[1,1,1]` 的中性版本对照，再在新目录做选定调色。上述暖色只是美术起点；保留蛋白/原子身份颜色，查看高对比主体和背景是否仍可区分。改变背景、材质或局部深度应回到相应生成入口。

`project_format` 是 Natron 的实际格式 ID：1280×720 用 `HD_720`，1920×1080 用 `HD`。新尺寸可以省略此项，让脚本添加 `CellBiology`；若该尺寸已存在于软件预设，脚本会列出格式 ID，需选择对应项。最终像素尺寸还会独立检查，不能用格式名称替代验证。范围内 gamma 要求正数；负 RGB 的非单位 gamma 需要另行核定。

准备脚本保存 `scene.ntp`、`natron-graph.json` 后关闭工程。该版本新建 ReadOIIO 节点可能尚未初始化图层菜单，因此**必须重开工程再渲染**；不要把准备阶段的空输出当作成功。Natron 某些错误仍返回退出码 0，须继续检查工程、审计文件和实际像素。

## 3. 重开检查、合成与显示

```bash
python3 - "$PROJECT/natron-job.json" "$PROJECT/natron-audit.json" <<'PY'
import json, sys
from pathlib import Path
job = json.loads(Path(sys.argv[1]).read_text()); job['action'] = 'audit'
Path(sys.argv[2]).write_text(json.dumps(job, indent=2)+'\n')
PY

for FRAME in 1 194 613; do
  CELLBIO_COMPOSITE_JOB="$PROJECT/natron-audit.json" \
    "$NATRON" "$PROJECT/natron/scene.ntp" \
    -l "$SKILL/scripts/natron_composite.py" -w LinearOutput "$FRAME"
done

"$BLENDER" -b --factory-startup --python-exit-code 1 \
  --python "$SKILL/scripts/finish_natron_composite.py" -- \
  --job "$PROJECT/natron-job.json"
```

重开检查核对读写器、图层连接、调色参数、帧范围及色彩设置。完成脚本要求所有选定帧存在，检查全尺寸浮点像素、原图的图层重组、实际 Grade/Over 结果和文字贡献，再输出 `display_####.png` 与 `composite-check.json`。RGBA 验证容差为绝对/相对各 2e-6；显示用原 Blender 版本与相同色彩资源，不再额外套 LUT。检查文件只覆盖列出的帧，不能替代整片或运动审查。

若用户在 Natron 里编辑了工程，旧审计自然失效。先评估新接法、更新与之对应的配方和数值检查，再生成新记录；不要仅改哈希绕过错误。检查完成后按 [保留规则](review-and-retention.md) 保留当前工程、输入、结果和证据，清理可重建的调试副本。

## 实测与来源

2026-11-06：PTP1B 和 AMSH-LP 各三张 1280×720 帧，分别执行中性与暖色配方；实际像素相对预期 Grade/Over 的最大误差 1.20e-7，受保护文字贡献误差同量级。中性显示 PNG 与原图最多相差一个 8-bit 色阶，未声称 PNG 位级相同。两个场景均保留完整双视图、原标签和自然景深背景；这次未重新渲染整片，也未增加生物学阶段。

接口依据：[Natron 脚本执行与重开](https://natron.readthedocs.io/en/rb-2.5/devel/natronexecution.html)、[App 的创建与保存接口](https://natron.readthedocs.io/en/rb-2.5/devel/PythonReference/NatronEngine/App.html)、[Grade](https://natron.readthedocs.io/en/rb-2.5/plugins/net.sf.openfx.GradePlugin.html)、[Merge](https://natron.readthedocs.io/en/rb-2.5/plugins/net.sf.openfx.MergePlugin.html)。本入口用实际软件参数及像素验证这些约定，不能仅凭文档判断兼容。
