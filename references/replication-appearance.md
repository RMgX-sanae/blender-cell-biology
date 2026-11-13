# 复制叉的参考外观迁移

用于已有 1BNA 核苷酸原型与 7PFO 完整 CMG 参考层的局部延伸示意。先按 [制作流程](production-workflow.md#复制连续与分段延伸) 准备输入、生成并检查工程。复制 [style-replication-reference.json](../assets/style-replication-reference.json) 到项目后按参考图修改参数。

| 设置 | 实际作用 |
|---|---|
| `surface_roles.protein` / `palette.heads` | 11 个完整已建模 CMG 亚基；保留共同移动与源表面 |
| `surface_roles.nucleic_acid` | 完整核苷酸单元及主链路径连接表示 |
| `surface_roles.primer_symbol` / `palette.accent` | RNA 引物的三单位符号，不表示核糖和碱基的原子化学 |
| `palette.dna2/protein/dna1` | 亲代 DNA、前导新链、后随新链；`protein` 是沿用的色板键名 |
| `palette.ink` | 宋体中文、Times New Roman 英文及数字的保护文字 |
| `lighting/world/exposure/view_transform/look` | 实际光源、世界填充及色彩显示 |
| `camera.scale_multiplier/center_offset_nm` | 正交缩放及相机 right/up 平移；不改变复制叉 |
| `context_blur_pixels_1280` | 按输出宽度缩放的 CMG 背景层虚化，0 为清晰 |
| `render_recipe` | 共用合成配方；背景、主体和文字分别处理 |

新外观使用 `--font-latin` 指定 Times New Roman；`--font` 指定宋体并选择中文，省略它使用英文。字体随项目打包，不分发进 skill。方向箭头和图例色条采用对应链的无光照语义色，避免相机前方的标记因场景光源而变暗。来源说明与简化边界放在 manifest 和文档，画面不加微小页脚或文字描边。

完整核苷酸共享原型网格，各链通过 **OBJECT 材质槽** 保留自己的颜色。必须读取 `obj.material_slots` 的有效材质；只检查 `obj.data.materials` 会把前导链、后随链误认成亲代链。不要为改色复制或拉伸每个核苷酸。

翻译、转录与复制共用 `layered_molecular_appearance.py`，各入口分别保留源几何与机制规则。新外观不接受旧 `compositor`、`representation` 或任意镜头方向，也不与 `--compositor` 同用。旧无角色调用保留兼容。

`audit_replication.py` 检查核苷酸完整性、亲代身份、互补配对、增长顺序、引物先决条件和 CMG 深度。分层工程按各对象实际所在的渲染层读取求值状态，不能依赖被排除对象的陈旧位置。新工程还调用 `audit_replication_appearance.py` 检查有效材质槽、颜色、光源、镜头、无光照箭头、背景路由、字体、文字边界和片段极性标签的出现时间。`test_replication_appearance.py` 对保存工程执行有针对性的内存破坏测试，并重新打开原文件恢复；不保存测试修改。

对照旧工程的完整运动时，`hide_viewport` 对象即使位于当前层，也可能在视口依赖图中保留旧矩阵。先核对未经修改的对象及几何数据动画通道（包括可见性、关键点、插值和控制柄），再在不保存的内存副本中暂时关闭视口隐藏，让全部对象参与求值；保留渲染可见性和变换通道，检查整数帧及必要的子帧。不要把隐藏对象的缓存位置差异当成运动变化，也不要直接跳过这些对象或放宽误差阈值。

最后查看正式尺寸的起始、中间和末尾关键帧，判断双链与新链是否可分、引物符号是否清楚、箭头在移动背景上是否可读。数值通过不能替代实际画面；配方能够运行不代表所有配方都已完成审美验收。

范围仍是构造的 72 bp 局部延伸窗口：1BNA 提供局部原子形状，CMG 整组后移作为相关结构背景，二者不是实验复合物界面。没有新增聚合酶/引物酶、能量化学、校对、引物去除、连接、起始或终止。背景虚化属于艺术表现，不是实验光学数据。
