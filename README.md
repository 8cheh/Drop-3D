# Drop-3D

**重建液滴的三维形状，并从中测量界面性质。**

A 3D re-build of droplets from pictures, and the interfacial properties measured from them.

---

## 这个库是做什么的 / What this repository is

围绕「从图像重建液滴并测量其界面性质」的一系列工作。当前的主力任务是 **P1**。

| 编号 | 任务 | 状态 |
|---|---|---|
| **P1** | 悬滴法表面张力、表面自由能、不确定度、分布外拒绝（面向动态/滚动液滴） | 🚧 进行中 |
| P2 | 待定 | — |

任务分工与里程碑见 [`docs/plan/roadmap.md`](docs/plan/roadmap.md)。

---

## P1：我们具体要做什么

**目标**：给定液滴图像（含液滴滚动等动态场景），可靠地给出

1. **表面张力 / 界面张力** —— 悬滴法，Young-Laplace 拟合；
2. **固体表面自由能** —— 由多种探测液体的接触角反解；
3. **不确定度** —— 每个数字都必须带一个诚实的误差条；
4. **分布外拒绝** —— 当输入超出方法适用域时，**拒绝给出数字**，而不是给一个自信的错值。

### 为什么「分布外拒绝」是核心而不是附属

现有的大量接触角/表面张力代码会在测量完全失败时仍然返回一个看起来合理的数字。
我们的立场是：**一个带错误条的正确数字，胜过一个不带错误条的错数字；而拒绝回答，
胜过两者之中任何一个错的。**

具体地，系统必须能识别并拒绝：

| 情形 | 举例 |
|---|---|
| 根本不是液滴 | 空图、纯背景、异物 |
| 分割错误 | 把背景纹理/反光当成液滴 |
| 几何假设被破坏 | 液滴非轴对称、明显倾斜、正在剧烈变形 |
| 物理上有效但条件数不足 | 液滴太接近球形，形状里不携带表面张力信息 |
| 数值失败 | 拟合不收敛，或收敛到错误的解分支 |

### 为什么强调「动态」

静态悬滴法是成熟方法。难点在于**液滴在动**——滚动、滑移、振荡。这带来两个后果：

- 一部分测量（轴对称 Young-Laplace 拟合）在动态下**假设本身失效**，必须识别出来并降级处理；
- 另一部分测量（体积、速度、前进/后退接触角、接触线形状）在动态下**反而成为主要信息**。

哪些能用、哪些不能用，是本任务要回答清楚的第一件事。结论会写进
[`docs/research/`](docs/research/)。

---

## 快速开始

```bash
pip install -e ".[dev]"
pytest
```

从一张图像到表面张力（完整链路）：

```python
from drop3d import (segment_drop, extract_profile, young_laplace_fit,
                    surface_tension, pixel_scale_from_needle, assess, AIR_DENSITY)

seg = segment_drop(image)                    # 二维灰度 numpy 数组
profile = extract_profile(image, seg)         # (2, N) 亚像素轮廓
res = young_laplace_fit(profile)
px_mm = pixel_scale_from_needle(1.5, 60.6)    # 1.5 mm 针头跨 60.6 px
gamma = surface_tension(998.0 - AIR_DENSITY, res.radius_px, px_mm, res.bond)

report = assess(profile, res, px_size_mm=px_mm,
                delta_rho=998.0 - AIR_DENSITY, needle_diameter_mm=1.5)
print(gamma, report.verdict)                  # 闸门不过就不要报这个数
```

**五条完整工作流**（静态悬滴 / 表面自由能 / 不确定度 / 振荡悬滴 / 滑移液滴），
每一段代码都被执行过：见 [`docs/guide/workflows.md`](docs/guide/workflows.md)。

命令行：

```bash
drop3d-ps                                     # 形状参数 vs Bond 数
drop3d-ps --target 0.15 --gamma 72 --delta-rho 998
```

回答的是**实验设计**问题：液滴要多大，轮廓里才带得动表面张力信息。

---

## 模块地图

| 模块 | 干什么 |
|---|---|
| `segmentation` | 图像入口：Otsu 阈值 + 测量对称轴 + 亚像素轮廓 |
| `younglaplace` | 轴对称 Young-Laplace 求解器（`Bo = 0` 精确退化为单位球） |
| `fitting` | 圆拟合、椭圆拟合（几何距离）、切线夹角 |
| `tensiometry` | 五参数拟合、尺度标定、γ 换算、Worthington 数、形状参数 |
| `surface_energy` | OWRK / Fowkes / Wu / 酸碱 / Zisman，**强制报模型间离散度** |
| `uncertainty` | GUM 传播、Monte Carlo、bootstrap、HAC 稳健协方差、覆盖性检验 |
| `conformal` | **分布无关**区间，按 Wo 分箱，**按区间宽度拒绝** |
| `validity` | 分布外拒绝：五道闸门 + 命名错误码 + 可靠性分级 |
| `dynamics` | 滑移/滚动关系式：Furmidge、Dunlop 精确式、Cox–Voinov、足迹几何 |
| `oscillation` | 振荡悬滴 → 膨胀模量 E′、E″（相位滞后必须拟合） |
| `tracking` | 跨帧跟踪同一液滴 + 速度与不确定度 |
| `hazards` | 采集风险闸门：在实验设计阶段拦住已知的坑 |

该用哪个入口、被拒绝了怎么办：见 [`docs/guide/`](docs/guide/)。

---

## 目录结构

```
.github/       CI 工作流、issue / PR 模板
docs/
  guide/       使用说明（模块地图、五条工作流、拒绝码含义）
  plan/        研发计划、里程碑
  log/         研发日志（一天一文件，当天多轮按轮次追加，带轮次索引）
  research/    调研报告 + 结论索引（含已发现的报告内部矛盾）
  decisions/   技术决策记录 (ADR)
  validation/  验证记录（每个数值主张的实测证据）
src/drop3d/    核心代码
src/drop3d/tools/  命令行工具（drop3d-ps）
tests/         测试
pyproject.toml 包配置（运行时依赖只有 numpy + scipy）
```

---

## 协作方式

见 [`CONTRIBUTING.md`](CONTRIBUTING.md)。要点：

- 主分支受保护思路：**所有改动走 PR**，不直接推 `main`；
- 每次工作结束在 `docs/log/` 追加一条日志；
- 任何影响接口或方法的决定，写一条 ADR 到 `docs/decisions/`。

---

## 许可证：**尚未确定，需要团队拍板**

本库目前**没有 LICENSE 文件**，这意味着在法律上默认「保留所有权利」。
请勿在未确认前对外分发或引用。

候选方案与理由见 [`docs/decisions/0001-license.md`](docs/decisions/0001-license.md)。
**这是当前最高优先级的待决事项之一。**

---

## 致谢与来源

本项目在早期调研中参考了 OpenDrop（GPL-3.0）及其所依据的公开文献。
**我们的代码是独立实现，不包含 OpenDrop 源码**：Young-Laplace 求解器与拟合器
来自本项目上一阶段的工作（`8cheh/AngleDrop`），数学依据为公开发表的论文。

完整引用（含 DOI）见 [`docs/research/`](docs/research/) 各报告末尾的参考文献列表，
以及 [`docs/research/README.md`](docs/research/README.md) 里汇总的结论与出处。
代码中每个引用文献的阈值/常数都在注释里标了出处与来源类别
（`literature` / `verified` / `engineering`）。
