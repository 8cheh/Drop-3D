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

```python
from drop3d import young_laplace_fit, surface_tension, pixel_scale_from_needle

res = young_laplace_fit(profile_xy)          # (2, N) image pixels, y downward
px_mm = pixel_scale_from_needle(1.5, 300.0)  # 1.5 mm needle spanning 300 px
gamma = surface_tension(997.0 - 1.184, res.radius_px, px_mm, res.bond)
```

---

## 目录结构

```
.github/       CI 工作流、issue / PR 模板
docs/
  plan/        研发计划、里程碑
  log/         研发日志（每次工作一条）
  research/    调研报告 + 结论索引
  decisions/   技术决策记录 (ADR)
src/drop3d/    核心代码
tests/         测试
pyproject.toml 包配置（依赖：numpy + scipy，OpenCV 为可选 extra）
```

命令行工具尚未提供。在它存在之前，这里不列它——
文档承诺一个不存在的入口，比不写更糟。

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
