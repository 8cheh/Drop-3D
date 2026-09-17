# 使用说明 / Guide

这个库把**液滴图像**变成**界面性质**，并在做不到的时候**拒绝作答**。
本目录是它的说明书；计划与进度在 [`docs/plan/roadmap.md`](../plan/roadmap.md)，
为什么这么设计在 [`docs/decisions/`](../decisions/)，验证证据在
[`docs/validation/`](../validation/)。

- **[workflows.md](workflows.md)** —— 五条可跑的工作流，照着抄就能用
- 本页 —— 模块地图、该用哪个入口、被拒绝了怎么办

---

## 1. 装与跑

```bash
pip install -e ".[dev]"
pytest
```

运行时依赖**只有 `numpy` + `scipy`**。没有 OpenCV、没有 PyTorch——
每加一个依赖都是一份许可负担和编译负担，所以加之前要有理由。

```python
import drop3d
print(drop3d.__version__)
```

---

## 2. 模块地图

11 个模块，按「你要做什么」分组。括号里是行数，方便判断分量。

### 入口：从图像到轮廓

| 模块 | 干什么 |
|---|---|
| `segmentation` (279) | **灰度图进，亚像素轮廓出**。Otsu 阈值（从直方图算，不是写死）+ 测量对称轴（不是假设在画面中心）+ 亚像素边缘插值 |

### 求解器与拟合

| 模块 | 干什么 |
|---|---|
| `younglaplace` (260) | 轴对称 Young-Laplace 求解器。`Bo = 0` 时精确退化为单位球，这是锁住整个方法的恒等式 |
| `fitting` (471) | 圆拟合（顶点曲率）、椭圆拟合（足迹）、切线夹角。椭圆用**几何距离**而不是代数残差 |
| `tensiometry` (382) | 五参数最小二乘拟合、尺度标定、γ 换算、Worthington 数、形状参数 P_s |

### 四类交付物

| 模块 | 交付物 |
|---|---|
| `surface_energy` (404) | **表面自由能**：OWRK / Fowkes / Wu / 酸碱 / Zisman，**强制同时报模型间离散度** |
| `uncertainty` (535) | **不确定度**：GUM 传播、Monte Carlo、bootstrap、HAC 稳健协方差、覆盖性检验 |
| `conformal` (351) | **分布无关的不确定度**：分割共形区间，按 Worthington 数分箱，**按区间宽度拒绝** |
| `validity` (571) | **分布外拒绝**：五道独立闸门 + 命名错误码 + 计数式可靠性分级 |

### 动态

| 模块 | 干什么 |
|---|---|
| `dynamics` (476) | 滑移/滚动液滴的关系式：Furmidge（`k` 显式）、Dunlop 精确力平衡、Cox–Voinov、足迹几何 |
| `oscillation` (407) | 振荡悬滴 → 膨胀模量 E′, E″。**仪器相位滞后必须拟合** |
| `tracking` (394) | 跨帧跟踪同一液滴 + 速度与不确定度 |
| `hazards` (396) | **采集风险闸门**：把无症状的失效在实验设计阶段就摆出来 |

### 命令行

```bash
drop3d-ps                    # 形状参数 vs Bond 数扫描
drop3d-ps --target 0.15 --gamma 72 --delta-rho 998
```

回答的是**实验设计**问题：液滴要多大，轮廓里才带得动表面张力信息。

---

## 3. 该用哪个入口

| 你手上有什么 | 用这个 | 得到 |
|---|---|---|
| 一张悬滴照片 + 针头外径 | `segmentation` → `tensiometry` | γ 及其不确定度 |
| 已经是轮廓点（不用我们分割） | `tensiometry.young_laplace_fit` | 同上 |
| 多种探测液体的接触角 | `surface_energy` | 表面自由能 + **模型间离散度** |
| 一批已知真值的测量 | `conformal` | 分布无关区间 + 宽度拒绝 |
| 一帧帧的 γ(t)、A(t) | `oscillation` | E′, E″ |
| 逐帧检测到的液滴位置 | `tracking` | 速度 + 不确定度 |
| 一帧帧的两侧接触角 | `dynamics` | 保持力、动态接触角 |
| 还没做实验，想先规划 | `hazards` + `drop3d-ps` | 采样率门槛、所需液滴尺寸 |

**两条最重要的分界**：

1. **滑移/滚动液滴上不要用轴对称 YL 拟合。** 接触线不是圆、表面不是回转面，
   拟合出的「表面张力」是形状补偿参数。滑移侧分析里 **γ 是输入**，
   来自单独的悬滴测量。见 [ADR-0002](../decisions/0002-axisymmetric-vs-non-axisymmetric.md)。
2. **逐帧的东西全部吃数组，不碰图像。** `oscillation` 和 `tracking`
   的输入是逐帧测量结果，所以它们可以脱离分割算法单独验证。

---

## 4. 三条设计承诺

### 一、每个数字都带一个诚实的误差条

不是「附赠 ±」，而是**校准过的**：`uncertainty.coverage_test` 检验报出的
95% 是否真的覆盖 95%。共形区间更直接——覆盖率是它的定义本身。

**已知局限**：几何拟合的残差天然相关，朴素协方差会**低估**。
默认已改用 HAC 稳健估计，但强相关下**仍偏低**——
代码会在 `info['note']` 里自曝这一点并指明改用 bootstrap。
证据见 [`docs/validation/hac-covariance.md`](../validation/hac-covariance.md)。

### 二、拒绝作答是一等公民

`validity.assess()` 的判决是三态，不是两态：

| 判决 | 含义 |
|---|---|
| `accept` | 每道**硬**闸门都跑了并且通过了 |
| `accept_with_warning` | 可用，但某事是边缘的，或某道硬闸门**没跑成** |
| `reject` | 至少一道硬闸门失败，**不报数字**，只报原因 |

「没跑成」不允许伪装成 `accept`——否则 `accept` 的含义会退化
成「我碰巧测了的都过了」。

### 三、查不到出处就不编

代码里每个阈值都带来源标签：

- `literature` —— 文献写了这个值
- `verified` —— 我们自己数值验证过
- `engineering` —— 工程判断，**并写明依据**

有些参数**必填、无默认值**（Cox–Voinov 的 `ln(b/a)`、共形的宽度容忍度、
跟踪的 `max_speed`）——给默认值等于把我们的方便伪装成要求。

---

## 5. 被拒绝了怎么办

拒绝理由是可执行的，不是「置信度低」。`ValidityReport.report()` 会打印
**每道闸门、实测值、阈值、以及物理含义**：

```python
rep = assess(pts, res, px_size_mm=..., delta_rho=..., ...)
print(rep.report())
print(rep.error_codes)        # 机器可读，例如 ['REJECT_ILL_CONDITIONED']
print(rep.reliability_class)  # 'reject (1 hard, 0 soft)'
```

常见错误码与**该改什么**（这张表是 `validity.ERROR_CODES` 的完整内容，
不是举例）：

| 错误码 | 来自哪道闸门 | 该改什么 |
|---|---|---|
| `REJECT_TOO_FEW_POINTS` | `profile_points` | 轮廓点太少——分割或提取环节 |
| `REJECT_BAD_PROFILE` | `profile_finite` | 轮廓含 NaN/Inf |
| `REJECT_DEGENERATE_PROFILE` | `profile_extent` | 轮廓跨度太小，不是液滴 |
| `REJECT_NON_CONVERGENT` | `fit_converged` | 拟合没收敛，看 `res.error` |
| `REJECT_NO_PARAMETERS` | `parameters_present` | 拟合报 ok 但参数不全 |
| `REJECT_POOR_FIT` | `residual_rms` | 残差 RMS 太大——分割质量，或背景没处理干净 |
| `WARN_NOT_AXISYMMETRIC`（软） | `axisymmetry` | 左右拟合不一样好。**液滴在动？滑移侧不能用 YL** |
| `WARN_RESIDUAL_STRUCTURE`（软） | `residual_structure` | 残差符号有结构——模型与数据不符，即使 RMS 很小 |
| `REJECT_ILL_CONDITIONED` | `bond_range` / `shape_parameter` / `worthington` | 形状或条件数不足。**换更大的液滴**——`drop3d-ps` 能算出要多大 |
| `WARN_ILL_CONDITIONED`（软） | `neumann` | Neumann 数偏低；注意**其阈值是我们的占位值**，无文献依据 |
| `WARN_TILTED`（软） | `rotation` | 拟合出明显相机倾角 |
| `WARN_IMPLAUSIBLE_SHAPE`（软） | `radius_plausible` | 顶点半径相对轮廓跨度太大，通常不是液滴 |

**软/硬是有区分的**：硬闸门失败 → `reject`；软闸门失败 → `accept_with_warning`。
**硬闸门「没跑成」也只会降到 `accept_with_warning`**，不会变成 `accept`
——缺元数据是调用的局限，不是反对该次测量的证据，但也不能算通过。

**拒绝时不要绕过闸门。** 每个阈值都可能过严，但正确做法是
**用 `thresholds=` 参数覆盖并记录理由**，而不是删掉检查——
覆盖是显式的、可审计的，删掉不是：

```python
assess(pts, res, thresholds={'min_worthington': 0.05})   # 有意放宽，写进日志
```

---

## 6. 现在还不能做什么

**所有验证都是合成数据。** 分割在真实照片上能否工作**未经检验**：
真实光学会带来镜面反射、半透明液滴、遮挡顶点的针头。
合成图是测试夹具，不是数据的替代品。

需要真实实验数据才能推进的：阈值标定（标准液体）、
真实光学下的分割、滑移侧验证、共形标定集。
需要团队拍板的：**许可证**、动态场景形态、精度目标、数据来源。
详见 [`docs/plan/roadmap.md`](../plan/roadmap.md) 第 6 节。
