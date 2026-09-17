# 工作流 / Workflows

五条能直接跑的工作流。**每一段都被执行过**（见 `tests/test_guide_examples.py`），
所以复制粘贴不会因为 API 改名而失效——文档里的代码和测试里的是同一份。

单位约定统一写在这里，因为这是最常见的错源：

| 量 | 单位 |
|---|---|
| 坐标、半径 | 像素 |
| 像素尺度 | mm/px |
| 密度差 `delta_rho` | kg/m³ |
| 表面张力 | mN/m |
| 长度（针头、体积） | mm、m³ |
| 角度 | 度（弧度在函数内部转换） |

> `γ ∝ 尺度²`，**精确成立**：1% 的尺度误差给 2.01% 的 γ 误差，不是 2.00%。

---

## 1. 静态悬滴：图像 → 表面张力

```python
import drop3d
from drop3d import (segment_drop, extract_profile, young_laplace_fit,
                    surface_tension, assess, AIR_DENSITY)

DELTA_RHO = 998.0 - AIR_DENSITY          # 水 vs 空气，20 °C
PX_MM = 0.02475                          # 由针头标定得到，见下

seg = segment_drop(image)                # image 是二维灰度 numpy 数组
assert seg.ok, seg.error
profile = extract_profile(image, seg)     # (2, N) 亚像素轮廓，像素坐标

res = young_laplace_fit(profile)
assert res.ok, res.error
gamma = surface_tension(DELTA_RHO, res.radius_px, PX_MM, res.bond)
```

**尺度从哪来**：这是整个测量里唯一的长度参考，而 γ 走它的平方。

```python
from drop3d import pixel_scale_from_needle
px_mm = pixel_scale_from_needle(needle_diameter_mm=1.5,
                                needle_diameter_px=60.6)
```

文献明确警告：**两轴标定不一致时结果「不可靠」**，所以要独立标定 x 与 y
并验证一致，不要假设方形像素、无畸变。

**报之前过一遍闸门**：

```python
report = assess(profile, res, px_size_mm=PX_MM, delta_rho=DELTA_RHO,
                needle_diameter_mm=1.5, volume_m3=2.0e-8)
if not report.ok:
    print(report.report())               # 每道闸门 + 实测值 + 阈值 + 物理含义
    print(report.error_codes)            # 机器可读
else:
    print(f'{gamma:.2f} mN/m   [{report.reliability_class}]')
```

---

## 2. 表面自由能：接触角 → 表面能

**表面能的数值本身是模型依赖的**，所以这个 API 强制你面对这件事。

```python
from drop3d.surface_energy import owrk, wu, compare_models, recommend

# 接触角顺序必须与液体顺序一一对应
theta = [48.0, 72.0, 65.0]               # 度
liquids = ['water', 'diiodomethane', 'ethylene_glycol']

res = owrk(theta, liquids)
print(res.total, res.dispersive, res.polar, res.rms_deg)
print(res.warnings)                       # 两种液体时无法检出模型不符

# 模型间差多少，比任何一个数字都重要
spread = compare_models(theta, liquids)
```

三条必须一起报的东西：

1. **模型间差异**在文献里是 **60–130%**，不是一个可以忽略的修正。
2. **探测液体参数本身分歧很大**——二碘甲烷的色散分量横跨 **44.1–50.8 mN/m**，
   而它是几乎所有拟合锚定 `γ_s^d` 的那个。**不要跨来源混用参数集。**
3. **扩散压常被忽略**。Young 方程里是 `γ_sv`，不是 `γ_s`；
   忽略它**单方向低估**表面能。对聚合物可忽略，对金属/氧化物/玻璃不可忽略。

---

## 3. 不确定度：两条路，一条假设少

### GUM 传播（快，但假设线性化成立）

```python
from drop3d import Budget, surface_tension_uncertainty

unc = surface_tension_uncertainty(delta_rho=DELTA_RHO,
                                  radius_px=res.radius_px,
                                  px_size_mm=PX_MM, bond=res.bond,
                                  u_bond=0.005 * res.bond,
                                  u_px_size_mm=0.005 * PX_MM)
print(f'{unc.value:.2f} +/- {unc.std:.2f} mN/m')
print(unc.contributions)                  # 哪一项主导——这决定你该改什么
```

自由度为 1–4 时要报 k > 2（Welch–Satterthwaite：3 种液体配 2 个参数时
ν_eff 只有 1–4，t 因子是 4.30/3.18，直接报 ±2SE 会差 1.5–2 倍）。

### 共形区间（不假设模型正确，但要求可交换）

```python
from drop3d.conformal import (calibrate, predict_interval, assess_width,
                              exchangeability_check, required_calibration_size)

# scores = 标定集上的 |测量值 − 真值|，单位同 γ
scores = [...]                            # 已知真值的测量残差
assert len(scores) >= required_calibration_size(0.05)   # 95% 需要 >= 19 个

cal = calibrate(scores, alpha=0.05, wo=worthington_numbers)  # 按 Wo 分箱
print(cal.half_width, cal.level_achieved)  # 注意是 achieved，不是名义 0.95

iv = predict_interval(cal, gamma, wo=wo)
decision = assess_width(iv, max_half_width=1.0)   # 容忍度必填——这是应用决策
if not decision['reportable']:
    print(decision['reason'])              # 区间是「正确但没用」，仍然拒绝
```

**分箱不是优化，是必需**：实测边缘区间总体覆盖率 0.950，
但低 Wo 箱只有 **0.690**、高 Wo 箱过度覆盖到 0.997——
而低 Wo 正是拟合已经退化的区间。

**喂新数据前先验可交换性**，因为区间不会自己报警：

```python
chk = exchangeability_check(cal, new_scores)
if not chk['exchangeable']:
    print(chk['note'])                     # 覆盖率不再有保证
```

---

## 4. 振荡悬滴：γ(t) → 膨胀模量

逐帧的 YL 拟合**不变**（每帧仍是轴对称悬滴），新增的只是时域分析。
`t` 秒、`area` 任意一致单位（只有比值进入模量）、`gamma` mN/m。

```python
from drop3d.oscillation import fit_oscillation

res = fit_oscillation(t, area, gamma,
                      frequency_hz=1.0,        # 仪器给，别让它估
                      fit_instrument_lag=True) # 这句是重点
print(res.storage, res.loss, res.modulus, res.delta_deg)
print(res.instrument_lag_deg)
print(res.warnings)
```

**仪器相位滞后必须拟合。** 实测去掉滞后项，注入 90° 滞后时：
面积振幅从 0.50 塌到 **0.00005**，模量炸到 **1.8×10⁷**——
而 γ 振幅**保持正确**，所以这个失败极难察觉。

> 调研说「会得到负振幅」。实测**只在滞后 ≥150° 且拟合未加正性约束时**才出现。
> 我们发布的实现有正性约束，所以失败形态是「数值巨大但看起来正常」——
> **比负号更危险**，因此该开关会主动警告。

采样：≥100 帧/周期（工程规则，文献点 750 fps @1 Hz）。
实测 100 帧/周期对应 E 误差约 **0.68%**，25 帧/周期约 1.2%——规则保守。

---

## 5. 滑移/滚动液滴：先问该不该测

**第一步不是跑代码，是查采集参数合不合格。**

```python
from drop3d.hazards import assess_dynamics

rep = assess_dynamics('contact_line_dynamics', fps=100.0,
                      um_per_px=0.7, contact_line_resolved=True,
                      temperature_c=22.0, relative_humidity_pct=45.0,
                      drop_volume_ul=45.0, tilt_rate_deg_s=1.0)
print(rep.verdict)          # unusable (1 blocking)
print(rep.report())
```

它在**实验设计阶段**就拦住已知的坑：解钉 <10 ms 会被 100 fps **混叠**掉，
而症状是「前进接触角图上有空洞」——读起来像数据缺失，
不像「这个测量做不了」。

**γ 在这里是输入，不是输出**：

```python
from drop3d.dynamics import (furmidge_force, cox_voinov_angle,
                             capillary_number, footprint_from_contact_line)
from drop3d.tracking import link_detections, fit_velocity, assess_track

# 保持力——k 不是 1，所以每次都返回范围
f = furmidge_force(width_m=2e-3, gamma_mN_m=72.0,
                   theta_a_deg=100.0, theta_r_deg=80.0)
print(f.force_mN, f.force_min_mN, f.force_max_mN, f.k_used)

# 足迹——先看 aspect_is_informative
fp = footprint_from_contact_line(x, y)
if fp['aspect_is_informative']:
    print(fp['aspect'])
else:
    print(fp['warnings'])   # 跨度不足，L/W 没被确定

# 跨帧跟踪
tracks = link_detections(times, detections,
                         max_speed=2.0,       # 两个都必填
                         position_sigma=0.05)
v = fit_velocity(tracks[0])
print(v['speed'], v['std_speed'], v['n_missing_frames'])
print(assess_track(tracks[0], v)['verdict'])
```

**用 Dunlop 精确式可以完全绕开 `k`**：给定实测的 `θ(φ)`，
`dunlop_bo_sin_alpha(fourier_c1(phi_deg, cos_theta))` 直接给出重力驱动项，
不需要任何几何前因子。

> **网格分辨率会进入结果。** `fourier_c1` 用梯形法，而实测的 `θ(φ)` 可能有
> 接近间断的地方，收敛会变慢。实测：2° 网格往返误差 < 1e-4，
> 0.1° 网格 < 1e-6——**加密必须让它变好**，如果不变好说明采样有问题。
> 所以报出 C₁ 时要连同接触线的采样间隔一起报。

---

## 附：几个容易踩的点

| 坑 | 事实 |
|---|---|
| `Bo` 有两个不同定义 | 求解器用的是 `ΔρgR₀²/γ`（顶点曲率半径）；坡度 Bond 数用体积或等效球半径，**两者差约 2.6 倍**，必须具名传入 |
| 共形区间「95%」 | 实际是 `⌈(n+1)(1−α)⌉/(n+1)`，**只会偏宽**；比名义水平低报了自己 |
| 拟合出的 ℓ 不是 1 | Furmidge 的几何前因子物理容许区间是 **[0.5, 0.884]**；`k ≥ 1` 通常意味着黏性阻力污染了测量 |
| 足迹略椭圆 ≠ 表面不均一 | 均匀表面上 `L/W` 就有 **1.011–1.097** |
| 滚落角不是材料常数 | 高度依赖液滴尺寸与倾角速率，**不报这两个就没法比较** |
| 更大不一定更好 | 形状参数 P_s 在 **Bo ≈ 0.45** 达峰后**下降**，不是单调的 |
