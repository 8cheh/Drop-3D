# 桌面版使用说明 / Desktop application

Drop-3D 的 Windows 图形界面。左侧侧边栏八个页面，覆盖五条工作流加实验设计。
**中英双语可切换**（右上角）。

界面风格是苹果式的毛玻璃，用 CSS 的 `backdrop-filter` 实现，理由见
[ADR-0003](../decisions/0003-desktop-gui.md)。

---

## 1. 跑起来

### 已经打好包

解压 `Drop-3D-0.1.0-win64.zip`，双击 `Drop-3D.exe`。不需要装 Python。

需要 **Microsoft Edge WebView2 运行时**：Windows 11 和多数 Windows 10 自带。
缺失时会弹出一个说明对话框，并把日志写到
`%LOCALAPPDATA%\Drop-3D\logs\desktop.log`。

### 从源码

```bash
pip install -e ".[gui]"
drop3d-gui
```

或者 `python -m drop3d_desktop`。

命令行参数：

| 参数 | 作用 |
|---|---|
| `--framed` | 用系统自带窗口边框，而不是自绘标题栏。窗口管理器不配合时的退路。 |
| `--acrylic` | 额外请求 Windows 的原生亚克力背景。**默认关闭**：实测它会让整个窗口变透明，界面在花哨的壁纸前不可读。 |
| `--debug` | 打开 WebView2 的开发者工具。 |
| `--selftest` | 不开窗口，跑一遍合成悬滴链路，把结果写到 `%LOCALAPPDATA%\Drop-3D\logs\selftest.json`。成功退出码 0，失败 1。 |

### 为什么有 `--selftest`

因为「exe 能启动」不等于「exe 能用」。冻结产物可以正常开窗、正常渲染界面，
却因为缺一个 scipy 子模块而在第一次调用求解器时失败——而 `console=False`
意味着这个失败没有任何地方可以看到。所以打好的包要能被问那个唯一重要的问题：
**你还能测量吗？**

打完包后建议跑一次：

```powershell
& "desktop\dist\Drop-3D\Drop-3D.exe" --selftest
```

---

## 2. 打包

```powershell
powershell -ExecutionPolicy Bypass -File desktop\build\Build-Windows.ps1
```

产物在 `desktop\dist\`：一个 `Drop-3D\` 文件夹和一个 zip。

脚本默认先跑完整测试套件再打包（`-SkipTests` 跳过，`-NoZip` 不出压缩包）。

### 打包时的三个坑

这三个都是踩过的，不是假想的：

1. **必须用干净的虚拟环境。** 用带 `--system-site-packages` 的环境打包，
   PyInstaller 会把整机 site-packages 全收进去：第一次是 **932 MB**，
   含 torch / OpenCV / transformers。干净环境 **126 MB**。
   脚本自己建 `.build-venv`，不共享。
2. **版本必须钉死。** `requirements-desktop.txt` 用 `numpy==2.2.6` / `scipy==1.15.2`，
   因为库的 282 项测试是在这两个版本上验证的。换成 numpy 2.5.3 / scipy 1.18.1，
   同一个合成测量的 γ 从 **75.81 变成 75.53 mN/m**（差 0.37%）。
3. **构建脚本必须是纯 ASCII。** Windows PowerShell 5.1 在 .ps1 没有 BOM 时按 ANSI 读，
   脚本里写死的中文路径会被吞掉。脚本里所有路径都从 `$PSScriptRoot` 推出来。
   同理，不要把 `$ErrorActionPreference` 设成 `Stop`：PyInstaller 往 stderr 写一条
   弃用提示就会让构建直接失败。

---

## 3. 页面

| 页面 | 做什么 | 来自哪个模块 |
|---|---|---|
| 概览 | 六条工作流入口、模块地图 | — |
| 静态悬滴 | 图像 → 轮廓 → Young-Laplace 拟合 → γ ± U，逐道闸门 | `segmentation` `tensiometry` `validity` |
| 表面自由能 | 接触角 → 五种模型 + 模型间离散度 | `surface_energy` |
| 不确定度 | GUM 贡献分解、Monte Carlo 交叉验证、共形区间 | `uncertainty` `conformal` |
| 振荡悬滴 | γ(t)、A(t) → E′、E″、相位滞后 | `oscillation` |
| 滑移 / 滚动 | 采集风险闸门 → 保持力、足迹、跨帧跟踪 | `hazards` `dynamics` `tracking` |
| 实验设计 | P_s vs Bond 数，反推所需液滴尺寸 | `tools.shape_parameter_scan` |
| 系统信息 | 版本、许可证状态、判决的含义 | — |

### 每个页面都有「生成合成数据」

因为**当前所有验证都基于合成数据**——真实图片上的分割尚未检验。
合成图是测试夹具，不是真实数据的替代品。界面里每一处合成数据都带标记。

### 拒绝作答是怎么呈现的

这是整个界面里唯一一处贯穿所有页面的设计决定：

硬闸门失败时，**头条数字的位置被「已拒绝作答」取代**，而不是把数字留着、
在旁边放一个警告。屏幕上的数字会被引用；带原因的拒绝会被修掉。
这个不对称就是这个界面存在的意义。

拒绝理由先给界面自己的语言（错误码 + 该改什么），
库自己的英文原话放在下面并标注「库原始输出」——
库的消息是写给看 traceback 的人看的，改写成中文等于替库说话。

---

## 4. 开发

```
desktop/
  src/drop3d_desktop/
    __main__.py      窗口、DPI、放置、启动失败处理
    api.py           JS↔Python 桥。纯函数，不 import webview
    serialize.py     dataclass/dict → JSON 的唯一转换点
    images.py        图像文件 → 二维灰度数组（Pillow）
    assets/          index.html / style.css / app.js / pages.js / i18n.js
  build/             PyInstaller spec 与构建脚本
```

### 三条约定

1. **`api.py` 不许 import `webview`。** 这是整套桥接能在无显示器 CI 上测试的唯一原因
   （见 `.github/workflows/ci.yml` 的 `desktop` job）。
2. **桥接方法不抛异常。** 一律返回 `{'ok': False, 'error': ...}`。
   抛到 JavaScript 是「promise 被拒 + 一堆乱码栈」，返回错误字典才能把消息原样显示。
3. **序列化只走 `serialize.py`。** 它会把 dataclass 的**属性**也带上——
   只走字段会丢掉 `ValidityReport.error_codes` 和 `DynamicsReport.verdict`，
   而界面上会表现为「这里本该有个答案」。

### 界面调试

没有自动化点击工具时，可以照 `_gui_spike/drive_gui.py` 的思路写一个驱动脚本：
用 `window.evaluate_js` 点按钮、翻页，再用 Pillow 的 `ImageGrab` 截图。
这比肉眼看快得多，也能抓到「某个页面渲染时抛异常」这类只在运行时出现的问题。

---

## 5. 测试

```bash
pytest tests/test_desktop_api.py tests/test_public_api.py
```

无窗口、不需要 pywebview。覆盖：序列化（含 NaN、numpy 类型、dataclass 属性）、
零参数方法可调用性、完整悬滴链路、标定缺失时的拒绝、相对不确定度的换算、
共形区间的宽度拒绝、跟踪速度、以及「桥接对垃圾输入不抛异常」。

`tests/test_public_api.py` 守 `__all__` 与实际可导入名的一致性——
做桌面端时发现 `footprint_from_contact_line` 被列进 `__all__` 却从未 import，
`from drop3d import *` 直接失败。见 ADR-0003。
