# ADR-0003：桌面端的技术栈与依赖

- 状态：**已接受**
- 日期：2026-09-18
- 相关：[ADR-0001 许可证](0001-license.md)、[ADR-0002 轴对称 vs 非轴对称](0002-axisymmetric-vs-non-axisymmetric.md)

---

## 背景

库的数值能力已经具备（`docs/plan/roadmap.md` 第 4 节），但只有 Python API 和一条命令行。
要把它交给不做数值的人试用——尤其是要让实验组在没有真实数据的情况下先跑通流程——
需要一个能在 Windows 上双击运行的界面。

这个决定必须回答三个问题：

1. 用什么画界面；
2. 为此引入哪些依赖，以及它们如何不污染核心库；
3. 冻结成可执行文件时，哪些约束会变成正确性问题。

---

## 决定

### 1. 界面用 pywebview + Edge WebView2，界面本身用原生 HTML/CSS/JS

界面进程是一个 WebView2 控件，Python 侧通过一个进程内对象暴露方法。

**毛玻璃是 CSS 的 `backdrop-filter`，不是操作系统的合成效果。** 界面里有一层
渐变壁纸，玻璃面板是半透明白色叠在它上面。壁纸不是装饰：`backdrop-filter`
只有在背后有结构时才是可见的，所以壁纸是整个效果的前提。

### 2. 依赖放进可选 extra `[gui]`，核心库的依赖不变

```toml
gui = ["pywebview>=5.0", "pillow>=10.0"]
```

核心运行时的依赖仍然只有 `numpy` + `scipy`。只想要数值的人不会因此下载一套 GUI。

`Pillow` 是**开发依赖之外**唯一新增的东西，理由是：
`drop3d.segment_drop` 接收的是二维灰度数组，库特意不做图像解码
（见 `pyproject.toml` 里关于 OpenCV 的注释）。相机拍出来的是 PNG，
总得有人把它变成数组，而在桌面端这个人是 `drop3d_desktop/images.py`。

许可证按 `CONTRIBUTING.md` 的要求逐个打开仓库核对：
pywebview 为 BSD-3-Clause，Pillow 为 HPND/MIT-CMU 型，pythonnet 为 MIT，
三者都是宽松许可、纯 wheel、无需编译。

### 3. 冻结用 PyInstaller，产物是文件夹 + zip

不用单文件 exe：单文件每次启动都要把 numpy/scipy 解压到临时目录，实测启动慢一个量级。
产物约 **126 MB / 395 个文件**。

---

## 被否决的方案

| 方案 | 否决理由 |
|---|---|
| **PySide6 / Qt** | 毛玻璃没有 CSS 的 `backdrop-filter` 可用，要靠 DWM 的 ctypes 手工调；再增加 200 MB 量级的体积，且要引入第二个 GUI 框架。 |
| **Electron / Tauri** | 需要 Node 工具链。本仓库的构建前提是「一个 Python 环境加一个 PyInstaller」，引入 Node 会让 `Build-Windows.ps1` 变成两套工具链。 |
| **Tkinter / CustomTkinter** | 没有可用的模糊，做不出要求的视觉语言。 |
| **本地 HTTP 服务 + 浏览器** | 会打开一个网络端口。进程内桥接没有任何网络面，这在只做内部测试、许可证尚未决定的阶段是明显更安全的默认。 |

---

## 后果

### 好的

- 界面与库彻底解耦：`api.py` 不 import `webview`，所以整套桥接可以在无显示器的 CI 上测试。
- 桥接的所有返回值都是纯 JSON：库的公开 API 返回 dataclass 和 dict，`serialize.py` 是唯一的转换点。
- 界面文件是数据，可以直接改，不需要重新打包 Python。

### 代价与约束

- **冻结必须用干净的虚拟环境。** 第一次打包用的是带 `--system-site-packages` 的 venv，
  PyInstaller 把整机 site-packages 都收了进去：**932 MB**，包含 torch、OpenCV、transformers
  ——三者本应用都不 import。干净环境是 126 MB。见 `desktop/build/Build-Windows.ps1`。
- **依赖版本必须钉死。** 库的 282 项测试是在 numpy 2.2.6 / scipy 1.15.2 上验证的；
  用 numpy 2.5.3 / scipy 1.18.1 冻结同一个测量，γ 从 75.81 变成 75.53 mN/m，
  差 0.37%。会悄悄改变产品数字的构建不算可复现，所以
  `desktop/requirements-desktop.txt` 用 `==` 而不是范围。
- **WebView2 运行时是前提。** Windows 11 和多数 Windows 10 自带；
  缺失时会弹原生对话框并写入日志，而不是静默退出。
- **`__all__` 的完整性变成了硬要求。** 见下。

### 顺带修掉的一个真 bug

做桌面端时发现 `drop3d.__all__` 里列着 `footprint_from_contact_line`，
但 `src/drop3d/__init__.py` 从来没有 import 它。后果是
`from drop3d import *` 直接抛 `AttributeError`。函数本身没问题，
从 `drop3d.dynamics` 单独 import 也能用——这正是能躲过整个绿色测试套件的那类 bug。
已修复，并加了 `tests/test_public_api.py` 守住 `__all__` 与包属性的一致性。

---

## 备选路径（如果以后要换）

界面与库之间只有一个 `api.py`。换掉前端（Qt、终端 TUI、Web 服务）
不需要动 `drop3d` 一行代码，`tests/test_desktop_api.py` 也继续有效。
这是把桥接设计成纯函数的直接收益。
