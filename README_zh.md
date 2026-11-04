# Dy、Nd 水解：box-PBC 配位数、三步表观活化能与 Shiery 2021 对比

[English](README.md) | **中文**

- 在线网页：英文首页 https://jiaaowang-ut.github.io/dy-nd-box-cn/ ，中文页面 https://jiaaowang-ut.github.io/dy-nd-box-cn/index_zh.html
- 独立的 pKa₁/pKa₂/pKa₃ 对比表（Shiery 2021 vs 本工作，含反应方程式与公式）：中文 https://jiaaowang-ut.github.io/dy-nd-box-cn/pka_table_zh.html ，英文 https://jiaaowang-ut.github.io/dy-nd-box-cn/pka_table.html ；由 `python build_pka_table.py` 从 `comparison/compare_pKa_stepwise.csv` 生成
- 打包下载（GitHub release）：https://github.com/JiaaoWANG-ut/dy-nd-box-cn/releases/latest
- 本地：打开 `index.html`（英文）或 `index_zh.html`（中文），或解压 `dy_nd_box_cn_data.zip` 后打开 `dy_nd_box_cn/index.html`。页面右上角的 English / 中文 按钮可切换语言。所有链接都是相对路径，可离线浏览。

网页中每个表、每张图下面都有对应原始数据的下载链接，第 8 节列出全部文件。网页数字全部来自本目录的 CSV/JSON，页面只做取范围、作差之类的简单汇总。

## 内容

1. 任务 (a)，REE 中心：REE–O(水)、REE–O(OH)、REE–O(总)。截断取总 REE–O g(r) 第一壳层之后的极小。
2. 任务 (b)，OH 中心：O–REE、O–H(水)、O–X(总)。分严格截断和壳层截断两种处理，见下文。
3. 三步表观活化能：先求 committed 跃迁速率，再得 Eyring ΔG‡（300 K）和 Arrhenius Ea（溶液，300–500 K）。
4. 从占据数得到逐步 pKa 等效值和 ΔG，与 Shiery, Cooper, Cantu, *Inorg. Chem.* 2021, 60, 10257 对比：Table 1、Table 2、SI Table S1/S2、SI 优化结构、SI 第一步 PMF。

只用三维周期 box（box_pbc）轨迹。体系为 Ln(OH)3 + 水，石墨体系另含一层石墨。

- 石墨：DPA4C-Plus、UMA-s-1.2.1；300 K；4 种 OH 起始构型（3near / 2near1far / 1near2far / 3far）× 50 条 × 150 ps。
- 溶液：DPA4-Pro、DPA4C-Plus、MACE-POLAR-1-M、UMA-s-1.2.1；300/350/400/450/500 K；每个温度 1 条 × 150 ps。

## 运行顺序

在本目录依次运行。数字参数 8 是并行进程数。

```
module load python/3.12-26.1.0
python kinetics.py extract 8    # 原始轨迹 -> kinetics_cache/（20 fs 时间序列，约 372 MB）
python kinetics.py screen       # 坍塌筛查 -> collapse_screen.csv
python analyze.py               # 任务 (a)(b)：g(r)、截断、CN、图
python kinetics.py analyze 8    # 速率、ΔG‡、Arrhenius、占据数、pKa 等效值 -> kinetics/
python kinetics.py figures      # 动力学图 -> kinetics/*.svg
python literature.py            # 论文表格、SI 结构与 PMF -> literature/（会自动下载 SI PDF）
python compare.py 8             # 与论文对比 -> comparison/
python build_page.py            # 网页 index.html（英文）、index_zh.html（中文）+ dy_nd_box_cn_data.zip，并检查链接和未翻译的说明
```

`analyze.py` 的选项：

- `--reuse-rdf`：复用 `graphite_rdf_accumulated.npz`，跳过石墨 RDF 累加。
- `--rdf-only`：只输出 g(r) 和截断。
- `--reuse-cn`：复用已有的逐轨迹 CN 表。

原始数据（路径相对于 `genesis/`）：

- 石墨：`Graphite-dy-nd-50-denser-dumping/`
- 溶液轨迹：
  - `md/ln-150ps-denser-dumping/`：DPA4C-Plus、UMA，每 2 fs 一帧；
  - `md/ln-150ps-every-step-dumping/`：DPA4-Pro、MACE，每 1 fs 一帧。
- 溶液 CN 和 g(r)：`md/dense-gr-cn/gr_cn.py` 已从上面的溶液轨迹算好（`cn_table.csv`、`results/*_gr.npz`、`results/*.json`），`analyze.py` 直接读取。

## 方法

### 坍塌筛查（`kinetics.py screen` → `collapse_screen.csv`）

部分 UMA 轨迹会坍塌：势能骤降，H/O 被吸到 Ln 上。判据有两条，满足任一即可：

- 势能：每 1 ps 的块平均比前 5 块的最低值低 25 eV（石墨）或 15 eV（溶液）；起点在 −10 eV 处细化。
- 几何：Ln–O < 1.85 Å 或 Ln–H < 1.90 Å，连续 10 个采样点。

轨迹截到坍塌起点前 1 ps。剩余 ≥ 20 ps 的记为“截断”，继续使用；否则记为“排除”。

- 共排除 91 条。
- 石墨 UMA：Dy 排除 79 条、截断 116 条；Nd 排除 7 条、截断 116 条。
- 溶液：UMA Dy 的 5 条全部排除，UMA Nd 500 K 截断。
- DPA4-Pro、DPA4C-Plus、MACE 都没有坍塌。

石墨的 CN/g(r) 和全部动力学只用可用时间段。溶液 CN/g(r) 取自 `md/dense-gr-cn`，它对坍塌轨迹不给 CN，所以 UMA Nd 500 K 溶液只有动力学结果（用前 66 ps）。

### 配位数（`analyze.py`）

- **O 物种：** 每帧把每个 H 按最小镜像距离归属到最近的 O。有 2 个 H 的 O 是水，1 个是 OH⁻，0 个或 ≥ 3 个是“其它”（O²⁻/H₃O⁺）。质子转移后 OH 和水的身份会随之交换。
- **REE–O(其它)：** 石墨中为 0；溶液未坍塌轨迹中 |总 − 水 − OH| ≤ 0.002。
- **(a) REE 中心：** r_c 取总 REE–O g(r) 第一壳层之后的极小。g(r) 网格 0.02 Å，5 点滑动平均。
  - 做法：先找主峰（石墨在 1.6–3.4 Å、溶液在 1.8–3.0 Å 内的最高点），再取主峰之后 1.6 Å（溶液 1.5 Å）内的最低点。
  - 这样可以跳过 OH（约 2.25 Å）与水（约 2.5 Å）两个子峰之间的小凹陷。
  - 水、OH、总三项都用这一个 r_c 计数。
- **(b) OH 中心：** 中心是每帧所有 OH 的 O。X 是除中心 O 和它自己的 H 之外的所有原子：REE、水的 O/H、其它 OH 的 O/H、“其它”O/H，石墨体系还包括 C。
- **采样与误差：**
  - 石墨：t ≥ 10 ps 起每 1 ps 一帧。误差是轨迹间的标准误：REE 中心的量按帧数加权，OH 中心的量按 OH 中心数加权。`ALL` 行合并 4 种起始构型。
  - 溶液：`md/dense-gr-cn` 的逐帧计数（t ≥ 10 ps，每 2 fs），误差是 5 段 block 标准误。例外是严格截断下的 O–REE 和 O–H(水)：它们由部分 g(r) 积分得到，没有误差棒；与直接计数最多相差 0.05。

### OH 溶剂化壳层不清晰时的处理

O–X(总) g(r) 由几个重叠的峰组成：氢键水 H（1.6–1.8 Å；DPA4C-Plus 还有 1.23–1.25 Å 的共享质子）、REE（2.2–2.3 Å）、水 O（2.4–2.7 Å）。所以“总 g(r) 的第一极小”不唯一。处理如下：

1. **严格截断：** O–X(总) 第一个真实峰（g > 0.5，且之后下降 ≥ 10%）后的第一个极小。它总在 O–REE 峰之前，所以 O–REE ≈ 0。
2. **壳层截断：** 从 O–REE 峰位置减 0.1 Å 出发（O–REE 峰 g ≤ 0.5 时改用 2.25 Å），沿 O–X(总) 先爬到局部极大，再下降到下一个极小（不超过 3.0 Å），这样 REE 也被计入。
3. **清晰度：** 峰到谷降幅 1 − g_min/g_peak ≥ 10% 记为清晰（`*_min_clear`）。只有石墨 DPA4C-Plus Nd 的壳层极小不清晰。
4. **辅助检查：**
   - 余量：壳层内 X − REE − H(水) 应 ≈ 0。
   - 对称参考：REE–O(OH)/3 不依赖 OH 中心截断，可以和壳层 O–REE 对照。

溶液的严格截断和壳层截断由 `md/dense-gr-cn` 按相同规则给出；清晰度由 `analyze.py` 用同样的 10% 标准补算。

### 动力学（`kinetics.py analyze`）

- **状态：** n(t) = Ln 周围 r_c 内 OH⁻ 的 O 的个数，上限 3，每 20 fs 一帧。r_c 取 Ln–O(OH) 径向数密度第一峰之后的极小（3.03–3.47 Å，见 `kinetics_cutoffs.csv`）。
- **步骤对应：**
  - MD 正向 n−1 → n（Ln 得到 OH⁻）对应论文第 n 步酸解离 Ln(OH)ₙ₋₁ + H₂O → Ln(OH)ₙ + H⁺。
  - MD 逆向 n → n−1（Ln 失去 OH⁻）是它的逆过程。
- **committed 转变：** 宽 W = 0.5 ps 的居中众数滤波，再把短于 W/2 的片段吸收掉。速率从 t ≥ 10 ps 开始统计。
- **速率与能垒：**
  - k = N / T_origin，ΔG‡ = RT ln(k_B T / h k)。
  - 95% CI 取 bootstrap（2000 次）与 Garwood 泊松区间的包络。石墨按轨迹重采样；溶液每个温度只有 1 条轨迹，按 10 ps 块重采样。
  - N = 0 时用泊松上限 3.689 / T_origin 给 ΔG‡ 的下限。
- **Arrhenius/Eyring：** 按 N 加权拟合。拟合质量分级：
  - 好：≥ 4 个温度、≥ 20 个事件，且每个温度 ≥ 3 个事件；
  - 尚可：≥ 3 个温度且 ≥ 8 个事件；
  - 差：其余情况；
  - 无：有事件的温度少于 2 个。
- **敏感性：** W = 0.2 / 1.0 / 2.0 ps；r_c = 2.90 / 3.30 Å；起点 t ≥ 1 / 50 / 90 ps。见 `kinetics_rates.csv` 的 `scenario` 列。

### 占据数与 pKa 等效值

- 占据数 Pₙ 从 t ≥ 10 ps（溶液）或 t ≥ 50 ps（石墨）开始统计。
- 解离常数：K_d = (Pₙ₋₁/Pₙ)·(4 − n)·c₁。c₁ 是一个分子在 box 中的浓度：溶液 0.4224 M，石墨 0.3594 M。
- pKaₙ(等效) = pK_w + log₁₀ K_d，取 pK_w = 14，按理想溶液处理，没有活度校正。
- 状态分类（`compare_pKa_stepwise.csv` 的 `status` 列，括号里是网页上的显示）：
  - `measured`（加粗数值）：两个状态的 committed 停留都 ≥ 1 ps；
  - `transient`（“瞬态”）：两个状态都出现过，但其中一个的 committed 停留 < 1 ps，只有量级意义；
  - `upper_bound` / `lower_bound`（“≤ / ≥ … 界限”）：只出现过一个状态。假设没出现的那个状态的占据数不超过出现的那个，取 Pₙ₋₁ = Pₙ 处的 pKa 作为界限；
  - `undetermined`（“无法确定”）：两个状态都没出现。
- 石墨 DPA4C-Plus 的占据数对起始构型有记忆，到约 90 ps 才收敛，见 `start90` 设定。

### 与论文对比（`literature.py`、`compare.py`）

- **论文方法**（取自 Shiery 博士论文第 5 章，即这篇文章）：
  - pKa₁：约束 AIMD 得到 PMF，取其终点。
    - 设置：Bluemoon，CP2K/PBE，64 个水 + 3 个 Cl⁻。
    - 9 个 r(O–H) 约束点，每点 15 ps NVT，约束力取 5–15 ps 的平均。
  - pKa₂、pKa₃：M06 团簇（第一层水和 OH 显式）+ CPCM 隐式溶剂（ORCA）的反应自由能，以 pKa₁ 为基准取相对值。
    - 团簇配位数固定：La–Sm 为 9，Eu–Lu 为 8，即 Nd 9、Dy 8。
  - 论文只给 298.15 K 的热力学，没有能垒或活化能。第一步 PMF 单调上升，没有过渡态。
- **Table 1/2：** 转录自博士论文 Table 5-1/5-2，与文章 Table 1/2 相同。用 Table 2 重算 Table 1 的 MAD 来核对：38/39 项一致（±0.05），只有 Yb 第 2 步不一致；另外 6 项没有实验值。
- **SI：** Table S1/S2、60 个 M06 优化结构，以及图 S1–S15 的 PMF（从 PDF 矢量路径读取）。
- **单位换算：** ΔG = RT ln10 · pKa。
- **第 2、3 步实验值的读法：**
  - Nd 的 Kragten 1984（Table 2 ref 99）为 8.1/16.2/24.3 = 8.1 × n，看起来是累积常数 −log *βₙ，而不是逐步常数。
  - 第 2、3 步的其它实验值在数值上也更像累积值。
  - 按哪种读法会决定结论（预测的 ⟨n⟩ ≈ 1 或 3）。`compare_pKa_cumulative.csv` 给出按累积常数读的对比。

## 输出文件

- 本目录：
  - 网页：`index.html`（英文首页）、`index_zh.html`（中文）；说明文档：`README.md`（英文）、`README_zh.md`（中文）
  - CN 汇总：`cn_summary_300K.csv`
  - 石墨：`graphite_cn_summary_box_pbc.csv`、`graphite_cn_by_replica_box_pbc.csv`、`graphite_cutoffs_box_pbc.csv`、`graphite_gr_box_pbc.csv`、`graphite_rdf_accumulated.npz/.json`
  - 溶液：`solution_cn_box_pbc.csv`、`solution_cutoffs_box_pbc_300K.csv`、`solution_gr_box_pbc_300K.csv`
  - 坍塌筛查：`collapse_screen.csv`
  - 图：`*_box_pbc_300K.svg`
- `kinetics/`：
  - 速率与能垒：`kinetics_rates.csv`、`kinetics_arrhenius.csv`
  - 热力学：`kinetics_occupancy.csv`、`kinetics_thermo.csv`、`kinetics_vanthoff.csv`
  - 截断与诊断：`kinetics_cutoffs.csv`、`kinetics_oh_density.csv`、`kinetics_events_by_traj.csv`、`kinetics_occupancy_vs_time.csv`、`kinetics_traces_example.csv`
  - 设定：`kinetics_metadata.json`
  - 图：`fig_*.svg`
- `comparison/`：
  - `compare_activation.csv`（三步 ΔG‡ 与 Ea 总表）
  - `compare_pKa_stepwise.csv`（逐步 pKa、ΔG、ΔΔG）
  - `compare_pKa_cumulative.csv`、`compare_box_occupancy.csv`、`compare_structure_by_nOH.csv`、`compare_PMF_step1.csv`
  - 图：`fig_compare_*.svg`
- `literature/`：论文 Table 1/2、S1/S2 的转录，SI 结构与几何指标，PMF 数字化与终点核对，`lit_metadata.json`。

每个文件的一句话说明见网页第 8 节。

## 没有打包的内容

- `kinetics_cache/`（约 372 MB）：用 `python kinetics.py extract 8` 从原始轨迹重新生成。
- 出版社的 SI 原文件 `literature/source/ic1c00662_si_001.pdf`：`literature.py` 会自动下载。
- 日志和 Python 缓存。

## 主要局限

详见网页第 7 节。

- **表观量：** ΔG‡ 和 Ea 依赖滤波窗口和截断。
  - 例如石墨 DPA4C-Plus Dy 3→2 的 ΔG‡ 在 W = 0.2 / 2.0 ps 时为 6.2 / 10.1 kJ/mol。
  - 三步的相对顺序（正向第 1 < 2 < 3 步，逆向相反）在 9 种设定 × 2 种金属下都成立，比绝对值可靠。
- **溶液统计量少：** 每个温度只有 1 条 150 ps 轨迹，Ea 的 95% CI 很宽。部分 Ea 为负，说明这是复合的表观速率常数，不是单一基元步骤。
- **pKa 等效值的假设：** pK_w = 14、理想溶液、标准态用 c₁ = 1 个分子 / box。游离 OH⁻ 浓度 (4 − n)c₁ 可达约 1 M，没有活度与离子对校正。
- **起始构型记忆：** 石墨 DPA4C-Plus 的占据数到约 90 ps 才收敛。
- **UMA 坍塌：** UMA 轨迹常坍塌，UMA 的结果只作参考。
- **势函数之间定性不同：** DPA4-Pro、MACE、UMA 中 3 个 OH 基本始终在 Ln 上，DPA4C-Plus 中 OH 会频繁离开。哪个更接近真实，需要 AIMD 或实验结构独立验证。
- **论文的可比性：** 论文只有 298 K 热力学，没有能垒，所以我们的 Ea 只能和论文的 ΔG 比量级。
- **OH 溶剂化壳层：** 石墨 DPA4C-Plus Nd 的 O–X(总) 壳层极小不清晰，不宜下结论，建议改看 O–REE(壳层) 与 REE–O(OH)/3。
