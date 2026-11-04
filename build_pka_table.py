"""Standalone page: stepwise pKa1/2/3 of Shiery 2021 (Table 1 DFT, Table 2 experiment) vs this work (MD).

Reads comparison/compare_pKa_stepwise.csv and writes pka_table.html (English) and pka_table_zh.html (Chinese).
"""
import csv
import html
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "comparison" / "compare_pKa_stepwise.csv"
OUT = {"en": HERE / "pka_table.html", "zh": HERE / "pka_table_zh.html"}

METALS = ("Nd", "Dy")
STEPS = (1, 2, 3)
PAPER_ROWS = (
    ("Shiery 2021 Table 1", "calc", "298.15", "Shiery2021 Table 1"),
    ("Shiery 2021 Table 2", "exp", "298.15", "Shiery2021 Table 2"),
)
MD_ROWS = (
    ("graphite", "DPA4C-Plus"),
    ("graphite", "uma-s-1.2.1"),
    ("solution", "DPA4C-Plus"),
    ("solution", "DPA4-Pro"),
    ("solution", "MACE-POLAR-1-M"),
    ("solution", "uma-s-1.2.1"),
)
MODEL_LABEL = {"uma-s-1.2.1": "UMA-s-1.2.1"}

T = {
    "zh": {
        "calc": "计算（DFT）", "exp": "实验（中位数）", "no_exp": "无实验值",
        "graphite": "石墨/水界面", "solution": "体相溶液",
        "this_work": "本工作 MD", "traj": lambda n: f"（{n} 条轨迹）",
    },
    "en": {
        "calc": "calculated (DFT)", "exp": "experiment (median)", "no_exp": "no exp. value",
        "graphite": "graphite/water", "solution": "bulk solution",
        "this_work": "This work, MD", "traj": lambda n: f" ({n} traj.)",
    },
}


def load():
    rows = list(csv.DictReader(SRC.open()))
    return rows


def find(rows, metal, step, **kw):
    for r in rows:
        if r["metal"] == metal and int(r["hydrolysis_step"]) == step and all(
            (r[k].startswith(v) if k == "source" else r[k] == v) for k, v in kw.items()
        ):
            return r
    return None


def small(text):
    return f'<span class="ci">{text}</span>'


def paper_cell(r, kind, lang):
    if r is None:
        return f'<td class="na">{T[lang]["no_exp"]}</td>'
    val = float(r["pKa"])
    if kind == "Shiery2021 Table 1":
        if r["pKa_CI95_lo"]:
            err = (float(r["pKa_CI95_hi"]) - float(r["pKa_CI95_lo"])) / 2
            return f"<td>{val:.1f} ± {err:.1f}</td>"
        return f"<td>{val:.1f}</td>"
    n = int(r["n_values"])
    lo, hi = float(r["pKa_CI95_lo"]), float(r["pKa_CI95_hi"])
    rng = small(f"[{lo:.2f}, {hi:.2f}]") if n > 1 else ""
    return f'<td>{val:.2f}{rng}<span class="n">n = {n}</span></td>'


def md_cell(r):
    if r is None:
        return '<td class="na">—</td>'
    st = r["status"]
    if st == "measured":
        return f'<td><b>{float(r["pKa"]):.2f}</b>{small(f"[{float(r["pKa_CI95_lo"]):.2f}, {float(r["pKa_CI95_hi"]):.2f}]")}</td>'
    if st == "transient":
        return (f'<td><b>{float(r["pKa"]):.2f}</b><sup>†</sup>'
                f'{small(f"[{float(r["pKa_CI95_lo"]):.2f}, {float(r["pKa_CI95_hi"]):.2f}]")}</td>')
    if st == "upper_bound":
        return f'<td class="bound">≤ {float(r["pKa_bound"]):.2f}<sup>‡</sup></td>'
    if st == "lower_bound":
        return f'<td class="bound">≥ {float(r["pKa_bound"]):.2f}<sup>‡</sup></td>'
    return '<td class="na">n.d.</td>'


def build_rows(rows, lang):
    t = T[lang]
    out = []
    for metal in METALS:
        md_present = [m for m in MD_ROWS
                      if any(find(rows, metal, s, source="this work", system=m[0], model=m[1]) for s in STEPS)]
        span = len(PAPER_ROWS) + len(md_present)
        first = True
        for label, method, temp, key in PAPER_ROWS:
            cells = "".join(paper_cell(find(rows, metal, s, source=key), key, lang) for s in STEPS)
            head = f'<th rowspan="{span}" class="metal">{metal}</th>' if first else ""
            first = False
            out.append(f'<tr class="paper">{head}<td>{label}</td><td>{t[method]}</td><td>{temp}</td>{cells}</tr>')
        for i, (system, model) in enumerate(md_present):
            recs = [find(rows, metal, s, source="this work", system=system, model=model) for s in STEPS]
            cells = "".join(md_cell(r) for r in recs)
            n_traj = max(int(r["n_values"]) for r in recs if r)
            cls = ' class="grp"' if i == 0 else ""
            out.append(f'<tr{cls}><td>{t["this_work"]} · {html.escape(MODEL_LABEL.get(model, model))}</td>'
                       f'<td>{t[system]}{t["traj"](n_traj)}</td><td>300</td>{cells}</tr>')
    return "\n".join(out)


# MathML snippets (rendered natively by current browsers; no external scripts).
M = {
    "ka": """<math display="block"><msub><mi>K</mi><mrow><mi mathvariant="normal">a</mi><mo>,</mo><mi>n</mi></mrow></msub><mo>=</mo>
<mfrac><mrow><mo>[</mo><msup><mrow><mi>Ln</mi><msub><mrow><mo>(</mo><mi>OH</mi><mo>)</mo></mrow><mi>n</mi></msub></mrow><mrow><mo>(</mo><mn>3</mn><mo>−</mo><mi>n</mi><mo>)</mo><mo>+</mo></mrow></msup><mo>]</mo><mo>[</mo><msup><mi mathvariant="normal">H</mi><mo>+</mo></msup><mo>]</mo></mrow>
<mrow><mo>[</mo><msup><mrow><mi>Ln</mi><msub><mrow><mo>(</mo><mi>OH</mi><mo>)</mo></mrow><mrow><mi>n</mi><mo>−</mo><mn>1</mn></mrow></msub></mrow><mrow><mo>(</mo><mn>4</mn><mo>−</mo><mi>n</mi><mo>)</mo><mo>+</mo></mrow></msup><mo>]</mo></mrow></mfrac>
<mo>,</mo><mspace width="1em"/><msub><mi>pK</mi><mrow><mi mathvariant="normal">a</mi><mo>,</mo><mi>n</mi></mrow></msub><mo>=</mo><mo>−</mo><msub><mi>log</mi><mn>10</mn></msub><msub><mi>K</mi><mrow><mi mathvariant="normal">a</mi><mo>,</mo><mi>n</mi></mrow></msub>
<mo>,</mo><mspace width="1em"/><mi>Δ</mi><msub><mi>G</mi><mi>n</mi></msub><mo>=</mo><mi>R</mi><mi>T</mi><mi>ln</mi><mn>10</mn><mo>·</mo><msub><mi>pK</mi><mrow><mi mathvariant="normal">a</mi><mo>,</mo><mi>n</mi></mrow></msub></math>""",
    "pmf": """<math display="block"><mi>A</mi><mo>(</mo><mi>r</mi><mo>)</mo><mo>=</mo><mo>−</mo><msubsup><mo>∫</mo><msub><mi>r</mi><mn>0</mn></msub><mi>r</mi></msubsup>
<msub><mrow><mo>⟨</mo><mi>f</mi><mo>⟩</mo></mrow><msup><mi>r</mi><mo>′</mo></msup></msub><mspace width="0.2em"/><mi>d</mi><msup><mi>r</mi><mo>′</mo></msup>
<mo>,</mo><mspace width="1em"/><msub><mi>pK</mi><mrow><mi mathvariant="normal">a</mi><mo>,</mo><mn>1</mn></mrow></msub><mo>=</mo>
<mfrac><mrow><mi>A</mi><mo>(</mo><msub><mi>r</mi><mi>end</mi></msub><mo>)</mo><mo>−</mo><mi>A</mi><mo>(</mo><msub><mi>r</mi><mn>0</mn></msub><mo>)</mo></mrow><mrow><mi>R</mi><mi>T</mi><mi>ln</mi><mn>10</mn></mrow></mfrac></math>""",
    "clus": """<math display="block"><msub><mi>pK</mi><mrow><mi mathvariant="normal">a</mi><mo>,</mo><mi>n</mi></mrow></msub><mo>=</mo><msub><mi>pK</mi><mrow><mi mathvariant="normal">a</mi><mo>,</mo><mn>1</mn></mrow></msub><mo>+</mo>
<mfrac><mrow><mi>Δ</mi><msubsup><mi>G</mi><mi>n</mi><mtext>M06/CPCM</mtext></msubsup><mo>−</mo><mi>Δ</mi><msubsup><mi>G</mi><mn>1</mn><mtext>M06/CPCM</mtext></msubsup></mrow><mrow><mi>R</mi><mi>T</mi><mi>ln</mi><mn>10</mn></mrow></mfrac>
<mo>,</mo><mspace width="1em"/><mi>n</mi><mo>=</mo><mn>2</mn><mo>,</mo><mn>3</mn></math>""",
    "kd": """<math display="block"><msub><mi>K</mi><mrow><mi mathvariant="normal">d</mi><mo>,</mo><mi>n</mi></mrow></msub><mo>=</mo>
<mfrac><mrow><mo>[</mo><mi>Ln</mi><msub><mrow><mo>(</mo><mi>OH</mi><mo>)</mo></mrow><mrow><mi>n</mi><mo>−</mo><mn>1</mn></mrow></msub><mo>]</mo><mo>[</mo><msup><mi>OH</mi><mo>−</mo></msup><mo>]</mo></mrow><mrow><mo>[</mo><mi>Ln</mi><msub><mrow><mo>(</mo><mi>OH</mi><mo>)</mo></mrow><mi>n</mi></msub><mo>]</mo></mrow></mfrac>
<mo>=</mo><mfrac><msub><mi>P</mi><mrow><mi>n</mi><mo>−</mo><mn>1</mn></mrow></msub><msub><mi>P</mi><mi>n</mi></msub></mfrac><mo>·</mo>
<mfrac><mrow><mo>(</mo><mn>4</mn><mo>−</mo><mi>n</mi><mo>)</mo><mspace width="0.15em"/><msub><mi>c</mi><mn>1</mn></msub></mrow><msup><mi>c</mi><mo>°</mo></msup></mfrac>
<mo>,</mo><mspace width="1em"/><msub><mi>c</mi><mn>1</mn></msub><mo>=</mo><mfrac><mn>1</mn><mrow><msub><mi>N</mi><mi mathvariant="normal">A</mi></msub><mi>V</mi></mrow></mfrac></math>""",
    "pka": """<math display="block"><msub><mi>K</mi><mrow><mi mathvariant="normal">a</mi><mo>,</mo><mi>n</mi></mrow></msub><mo>=</mo><mfrac><msub><mi>K</mi><mi mathvariant="normal">w</mi></msub><msub><mi>K</mi><mrow><mi mathvariant="normal">d</mi><mo>,</mo><mi>n</mi></mrow></msub></mfrac>
<mspace width="0.6em"/><mo>⇒</mo><mspace width="0.6em"/>
<msub><mi>pK</mi><mrow><mi mathvariant="normal">a</mi><mo>,</mo><mi>n</mi></mrow></msub><mo>=</mo><msub><mi>pK</mi><mi mathvariant="normal">w</mi></msub><mo>+</mo><msub><mi>log</mi><mn>10</mn></msub><msub><mi>K</mi><mrow><mi mathvariant="normal">d</mi><mo>,</mo><mi>n</mi></mrow></msub>
<mo>=</mo><mn>14</mn><mo>+</mo><msub><mi>log</mi><mn>10</mn></msub><mrow><mo>[</mo><mfrac><msub><mi>P</mi><mrow><mi>n</mi><mo>−</mo><mn>1</mn></mrow></msub><msub><mi>P</mi><mi>n</mi></msub></mfrac><mo>(</mo><mn>4</mn><mo>−</mo><mi>n</mi><mo>)</mo><msub><mi>c</mi><mn>1</mn></msub><mo>]</mo></mrow></math>""",
    "bound": """<math display="block"><msub><mi>pK</mi><mrow><mi mathvariant="normal">a</mi><mo>,</mo><mi>n</mi></mrow></msub><mo>≤</mo><mn>14</mn><mo>+</mo><msub><mi>log</mi><mn>10</mn></msub><mrow><mo>[</mo><mo>(</mo><mn>4</mn><mo>−</mo><mi>n</mi><mo>)</mo><msub><mi>c</mi><mn>1</mn></msub><mo>]</mo></mrow>
<mspace width="1em"/><mo>(</mo><msub><mi>P</mi><mrow><mi>n</mi><mo>−</mo><mn>1</mn></mrow></msub><mo>:=</mo><msub><mi>P</mi><mi>n</mi></msub><mo>)</mo></math>""",
}

CSS = """
body{font-family:-apple-system,"Segoe UI","PingFang SC","Noto Sans CJK SC","Microsoft YaHei",sans-serif;
 color:#222;max-width:1100px;margin:32px auto;padding:0 20px;line-height:1.55;font-size:15px}
h1{font-size:20px;margin:0 0 4px}
.ref{color:#666;font-size:13px;margin:0 0 18px}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{padding:6px 10px;border-bottom:1px solid #e3e3e3;text-align:left;vertical-align:top}
thead th{border-bottom:2px solid #333;background:#fafafa}
tbody td:nth-last-child(-n+4),thead th:nth-last-child(-n+4){text-align:right;white-space:nowrap}
tr.paper td{background:#f4f7fb}
tr.grp td{border-top:1px solid #bbb}
th.metal{font-size:17px;vertical-align:middle;text-align:center;border-right:1px solid #e3e3e3;background:#fff}
tbody tr:last-child td,tbody tr:last-child th{border-bottom:2px solid #333}
.ci{display:block;color:#777;font-size:12px}
.n{display:block;color:#999;font-size:11px}
.na{color:#aaa}
.bound{color:#555}
h2{font-size:16px;margin:28px 0 6px;border-bottom:1px solid #ddd;padding-bottom:3px}
h3{font-size:14.5px;margin:16px 0 4px}
.rx{font-family:"Cambria Math","STIX Two Math","Times New Roman",serif;font-size:15.5px;margin:4px 0 4px 18px}
math{font-family:"Cambria Math","STIX Two Math","Latin Modern Math",math}
math[display=block]{margin:8px 0 8px 18px;font-size:16px}
ul{margin:4px 0;padding-left:22px}
li{margin:2px 0}
.note{font-size:13px;color:#555}
.lang{float:right;font-size:13px;margin-top:4px}
.lang a,.lang b{margin-left:8px}
"""


def lnoh(k, c):
    return f"Ln(OH)<sub>{k}</sub><sup>{c}</sup>" if k else f"Ln<sup>{c}</sup>"


RXN = [
    f'{lnoh(0, "3+")}(H<sub>2</sub>O) ⇌ {lnoh(1, "2+")} + H<sup>+</sup>',
    f'{lnoh(1, "2+")}(H<sub>2</sub>O) ⇌ {lnoh(2, "+")} + H<sup>+</sup>',
    f'{lnoh(2, "+")}(H<sub>2</sub>O) ⇌ {lnoh(3, "")} + H<sup>+</sup>',
]
RXN_HTML = "\n".join(f'<div class="rx">n = {i + 1}:&nbsp; {r}</div>' for i, r in enumerate(RXN))
RXN_MD = f'<div class="rx">{lnoh("n", "(3−n)+")} ⇌ {lnoh("n−1", "(4−n)+")} + OH<sup>−</sup></div>'
REF = ('Shiery, Cooper, Cantu, <i>Inorg. Chem.</i> <b>2021</b>, 60, 10257 · '
       '<a href="https://doi.org/10.1021/acs.inorgchem.1c00662">doi:10.1021/acs.inorgchem.1c00662</a>')

BODY = {
    "zh": f"""
<p class="note">
本工作数值为点估计，灰色方括号为 bootstrap 95% CI；实验行为 Table 2 所列全部实验值的中位数，方括号为 [最小, 最大]，n 为实验值个数。
† 两个态都出现，但其中一个态的 committed 停留 &lt; 1 ps，仅作数量级参考。
‡ 只出现一个态，给出上界（见下文）。n.d.：两个态都未出现，无法确定。
每条轨迹 150 ps。Dy 的 UMA-s-1.2.1 溶液轨迹全部坍塌，无数据；石墨 UMA 已剔除或截断坍塌轨迹。</p>

<h2>化学反应与定义</h2>
<p>第 n 步酸解离（n = 1, 2, 3；省略配位水）：</p>
{RXN_HTML}
{M["ka"]}

<h2>计算方法</h2>
<h3>Shiery 2021 — 计算值（Table 1，298.15 K）</h3>
<ul>
<li><b>pK<sub>a1</sub></b>：约束 AIMD（Bluemoon，CP2K/PBE，64 H<sub>2</sub>O + 3 Cl<sup>−</sup>）。以配位水的 r(O–H) 为反应坐标，0.99 → 1.40 Å 取 9 个约束点，每点 NVT 15 ps；对约束力 ⟨f⟩ 积分得 PMF，取端点自由能换算 pK<sub>a1</sub>（误差为 Table 1 所给）。</li>
</ul>
{M["pmf"]}
<ul>
<li><b>pK<sub>a2</sub>、pK<sub>a3</sub></b>：M06 团簇（第一壳层水与 OH 显式）+ CPCM 隐式溶剂计算三步反应自由能，以 pK<sub>a1</sub> 为锚点取相对值；团簇配位数固定：Nd 为 9，Dy 为 8。</li>
</ul>
{M["clus"]}
<h3>Shiery 2021 — 实验值（Table 2）</h3>
<ul><li>文献实验值汇编，表中取中位数。注意：第 2、3 步的实验值在数值上更接近累积常数 −log*β<sub>n</sub>（Ln<sup>3+</sup> + n H<sub>2</sub>O ⇌ Ln(OH)<sub>n</sub><sup>(3−n)+</sup> + n H<sup>+</sup>），而论文按逐步常数使用。</li></ul>

<h3>本工作 — MD 占据数 → pK<sub>a</sub> 等效值（300 K）</h3>
<ul>
<li>体系：Ln<sup>3+</sup> + 3 OH<sup>−</sup> + 水（溶液 128 H<sub>2</sub>O；石墨/水界面 4 种 OH 起始构型 × 50 条），三维周期盒子。</li>
<li>状态 n = Ln 第一壳层内 OH<sup>−</sup> 的个数（OH = 恰有 1 个最近邻 H 的 O；截断 r<sub>c</sub> = Ln–O(OH) 径向密度第一极小，3.03–3.47 Å），每 20 fs 取一帧。</li>
<li>P<sub>n</sub>：n 态的平衡占据概率（溶液取 t ≥ 10 ps，石墨取 t ≥ 50 ps）；未结合的 OH<sup>−</sup> 有 (4 − n) 个。</li>
<li>c<sub>1</sub>：盒子中 1 个分子对应的浓度，溶液 0.4224 M，石墨 0.3594 M；c° = 1 M。理想溶液，无活度修正。</li>
</ul>
<p>MD 中直接观测的是 OH<sup>−</sup> 解离：</p>
{RXN_MD}
{M["kd"]}
<p>与水的自电离 H<sub>2</sub>O ⇌ H<sup>+</sup> + OH<sup>−</sup>（K<sub>w</sub>，pK<sub>w</sub> = 14）相加即得上面的酸解离反应：</p>
{M["pka"]}
<p>若只出现 n 态（n − 1 态从未出现），假设未出现的态不比出现的态更多，取 P<sub>n−1</sub> = P<sub>n</sub> 得上界：</p>
{M["bound"]}
<ul>
<li>95% CI：2000 次 bootstrap（石墨按轨迹、按起始构型分层重采样；溶液按 10 ps 块重采样）。</li>
<li>换算：ΔG = RT ln10 · pK<sub>a</sub>，298.15 K 时 RT ln10 = 5.708 kJ/mol，300 K 时 5.743 kJ/mol。</li>
</ul>
""",
    "en": f"""
<p class="note">
This-work values are point estimates with bootstrap 95% CIs in grey brackets. Experiment rows give the median of all Table 2 experimental values,
with [min, max] in brackets and n = number of values.
† Both states occur, but one has committed dwells &lt; 1 ps; order of magnitude only.
‡ Only one state occurs; an upper bound is given (see below). n.d.: neither state occurs, not determinable.
Each trajectory is 150 ps. All UMA-s-1.2.1 Dy solution trajectories collapsed (no data); collapsed graphite UMA trajectories were excluded or truncated.</p>

<h2>Reactions and definitions</h2>
<p>Step-n acid dissociation (n = 1, 2, 3; spectator waters omitted):</p>
{RXN_HTML}
{M["ka"]}

<h2>Methods</h2>
<h3>Shiery 2021, calculated (Table 1, 298.15 K)</h3>
<ul>
<li><b>pK<sub>a1</sub></b>: constrained AIMD (Bluemoon, CP2K/PBE, 64 H<sub>2</sub>O + 3 Cl<sup>−</sup>). The r(O–H) of a coordinated water is the reaction coordinate, with 9 constraint points from 0.99 to 1.40 Å and 15 ps NVT each; integrating the mean constraint force ⟨f⟩ gives the PMF, whose end-point free energy is converted to pK<sub>a1</sub> (error as given in Table 1).</li>
</ul>
{M["pmf"]}
<ul>
<li><b>pK<sub>a2</sub>, pK<sub>a3</sub></b>: reaction free energies of the three steps from M06 clusters (explicit first-shell water and OH) + CPCM implicit solvent, taken relative to pK<sub>a1</sub>; the cluster coordination number is fixed at 9 for Nd and 8 for Dy.</li>
</ul>
{M["clus"]}
<h3>Shiery 2021, experiment (Table 2)</h3>
<ul><li>Compilation of literature values; the median is listed. Note: the step-2 and step-3 values are numerically closer to cumulative constants −log*β<sub>n</sub> (Ln<sup>3+</sup> + n H<sub>2</sub>O ⇌ Ln(OH)<sub>n</sub><sup>(3−n)+</sup> + n H<sup>+</sup>), whereas the paper uses them as stepwise constants.</li></ul>

<h3>This work: MD occupancy → pK<sub>a</sub> equivalent (300 K)</h3>
<ul>
<li>Systems: Ln<sup>3+</sup> + 3 OH<sup>−</sup> + water (solution: 128 H<sub>2</sub>O; graphite/water interface: 4 initial OH configurations × 50 trajectories), 3D-periodic box.</li>
<li>State n = number of OH<sup>−</sup> in the first shell of Ln (OH = an O with exactly one nearest H; cutoff r<sub>c</sub> = first minimum of the Ln–O(OH) radial density, 3.03–3.47 Å), sampled every 20 fs.</li>
<li>P<sub>n</sub>: equilibrium probability of state n (t ≥ 10 ps for solution, t ≥ 50 ps for graphite); there are (4 − n) unbound OH<sup>−</sup>.</li>
<li>c<sub>1</sub>: concentration of one molecule in the box, 0.4224 M (solution) and 0.3594 M (graphite); c° = 1 M. Ideal solution, no activity corrections.</li>
</ul>
<p>What MD observes directly is OH<sup>−</sup> dissociation:</p>
{RXN_MD}
{M["kd"]}
<p>Adding water autoionization H<sub>2</sub>O ⇌ H<sup>+</sup> + OH<sup>−</sup> (K<sub>w</sub>, pK<sub>w</sub> = 14) gives the acid dissociation above:</p>
{M["pka"]}
<p>If only state n occurs (state n − 1 never seen), assuming the unseen state is not more populated than the seen one and setting P<sub>n−1</sub> = P<sub>n</sub> gives an upper bound:</p>
{M["bound"]}
<ul>
<li>95% CI: 2000 bootstrap resamples (graphite: by trajectory, stratified by initial configuration; solution: 10 ps blocks).</li>
<li>Conversion: ΔG = RT ln10 · pK<sub>a</sub>; RT ln10 = 5.708 kJ/mol at 298.15 K and 5.743 kJ/mol at 300 K.</li>
</ul>
""",
}

HEAD = {
    "zh": {"html_lang": "zh-CN", "title": "Nd/Dy 逐步水解 pKa：Shiery 2021 vs 本工作",
           "h1": "Nd、Dy 逐步水解 pKa：Shiery 2021 与本工作对比",
           "cols": ("Ln", "来源", "方法 / 体系", "T (K)"),
           "switch": '<a href="pka_table.html">English</a><b>中文</b>'},
    "en": {"html_lang": "en", "title": "Nd/Dy stepwise hydrolysis pKa: Shiery 2021 vs this work",
           "h1": "Nd and Dy stepwise hydrolysis pKa: Shiery 2021 vs this work",
           "cols": ("Ln", "Source", "Method / system", "T (K)"),
           "switch": '<b>English</b><a href="pka_table_zh.html">中文</a>'},
}


def page(rows, lang):
    h = HEAD[lang]
    cols = "".join(f"<th>{c}</th>" for c in h["cols"])
    return f"""<!DOCTYPE html>
<html lang="{h['html_lang']}"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{h['title']}</title>
<style>{CSS}</style></head><body>
<div class="lang">{h['switch']}</div>
<h1>{h['h1']}</h1>
<p class="ref">{REF}</p>
<table>
<thead><tr>{cols}<th>pK<sub>a1</sub></th><th>pK<sub>a2</sub></th><th>pK<sub>a3</sub></th></tr></thead>
<tbody>
{build_rows(rows, lang)}
</tbody></table>
{BODY[lang]}
</body></html>
"""


def main():
    rows = load()
    for lang, path in OUT.items():
        path.write_text(page(rows, lang), encoding="utf-8")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
