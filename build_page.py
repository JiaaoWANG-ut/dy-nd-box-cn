#!/usr/bin/env python3
"""Static report page + data bundle for the Dy/Nd box-PBC analysis.

Reads the CSV/JSON/SVG outputs of analyze.py, kinetics.py, literature.py and
compare.py (nothing is recomputed) and writes

  index.html               Chinese report; every table/figure links its raw data
  dy_nd_box_cn_data.zip    all outputs + scripts + index.html.  kinetics_cache/
                           (regenerable, 380 MB) and the publisher's SI PDF
                           (literature.py downloads it again) are not bundled.

Run:  python build_page.py
"""

from __future__ import annotations

import csv
import html
import json
import math
import os
import re
import zipfile
from collections import Counter, defaultdict
from datetime import datetime
from html.parser import HTMLParser
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
OUT_HTML = HERE / "index.html"
OUT_ZIP = HERE / "dy_nd_box_cn_data.zip"
ZIP_ROOT = "dy_nd_box_cn"
SKIP_DIR_NAMES = {"kinetics_cache", "__pycache__"}
SKIP_RELS = {"literature/source/ic1c00662_si_001.pdf"}
SKIP_FILES = {"kinetics_extract.log", OUT_ZIP.name}

ARTICLE_URL = "https://pubs.acs.org/doi/10.1021/acs.inorgchem.1c00662"
SI_URL = "https://ndownloader.figshare.com/files/28696985"
THESIS_URL = "https://scholarwolf.unr.edu/bitstreams/ae32a3d1-8c2e-40ab-a18f-9499bf7f21e1/download"
EXTERNAL_URLS = {ARTICLE_URL, SI_URL, THESIS_URL}

MINUS = "\u2212"
METALS = ("Dy", "Nd")
MODEL_ORDER = ("DPA4C-Plus", "DPA4-Pro", "MACE-POLAR-1-M", "uma-s-1.2.1")
MODEL_LABEL = {"DPA4C-Plus": "DPA4C-Plus", "DPA4-Pro": "DPA4-Pro",
               "MACE-POLAR-1-M": "MACE-POLAR-1-M", "uma-s-1.2.1": "UMA-s-1.2.1"}
SYSTEM_LABEL = {"graphite": "石墨表面", "solution": "溶液", "cluster": "DFT 团簇（M06 + CPCM）"}
SERIES = (("graphite", "DPA4C-Plus"), ("graphite", "uma-s-1.2.1"),
          ("solution", "DPA4C-Plus"), ("solution", "DPA4-Pro"),
          ("solution", "MACE-POLAR-1-M"), ("solution", "uma-s-1.2.1"))
GRAPHITE_STARTS = ("3near", "2near1far", "1near2far", "3far")
SCENARIOS = ("primary", "W0.2", "W1.0", "W2.0", "rc2.90", "rc3.30", "start1", "start50", "start90")
SCENARIO_LABEL = {
    "primary": "主分析：W = 0.5 ps，r<sub>c</sub> = 密度极小，t ≥ 10 ps",
    "W0.2": "过滤窗口 W = 0.2 ps", "W1.0": "过滤窗口 W = 1.0 ps", "W2.0": "过滤窗口 W = 2.0 ps",
    "rc2.90": "固定 r<sub>c</sub> = 2.90 Å", "rc3.30": "固定 r<sub>c</sub> = 3.30 Å",
    "start1": "分析起点 t ≥ 1 ps", "start50": "分析起点 t ≥ 50 ps", "start90": "分析起点 t ≥ 90 ps",
}
DIRECTIONS = (("association", "正向"), ("dissociation", "逆向"))

esc = html.escape


# ------------------------------------------------------------------ formatting
def num(value) -> float:
    if value is None:
        return math.nan
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if text.lower() in ("", "nan", "none", "null"):
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def finite(value) -> bool:
    return math.isfinite(num(value))


def inum(value) -> int:
    x = num(value)
    return int(round(x)) if math.isfinite(x) else 0


def fmt(value, digits: int = 2, sign: bool = False) -> str:
    x = num(value)
    if math.isnan(x):
        return "—"
    if math.isinf(x):
        return "∞" if x > 0 else MINUS + "∞"
    text = f"{x:+.{digits}f}" if sign else f"{x:.{digits}f}"
    if text.lstrip("+-") == f"{0:.{digits}f}":
        text = text.lstrip("+-")
    return text.replace("-", MINUS)


def fmt_p(value) -> str:
    """Probability: 3 decimals, scientific notation below 1e-3."""
    x = num(value)
    if math.isnan(x):
        return "—"
    if 0 < x < 1e-3:
        mant, expo = f"{x:.1e}".split("e")
        return f"{mant}×10<sup>{MINUS if int(expo) < 0 else ''}{abs(int(expo))}</sup>"
    return fmt(x, 3)


def ci(lo, hi, digits: int = 2) -> str:
    a, b = num(lo), num(hi)
    if math.isnan(a) and math.isnan(b):
        return ""
    return f"[{fmt(a, digits)}, {fmt(b, digits)}]"


def pm(mean, sem, digits: int = 2) -> str:
    m = num(mean)
    if math.isnan(m):
        return "—"
    s = num(sem)
    if math.isnan(s) or s == 0:
        return fmt(m, digits)
    tiny = 0.5 * 10 ** -digits
    err = f"&lt;{fmt(2 * tiny, digits)}" if s < tiny else fmt(s, digits)
    return f'{fmt(m, digits)} <span class="pm">± {err}</span>'


def rng(values, digits: int = 2) -> str:
    vals = [x for x in (num(v) for v in values) if math.isfinite(x)]
    if not vals:
        return "—"
    lo, hi = min(vals), max(vals)
    if fmt(lo, digits) == fmt(hi, digits):
        return fmt(lo, digits)
    sep = "–" if lo >= 0 else " … "
    return f"{fmt(lo, digits)}{sep}{fmt(hi, digits)}"


def plain(fragment: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", "", fragment))


def metal_html(metal: str) -> str:
    return f'<span class="{esc(metal)}">{esc(metal)}</span>'


def pill(text: str, kind: str = "") -> str:
    return f'<span class="pill {kind}">{text}</span>' if kind else f'<span class="pill">{text}</span>'


QUALITY_ZH = {"good": "好", "fair": "尚可", "poor": "差", "none": "无"}


def quality_pill(quality: str) -> str:
    kind = {"good": "ok", "fair": "warn", "poor": "bad"}.get(quality, "")
    return pill(QUALITY_ZH.get(quality or "none", esc(quality)), kind)


def size_label(nbytes: int) -> str:
    if nbytes < 1024:
        return f"{nbytes} B"
    if nbytes < 1024 ** 2:
        kb = nbytes / 1024
        return f"{kb:.1f} KB" if kb < 10 else f"{kb:.0f} KB"
    return f"{nbytes / 1024 ** 2:.1f} MB"


# ------------------------------------------------------------------ links
LINKED: set[str] = set()


def _register(rel: str) -> Path:
    path = HERE / rel
    if not path.is_file():
        raise FileNotFoundError(f"page links a missing file: {rel}")
    LINKED.add(rel)
    return path


def dl(rel: str, text: str | None = None) -> str:
    path = _register(rel)
    return (f'<a class="dl" href="{esc(rel)}" download>{esc(text or path.name)}</a>'
            f'<span class="size">{size_label(path.stat().st_size)}</span>')


def data_line(*rels: str, label: str = "原始数据") -> str:
    return f'<p class="files">{label}：' + " · ".join(dl(r) for r in rels) + "</p>"


def figure(rel: str, caption: str, data: tuple[str, ...] = ()) -> str:
    _register(rel)
    links = [f'<a href="{esc(rel)}" target="_blank" rel="noopener">新窗口打开</a>', dl(rel, "下载 SVG")]
    if data:
        links.append("作图数据：" + " · ".join(dl(r) for r in data))
    return (f'<figure class="fig"><a href="{esc(rel)}" target="_blank" rel="noopener">'
            f'<img src="{esc(rel)}" alt="{esc(plain(caption))}" loading="lazy"></a>'
            f'<figcaption>{caption}<div class="files">{" · ".join(links)}</div></figcaption></figure>')


# ------------------------------------------------------------------ tables
def table(head, rows, num_cols=(), cls: str = "", thead: str | None = None) -> str:
    """rows: list of cell lists or {"cells": [...], "cls": "..."}; a cell may be (html, attrs)."""
    parts = ['<div class="table-wrap">', f'<table class="{cls}">' if cls else "<table>"]
    parts.append(thead if thead is not None else
                 "<thead><tr>" + "".join(f"<th>{h}</th>" for h in head) + "</tr></thead>")
    parts.append("<tbody>")
    for row in rows:
        cells, rcls = (row["cells"], row.get("cls", "")) if isinstance(row, dict) else (row, "")
        tds = []
        for j, cell in enumerate(cells):
            if isinstance(cell, tuple):
                tds.append(f"<td {cell[1]}>{cell[0]}</td>")
            elif j in num_cols:
                tds.append(f'<td class="num">{cell}</td>')
            else:
                tds.append(f"<td>{cell}</td>")
        parts.append((f'<tr class="{rcls}">' if rcls else "<tr>") + "".join(tds) + "</tr>")
    parts.append("</tbody></table></div>")
    return "\n".join(parts)


def thead2(fixed, groups) -> str:
    """Two-row header: `fixed` labels span both rows; groups = [(label, [sub labels])]."""
    top = "".join(f'<th rowspan="2">{h}</th>' for h in fixed)
    top += "".join(f'<th class="grp" colspan="{len(subs)}">{g}</th>' for g, subs in groups)
    low = "".join(f"<th>{s}</th>" for _, subs in groups for s in subs)
    return f"<thead><tr>{top}</tr><tr>{low}</tr></thead>"


def details(summary: str, body: str) -> str:
    return f"<details><summary>{summary}</summary>{body}</details>"


def card(anchor: str, title: str, sub: str, body: str) -> str:
    return f'<section class="card" id="{anchor}"><h2>{title}</h2><p class="sub">{sub}</p>{body}</section>'


# ------------------------------------------------------------------ data
def read(rel: str) -> list[dict]:
    with open(HERE / rel, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def where(rows: list[dict], **conds) -> list[dict]:
    out = []
    for row in rows:
        for key, want in conds.items():
            have = row.get(key)
            if callable(want):
                if not want(have):
                    break
            elif isinstance(want, (int, float)) and not isinstance(want, bool):
                if num(have) != want:
                    break
            elif have != want:
                break
        else:
            out.append(row)
    return out


def first(rows: list[dict], **conds) -> dict | None:
    hits = where(rows, **conds)
    return hits[0] if hits else None


def load() -> SimpleNamespace:
    return SimpleNamespace(
        cn=read("cn_summary_300K.csv"),
        sol_all=read("solution_cn_box_pbc.csv"),
        scut=read("solution_cutoffs_box_pbc_300K.csv"),
        gcut=read("graphite_cutoffs_box_pbc.csv"),
        screen=read("collapse_screen.csv"),
        rates=read("kinetics/kinetics_rates.csv"),
        arr=read("kinetics/kinetics_arrhenius.csv"),
        occ=read("kinetics/kinetics_occupancy.csv"),
        thermo=read("kinetics/kinetics_thermo.csv"),
        vth=read("kinetics/kinetics_vanthoff.csv"),
        kcut=read("kinetics/kinetics_cutoffs.csv"),
        kmeta=json.loads((HERE / "kinetics/kinetics_metadata.json").read_text()),
        act=read("comparison/compare_activation.csv"),
        stepw=read("comparison/compare_pKa_stepwise.csv"),
        cum=read("comparison/compare_pKa_cumulative.csv"),
        boxocc=read("comparison/compare_box_occupancy.csv"),
        struct=read("comparison/compare_structure_by_nOH.csv"),
        pmf=read("comparison/compare_PMF_step1.csv"),
        t1=read("literature/lit_table1_calc_pKa.csv"),
        t2=read("literature/lit_table2_pKa_all.csv"),
        madchk=read("literature/lit_table2_mad_check.csv"),
        s1=read("literature/lit_tableS1_multiplicity.csv"),
        s2=read("literature/lit_tableS2_basis_set_MAD.csv"),
        si=read("literature/lit_SI_structures_metrics.csv"),
        pmf_end=read("literature/lit_SI_PMF_endpoint_check.csv"),
        lmeta=json.loads((HERE / "literature/lit_metadata.json").read_text()),
    )


def cn_get(d, model: str, metal: str, variant: str) -> dict | None:
    return first(d.cn, model=model, metal=metal, variant=variant)


def cut_get(d, system: str, model: str, metal: str) -> dict | None:
    if system == "graphite":
        return first(d.gcut, model=model, metal=metal, variant="ALL")
    return first(d.scut, model=model, metal=metal)


def act_get(d, system: str, model: str, metal: str, step: int, direction: str) -> dict | None:
    return first(d.act, system=system, model=model, metal=metal, hydrolysis_step=str(step), direction=direction)


def series_variant(system: str) -> str:
    return "ALL" if system == "graphite" else "box_pbc"


# ------------------------------------------------------------------ page chrome
CSS = """
:root{--ink:#17263b;--muted:#607086;--line:#dce4ee;--paper:#f5f8fc;--card:#fff;--nd:#d95f02;--dy:#0072b2;--accent:#2459a6;--warn:#fff5df;--warnline:#e6b950}
*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.58 system-ui,-apple-system,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif}
header{background:linear-gradient(120deg,#132b49,#244b78);color:#fff;padding:30px max(24px,calc((100vw - 1200px)/2));}
.head{max-width:1200px;margin:auto}.eyebrow{font-size:12px;letter-spacing:.13em;opacity:.78}
h1{font-size:clamp(24px,3.6vw,34px);line-height:1.25;margin:13px 0 8px}header p{max-width:900px;margin:0;color:#e3edf8}
.btns{display:flex;gap:10px;flex-wrap:wrap;margin-top:16px}.btn{display:inline-block;background:#fff;color:#183558;padding:9px 15px;border-radius:9px;font-weight:650;text-decoration:none}
.btn.ghost{background:transparent;color:#fff;border:1px solid #ffffff88}
main{max-width:1200px;margin:24px auto 60px;padding:0 20px}
nav.toc{display:flex;flex-wrap:wrap;gap:8px;margin:4px 0 6px}nav.toc a{background:#fff;border:1px solid var(--line);border-radius:99px;padding:4px 12px;color:var(--accent);text-decoration:none;font-size:13.5px}
.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:16px 0}
.metric,.card{background:var(--card);border:1px solid var(--line);border-radius:13px;box-shadow:0 3px 14px #1a365008}.metric{padding:15px 17px}.metric b{display:block;font-size:21px;line-height:1.25}.metric span{color:var(--muted);font-size:13px}
.notice{border:1px solid var(--warnline);background:var(--warn);padding:16px 18px;border-radius:12px;margin:18px 0 22px}.notice strong{display:block;margin-bottom:4px}.notice ul{margin:6px 0 0;padding-left:20px}
.card{padding:20px;margin:18px 0}.card h2{margin:0 0 5px;font-size:21px}.card h3{font-size:17px;margin:24px 0 6px}.card h4{font-size:15px;margin:16px 0 6px}.sub{color:var(--muted);margin:0 0 15px;font-size:13px}
.table-wrap{overflow:auto;margin:8px 0}table{width:100%;border-collapse:collapse;font-size:14px}th,td{text-align:left;padding:8px 9px;border-bottom:1px solid var(--line);vertical-align:top}th{background:#f4f7fb;font-weight:650}
th.grp{text-align:center;border-left:2px solid #fff}td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
table.compact{font-size:13px}table.compact th,table.compact td{padding:5px 7px}
td small{display:block;color:var(--muted);font-size:11.5px;line-height:1.3;white-space:nowrap}.pm{color:var(--muted);font-size:12px}
tr.hl td{background:#fff8ec}tr.paper td{background:#eef3fa}tr.weak td{color:var(--muted)}td.weak,td.weak b{color:var(--muted)}
.pill{display:inline-block;border-radius:99px;padding:1px 8px;background:#edf3fa;color:#345372;font-size:12px;white-space:nowrap}
.pill.ok{background:#e7f6ec;color:#1d6b3a}.pill.warn{background:#fff1d6;color:#8a5a00}.pill.bad{background:#fde8e8;color:#9b1c1c}
.Dy{color:var(--dy);font-weight:650}.Nd{color:var(--nd);font-weight:650}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:14px}.mini{background:#f7f9fc;padding:14px;border-radius:10px}.mini h3,.mini h4{margin:0 0 8px;font-size:16px}
figure.fig{margin:16px 0;background:#fff;border:1px solid var(--line);border-radius:10px;padding:10px}figure.fig img{width:100%;height:auto;display:block}
figcaption{font-size:13px;color:var(--muted);margin-top:6px}.files{font-size:13px;margin:6px 0}a.dl{font-weight:600}.size{color:var(--muted);font-size:12px;margin-left:4px}
details{margin:12px 0}details>summary{cursor:pointer;font-weight:650;color:var(--accent)}
.formula{font-family:"Cambria Math","STIX Two Math","Times New Roman",serif;font-size:15.5px;background:#f7f9fc;padding:8px 12px;border-radius:8px;margin:6px 0;overflow-x:auto}
ul.tight li{margin:4px 0}pre{background:#0f1f33;color:#e3edf8;padding:12px 14px;border-radius:10px;overflow:auto;font-size:13px}
.foot{color:var(--muted);font-size:13px;padding:10px 4px 30px}.muted{color:var(--muted)}code{background:#f0f3f8;padding:1px 5px;border-radius:4px}
@media(max-width:760px){.metrics{grid-template-columns:repeat(2,1fr)}.grid2{grid-template-columns:1fr}.card{padding:15px}header{padding:22px 18px}}
@media print{body{background:#fff}.card,.metric{box-shadow:none;break-inside:avoid}header{background:#fff;color:var(--ink);padding:10px}header p{color:var(--muted)}.btns{display:none}main{margin-top:8px}details{display:block}}
"""

TOC = (("summary", "1 结论摘要"), ("methods", "2 方法"), ("cn", "3 配位数"), ("kinetics", "4 三步表观活化能"),
       ("occupancy", "5 占据数与 pKa 等效"), ("paper", "6 与 Shiery 2021 对比"), ("limits", "7 局限性"),
       ("downloads", "8 数据下载"), ("sources", "9 来源与复现"))


def page_header(zip_label: str) -> str:
    _register("README.md")
    return f"""<header><div class="head">
<div class="eyebrow">Dy / Nd · Ln(OH)<sub>3</sub> + H<sub>2</sub>O · 3D 周期 box (PBC) · 机器学习势函数</div>
<h1>Dy、Nd 水解：配位数、三步表观活化能，以及与 Shiery 2021 (Inorg. Chem. 60, 10257) 的对比</h1>
<p>石墨表面与溶液两种体系（都只用 box PBC 轨迹），比较 DPA4C-Plus、DPA4-Pro、MACE-POLAR-1-M、UMA-s-1.2.1 四个势函数。每张表和图下面都有对应原始数据 (CSV/SVG/JSON) 的下载链接，也可以一次下载全部。</p>
<div class="btns"><a class="btn" href="{OUT_ZIP.name}" download>⬇ 下载全部数据与脚本 (zip, {zip_label})</a><a class="btn ghost" href="#downloads">逐个文件下载</a><a class="btn ghost" href="README.md" target="_blank" rel="noopener">README</a></div>
</div></header>"""


def toc() -> str:
    return '<nav class="toc">' + "".join(f'<a href="#{a}">{t}</a>' for a, t in TOC) + "</nav>"


def metrics(d) -> str:
    n_tot = len(d.screen)
    n_gra = sum(1 for r in d.screen if r["system"] == "graphite")
    status = Counter(r["status"] for r in d.screen)
    gd, gn = cn_get(d, "DPA4C-Plus", "Dy", "ALL"), cn_get(d, "DPA4C-Plus", "Nd", "ALL")
    sol_tot = {m: [r["REE_O_total_mean"] for r in where(d.cn, variant="box_pbc", metal=m)] for m in METALS}
    main_rows = [r for r in d.cn if r["variant"] in ("ALL", "box_pbc")]
    oree_dpa = [r["OH_REE_shell_mean"] for r in main_rows if r["model"] == "DPA4C-Plus"]
    oree_oth = [r["OH_REE_shell_mean"] for r in main_rows if r["model"] != "DPA4C-Plus"]
    dg = {m: [r["dG_act_300K_kJmol"] for r in where(d.act, system="graphite", model="DPA4C-Plus", metal=m)]
          for m in METALS}
    ea = [r["Ea_kJmol"] for r in d.act if r["Ea_quality"] in ("good", "fair")]
    pka1 = {m: first(d.stepw, metal=m, hydrolysis_step="1", system="graphite", model="DPA4C-Plus") for m in METALS}
    t1 = {m: first(d.t1, Ln=m) for m in METALS}
    pmf = {m: first(d.pmf, metal=m, source=lambda s: s.startswith("Shiery2021 SI PMF")) for m in METALS}
    cards = [
        (f"{n_tot} 条轨迹", f"石墨 {n_gra} + 溶液 {n_tot - n_gra}；坍塌筛查排除 {status['excluded']} 条、截断 {status['truncated']} 条"),
        (f"Dy {fmt(gd['REE_O_total_mean'])} · Nd {fmt(gn['REE_O_total_mean'])}", "石墨 DPA4C-Plus：CN(REE–O 总)，300 K"),
        (f"Dy {rng(sol_tot['Dy'])} · Nd {rng(sol_tot['Nd'])}", "溶液 4 个势函数的 CN(REE–O 总) 范围，300 K"),
        (f"{rng(oree_dpa)} | {rng(oree_oth)}", "OH 中心 O–REE（壳层截断）：DPA4C-Plus | 其它势函数"),
        (f"Dy {rng(dg['Dy'], 1)} · Nd {rng(dg['Nd'], 1)}", "石墨 DPA4C-Plus 三步表观 ΔG<sup>‡</sup> (kJ/mol, 300 K, 6 个方向)"),
        (f"{rng(ea, 1)}", "溶液 Arrhenius 表观 E<sub>a</sub> 范围 (kJ/mol；拟合质量“好/尚可”)"),
        (f"Nd {fmt(pka1['Nd']['pKa'], 1)} vs {fmt(t1['Nd']['pKa1_calc'], 1)} · Dy {fmt(pka1['Dy']['pKa'], 1)} vs {fmt(t1['Dy']['pKa1_calc'], 1)}",
         "pKa<sub>1</sub>：石墨 DPA4C-Plus (MD 等效) vs 论文 Table 1"),
        (f"Nd {fmt(pmf['Nd']['value_kJmol'], 1)} · Dy {fmt(pmf['Dy']['value_kJmol'], 1)}",
         "论文第一步 PMF 终点 (kJ/mol)：单调、无能垒 → 论文没有可比的活化能"),
    ]
    return '<section class="metrics">' + "".join(
        f'<div class="metric"><b>{v}</b><span>{t}</span></div>' for v, t in cards) + "</section>"


def notice(d) -> str:
    c1 = d.kmeta["c_single_M"]
    return f"""<div class="notice"><strong>阅读前请注意</strong><ul>
<li>速率与 ΔG<sup>‡</sup> 是 committed-transition 的<b>表观</b>值，依赖过滤窗口 W 和截断 r<sub>c</sub>（见 4.4 敏感性）；三步之间的相对顺序稳定，绝对值会移动几个 kJ/mol。</li>
<li>pKa“等效值”假设 pK<sub>w</sub> = {fmt(d.kmeta['pKw'], 0)}，并把 box 中的 OH<sup>−</sup> 当理想溶液（单分子浓度 c<sub>1</sub> = {fmt(c1['solution'], 3)} M 溶液 / {fmt(c1['graphite'], 3)} M 石墨）；实际是高浓度、有离子对，不是无限稀释。</li>
<li>溶液每个温度只有 1 条 150 ps 轨迹；石墨 DPA4C-Plus 在 ~90 ps 前仍有起始构型记忆；UMA 大量轨迹坍塌（截断/排除），剩余轨迹中 OH 几乎不离开 Ln（非遍历）。</li>
<li>论文的 pKa<sub>1</sub> 来自约束 AIMD（显式水）的 PMF 终点，pKa<sub>2</sub>、pKa<sub>3</sub> 来自团簇 + 隐式溶剂（M06/CPCM）相对第 1 步的自由能差，都是 298 K 热力学量，<b>没有</b>给出能垒/活化能；与论文能比较的是 ΔG、pKa、结构，而不是 E<sub>a</sub>。</li>
</ul></div>"""


def sec_summary(d) -> str:
    gd, gn = cn_get(d, "DPA4C-Plus", "Dy", "ALL"), cn_get(d, "DPA4C-Plus", "Nd", "ALL")
    sd, sn = cn_get(d, "DPA4C-Plus", "Dy", "box_pbc"), cn_get(d, "DPA4C-Plus", "Nd", "box_pbc")
    ud, un = cn_get(d, "uma-s-1.2.1", "Dy", "ALL"), cn_get(d, "uma-s-1.2.1", "Nd", "ALL")
    main_rows = [r for r in d.cn if r["variant"] in ("ALL", "box_pbc")]
    oth_sol = [r for r in where(d.cn, variant="box_pbc") if r["model"] != "DPA4C-Plus"]
    rc_ree = [r["rc_REE_O_A"] for r in main_rows if finite(r["REE_O_total_mean"])]
    cuts = [cut_get(d, s, m, metal) for s, m in SERIES for metal in METALS]
    cuts = [c for c in cuts if c]
    odd = cn_get(d, "DPA4C-Plus", "Nd", "box_pbc")
    odd_rest = num(odd["OH_Xtotal_shell_mean"]) - num(odd["OH_REE_shell_mean"]) - num(odd["OH_Hwater_shell_mean"])
    ox_norm = [r for r in main_rows if r is not odd and finite(r["OH_Xtotal_shell_mean"])]
    ox_dpa = [r for r in ox_norm if r["model"] == "DPA4C-Plus"]
    ox_oth = [r for r in ox_norm if r["model"] != "DPA4C-Plus"]

    def shell(rows, key):
        return rng([r[key] for r in rows])
    nd_drop = num(cut_get(d, "graphite", "DPA4C-Plus", "Nd")["OH_X_shell_drop_fraction"]) * 100

    def a(metal, step, direction):
        return num(act_get(d, "graphite", "DPA4C-Plus", metal, step, direction)["dG_act_300K_kJmol"])

    fwd = {m: [a(m, k, "association") for k in (1, 2, 3)] for m in METALS}
    rev = {m: [a(m, k, "dissociation") for k in (1, 2, 3)] for m in METALS}
    gact = where(d.act, system="graphite", model="DPA4C-Plus")
    n_ev = [inum(r["N_events_300K"]) for r in gact]
    ptl = [num(r["frac_proton_transfer_like"]) for r in gact if finite(r["frac_proton_transfer_like"])]
    diff = [fwd["Dy"][i] - fwd["Nd"][i] for i in range(3)] + [rev["Dy"][i] - rev["Nd"][i] for i in range(3)]
    ea_rows = [r for r in d.act if r["Ea_quality"] in ("good", "fair")]
    nd21 = act_get(d, "solution", "DPA4C-Plus", "Nd", 2, "dissociation")
    vh1 = first(d.vth, system="solution", model="DPA4C-Plus", metal="Nd", hydrolysis_step="1")
    pmf = {m: first(d.pmf, metal=m, source=lambda s: s.startswith("Shiery2021 SI PMF")) for m in METALS}
    t1 = {m: first(d.t1, Ln=m) for m in METALS}
    md = [r for r in d.stepw if r["source"].startswith("this work") and r["status"] in ("measured", "transient")]
    md_dpa = [r for r in md if r["model"] == "DPA4C-Plus"]
    ddg1 = [r["ddG_vs_calc_kJmol_298K"] for r in md_dpa if r["hydrolysis_step"] == "1"]
    ddg3 = [r["ddG_vs_calc_kJmol_298K"] for r in md_dpa if r["hydrolysis_step"] == "3"]
    ddg3_oth = [r["ddG_vs_calc_kJmol_298K"] for r in md if r["hydrolysis_step"] == "3" and r["model"] != "DPA4C-Plus"]
    agree = Counter(r["agrees_within_0.05"] for r in d.madchk)
    paper = {m: "/".join(fmt(t1[m][f"pKa{k}_calc"], 1) for k in (1, 2, 3)) for m in METALS}
    items = [
        f"""<b>(a) REE–O 配位（300 K，box PBC）：</b>DPA4C-Plus 在石墨上 Dy/Nd 的 REE–O(总) = {fmt(gd['REE_O_total_mean'])}/{fmt(gn['REE_O_total_mean'])}
（水 {fmt(gd['REE_O_water_mean'])}/{fmt(gn['REE_O_water_mean'])}，OH {fmt(gd['REE_O_hydroxide_mean'])}/{fmt(gn['REE_O_hydroxide_mean'])}）；
在溶液中 Dy = {fmt(sd['REE_O_total_mean'])}（水 {fmt(sd['REE_O_water_mean'])}，OH {fmt(sd['REE_O_hydroxide_mean'])}），Nd = {fmt(sn['REE_O_total_mean'])}（水 {fmt(sn['REE_O_water_mean'])}，OH {fmt(sn['REE_O_hydroxide_mean'])}）。
DPA4-Pro、MACE-POLAR-1-M、UMA 在溶液中 REE–O(OH) = {rng([r['REE_O_hydroxide_mean'] for r in oth_sol])}，即 3 个 OH 始终全部配位，水相应减少，总 CN 更低
（Dy {rng([r['REE_O_total_mean'] for r in oth_sol if r['metal'] == 'Dy'])}，Nd {rng([r['REE_O_total_mean'] for r in oth_sol if r['metal'] == 'Nd'])}；UMA 石墨 Dy/Nd = {fmt(ud['REE_O_total_mean'])}/{fmt(un['REE_O_total_mean'])}）。截断取总 REE–O g(r) 第一极小（{rng(rc_ree)} Å）。""",
        f"""<b>(b) OH 中心配位：</b>O–X(总) g(r) 的第一个极小（{rng([c['OH_X_first_min_A'] for c in cuts])} Å）在 O–REE 峰（{rng([c['OH_REE_peak_A'] for c in cuts])} Å）之前，
所以“严格第一极小”下 O–REE 恒为 0，只数到氢键给体 H（DPA4C-Plus 在 1.23–1.25 Å 还有一个共享质子峰）。因此同时给出“含 REE 的壳层”截断（O–REE 峰之后的极小，{rng([c['OH_X_shell_min_A'] for c in cuts])} Å）：
DPA4-Pro、MACE、UMA 的 O–X(总) = {shell(ox_oth, 'OH_Xtotal_shell_mean')}，即 O–REE {shell(ox_oth, 'OH_REE_shell_mean')} + O–H(水) {shell(ox_oth, 'OH_Hwater_shell_mean')}（≈ 1 个 REE + 2 个多水 H）；
DPA4C-Plus 的 O–X(总) = {shell(ox_dpa, 'OH_Xtotal_shell_mean')} 相近，但组成不同：O–REE 只有 {shell(ox_dpa, 'OH_REE_shell_mean')}（OH 多数不在 REE 上），O–H(水) 更多（{shell(ox_dpa, 'OH_Hwater_shell_mean')}）。
DPA4C-Plus Nd 石墨的壳层极小几乎不存在（峰→谷降幅 {nd_drop:.1f}%），DPA4C-Plus Nd 溶液的壳层峰其实是 O···O(水)（O–X = {fmt(odd['OH_Xtotal_shell_mean'])}，其中 O–REE {fmt(odd['OH_REE_shell_mean'])}，含 ≈ {fmt(odd_rest, 1)} 个水 O）——这两项只作诊断。""",
        f"""<b>三步表观活化能（Eyring，300 K）：</b>统计最好的是石墨 DPA4C-Plus（每个方向 {min(n_ev)}–{max(n_ev)} 个 committed 事件）。
水解正向（Ln 得 OH<sup>−</sup>，(k−1)→k）第 1/2/3 步：Dy {'/'.join(fmt(x, 1) for x in fwd['Dy'])}，Nd {'/'.join(fmt(x, 1) for x in fwd['Nd'])} kJ/mol；
逆向（Ln 失 OH<sup>−</sup>，k→(k−1)）：Dy {'/'.join(fmt(x, 1) for x in rev['Dy'])}，Nd {'/'.join(fmt(x, 1) for x in rev['Nd'])} kJ/mol。
第 1 步正向最快、逆向最慢，第 3 步相反；Nd 比 Dy 低 {rng(diff, 1)} kJ/mol。≥ {math.floor(min(ptl) * 100)}% 的事件是质子转移型（OH 的 O 留在 Ln 上，只交换 H）。""",
        f"""<b>温度依赖（溶液，300–500 K）：</b>只有 DPA4C-Plus（和 MACE Nd 第 3 步）事件足够做 Arrhenius；拟合质量为“好/尚可”的表观 E<sub>a</sub> = {rng([r['Ea_kJmol'] for r in ea_rows], 1)} kJ/mol，95% CI 宽。
失 OH 方向出现负 E<sub>a</sub>（如 Nd 2→1：{fmt(nd21['Ea_kJmol'])} {ci(nd21['Ea_CI95_lo'], nd21['Ea_CI95_hi'])}），与 van 't Hoff ΔH°<sub>d</sub> &lt; 0（Nd 第 1 步 {fmt(vh1['dH0_dissociation_kJmol'])} kJ/mol）一致：升温使 OH 更稳定地留在 Ln 上。
这些是复合的表观速率常数，不是单一基元步骤的能垒；DPA4-Pro、MACE、UMA 在溶液中 OH 几乎从不离开 Ln，只能给 ΔG<sup>‡</sup> 下限。""",
        f"""<b>与 Shiery 2021 对比：</b>论文只有热力学（pKa、ΔG，298 K）。SI 第一步水解 PMF 单调上升（终点 Nd {fmt(pmf['Nd']['value_kJmol'], 1)}、Dy {fmt(pmf['Dy']['value_kJmol'], 1)} kJ/mol ≈ Table 1 ΔG<sub>1</sub> {fmt(t1['Nd']['dG1_kJmol_298K'], 1)}/{fmt(t1['Dy']['dG1_kJmol_298K'], 1)}），没有能垒，因此活化能无法与论文直接比较。
DPA4C-Plus 的 pKa 等效值（pK<sub>w</sub> = 14）三步都落在 {rng([r['pKa'] for r in md_dpa], 1)}，三步几乎没有逐级分离（论文 Nd {paper['Nd']}，Dy {paper['Dy']}）：
pKa<sub>1</sub> 比论文高（ΔΔG<sub>1</sub> = {rng(ddg1, 1)} kJ/mol），pKa<sub>3</sub> 比论文低（ΔΔG<sub>3</sub> = {rng(ddg3, 1)} kJ/mol；其它势函数第 3 步 {rng(ddg3_oth, 1)} kJ/mol）。
其它势函数在溶液中 OH 基本不离开 Ln，第 1、2 步多为“无法确定”或只有界限。""",
        f"""<b>Table 2 的读法：</b>Nd 的 Kragten 1984（ref 99）8.1/16.2/24.3 = 8.1 × n，看起来是累积常数 −log*β<sub>n</sub>；第 2、3 步的实验中位数在数值上也更像累积值。
论文 MAD 是把它们当逐步 pKa 计算的（我们用 Table 2 重算 MAD：{agree['yes']}/{agree['yes'] + agree['NO']} 与 Table 1 一致，另有 {agree['']} 项无实验值）。第 6 节同时给出两种读法。""",
    ]
    body = '<ul class="tight">' + "".join(f"<li>{x}</li>" for x in items) + "</ul>"
    return card("summary", "1 结论摘要", "所有数字都直接从下面各节的 CSV 读出", body)


# ------------------------------------------------------------------ 2 methods
def median(values) -> float:
    vals = sorted(x for x in (num(v) for v in values) if math.isfinite(x))
    if not vals:
        return math.nan
    mid = len(vals) // 2
    return vals[mid] if len(vals) % 2 else 0.5 * (vals[mid - 1] + vals[mid])


def sol300(d, model: str, metal: str) -> dict | None:
    return first(d.sol_all, model=model, metal=metal, T_K="300")


def sol_screen(d, s: dict) -> dict | None:
    return first(d.screen, system="solution", model=s["model"], metal=s["metal"], T_K=s["T_K"])


def series_rows():
    return [(system, model, metal) for system, model in SERIES for metal in METALS]


def lead_cells(system: str, model: str, metal: str) -> list[str]:
    return [SYSTEM_LABEL[system], MODEL_LABEL[model], metal_html(metal)]


def sec_methods(d) -> str:
    groups = defaultdict(list)
    for r in d.screen:
        groups[(r["system"], r["model"], r["metal"])].append(r)
    rows = []
    for system, model, metal in series_rows():
        g = groups.get((system, model, metal), [])
        st = Counter(x["status"] for x in g)
        temps = sorted({inum(x["T_K"]) for x in g})
        starts = sorted({x["variant"] for x in g})
        layout = (f"{len(starts)} 种起始 × {len(g) // max(len(starts), 1)}" if system == "graphite"
                  else "每个 T 1 条")
        usable = [x["usable_end_ps"] for x in g if x["status"] != "excluded"]
        rows.append({"cells": lead_cells(system, model, metal) + [
            str(len(g)), "/".join(map(str, temps)), layout, str(st["ok"]), str(st["truncated"]), str(st["excluded"]),
            fmt(median(usable), 0) if usable else "—"],
            "cls": "weak" if st["excluded"] == len(g) else ""})
    systems = table(["体系", "势函数", "金属", "轨迹", "T (K)", "布局", "正常", "截断", "排除", "可用时长中位数 (ps)"],
                    rows, num_cols=(3, 6, 7, 8, 9))
    n_used = sum(1 for r in d.screen if r["status"] != "excluded")
    g_other = max(num(r["REE_O_other_mean"]) for r in d.cn if r["variant"] != "box_pbc")
    s_other = max(abs(num(r["REE_Ototal_mean"]) - num(r["REE_Ow_mean"]) - num(r["REE_OOH_mean"]))
                  for r in d.sol_all if r["collapsed"] == "0")
    sol_trunc = [s for s in d.sol_all if s["collapsed"] != "0" and (sc := sol_screen(d, s)) and sc["status"] == "truncated"]
    trunc_txt = "；".join(f"{MODEL_LABEL[s['model']]} {s['metal']} {s['T_K']} K 只有动力学结果（用前 {fmt(sol_screen(d, s)['usable_end_ps'], 0)} ps）"
                          for s in sol_trunc)
    body = f"""
<h3>2.1 体系与轨迹</h3>
<p>所有轨迹都是三维周期 box（PBC），每条 150 ps：Ln(OH)<sub>3</sub> + 水，石墨体系另含一层石墨。石墨只有 DPA4C-Plus 和 UMA 的轨迹，300 K，4 种 OH 起始构型（3near / 2near1far / 1near2far / 3far，指 3 个 OH 相对石墨的远近）× 50 条；
溶液四个势函数各 1 条 / 温度（300–500 K）。共 {len(d.screen)} 条，坍塌筛查后可用 {n_used} 条。</p>
{systems}
{data_line("collapse_screen.csv")}
<h3>2.2 配位数定义（任务 a、b）</h3>
<ul class="tight">
<li><b>O 物种：</b>每个 H 按最小镜像距离归属最近的 O；有 2 个 H 的 O = 水，1 个 H = OH<sup>−</sup>，0 或 ≥ 3 个 = 其它（O<sup>2−</sup>/H<sub>3</sub>O<sup>+</sup>）。逐帧重新归属，所以质子转移后 OH/水 的身份会交换。REE–O(其它)：石墨体系为 {'0' if g_other == 0 else fmt(g_other)}；溶液（未坍塌轨迹）中 |REE–O(总) − 水 − OH| ≤ {fmt(s_other, 3)}，可以忽略，所以 REE–O(总) ≈ 水 + OH。</li>
<li><b>(a) REE 中心：</b>截断 r<sub>c</sub> = 总 REE–O g(r) 第一峰之后的极小（g(r) 网格 0.02 Å，5 点平滑）。REE–O(水)、REE–O(OH)、REE–O(总) 都用这一个 r<sub>c</sub> 直接计数。</li>
<li><b>(b) OH 中心：</b>中心 = 每帧所有 OH 的 O；X = 除中心 O 自身和它自己的 H 以外的所有原子（REE、水 O/H、其它 OH 的 O/H，石墨体系还包括 C）。给出 O–REE、O–H(水)、O–X(总)。</li>
<li><b>采样与误差：</b>石墨 t ≥ 10 ps 起每 1 ps 一帧，误差为轨迹间（按帧数加权）标准误；溶液用 <code>md/dense-gr-cn</code> 的逐帧计数（t ≥ 10 ps，每 2 fs，{inum(sol300(d, 'DPA4C-Plus', 'Nd')['n_frames']):,} 帧），误差为 5 段 block 标准误。</li>
</ul>
<h3>2.3 OH 溶剂化壳层不清晰时怎么处理</h3>
<p>OH<sup>−</sup> 的 O–X(总) g(r) 由好几个重叠的峰组成：给 OH 氢键的水 H（1.6–1.8 Å，DPA4C-Plus 还有 1.23–1.25 Å 的共享质子）、REE（2.2–2.3 Å）、水 O（2.4–2.7 Å）。所以“总 g(r) 的第一极小”不唯一，我们的做法是：</p>
<ol class="tight">
<li><b>严格截断</b>（按要求的定义）：O–X(总) 第一个真实峰（g &gt; 0.5 且之后下降 ≥ 10%）之后的第一个极小。它总在 O–REE 峰之前，所以只数到氢键 H，O–REE ≈ 0。</li>
<li><b>壳层截断</b>：从 O–REE g(r) 峰（1.9–2.8 Å 内最大值）出发，取 O–X(总) 随后的第一个极小，使 REE 包含在内。</li>
<li><b>清晰度：</b>峰→谷降幅 1 − g<sub>min</sub>/g<sub>peak</sub> ≥ 10% 记为“清晰”，否则只作诊断。</li>
<li><b>余量检查：</b>壳层内 X − REE − H(水) 应 ≈ 0；如果明显 &gt; 0，说明壳层截断已经把水 O 包进来了。</li>
<li><b>对称参考：</b>REE–O(OH)/3 不依赖 OH 中心截断（用 REE–O 的 r<sub>c</sub>），可以和壳层 O–REE 对照。</li>
</ol>
<h3>2.4 坍塌筛查</h3>
<p>部分 UMA 轨迹会坍塌（势能骤降、H/O 被吸到 Ln 上）。判据：(i) 势能：每 50 个采样点（1 ps）的块平均比前 5 块最低值低 25 eV（石墨）/ 15 eV（溶液），在 −10 eV 处细化起点；(ii) 几何：Ln–O &lt; 1.85 Å 或 Ln–H &lt; 1.90 Å 连续 10 个采样点。
轨迹截到坍塌前 1 ps；剩余 ≥ 20 ps 记为“截断”并继续使用，否则“排除”。共排除 {Counter(r['status'] for r in d.screen)['excluded']} 条（UMA Dy 溶液 5 条全部排除）。石墨 CN/g(r) 和全部动力学都只用可用时间段；溶液 CN/g(r) 取自 <code>md/dense-gr-cn</code>，它对坍塌的轨迹不给 CN{('，所以 ' + trunc_txt) if trunc_txt else ''}。</p>
"""
    return card("methods", "2 方法", "CN/g(r)：analyze.py · 动力学：kinetics.py · 文献：literature.py · 对比：compare.py", body)


# ------------------------------------------------------------------ 3 CN
def sample_label(d, system: str, model: str, metal: str, row: dict) -> str:
    if system == "graphite":
        extra = []
        if inum(row["n_truncated"]):
            extra.append(f"截断 {inum(row['n_truncated'])}")
        if inum(row["n_excluded"]):
            extra.append(f"排除 {inum(row['n_excluded'])}")
        tail = f"（{'，'.join(extra)}）" if extra else ""
        return f"{inum(row['n_replicas'])} 条{tail}<small>{inum(row['n_samples']):,} 帧 × 1 ps</small>"
    s = sol300(d, model, metal)
    return f"1 条<small>{inum(s['n_frames']):,} 帧 × 2 fs</small>"


def collapsed_cells(n_cols: int) -> tuple[str, str]:
    return ("轨迹坍塌（5/5 条被排除），无可用数据", f'colspan="{n_cols}" class="weak"')


def shell_flags(c: dict) -> list[str]:
    flags = []
    if c["OH_X_shell_min_clear"] != "True":
        flags.append(f"壳层极小不清晰（降幅 {fmt(num(c['OH_X_shell_drop_fraction']) * 100, 1)}%）")
    if abs(num(c["OH_X_shell_peak_A"]) - num(c["OH_REE_peak_A"])) > 0.05:
        flags.append(f"壳层峰 {fmt(c['OH_X_shell_peak_A'])} Å ≠ O–REE 峰 {fmt(c['OH_REE_peak_A'])} Å（是 O···O 水峰；O–REE 峰高仅 g = {fmt(c['OH_REE_peak_g'], 1)}）")
    return flags


def sec_cn(d) -> str:
    # (a) REE-centred
    rows_a = []
    for system, model, metal in series_rows():
        r = cn_get(d, model, metal, series_variant(system))
        lead = lead_cells(system, model, metal)
        if r is None or r["collapsed"] == "1" or not finite(r["REE_O_total_mean"]):
            rows_a.append({"cells": lead + [collapsed_cells(5)], "cls": "weak"})
            continue
        rows_a.append(lead + [fmt(r["rc_REE_O_A"]), pm(r["REE_O_water_mean"], r["REE_O_water_sem"]),
                              pm(r["REE_O_hydroxide_mean"], r["REE_O_hydroxide_sem"]),
                              "<b>" + pm(r["REE_O_total_mean"], r["REE_O_total_sem"]) + "</b>",
                              sample_label(d, system, model, metal, r)])
    tab_a = table(["体系", "势函数", "金属", "r<sub>c</sub> (Å)", "REE–O (水)", "REE–O (OH)", "REE–O (总)", "样本"],
                  rows_a, num_cols=(3, 4, 5, 6))

    # (b) OH-centred: cutoffs and clarity
    def drop(c, key):
        ok = c[f"OH_X_{key}_min_clear"] == "True"
        return f"{fmt(num(c[f'OH_X_{key}_drop_fraction']) * 100, 1)}% " + (pill("清晰", "ok") if ok else pill("不清晰", "bad"))

    rows_c, cuts = [], []
    for system, model, metal in series_rows():
        c = cut_get(d, system, model, metal)
        lead = lead_cells(system, model, metal)
        if c is None:
            rows_c.append({"cells": lead + [collapsed_cells(8)], "cls": "weak"})
            continue
        cuts.append(c)
        flags = shell_flags(c)
        rows_c.append({"cells": lead + [
            fmt(c["OH_X_first_peak_A"]), fmt(c["OH_X_first_min_A"]), (drop(c, "first"), 'class="num"'),
            f"{fmt(c['OH_REE_peak_A'])}<small>g = {fmt(c['OH_REE_peak_g'], 0)}</small>",
            fmt(c["OH_X_shell_peak_A"]), fmt(c["OH_X_shell_min_A"]), (drop(c, "shell"), 'class="num"'),
            "；".join(flags) + " ‡" if flags else pill("两种截断都可用", "ok")], "cls": "hl" if flags else ""})
    head_c = thead2(["体系", "势函数", "金属"],
                    [("严格：O–X(总) 第一极小", ["峰 (Å)", "极小 (Å)", "降幅"]),
                     ("O–REE g(r)", ["峰 (Å)"]),
                     ("壳层：O–REE 峰之后的 O–X 极小", ["峰 (Å)", "极小 (Å)", "降幅"]),
                     ("判定", [""])])
    tab_c = table([], rows_c, num_cols=(3, 4, 6, 7), thead=head_c)

    # (b) OH-centred: counts
    rows_b, oree_first, gr_diff = [], [], []
    for system, model, metal in series_rows():
        r = cn_get(d, model, metal, series_variant(system))
        lead = lead_cells(system, model, metal)
        if r is None or r["collapsed"] == "1" or not finite(r["OH_Xtotal_shell_mean"]):
            rows_b.append({"cells": lead + [collapsed_cells(10)], "cls": "weak"})
            continue
        dag = " †" if system == "solution" else ""
        rest = num(r["OH_Xtotal_shell_mean"]) - num(r["OH_REE_shell_mean"]) - num(r["OH_Hwater_shell_mean"])
        oree_first.append(num(r["OH_REE_first_mean"]))
        flags = shell_flags(cut_get(d, system, model, metal))
        rows_b.append({"cells": lead + [
            fmt(r["rc_OX_first_A"]), pm(r["OH_REE_first_mean"], r["OH_REE_first_sem"], 3) + dag,
            pm(r["OH_Hwater_first_mean"], r["OH_Hwater_first_sem"]) + dag,
            pm(r["OH_Xtotal_first_mean"], r["OH_Xtotal_first_sem"]),
            fmt(r["rc_OX_shell_A"]), "<b>" + pm(r["OH_REE_shell_mean"], r["OH_REE_shell_sem"]) + "</b>",
            pm(r["OH_Hwater_shell_mean"], r["OH_Hwater_shell_sem"]),
            "<b>" + pm(r["OH_Xtotal_shell_mean"], r["OH_Xtotal_shell_sem"]) + "</b>" + (" ‡" if flags else ""),
            fmt(rest), fmt(num(r["REE_O_hydroxide_mean"]) / 3)], "cls": "hl" if flags else ""})
    for s in d.sol_all:
        if s["collapsed"] != "0":
            continue
        for direct, via_gr in (("OH_Xtotal_mean", "OH_Xtotal_first_from_gr"), ("OH_REE_shell_mean", "OH_REE_shell_from_gr"),
                               ("OH_Hw_shell_mean", "OH_Hwater_shell_from_gr"), ("OH_Xtotal_shell_mean", "OH_Xtotal_shell_from_gr")):
            if finite(s[direct]) and finite(s[via_gr]):
                gr_diff.append(abs(num(s[direct]) - num(s[via_gr])))
    head_b = thead2(["体系", "势函数", "金属"],
                    [("严格截断（第一极小）", ["r<sub>c</sub> (Å)", "O–REE", "O–H(水)", "O–X(总)"]),
                     ("壳层截断（含 REE）", ["r<sub>c</sub> (Å)", "O–REE", "O–H(水)", "O–X(总)", "余量 X−REE−H<sub>w</sub>"]),
                     ("对称参考", ["REE–O(OH)/3"])])
    tab_b = table([], rows_b, num_cols=tuple(range(3, 13)), thead=head_b)
    odd = cn_get(d, "DPA4C-Plus", "Nd", "box_pbc")
    dpa_first_peak = [c["OH_X_first_peak_A"] for c in cuts if c["model"] == "DPA4C-Plus"]
    notes = f"""<ul class="tight muted">
<li>严格第一极小（{rng([c['OH_X_first_min_A'] for c in cuts])} Å）都在 O–REE 峰（{rng([c['OH_REE_peak_A'] for c in cuts])} Å）之前，所以严格截断下 O–REE ≤ {fmt(max(oree_first), 3)}：OH 的“第一层”只有给它氢键的水 H（DPA4C-Plus 的第一峰在 {rng(dpa_first_peak)} Å，是共享质子）。</li>
<li>壳层截断把 REE 包进来：DPA4-Pro / MACE / UMA 的 O–REE ≈ 1（每个 OH 都在 Ln 上），DPA4C-Plus 只有 {rng([r['OH_REE_shell_mean'] for r in d.cn if r['model'] == 'DPA4C-Plus' and r['variant'] in ('ALL', 'box_pbc')])}（OH 大部分时间不在 Ln 上）。
壳层截断略小于 REE–O 的 r<sub>c</sub>，所以 O–REE(壳层) ≤ REE–O(OH)/3，差值就是落在两者之间的 OH。</li>
<li>† 溶液严格截断下的 O–REE、O–H(水) 由 g(r) 积分得到（无误差棒）；其余溶液数值是直接计数（5 段 block SEM），直接计数与 g(r) 积分最多相差 {fmt(max(gr_diff), 3)}。</li>
<li>‡ 标黄的两行只作诊断：DPA4C-Plus Nd 石墨的壳层极小几乎不存在；DPA4C-Plus Nd 溶液的“壳层”其实到了 O···O(水) 峰之后，余量 ≈ {fmt(num(odd['OH_Xtotal_shell_mean']) - num(odd['OH_REE_shell_mean']) - num(odd['OH_Hwater_shell_mean']))} 个水 O。对这两行建议用 O–REE(壳层) 与 REE–O(OH)/3 来描述 OH 是否在 Ln 上，不用 O–X(总)。</li>
</ul>"""
    body = f"""
<h3>3.1 (a) REE 中心：REE–O(水)、REE–O(OH)、REE–O(总)，300 K</h3>
{tab_a}
{data_line("cn_summary_300K.csv", "graphite_cn_summary_box_pbc.csv", "graphite_cn_by_replica_box_pbc.csv", "solution_cn_box_pbc.csv")}
<h3>3.2 (b) OH 中心：截断位置与清晰度</h3>
{tab_c}
{data_line("graphite_cutoffs_box_pbc.csv", "solution_cutoffs_box_pbc_300K.csv")}
<h3>3.3 (b) OH 中心：O–REE、O–H(水)、O–X(总)，300 K</h3>
{tab_b}
{notes}
{data_line("cn_summary_300K.csv", "solution_cn_box_pbc.csv")}
{details("石墨：按 4 种起始构型分开的 CN", graphite_start_table(d))}
{details("溶液：300–500 K 全部温度的 CN", solution_all_t_table(d))}
<h3>3.4 图</h3>
<div class="grid2">
{figure("graphite_REE_O_total_box_pbc_300K.svg", "石墨：REE–O g(r)（水 / OH / 总）与 r<sub>c</sub>", ("graphite_gr_box_pbc.csv", "graphite_cutoffs_box_pbc.csv"))}
{figure("solution_REE_O_total_box_pbc_300K.svg", "溶液：REE–O g(r)（水 / OH / 总）与 r<sub>c</sub>", ("solution_gr_box_pbc_300K.csv", "solution_cutoffs_box_pbc_300K.csv"))}
{figure("graphite_OH_X_total_box_pbc_300K.svg", "石墨：OH 中心 O–X(总)、O–REE、O–H(水) g(r) 与两种截断", ("graphite_gr_box_pbc.csv", "graphite_cutoffs_box_pbc.csv"))}
{figure("solution_OH_X_total_box_pbc_300K.svg", "溶液：OH 中心 O–X(总)、O–REE、O–H(水) g(r) 与两种截断", ("solution_gr_box_pbc_300K.csv", "solution_cutoffs_box_pbc_300K.csv"))}
{figure("graphite_CN_box_pbc_300K.svg", "石墨：各配位数汇总", ("graphite_cn_summary_box_pbc.csv",))}
{figure("solution_CN_box_pbc_300K.svg", "溶液：各配位数汇总", ("solution_cn_box_pbc.csv",))}
</div>
<p class="files">石墨累积 g(r) 直方图（可重新归一化）：{dl("graphite_rdf_accumulated.npz")} · {dl("graphite_rdf_accumulated.json")}</p>
"""
    return card("cn", "3 配位数（box PBC，300 K）", "任务 (a) REE 中心与 (b) OH 中心；石墨 = 全部起始构型合并 (ALL)，溶液 = 300 K 轨迹", body)


def graphite_start_table(d) -> str:
    rows = []
    for model in ("DPA4C-Plus", "uma-s-1.2.1"):
        for metal in METALS:
            for variant in ("ALL",) + GRAPHITE_STARTS:
                r = cn_get(d, model, metal, variant)
                if r is None:
                    continue
                rows.append({"cells": [MODEL_LABEL[model], metal_html(metal), "全部合并" if variant == "ALL" else variant,
                                       sample_label(d, "graphite", model, metal, r),
                                       pm(r["REE_O_water_mean"], r["REE_O_water_sem"]),
                                       pm(r["REE_O_hydroxide_mean"], r["REE_O_hydroxide_sem"]),
                                       pm(r["REE_O_total_mean"], r["REE_O_total_sem"]),
                                       pm(r["OH_REE_shell_mean"], r["OH_REE_shell_sem"]),
                                       pm(r["OH_Hwater_shell_mean"], r["OH_Hwater_shell_sem"]),
                                       pm(r["OH_Xtotal_shell_mean"], r["OH_Xtotal_shell_sem"])],
                             "cls": "paper" if variant == "ALL" else ""})
    return table(["势函数", "金属", "起始构型", "样本", "REE–O(水)", "REE–O(OH)", "REE–O(总)",
                  "O–REE 壳层", "O–H(水) 壳层", "O–X 壳层"], rows, num_cols=tuple(range(4, 10)), cls="compact") + \
        data_line("cn_summary_300K.csv", "graphite_cn_by_replica_box_pbc.csv") + \
        '<p class="muted">截断用全部合并 (ALL) 的 g(r) 确定，各起始构型共用。</p>'


def solution_all_t_table(d) -> str:
    rows = []
    for model in MODEL_ORDER:
        for metal in METALS:
            for s in sorted(where(d.sol_all, model=model, metal=metal), key=lambda x: num(x["T_K"])):
                lead = [MODEL_LABEL[model], metal_html(metal), s["T_K"]]
                if s["collapsed"] != "0":
                    sc = sol_screen(d, s)
                    msg = (f"{fmt(sc['usable_end_ps'], 0)} ps 后坍塌：无 CN（动力学只用前 {fmt(sc['usable_end_ps'], 0)} ps）"
                           if sc and sc["status"] == "truncated" else "坍塌（排除）")
                    rows.append({"cells": lead + [(msg, 'colspan="9" class="weak"')], "cls": "weak"})
                    continue
                rows.append(lead + [
                    f"{fmt(s['rc_Ln_O_total_A'])} / {fmt(s['rc_OH_X_first_A'])} / {fmt(s['rc_OH_X_shell_A'])}",
                    pm(s["REE_Ow_mean"], s["REE_Ow_sem"]), pm(s["REE_OOH_mean"], s["REE_OOH_sem"]),
                    pm(s["REE_Ototal_mean"], s["REE_Ototal_sem"]), pm(s["OH_Xtotal_mean"], s["OH_Xtotal_sem"]),
                    pm(s["OH_REE_shell_mean"], s["OH_REE_shell_sem"]), pm(s["OH_Hw_shell_mean"], s["OH_Hw_shell_sem"]),
                    pm(s["OH_Xtotal_shell_mean"], s["OH_Xtotal_shell_sem"]), fmt(s["n_OH_mean"])])
    return table(["势函数", "金属", "T (K)", "r<sub>c</sub>: REE–O / 严格 / 壳层 (Å)", "REE–O(水)", "REE–O(OH)", "REE–O(总)",
                  "O–X 严格", "O–REE 壳层", "O–H(水) 壳层", "O–X 壳层", "n<sub>OH</sub>"],
                 rows, num_cols=tuple(range(3, 12)), cls="compact") + data_line("solution_cn_box_pbc.csv") + \
        '<p class="muted">每个温度单独确定截断；UMA-s-1.2.1 Dy 溶液所有温度都坍塌。</p>'


# ------------------------------------------------------------------ 4 kinetics
def transition(step: str, direction: str) -> str:
    k = int(step)
    return f"{k - 1}→{k}" if direction == "association" else f"{k}→{k - 1}"


def dg_cell(n, dwell, dg, lo, hi, bound, compact: bool = False) -> tuple[str, str]:
    n, dwell = inum(n), num(dwell)
    if not dwell > 0:
        return ("—<small>起始态未出现</small>" if not compact else "—", 'class="num weak"')
    if n == 0:
        tail = f"<small>N = 0，T<sub>origin</sub> = {fmt(dwell, 0)} ps</small>" if not compact else ""
        return (f"&gt; {fmt(bound, 1)}{tail}", 'class="num weak"')
    tail = f"<small>{ci(lo, hi, 1)} · N = {n}</small>" if not compact else f"<small>N = {n}</small>"
    return (f"<b>{fmt(dg)}</b>{tail}", 'class="num weak"' if n < 5 else 'class="num"')


def rate_index(d) -> dict:
    return {(r["scenario"], r["system"], r["model"], r["metal"], r["variant"], inum(r["T_K"]), r["hydrolysis_step"],
             r["direction"]): r for r in d.rates}


def rate_cell(r: dict | None, compact: bool = False) -> tuple[str, str]:
    if r is None:
        return ("—", 'class="num weak"')
    return dg_cell(r["N_events"], r["T_origin_state_ps"], r["dG_act_kJmol"], r["dG_act_CI95_lo"], r["dG_act_CI95_hi"],
                   r["dG_act_lower_bound_if_no_events"], compact)


STEP_GROUPS = [(f"第 {k} 步", [f"正向 {k - 1}→{k}", f"逆向 {k}→{k - 1}"]) for k in (1, 2, 3)]


def sec_kinetics(d) -> str:
    mapping = table(["论文第 k 步", "论文：酸解离（298 K；第 1 步约束 AIMD，第 2、3 步 DFT 团簇 + CPCM）", "MD 正向：Ln 得 OH<sup>−</sup>", "MD 逆向：Ln 失 OH<sup>−</sup>"], [
        ["1", "Ln<sup>3+</sup> + H<sub>2</sub>O → LnOH<sup>2+</sup> + H<sup>+</sup>", "n：0 → 1", "n：1 → 0"],
        ["2", "LnOH<sup>2+</sup> + H<sub>2</sub>O → Ln(OH)<sub>2</sub><sup>+</sup> + H<sup>+</sup>", "n：1 → 2", "n：2 → 1"],
        ["3", "Ln(OH)<sub>2</sub><sup>+</sup> + H<sub>2</sub>O → Ln(OH)<sub>3</sub> + H<sup>+</sup>", "n：2 → 3", "n：3 → 2"]])
    kcut_rows = [lead_cells(r["system"], r["model"], r["metal"]) + [r["n_trajectories"], fmt(r["rc_A"]), fmt(r["first_peak_A"]),
                                                                    fmt(r["second_shell_peak_A"])] for r in d.kcut]
    kcut = table(["体系", "势函数", "金属", "轨迹", "r<sub>c</sub> (Å)", "Ln–O(OH) 第一峰 (Å)", "第二层峰 (Å)"], kcut_rows,
                 num_cols=(3, 4, 5, 6), cls="compact")
    gact = where(d.act, system="graphite", model="DPA4C-Plus")
    ptl = [num(r["frac_proton_transfer_like"]) for r in gact if finite(r["frac_proton_transfer_like"])]
    meta = d.kmeta
    method = f"""
<div class="grid2"><div class="mini"><h4>状态与事件</h4><ul class="tight">
<li>n(t) = Ln 周围 r<sub>c</sub> 内 OH<sup>−</sup> 的 O 的个数（OH = 恰有 1 个最近 H 的 O），上限 3，每 20 fs 一帧。</li>
<li>r<sub>c</sub> = Ln–O(OH) 径向密度第一峰之后的极小（{rng([r['rc_A'] for r in d.kcut])} Å，见下表）。</li>
<li>committed 转变：宽 W = 0.5 ps 的居中众数滤波，再吸收短于 W/2 的片段，去掉 O–H 振动造成的来回跳。</li>
<li>速率从 t ≥ 10 ps 开始统计，占据数从 t ≥ {fmt(meta['thermo_start_ps']['solution'], 0)} ps（溶液）/ {fmt(meta['thermo_start_ps']['graphite'], 0)} ps（石墨）开始。</li>
<li>石墨 DPA4C-Plus 中 {math.floor(min(ptl) * 100)}–{math.ceil(max(ptl) * 100)}% 的事件是“质子转移型”：O 留在 Ln 上，只交换 H（Ln–OH<sub>2</sub> + OH<sup>−</sup> ⇌ Ln–OH + H<sub>2</sub>O），这正是论文酸解离那一步，只是接受质子的是 OH<sup>−</sup> 而不是水。</li>
</ul></div><div class="mini"><h4>速率与表观能垒</h4>
<div class="formula">k = N / T<sub>origin</sub></div>
<div class="formula">ΔG<sup>‡</sup> = RT ln(k<sub>B</sub>T / hk)</div>
<div class="formula">ln k = ln A − E<sub>a</sub>/RT；ln(k/T) = ln(k<sub>B</sub>/h) + ΔS<sup>‡</sup>/R − ΔH<sup>‡</sup>/RT</div>
<ul class="tight"><li>N = 转变次数，T<sub>origin</sub> = 起始态的总停留时间。95% CI = bootstrap（{meta['n_bootstrap']} 次；石墨按轨迹重采样，溶液每个 T 只有 1 条轨迹，按 {fmt(meta['block_ps'], 0)} ps 块重采样）与 Garwood 泊松区间的包络。N 很小时 CI 可能延伸到 ΔG<sup>‡</sup> &lt; 0 或 ∞（重采样里出现 0 个事件），这类数值只有量级意义。</li>
<li>N = 0 时给下限：用泊松上限速率 3.689/T<sub>origin</sub> 代入 Eyring。</li>
<li>Arrhenius 按 N 加权拟合。拟合质量：好 = ≥ 4 个 T、≥ 20 个事件且每个 T ≥ 3 个；尚可 = ≥ 3 个 T 且 ≥ 8 个事件；其余为差；有事件的 T 少于 2 个为“无”。</li></ul>
</div></div>
{details("动力学截断 r<sub>c</sub>（Ln–O(OH) 密度极小）", kcut + data_line("kinetics/kinetics_cutoffs.csv", "kinetics/kinetics_oh_density.csv"))}"""

    # 4.2 dG table
    rows = []
    for metal in METALS:
        for system, model in SERIES:
            lead = [metal_html(metal), SYSTEM_LABEL[system], MODEL_LABEL[model]]
            if not where(d.act, system=system, model=model, metal=metal):
                rows.append({"cells": lead + [collapsed_cells(6)], "cls": "weak"})
                continue
            cells = []
            for k in (1, 2, 3):
                for direction, _ in DIRECTIONS:
                    r = act_get(d, system, model, metal, k, direction)
                    cells.append(dg_cell(r["N_events_300K"], r["dwell_origin_state_ps"], r["dG_act_300K_kJmol"],
                                         r["dG_act_CI95_lo"], r["dG_act_CI95_hi"], r["dG_act_lower_bound_kJmol"]))
            rows.append({"cells": lead + cells, "cls": "hl" if (system, model) == ("graphite", "DPA4C-Plus") else ""})
        paper = [(f"ΔG<sub>{k}</sub> = {fmt(act_get(d, 'graphite', 'DPA4C-Plus', metal, k, 'association')['paper_dG_acid_step_kJmol_298K'])}"
                  "<small>热力学（酸解离），无能垒</small>", 'colspan="2" class="num"') for k in (1, 2, 3)]
        rows.append({"cells": [metal_html(metal), ("Shiery 2021 Table 1（298 K）", 'colspan="2"')] + paper, "cls": "paper"})
    dg_tab = table([], rows, thead=thead2(["金属", "体系", "势函数"], STEP_GROUPS))

    # 4.3 Arrhenius
    arr_rows = []
    for r in d.arr:
        if r["scenario"] != "primary" or r["quality"] == "none":
            continue
        arr_rows.append([MODEL_LABEL[r["model"]], metal_html(r["metal"]), f"第 {r['hydrolysis_step']} 步",
                         ("正向 " if r["direction"] == "association" else "逆向 ") + transition(r["hydrolysis_step"], r["direction"]),
                         r["T_used"].replace(";", "/"), r["N_per_T"].replace(";", "/"), quality_pill(r["quality"]),
                         f"<b>{fmt(r['Ea_kJmol'])}</b><small>{ci(r['Ea_CI95_lo'], r['Ea_CI95_hi'])}</small>",
                         f"{fmt(r['dH_act_kJmol'])}<small>{ci(r['dH_act_CI95_lo'], r['dH_act_CI95_hi'])}</small>",
                         f"{fmt(r['dS_act_JmolK'], 1)}<small>{ci(r['dS_act_CI95_lo'], r['dS_act_CI95_hi'], 1)}</small>",
                         fmt(r["dG_act_300K_from_fit_kJmol"])])
    arr_tab = table(["势函数", "金属", "步", "方向", "T (K)", "每个 T 的 N", "质量", "E<sub>a</sub> (kJ/mol)",
                     "ΔH<sup>‡</sup> (kJ/mol)", "ΔS<sup>‡</sup> (J/mol/K)", "ΔG<sup>‡</sup>(300 K, 拟合)"],
                    arr_rows, num_cols=(7, 8, 9, 10))
    n_none = sum(1 for r in d.arr if r["scenario"] == "primary" and r["quality"] == "none")

    # 4.4 sensitivity (graphite DPA4C)
    idx = rate_index(d)
    sens_html, robust = [], {}
    for metal in METALS:
        def vals(sc):
            return [num(idx[(sc, "graphite", "DPA4C-Plus", metal, "ALL", 300, str(k), dr)]["dG_act_kJmol"])
                    for k in (1, 2, 3) for dr, _ in DIRECTIONS]

        def order(v):
            fwd, rev = v[0::2], v[1::2]
            return (sorted(range(3), key=lambda i: fwd[i]), sorted(range(3), key=lambda i: rev[i]))

        ref = order(vals("primary"))
        srows, n_same = [], 0
        for sc in SCENARIOS:
            same = order(vals(sc)) == ref
            n_same += same
            cells = [rate_cell(idx.get((sc, "graphite", "DPA4C-Plus", metal, "ALL", 300, str(k), dr)), compact=True)
                     for k in (1, 2, 3) for dr, _ in DIRECTIONS]
            srows.append({"cells": [SCENARIO_LABEL[sc]] + cells + [pill("一致", "ok") if same else pill("改变", "bad")],
                          "cls": "hl" if sc == "primary" else ""})
        robust[metal] = n_same
        sens_html.append(f"<h4>{metal_html(metal)} · 石墨 DPA4C-Plus</h4>" +
                         table([], srows, cls="compact", thead=thead2(["分析设定"], STEP_GROUPS + [("三步顺序", ["vs 主分析"])])))
    shifts = []
    for metal in METALS:
        for k in (1, 2, 3):
            for dr, _ in DIRECTIONS:
                v = [num(idx[(sc, "graphite", "DPA4C-Plus", metal, "ALL", 300, str(k), dr)]["dG_act_kJmol"]) for sc in SCENARIOS]
                v = [x for x in v if math.isfinite(x)]
                shifts.append(max(v) - min(v))

    # 4.5 collapsibles: solution Ea vs scenario, DPA4C solution dG per T
    arr_idx = {(r["scenario"], r["model"], r["metal"], r["hydrolysis_step"], r["direction"]): r for r in d.arr}
    keys = sorted({(r["model"], r["metal"], r["hydrolysis_step"], r["direction"]) for r in d.arr if r["quality"] in ("good", "fair")},
                  key=lambda t: (MODEL_ORDER.index(t[0]), t[1], t[2], t[3]))
    ea_rows = []
    for model, metal, k, dr in keys:
        cells = []
        for sc in SCENARIOS:
            r = arr_idx.get((sc, model, metal, k, dr))
            if r is None or r["quality"] == "none":
                cells.append(("—", 'class="num weak"'))
            else:
                cells.append((f"{fmt(r['Ea_kJmol'], 1)}<small>{QUALITY_ZH.get(r['quality'], esc(r['quality']))}</small>",
                              'class="num"' if r["quality"] in ("good", "fair") else 'class="num weak"'))
        ea_rows.append([MODEL_LABEL[model], metal_html(metal), ("正向 " if dr == "association" else "逆向 ") + transition(k, dr)] + cells)
    short = {"primary": "主分析", "W0.2": "W = 0.2 ps", "W1.0": "W = 1.0 ps", "W2.0": "W = 2.0 ps",
             "rc2.90": "r<sub>c</sub> = 2.90 Å", "rc3.30": "r<sub>c</sub> = 3.30 Å",
             "start1": "t ≥ 1 ps", "start50": "t ≥ 50 ps", "start90": "t ≥ 90 ps"}
    ea_tab = table(["势函数", "金属", "过程"] + [short[sc] for sc in SCENARIOS], ea_rows, cls="compact")
    temps = (300, 350, 400, 450, 500)
    pt_rows = []
    for metal in METALS:
        for k in (1, 2, 3):
            for dr, label in DIRECTIONS:
                pt_rows.append([metal_html(metal), f"{label} {transition(str(k), dr)}"] +
                               [rate_cell(idx.get(("primary", "solution", "DPA4C-Plus", metal, "box_pbc", t, str(k), dr)), compact=True)
                                for t in temps])
    pt_tab = table(["金属", "过程"] + [f"{t} K" for t in temps], pt_rows, cls="compact")

    body = f"""
<h3>4.1 步骤对应与方法</h3>
{mapping}
{method}
<h3>4.2 300 K 三步表观 ΔG<sup>‡</sup> (kJ/mol)</h3>
<p class="muted">粗体 = 点估计，小字 = 95% CI 与事件数 N；灰色 = N &lt; 5 或只有下限；“&gt; x”= 没有事件，只给下限。标黄行（石墨 DPA4C-Plus）统计最好。论文行是热力学 ΔG（酸解离自由能），不是能垒，放在这里只为对照量级。</p>
{dg_tab}
{data_line("comparison/compare_activation.csv", "kinetics/kinetics_rates.csv", "kinetics/kinetics_events_by_traj.csv")}
<h3>4.3 表观活化能 E<sub>a</sub>（溶液 300–500 K，Arrhenius）</h3>
<p class="muted">只列出拟合质量不是“无”的 {len(arr_rows)} 项；另外 {n_none} 项事件太少（多为 DPA4-Pro / MACE / UMA，OH 几乎不离开 Ln）。石墨只有 300 K，没有 E<sub>a</sub>。负 E<sub>a</sub> 表示升温反而变慢，说明这是复合的表观速率常数（包含 OH<sup>−</sup> 扩散、氢键网络重排等），不是单一基元步骤。</p>
{arr_tab}
{data_line("kinetics/kinetics_arrhenius.csv", "kinetics/kinetics_rates.csv")}
<h3>4.4 敏感性：石墨 DPA4C-Plus 在 9 种分析设定下的 ΔG<sup>‡</sup></h3>
<p class="muted">三步的相对顺序（正向 1 &lt; 2 &lt; 3，逆向 1 &gt; 2 &gt; 3）在 Dy {robust['Dy']}/{len(SCENARIOS)}、Nd {robust['Nd']}/{len(SCENARIOS)} 种设定下与主分析一致；单个 ΔG<sup>‡</sup> 随设定变化 {fmt(min(shifts), 1)}–{fmt(max(shifts), 1)} kJ/mol（W 越小，记到的快速来回事件越多，ΔG<sup>‡</sup> 越低）。</p>
{''.join(sens_html)}
{data_line("kinetics/kinetics_rates.csv")}
{details("溶液：E<sub>a</sub> 随分析设定的变化（kJ/mol，小字 = 拟合质量）", ea_tab + data_line("kinetics/kinetics_arrhenius.csv"))}
{details("溶液 DPA4C-Plus：各温度的 ΔG<sup>‡</sup>（主分析，kJ/mol）", pt_tab + data_line("kinetics/kinetics_rates.csv"))}
<h3>4.5 图</h3>
{figure("kinetics/fig_kinetics_dGact_300K.svg", "300 K 三步表观 ΔG<sup>‡</sup>：所有势函数、两种体系", ("comparison/compare_activation.csv", "kinetics/kinetics_rates.csv"))}
{figure("kinetics/fig_kinetics_arrhenius_DPA4C_solution.svg", "DPA4C-Plus 溶液：ln k vs 1/T（Arrhenius）", ("kinetics/kinetics_rates.csv", "kinetics/kinetics_arrhenius.csv"))}
{figure("kinetics/fig_kinetics_sensitivity.svg", "ΔG<sup>‡</sup> 对过滤窗口、r<sub>c</sub>、分析起点的敏感性", ("kinetics/kinetics_rates.csv",))}
<div class="grid2">
{figure("kinetics/fig_kinetics_traces.svg", "n(t) 示例：原始（灰）与 committed（彩色）", ("kinetics/kinetics_traces_example.csv",))}
{figure("kinetics/fig_kinetics_OH_density_cutoff.svg", "Ln–O(OH) 径向密度与动力学截断 r<sub>c</sub>", ("kinetics/kinetics_oh_density.csv", "kinetics/kinetics_cutoffs.csv"))}
</div>
"""
    return card("kinetics", "4 三步表观活化能", "committed-transition 速率 → Eyring ΔG<sup>‡</sup>（300 K）与 Arrhenius E<sub>a</sub>（溶液 300–500 K）", body)


# ------------------------------------------------------------------ 5 occupancy
def occ_get(d, scenario: str, system: str, model: str, metal: str, variant: str, temp: int) -> dict | None:
    return first(d.occ, scenario=scenario, system=system, model=model, metal=metal, variant=variant,
                 T_K=lambda t: inum(t) == temp)


def sec_occupancy(d) -> str:
    rows = []
    for system, model, metal in series_rows():
        r = occ_get(d, "primary", system, model, metal, series_variant(system), 300)
        lead = lead_cells(system, model, metal)
        if r is None:
            rows.append({"cells": lead + [collapsed_cells(6)], "cls": "weak"})
            continue
        cells = [f"{fmt_p(r[f'P{n}'])}<small>[{fmt_p(r[f'P{n}_CI95_lo'])}, {fmt_p(r[f'P{n}_CI95_hi'])}]</small>" for n in range(4)]
        rows.append(lead + [r["n_trajectories"]] + cells + [fmt(r["mean_n_OH_bound"])])
    occ_tab = table(["体系", "势函数", "金属", "轨迹", "P<sub>0</sub>", "P<sub>1</sub>", "P<sub>2</sub>", "P<sub>3</sub>", "⟨n⟩"],
                    rows, num_cols=(3, 4, 5, 6, 7, 8))
    c1 = d.kmeta["c_single_M"]
    vrows = []
    for r in d.vth:
        if inum(r["n_T"]) < 2:
            continue
        vrows.append([MODEL_LABEL[r["model"]], metal_html(r["metal"]), f"第 {r['hydrolysis_step']} 步（{r['step']}）",
                      r["T_used"].replace(";", "/"),
                      f"{fmt(r['dH0_dissociation_kJmol'])}<small>{ci(r['dH0_CI95_lo'], r['dH0_CI95_hi'])}</small>",
                      f"{fmt(r['dS0_dissociation_JmolK'], 1)}<small>{ci(r['dS0_CI95_lo'], r['dS0_CI95_hi'], 1)}</small>",
                      fmt(r["dG0_dissociation_300K_fit_kJmol"])])
    vh_tab = table(["势函数", "金属", "步（失 OH 方向）", "T (K)", "ΔH°<sub>d</sub> (kJ/mol)", "ΔS°<sub>d</sub> (J/mol/K)",
                    "ΔG°<sub>d</sub>(300 K, 拟合)"], vrows, num_cols=(4, 5, 6))
    th = {(r["scenario"], r["system"], r["model"], r["metal"], r["hydrolysis_step"]): r for r in d.thermo
          if r["variant"] in ("ALL", "box_pbc") and inum(r["T_K"]) == 300}

    def s90_cell(r, k: str) -> str:
        if r is None:
            return "—"
        p = num(r["pKa_equivalent"])
        if math.isfinite(p):
            return fmt(p) + (" " + pill("瞬态", "warn") if r["both_states_committed"] != "True" else "")
        missing = [str(n) for n, key in ((int(k) - 1, "P_n_minus_1"), (int(k), "P_n")) if not num(r[key]) > 0]
        return f"n = {'、'.join(missing)} 未出现"

    srows = []
    for system in ("graphite", "solution"):
        for metal in METALS:
            for k in ("1", "2", "3"):
                a, b = th.get(("primary", system, "DPA4C-Plus", metal, k)), th.get(("start90", system, "DPA4C-Plus", metal, k))
                pa, pb = num(a["pKa_equivalent"]) if a else math.nan, num(b["pKa_equivalent"]) if b else math.nan
                srows.append([SYSTEM_LABEL[system], metal_html(metal), f"第 {k} 步", s90_cell(a, k), s90_cell(b, k),
                              fmt(pb - pa, 2, sign=True) if math.isfinite(pa) and math.isfinite(pb) else "—"])
    s90 = table(["体系", "金属", "步", "pKa 等效（t ≥ 10/50 ps）", "pKa 等效（t ≥ 90 ps）", "差"], srows, num_cols=(3, 4, 5), cls="compact")
    o90 = occ_get(d, "start90", "solution", "DPA4C-Plus", "Nd", "box_pbc", 300)
    seen = [str(n) for n in range(4) if o90 and num(o90[f"P{n}"]) > 0]
    s90_note = f"（例如溶液 Nd 在 t ≥ 90 ps 只出现 n = {'、'.join(seen)}）" if o90 and len(seen) < 4 else ""
    vsig = sorted((r for r in d.vth if inum(r["n_T"]) >= 2
                   and (num(r["dH0_CI95_hi"]) < 0 or num(r["dH0_CI95_lo"]) > 0)),
                  key=lambda r: (r["model"], METALS.index(r["metal"]), r["hydrolysis_step"]))
    vsig_txt = "、".join(f"{MODEL_LABEL[r['model']]} {r['metal']} 第 {r['hydrolysis_step']} 步" for r in vsig) or "无"
    body = f"""
<p>把 box 里 3 个 OH<sup>−</sup> 在 Ln 上的占据数分布 P<sub>n</sub> 当作平衡分布，就能得到每一步的“pKa 等效值”，用来和论文的 pKa 对比（第 6 节）。</p>
<h3>5.1 300 K 占据数 P<sub>n</sub>（n = Ln 上的 OH 数）</h3>
{occ_tab}
{data_line("kinetics/kinetics_occupancy.csv", "kinetics/kinetics_occupancy_vs_time.csv")}
<h3>5.2 从占据数到 pKa</h3>
<div class="formula">Ln(OH)<sub>n</sub> ⇌ Ln(OH)<sub>n−1</sub> + OH<sup>−</sup>：K<sub>d</sub> = (P<sub>n−1</sub>/P<sub>n</sub>) · (4 − n) · c<sub>1</sub></div>
<div class="formula">pKa<sub>n</sub>(等效) = pK<sub>w</sub> + log<sub>10</sub> K<sub>d</sub> = 14 + log<sub>10</sub> K<sub>d</sub>；ΔG<sub>acid</sub> = RT ln10 · pKa（300 K：{fmt(8.314462618e-3 * 300 * math.log(10), 3)} kJ/mol 每单位）</div>
<p class="muted">c<sub>1</sub> = 一个分子在 box 中的浓度（溶液 {fmt(c1['solution'], 4)} M，石墨 {fmt(c1['graphite'], 4)} M）；(4 − n) = Ln 有 n−1 个 OH 时游离 OH<sup>−</sup> 的个数。
两个状态的 committed 停留都 ≥ 1 ps 时直接给出数值（加粗）；两个状态都出现过、但其中一个 committed 停留 &lt; 1 ps 记为“瞬态”（只有量级意义）；只有一个状态出现过时给界限（假设没出现的状态占据数不超过出现的那个，即 P<sub>n−1</sub> = P<sub>n</sub> 处的 pKa）；两个都没出现记为“无法确定”。</p>
{data_line("kinetics/kinetics_thermo.csv", "comparison/compare_pKa_stepwise.csv")}
<h3>5.3 van 't Hoff：失 OH 方向的 ΔH°、ΔS°（溶液，≥ 2 个温度）</h3>
{vh_tab}
<p class="muted">ΔH°<sub>d</sub> &lt; 0 表示升温使 OH 更倾向于留在 Ln 上（与 4.3 中负的失 OH E<sub>a</sub> 一致）。CI 都很宽，95% CI 不含 0 的只有 {vsig_txt}。</p>
{data_line("kinetics/kinetics_vanthoff.csv")}
<h3>5.4 起始构型记忆：只用 t ≥ 90 ps 的结果（DPA4C-Plus）</h3>
{s90}
<p class="muted">石墨 4 种起始构型的 P<sub>n</sub> 到 ~90 ps 才基本一致，所以用 t ≥ 90 ps 重算一遍。“瞬态”= 其中一个状态的 committed 停留 &lt; 1 ps，数值只有量级意义；“n = … 未出现”= 该设定下这个状态从未出现，这一步的 pKa 无法确定{s90_note}。</p>
{figure("kinetics/fig_kinetics_occupancy.svg", "占据数 P<sub>n</sub>：各势函数、各起始构型、随时间", ("kinetics/kinetics_occupancy.csv", "kinetics/kinetics_occupancy_vs_time.csv"))}
"""
    return card("occupancy", "5 占据数与 pKa 等效值", "平衡占据数 → K<sub>d</sub> → pKa（pK<sub>w</sub> = 14）", body)


# ------------------------------------------------------------------ 6 paper comparison
NOTE_ZH = (
    (r"median of (\d+) Table 2 step-(\d) entries, read as cumulative", r"Table 2 第 \2 步 \1 个值的中位数，按累积常数读"),
    (r"8\.1/16\.2/24\.3 = 8\.1 x n, i\.e\. cumulative constants", "8.1/16.2/24.3 = 8.1 × n，即累积常数"),
    (r"sum of stepwise pKa1\.\.pKa_n", "逐步 pKa<sub>1</sub> … pKa<sub>n</sub> 求和"),
    (r"sum of stepwise pKa_eq \(measured\); CI combined in quadrature \(approximate\)",
     "逐步 pKa 等效值求和；CI 按平方和合成（近似）"),
    (r"sum of stepwise pKa_eq \(transient\); CI combined in quadrature \(approximate\)",
     "逐步 pKa 等效值求和（含瞬态值）；CI 按平方和合成（近似）"),
    (r"monotonic \((\d+) points, smallest increment ([\d.]+) kJ/mol\): no barrier, so no activation energy is defined; "
     r"pKa from endpoint ([\d.]+)", r"单调上升（\1 个点，最小增量 \2 kJ/mol）：没有能垒，无法定义活化能；由终点得 pKa = \3"),
    (r"Table 1", "Table 1"),
    (r"both states visited with committed dwell >= 1 ps", "两个状态的 committed 停留都 ≥ 1 ps"),
    (r"one state only visited transiently \(committed dwell < 1 ps\); low confidence",
     "其中一个状态 committed 停留 &lt; 1 ps（瞬态），可信度低"),
    (r"apparent Eyring barrier, N = (\d+) events; no paper counterpart", r"表观 Eyring 能垒，N = \1；论文无对应量"),
)
UNTRANSLATED: set[str] = set()


def note_zh(note: str) -> str:
    for pat, rep in NOTE_ZH:
        if re.fullmatch(pat, note):
            return re.sub(pat, rep, note)
    UNTRANSLATED.add(note)
    return esc(note)


def temp_label(value) -> str:
    t = num(value)
    return f"{t:g}" if math.isfinite(t) else "—"


def src_label(source: str, system: str = "", model: str = "") -> str:
    s = source
    if s.startswith("Shiery2021 Table 1"):
        return "Shiery 2021 Table 1（DFT）" + ("，逐步求和" if "summed" in s else "")
    if s.startswith("Shiery2021 Table 2"):
        return "Shiery 2021 Table 2 实验（中位数）"
    if s.startswith("Shiery2021 SI PMF"):
        return "Shiery 2021 SI 第一步 PMF"
    if s.startswith("Shiery2021 SI"):
        return "Shiery 2021 SI 优化团簇"
    if s.startswith("Table 2 ref 99") or s.startswith("Kragten 1984 (Table 2"):
        return "Kragten 1984（Table 2 ref 99），按原值"
    if s.startswith("Kragten 1984 read as cumulative"):
        return "Kragten 1984，按累积常数读（逐步 8.1/8.1/8.1）"
    if s.startswith("this work"):
        m = re.match(r"this work: (.+?) \((\w+)\)", s)
        if m:
            mdl = m.group(1)
            sysname = SYSTEM_LABEL.get(m.group(2), m.group(2))
            return f"本工作：{sysname} {MODEL_LABEL.get(mdl, mdl)}"
        if system:
            return f"本工作：{SYSTEM_LABEL.get(system, system)} {MODEL_LABEL.get(model, model)}"
        return "本工作 (MD)"
    return esc(s)


def pka_cell(r: dict | None) -> str:
    if r is None:
        return '<span class="muted">无数据</span>'
    st = r["status"]
    if st == "undetermined":
        return pill("无法确定", "bad")
    if st in ("upper_bound", "lower_bound"):
        return f"{'≤' if st == 'upper_bound' else '≥'} {fmt(r['pKa_bound'])} {pill('界限', 'warn')}"
    small = f"<small>{ci(r['pKa_CI95_lo'], r['pKa_CI95_hi'])}</small>" if finite(r["pKa_CI95_lo"]) else ""
    return f"<b>{fmt(r['pKa'])}</b>{' ' + pill('瞬态', 'warn') if st == 'transient' else ''}{small}"


def ddg_cell(r: dict | None) -> str:
    if r is None or r["status"] == "undetermined" or not finite(r["ddG_vs_calc_kJmol_298K"]):
        return "—"
    sym = {"upper_bound": "≤ ", "lower_bound": "≥ "}.get(r["status"], "")
    return sym + fmt(r["ddG_vs_calc_kJmol_298K"], 1, sign=True)


def stepw_get(d, metal: str, k: str, kind: str, system: str = "", model: str = "") -> dict | None:
    if kind in ("calc", "experiment"):
        return first(d.stepw, metal=metal, hydrolysis_step=k, status=kind)
    return first(d.stepw, metal=metal, hydrolysis_step=k, system=system, model=model,
                 source=lambda s: s.startswith("this work"))


def sec_paper_a(d) -> str:
    # 6.1 Table 1
    t1_rows = [{"cells": [f"<b>{r['Ln']}</b>" if r["Ln"] in METALS else r["Ln"], r["assumed_CN_aqua_ion"],
                          f"{fmt(r['pKa1_calc'], 1)}<small>± {fmt(r['pKa1_error'], 1)}</small>", fmt(r["pKa1_MAD"], 1),
                          fmt(r["pKa2_calc"], 1), fmt(r["pKa2_MAD"], 1), fmt(r["pKa3_calc"], 1), fmt(r["pKa3_MAD"], 1),
                          fmt(r["dG1_kJmol_298K"]), fmt(r["dG2_kJmol_298K"]), fmt(r["dG3_kJmol_298K"])],
                "cls": "hl" if r["Ln"] in METALS else ""} for r in d.t1]
    t1_tab = table([], t1_rows, num_cols=tuple(range(1, 11)), cls="compact",
                   thead=thead2(["Ln", "水合 CN"], [("pKa<sub>1</sub>", ["计算", "MAD"]), ("pKa<sub>2</sub>", ["计算", "MAD"]),
                                                    ("pKa<sub>3</sub>", ["计算", "MAD"]), ("ΔG (kJ/mol, 298 K)", ["第 1 步", "第 2 步", "第 3 步"])]))

    # 6.2 stepwise pivot
    piv = []
    for metal in METALS:
        rows = []
        for kind, system, model in [("calc", "", ""), ("experiment", "", "")] + [("md", s, m) for s, m in SERIES]:
            recs = [stepw_get(d, metal, k, kind, system, model) for k in ("1", "2", "3")]
            if kind == "calc":
                label = "Shiery 2021 Table 1（DFT，298 K）"
            elif kind == "experiment":
                label = "Table 2 实验（中位数 [最小, 最大]）"
            else:
                label = f"本工作：{SYSTEM_LABEL[system]} {MODEL_LABEL[model]}"
            if kind == "md" and all(r is None for r in recs):
                rows.append({"cells": [label, collapsed_cells(9)], "cls": "weak"})
                continue
            cells = []
            for r in recs:
                if r is None:
                    cells += [('<span class="muted">无实验值</span>', 'class="num"'), "—", "—"]
                    continue
                dg = fmt(r["dG_acid_kJmol"], 1) if finite(r["dG_acid_kJmol"]) else "—"
                cells += [(pka_cell(r), 'class="num"'), dg, ddg_cell(r)]
            rows.append({"cells": [label] + cells,
                         "cls": "paper" if kind != "md" else ("hl" if model == "DPA4C-Plus" else "")})
        piv.append(f"<h4>{metal_html(metal)}</h4>" + table(
            [], rows, num_cols=(2, 3, 5, 6, 8, 9), cls="compact",
            thead=thead2(["来源"], [(f"第 {k} 步", ["pKa", "ΔG<sub>acid</sub>", "ΔΔG vs DFT"]) for k in (1, 2, 3)])))
    mad_rows = []
    for metal in METALS:
        for kind, system, model in [("calc", "", "")] + [("md", s, m) for s, m in SERIES]:
            recs = [stepw_get(d, metal, k, kind, system, model) for k in ("1", "2", "3")]
            if not any(r and finite(r["MAD_vs_table2_as_listed"]) for r in recs):
                continue
            label = "Shiery 2021 Table 1（DFT）" if kind == "calc" else f"本工作：{SYSTEM_LABEL[system]} {MODEL_LABEL[model]}"
            cells = []
            for r in recs:
                st = r["status"] if r else ""
                mad = fmt(r["MAD_vs_table2_as_listed"], 2) if r else "—"
                if st == "transient" and mad != "—":
                    mad += " " + pill("瞬态", "warn")
                dmed = "—"
                if r and finite(r["dpKa_vs_exp_median"]):
                    dmed = ("≤ " if st == "upper_bound" else "") + fmt(r["dpKa_vs_exp_median"], 2, sign=True)
                cells += [mad, dmed]
            mad_rows.append({"cells": [metal_html(metal), label] + cells, "cls": "paper" if kind == "calc" else ""})
    mad_tab = table([], mad_rows, num_cols=tuple(range(2, 8)), cls="compact",
                    thead=thead2(["金属", "来源"], [(f"第 {k} 步", ["MAD", "Δ 中位数"]) for k in (1, 2, 3)]))
    md_dpa = [r for r in d.stepw if r["model"] == "DPA4C-Plus" and r["status"] in ("measured", "transient")]
    spread = {m: max(num(r["pKa"]) for r in md_dpa if r["metal"] == m) - min(num(r["pKa"]) for r in md_dpa if r["metal"] == m)
              for m in METALS}
    paper_spread = {m: num(first(d.t1, Ln=m)["pKa3_calc"]) - num(first(d.t1, Ln=m)["pKa1_calc"]) for m in METALS}

    # 6.3 Table 2
    def t2_row(r):
        cit = r["citation"]
        calc = r["source"] == "this_work_calc"
        short = "Shiery 2021 自己的计算值（= Table 1）" if calc else (cit if len(cit) <= 70 else cit[:68] + "…")
        return {"cells": [r["Ln"], r["step"], f"{fmt(r['pKa'], 2)}" + (f"<small>± {fmt(r['pKa_error'], 1)}</small>" if finite(r["pKa_error"]) else ""),
                          "论文计算" if calc else "实验", r["ref_number"] or "—", f'<span title="{esc(cit)}">{esc(short)}</span>',
                          fmt(r["dG_kJmol_298K"])], "cls": "paper" if calc else ""}

    t2_head = ["Ln", "步", "pKa", "类型", "文献号", "引文", "ΔG (kJ/mol)"]
    t2_sel = table(t2_head, [t2_row(r) for r in d.t2 if r["Ln"] in METALS], num_cols=(1, 2, 4, 6), cls="compact")
    t2_all = table(t2_head, [t2_row(r) for r in d.t2], num_cols=(1, 2, 4, 6), cls="compact")
    agree = Counter(r["agrees_within_0.05"] for r in d.madchk)
    bad = [r for r in d.madchk if r["agrees_within_0.05"] == "NO"]
    bad_txt = "；".join(f"{r['Ln']} 第 {r['step']} 步：重算 {fmt(r['MAD_recomputed'], 3)} vs Table 1 {fmt(r['MAD_table1'], 2)}" for r in bad)
    chk_rows = [{"cells": [r["Ln"], r["step"], fmt(r["pKa_calc"], 1), r["n_exp"], fmt(r["exp_mean"]), fmt(r["exp_median"]),
                           fmt(r["exp_min"]), fmt(r["exp_max"]), fmt(r["MAD_recomputed"], 3), fmt(r["MAD_table1"], 2),
                           {"yes": pill("一致", "ok"), "NO": pill("不一致", "bad")}.get(r["agrees_within_0.05"], pill("无实验值"))],
                 "cls": "hl" if r["agrees_within_0.05"] == "NO" else ""} for r in d.madchk]
    chk_tab = table(["Ln", "步", "pKa 计算", "实验个数", "均值", "中位数", "最小", "最大", "MAD 重算", "MAD Table 1", "一致？"],
                    chk_rows, num_cols=tuple(range(1, 10)), cls="compact")

    # 6.4 cumulative
    cum_rows = []
    for r in d.cum:
        cum_rows.append({"cells": [metal_html(r["metal"]), r["n_OH"], src_label(r["source"]),
                                   f"<b>{fmt(r['value'])}</b>" + (f"<small>{ci(r['lo'], r['hi'])}</small>" if finite(r["lo"]) else ""),
                                   fmt(r["diff_vs_calc_cumulative"], 2, sign=True), f'<span class="muted">{note_zh(r["note"])}</span>'],
                         "cls": "paper" if r["source"].startswith(("Shiery", "Table 2", "Kragten")) else ""})
    cum_tab = table(["金属", "n", "来源", "−log *β<sub>n</sub>", "与 DFT 累积值之差", "说明"], cum_rows, num_cols=(1, 3, 4), cls="compact")
    kr = [r for r in d.t2 if r["Ln"] == "Nd" and r["ref_number"] == "99"]
    kr_txt = "/".join(fmt(r["pKa"], 1) for r in sorted(kr, key=lambda x: x["step"]))
    exp_cum = {(r["metal"], r["n_OH"]): num(r["value"]) for r in d.cum if r["source"].startswith("Shiery2021 Table 2")}
    dft_cum = {(r["metal"], r["n_OH"]): num(r["value"]) for r in d.cum if r["source"].startswith("Shiery2021 Table 1")}
    dft_minus_exp = [dft_cum[k] - v for k, v in exp_cum.items() if k in dft_cum and k[1] != "1"]
    md_cum = [r for r in d.cum if r["source"].startswith("this work") and (r["metal"], r["n_OH"]) in exp_cum]
    md_minus_exp = {grp: [num(r["value"]) - exp_cum[(r["metal"], r["n_OH"])] for r in md_cum if (r["n_OH"] == "1") == (grp == 1)]
                    for grp in (1, 2)}

    body = f"""
<p>文献：Shiery, Cooper, Cantu, <i>Inorg. Chem.</i> 2021, 60, 10257（<a href="{ARTICLE_URL}" target="_blank" rel="noopener">文章</a> · <a href="{SI_URL}" target="_blank" rel="noopener">SI</a>）。
论文方法（据作者学位论文第 5 章）：pKa<sub>1</sub> 用约束 AIMD（Bluemoon，CP2K/PBE，64 个显式水）在 9 个固定 r(O–H) 上各跑 15 ps，积分约束力得到 PMF，取 PMF 终点；pKa<sub>2</sub>、pKa<sub>3</sub> 用 M06 团簇（第一层水和 OH 显式）+ CPCM 隐式溶剂计算三步反应自由能，再以 pKa<sub>1</sub> 为基准取相对值；La–Sm 按 9 配位、Eu–Lu 按 8 配位。所有值都对应 298.15 K，ΔG = RT ln10 · pKa = 5.708 kJ/mol × pKa。
<b>论文没有给出任何能垒或活化能</b>：第一步的 PMF 单调上升（6.7 节），第 2、3 步只有反应自由能。所以下面对比的是热力学量（pKa、ΔG）和结构，我们的活化能只能和自己的热力学比较。</p>
<h3>6.1 论文 Table 1（全部 15 个 Ln）</h3>
{t1_tab}
{data_line("literature/lit_table1_calc_pKa.csv")}
<h3>6.2 逐步 pKa、ΔG 与论文的差 ΔΔG</h3>
<p class="muted">MD 的 pKa 是第 5 节的等效值（300 K）；ΔG<sub>acid</sub> = RT ln10 · pKa（MD 用 300 K，论文 298.15 K）；ΔΔG = 5.708 × (pKa − pKa<sub>DFT</sub>) kJ/mol。“≤”为界限（另一个状态从未出现）；标黄 = DPA4C-Plus。实验行的方括号是 Table 2 中该步所有实验值的最小、最大值。</p>
{''.join(piv)}
<p>DPA4C-Plus 的三步 pKa 等效值之间只差 Dy {fmt(spread['Dy'], 1)}、Nd {fmt(spread['Nd'], 1)} 个单位，而论文第 1 → 第 3 步差 Dy {fmt(paper_spread['Dy'], 1)}、Nd {fmt(paper_spread['Nd'], 1)}：MD 中第 1 步明显比论文难（ΔΔG<sub>1</sub> &gt; 0），第 3 步明显比论文容易（ΔΔG<sub>3</sub> &lt; 0）。</p>
<h4>与 Table 2 实验值的 MAD（按论文把第 2、3 步当逐步常数的读法）</h4>
{mad_tab}
<p class="muted">MAD = |pKa − 每个实验值| 的平均；Δ 中位数 = pKa − 实验中位数。“≤”= MD 只给出上界（另一个状态从未出现），此时只给 Δ 中位数的上界、不算 MAD。</p>
{data_line("comparison/compare_pKa_stepwise.csv")}
{figure("comparison/fig_compare_pKa_stepwise.svg", "逐步 pKa：论文 DFT、Table 2 实验、本工作 MD 等效值", ("comparison/compare_pKa_stepwise.csv",))}
<h3>6.3 论文 Table 2（实验值汇总）与 MAD 复核</h3>
<p>用 Table 2 重算论文 Table 1 的 MAD：{agree['yes']}/{agree['yes'] + agree['NO']} 项一致（±0.05），{agree['']} 项无实验值；不一致：{bad_txt or '无'}。说明转录的 Table 1/2 可靠。下表是 Nd、Dy 的全部条目（蓝底 = 论文自己的计算值）。</p>
{t2_sel}
{details(f"全部 15 个 Ln 的 Table 2（{len(d.t2)} 条）", t2_all)}
{details("MAD 复核表（全部 Ln）", chk_tab + data_line("literature/lit_table2_mad_check.csv"))}
{data_line("literature/lit_table2_pKa_all.csv", "literature/lit_references.csv", "literature/lit_table2_mad_check.csv")}
<h3>6.4 累积常数读法</h3>
<p>Nd 的 Kragten 1984（ref 99）在 Table 2 中三步为 {kr_txt} = 8.1 × n，这是累积常数 −log *β<sub>n</sub>（Ln<sup>3+</sup> + n H<sub>2</sub>O ⇌ Ln(OH)<sub>n</sub> + n H<sup>+</sup>）的特征；
第 2、3 步的其它实验值在数值上也更接近累积值。按累积常数读时，n ≥ 2 的论文 DFT 累积值（逐步求和）比实验中位数高 {rng(dft_minus_exp, 1)}；DPA4C-Plus MD 的累积值（逐步等效值求和）在 n ≥ 2 时比实验高 {rng(md_minus_exp[2], 1)}，n = 1 时高 {rng(md_minus_exp[1], 1)}。</p>
{cum_tab}
{data_line("comparison/compare_pKa_cumulative.csv")}
{figure("comparison/fig_compare_pKa_cumulative.svg", "累积 −log *β<sub>n</sub>：论文 DFT（求和）、Table 2 实验、本工作 MD", ("comparison/compare_pKa_cumulative.csv",))}
"""
    return body


def chem(label: str) -> str:
    """'[Nd(OH)2(H2O)7]+' -> HTML with sub/superscripts."""
    m = re.match(r"^\[(.*)\](\d*[+-]?|0)$", label)
    if not m:
        return esc(label)
    inner = esc(m.group(1)).replace("H2O", "H<sub>2</sub>O")
    inner = re.sub(r"\)(\d+)", r")<sub>\1</sub>", inner)
    return f"[{inner}]<sup>{m.group(2).replace('-', MINUS)}</sup>"


def sec_paper_b(d) -> str:
    # 6.5 box occupancy
    box_rows = []
    for r in d.boxocc:
        md = r["kind"].startswith("MD")
        box_rows.append({"cells": [metal_html(r["metal"]), SYSTEM_LABEL[r["system"]], src_label(r["source"]),
                                   "MD 占据数（300 K）" if md else "由 pKa 预测",
                                   fmt(r["pKa1"], 1), fmt(r["pKa2"], 1), fmt(r["pKa3"], 1)]
                                  + [fmt_p(r[f"P{n}"]) for n in range(4)] + [f"<b>{fmt(r['mean_n_OH'])}</b>"],
                         "cls": "" if md else "paper"})
    box_tab = table(["金属", "体系", "来源", "类型", "pKa<sub>1</sub>", "pKa<sub>2</sub>", "pKa<sub>3</sub>",
                     "P<sub>0</sub>", "P<sub>1</sub>", "P<sub>2</sub>", "P<sub>3</sub>", "⟨n⟩"], box_rows,
                    num_cols=tuple(range(4, 12)), cls="compact")

    def mean_n(metal, system, prefix):
        r = first(d.boxocc, metal=metal, system=system, source=lambda s: s.startswith(prefix))
        return fmt(r["mean_n_OH"]) if r else "—"

    pred = "；".join(f"{m} {mean_n(m, 'solution', 'Shiery2021 Table 1')}（溶液）/ {mean_n(m, 'graphite', 'Shiery2021 Table 1')}（石墨）" for m in METALS)
    mdn = "；".join(f"{m} {mean_n(m, 'solution', 'this work: DPA4C-Plus')} / {mean_n(m, 'graphite', 'this work: DPA4C-Plus')}" for m in METALS)
    others = [num(r["mean_n_OH"]) for r in d.boxocc if r["kind"].startswith("MD") and "DPA4C-Plus" not in r["source"]]

    # 6.6 structure
    def srow(r, show_metal=True):
        si = r["system"] == "cluster"
        lab = chem(r["label"]) if si else src_label("this work: " + r["model"] + f" ({r['system']})")
        rw = f"{fmt(r['Ln_O_water_mean_A'], 3)}" + (f"<small>± {fmt(r['Ln_O_water_sd_A'], 3)}</small>" if finite(r["Ln_O_water_sd_A"]) else "")
        ro = f"{fmt(r['Ln_O_OH_mean_A'], 3)}" + (f"<small>± {fmt(r['Ln_O_OH_sd_A'], 3)}</small>" if finite(r["Ln_O_OH_sd_A"]) else "")
        return {"cells": ([metal_html(r["metal"])] if show_metal else []) + [r["n_OH"], lab,
                          fmt(r["fraction"], 3) if finite(r["fraction"]) else "—", fmt(r["CN_water"]), fmt(r["CN_OH"]),
                          f"<b>{fmt(r['CN_total'])}</b>", rw, ro, fmt(r["Ln_O_all_mean_A"], 3)],
                "cls": "paper" if si else ""}

    main_rows, d_cn, d_rw, d_ro = [], [], [], []
    cn_drop = []
    for metal in METALS:
        for s in ("graphite", "solution"):
            a = first(d.struct, metal=metal, n_OH="0", system=s, model="DPA4C-Plus", T_K="300")
            b = first(d.struct, metal=metal, n_OH="3", system=s, model="DPA4C-Plus", T_K="300")
            if a and b:
                cn_drop.append(num(a["CN_total"]) - num(b["CN_total"]))
        for n in ("0", "1", "2", "3"):
            si = first(d.struct, metal=metal, n_OH=n, system="cluster")
            mds = [first(d.struct, metal=metal, n_OH=n, system=s, model="DPA4C-Plus", T_K="300") for s in ("graphite", "solution")]
            for r in [si] + mds:
                if r:
                    main_rows.append(srow(r))
            for r in mds:
                if r and si:
                    d_cn.append(num(r["CN_total"]) - num(si["CN_total"]))
                    d_rw.append(num(r["Ln_O_water_mean_A"]) - num(si["Ln_O_water_mean_A"]))
                    if finite(r["Ln_O_OH_mean_A"]):
                        d_ro.append(num(r["Ln_O_OH_mean_A"]) - num(si["Ln_O_OH_mean_A"]))
    s_head = ["金属", "n<sub>OH</sub>", "结构 / 来源", "占比", "CN(水)", "CN(OH)", "CN(总)", "⟨Ln–O(水)⟩ (Å)", "⟨Ln–O(OH)⟩ (Å)", "⟨Ln–O⟩ (Å)"]
    s_tab = table(s_head, main_rows, num_cols=tuple(range(3, 10)), cls="compact")
    md_struct = [r for r in d.struct if r["system"] != "cluster"]
    other_big = sorted((r for r in md_struct if num(r["CN_other"]) >= 0.005), key=lambda r: -num(r["CN_other"]))
    other_small = [num(r["CN_other"]) for r in md_struct if 0 < num(r["CN_other"]) < 0.005]
    other_txt = "、".join(f"{MODEL_LABEL[r['model']]} {r['metal']} {SYSTEM_LABEL[r['system']]} {r['T_K']} K n = {r['n_OH']}"
                          f"（{inum(r['frames']):,} 帧，{fmt(r['CN_other'])}）" for r in other_big) or "无"
    if other_small:
        other_txt += f"；其余 ≤ {fmt(max(other_small), 3)}"
    shell_gap = max(num(r["CN_total"]) - num(r["CN_water"]) - num(r["CN_OH"]) - num(r["CN_other"]) for r in md_struct)
    s_note = (f'<p class="muted">MD 行：t ≥ 10 ps 的帧按当时 Ln 上的 OH 数 n 分组，n 用第 4 节动力学的 r<sub>c</sub>(OH) 计数，CN(OH) = n；'
              f'CN(水)、CN(总) 用第 3 节的 REE–O 总截断计数。CN(总) 还包含 CN(其它)（O<sup>2−</sup>/H<sub>3</sub>O<sup>+</sup>；≥ 0.005 的只有 {other_txt}）'
              f'和落在 r<sub>c</sub>(OH) 与总截断之间的 OH（≤ {fmt(shell_gap, 3)}），所以 水 + OH 可能略小于 总。'
              '占比 = 该 n 状态的帧占比；⟨Ln–O⟩ 后的 ± 是分布标准差，不是标准误。团簇行：CN(总) = 3 Å 内的 O 数。</p>')
    all_rows = []
    for r in d.struct:
        all_rows.append([src_label(r["source"]) if r["system"] == "cluster" else "MD", SYSTEM_LABEL[r["system"]],
                         MODEL_LABEL.get(r["model"], r["model"]), metal_html(r["metal"]), r["T_K"] or "—", r["n_OH"],
                         r["n_trajectories"] or "—", f"{inum(r['frames']):,}" if finite(r["frames"]) else "—",
                         fmt(r["fraction"], 3) if finite(r["fraction"]) else "—", fmt(r["CN_water"]), fmt(r["CN_OH"]),
                         fmt(r["CN_other"]), fmt(r["CN_total"]), fmt(r["Ln_O_water_mean_A"], 3), fmt(r["Ln_O_OH_mean_A"], 3)])
    s_all = table(["来源", "体系", "势函数", "金属", "T (K)", "n<sub>OH</sub>", "轨迹", "帧", "占比", "CN(水)", "CN(OH)", "CN(其它)", "CN(总)",
                   "⟨Ln–O(水)⟩", "⟨Ln–O(OH)⟩"], all_rows, num_cols=tuple(range(4, 15)), cls="compact")
    si_rows = [{"cells": [r["Ln"], chem(r["label"]), r["n_atoms"], f"{r['n_oh_label']} / {r['n_oh_geom']}", r["n_water"], r["CN_O_3A"],
                          f"{fmt(r['Ln_O_water_mean'], 3)}<small>{ci(r['Ln_O_water_min'], r['Ln_O_water_max'], 3)}</small>",
                          f"{fmt(r['Ln_O_OH_mean'], 3)}<small>{ci(r['Ln_O_OH_min'], r['Ln_O_OH_max'], 3)}</small>" if finite(r["Ln_O_OH_mean"]) else "—",
                          fmt(r["Ln_O_all_mean"], 3), fmt(r["O_H_covalent_max"], 3)],
                "cls": "hl" if r["Ln"] in METALS else ""} for r in d.si]
    si_tab = table(["Ln", "结构", "原子数", "n<sub>OH</sub>（标签 / 几何）", "水", "CN (≤ 3 Å)", "⟨Ln–O(水)⟩ [最小, 最大]", "⟨Ln–O(OH)⟩ [最小, 最大]",
                    "⟨Ln–O⟩", "最长 O–H"], si_rows, num_cols=tuple(range(2, 10)), cls="compact")
    odd_labels = [r["label"] for r in d.si if re.search(r"\(OH\d\)", r["label"])]
    if odd_labels:
        si_tab += ('<p class="muted">结构名称照抄 SI 原文（' + "、".join(chem(x) for x in odd_labels) +
                   " 是 SI 中的写法，即 n<sub>OH</sub> = 3）。</p>")

    # 6.7 PMF
    def qlabel(q):
        if q.startswith("dG(end)"):
            return "PMF 终点 − 起点（脱质子自由能）"
        if q.startswith("dG1 (acid"):
            return "ΔG<sub>1</sub>（酸解离，第 1 步）"
        if q.startswith("dG1 acid-equivalent"):
            return "ΔG<sub>1</sub> 等效 = RT ln10 · pKa<sub>1</sub>(等效)"
        if q.startswith("dG_act Ln-OH loss"):
            return "ΔG<sup>‡</sup>：失 OH 1→0（第 1 步逆向）"
        if q.startswith("dG_act OH uptake"):
            return "ΔG<sup>‡</sup>：得 OH 0→1（第 1 步正向）"
        return esc(q)

    pmf_rows = [{"cells": [metal_html(r["metal"]), src_label(r["source"]), qlabel(r["quantity"]),
                           f"<b>{fmt(r['value_kJmol'])}</b>" + (f"<small>{ci(r['lo'], r['hi'])}</small>" if finite(r["lo"]) else ""),
                           temp_label(r["T_K"]), f'<span class="muted">{note_zh(r["note"])}</span>'],
                 "cls": "paper" if r["source"].startswith("Shiery") else ""} for r in d.pmf]
    pmf_tab = table(["金属", "来源", "量", "kJ/mol", "T (K)", "说明"], pmf_rows, num_cols=(3, 4), cls="compact")
    pe_diff = [abs(num(r["pKa_from_PMF_endpoint"]) - num(r["pKa1_table1"])) for r in d.pmf_end]
    act_rows = [(num(r["value_kJmol"]), int(m.group(1)) if (m := re.search(r"N = (\d+)", r["note"])) else 0,
                 f"{r['metal']} {SYSTEM_LABEL['solution' if '(solution)' in r['source'] else 'graphite']}")
                for r in d.pmf if r["quantity"].startswith("dG_act")]
    md_act = [v for v, n, _ in act_rows if n >= 5]
    act_few = sorted((v, n, lab) for v, n, lab in act_rows if n < 5)
    act_few_txt = ("（不计入事件数 N &lt; 5 的：" + "、".join(f"{lab} {fmt(v, 1)}（N = {n}）" for v, n, lab in act_few) + "）"
                   if act_few else "")
    paper_dg1 = [num(r["value_kJmol"]) for r in d.pmf if r["quantity"].startswith("dG1 (acid")]
    md_dg1 = [num(r["value_kJmol"]) for r in d.pmf if r["quantity"].startswith("dG1 acid-equivalent")]

    def s2_count(k):
        pairs = [(num(r[f"pKa{k}_MAD_manuscript"]), num(r[f"pKa{k}_MAD_ma_def2_TZVPP_opt"])) for r in d.s2]
        pairs = [p for p in pairs if finite(p[0]) and finite(p[1])]
        return sum(b > a for a, b in pairs), len(pairs)

    (w2, n2), (w3, n3) = s2_count(2), s2_count(3)
    pe_rows = [{"cells": [r["Ln"], r["n_points"], fmt(r["r_first_A"], 3), fmt(r["r_last_A"], 3), fmt(r["dG_first"]), fmt(r["dG_last"]),
                          f"<b>{fmt(r['dG_last_minus_first'])}</b>", fmt(r["pKa_from_PMF_endpoint"]), fmt(r["pKa1_table1"], 1),
                          fmt(num(r["pKa_from_PMF_endpoint"]) - num(r["pKa1_table1"]), 2, sign=True)],
                "cls": "hl" if r["Ln"] in METALS else ""} for r in d.pmf_end]
    pe_tab = table(["Ln", "点数", "r(O–H) 起点 (Å)", "终点 (Å)", "ΔG 起点", "ΔG 终点", "终点 − 起点", "pKa（PMF 终点）", "pKa<sub>1</sub>（Table 1）", "差"],
                   pe_rows, num_cols=tuple(range(1, 10)), cls="compact")

    # 6.8 S1 / S2
    s1_tab = table(["Ln", "4f 电子组态", "自旋多重度"],
                   [[r["Ln"], re.sub(r"([spdf])(\d+)", r"\1<sup>\2</sup>", esc(r["electronic_configuration"])), r["multiplicity"]]
                    for r in d.s1], num_cols=(2,), cls="compact")
    s2_tab = table([], [[r["Ln"], fmt(r["pKa2_MAD_manuscript"], 1), fmt(r["pKa2_MAD_ma_def2_TZVPP_opt"], 1),
                         fmt(r["pKa3_MAD_manuscript"], 1), fmt(r["pKa3_MAD_ma_def2_TZVPP_opt"], 1)] for r in d.s2],
                   num_cols=(1, 2, 3, 4), cls="compact",
                   thead=thead2(["Ln"], [("pKa<sub>2</sub> MAD", ["正文", "ma-def2-TZVPP 优化"]), ("pKa<sub>3</sub> MAD", ["正文", "ma-def2-TZVPP 优化"])]))
    return f"""
<h3>6.5 把论文 pKa 放进我们的 box：预测的占据数 vs MD</h3>
<p>用论文 pKa 和 5.2 节同一公式（box 内 3 个 OH<sup>−</sup>，c<sub>1</sub> 同 MD）反推平衡占据数。论文 DFT 预测 ⟨n⟩：{pred}；DPA4C-Plus MD（溶液 / 石墨）：{mdn}；其它势函数 ⟨n⟩ = {rng(others)}（OH 基本不离开 Ln）。
Kragten 两种读法给出完全相反的预测（⟨n⟩ ≈ 1 或 3），说明 Table 2 的读法本身就会决定结论。</p>
{box_tab}
{data_line("comparison/compare_box_occupancy.csv")}
{figure("comparison/fig_compare_box_occupancy.svg", "box 中的占据数：论文 pKa 预测 vs MD", ("comparison/compare_box_occupancy.csv",))}
<h3>6.6 结构：SI 优化团簇 vs MD（按 Ln 上的 OH 数 n 分开）</h3>
<p>论文团簇固定 CN（Nd 9、Dy 8），MD 中 CN 由动力学决定。DPA4C-Plus（300 K）与团簇相比：CN(总) 差 {rng(d_cn)}，⟨Ln–O(水)⟩ 差 {rng(d_rw, 3)} Å，⟨Ln–O(OH)⟩ 差 {rng(d_ro, 3)} Å。两者都显示随 n 增加 Ln–O(水) 变长、CN(水) 减少；不同的是 MD 中 CN(总) 从 n = 0 到 3 也减少了 {rng(cn_drop)}，而团簇的 CN(总) 按构造固定不变。</p>
{s_tab}
{s_note}
{data_line("comparison/compare_structure_by_nOH.csv", "literature/lit_SI_structures_Nd.xyz", "literature/lit_SI_structures_Dy.xyz")}
{details("全部势函数、全部温度的结构表", s_all)}
{details("SI 中 15 个 Ln × 4 个优化结构的几何指标", si_tab + data_line("literature/lit_SI_structures_metrics.csv"))}
{figure("comparison/fig_compare_structure.svg", "CN 与 Ln–O 距离随 n<sub>OH</sub> 的变化：SI 团簇 vs MD", ("comparison/compare_structure_by_nOH.csv",))}
<h3>6.7 第一步 PMF：论文没有能垒</h3>
<p>SI 图 S1–S15 给出第一步水解沿 r(O–H) 0.99 → 1.40 Å 的 PMF：单调上升，没有过渡态。15 个 Ln 的 PMF 终点换算成 pKa 后与 Table 1 pKa<sub>1</sub> 最多相差 {fmt(max(pe_diff), 2)}，即 Table 1 的 ΔG<sub>1</sub> 就是 PMF 终点。
因此论文里不存在可与我们 ΔG<sup>‡</sup> / E<sub>a</sub> 对比的量。下表把第一步的论文热力学与我们的热力学和能垒并列，方便比较量级：MD 第一步的 ΔG<sup>‡</sup> 为 {rng(md_act, 1)} kJ/mol{act_few_txt}，论文 ΔG<sub>1</sub> 为 {rng(paper_dg1, 1)} kJ/mol，MD 的 ΔG<sub>1</sub> 等效值为 {rng(md_dg1, 1)} kJ/mol。</p>
{pmf_tab}
{data_line("comparison/compare_PMF_step1.csv")}
{details("PMF 终点复核（全部 15 个 Ln）", pe_tab + data_line("literature/lit_SI_PMF_endpoint_check.csv", "literature/lit_SI_PMF_first_hydrolysis_digitized.csv"))}
{figure("comparison/fig_compare_PMF_step1.svg", "第一步：论文 PMF（数字化）与 MD 的 ΔG<sup>‡</sup>、ΔG<sub>1</sub> 等效", ("comparison/compare_PMF_step1.csv", "literature/lit_SI_PMF_first_hydrolysis_digitized.csv"))}
<h3>6.8 SI Table S1、S2</h3>
<div class="grid2"><div>{s1_tab}{data_line("literature/lit_tableS1_multiplicity.csv")}</div>
<div>{s2_tab}<p class="muted">用 ma-def2-TZVPP 重新优化结构后，pKa<sub>2</sub> 的 MAD 在 {w2}/{n2} 个元素上变大，pKa<sub>3</sub> 在 {w3}/{n3} 个元素上变大，并没有系统性改善。</p>{data_line("literature/lit_tableS2_basis_set_MAD.csv")}</div></div>
"""


def sec_paper(d) -> str:
    return card("paper", "6 与 Shiery 2021（Inorg. Chem. 60, 10257）的对比",
                "论文所有表格（Table 1、2、S1、S2）+ SI 结构与 PMF；逐步 pKa、ΔG、ΔΔG", sec_paper_a(d) + sec_paper_b(d))


# ------------------------------------------------------------------ 7 limits
def sec_limits(d) -> str:
    idx = rate_index(d)
    sol_n = [inum(r["N_events_300K"]) for r in where(d.act, system="solution", model="DPA4C-Plus")]
    uma = {m: [r for r in d.screen if r["system"] == "graphite" and r["model"] == "uma-s-1.2.1" and r["metal"] == m] for m in METALS}
    uma_txt = "；".join(
        f"{m} 排除 {sum(r['status'] == 'excluded' for r in uma[m])}、截断 {sum(r['status'] == 'truncated' for r in uma[m])}，"
        f"可用时长中位数 {fmt(median([r['usable_end_ps'] for r in uma[m] if r['status'] != 'excluded']), 0)} ps" for m in METALS)
    w02 = num(idx[("W0.2", "graphite", "DPA4C-Plus", "Dy", "ALL", 300, "3", "dissociation")]["dG_act_kJmol"])
    w20 = num(idx[("W2.0", "graphite", "DPA4C-Plus", "Dy", "ALL", 300, "3", "dissociation")]["dG_act_kJmol"])

    def step_order_holds(sc: str, metal: str) -> bool:
        g = {(k, dr): num(idx[(sc, "graphite", "DPA4C-Plus", metal, "ALL", 300, str(k), dr)]["dG_act_kJmol"])
             for k in (1, 2, 3) for dr, _ in DIRECTIONS}
        return (g[(1, "association")] < g[(2, "association")] < g[(3, "association")]
                and g[(1, "dissociation")] > g[(2, "dissociation")] > g[(3, "dissociation")])

    n_order = sum(step_order_holds(sc, m) for sc in SCENARIOS for m in METALS)
    n_combo = len(SCENARIOS) * len(METALS)
    order_txt = (f"三步的相对顺序（正向 1 &lt; 2 &lt; 3，逆向 1 &gt; 2 &gt; 3）在全部 {n_combo} 个“设定 × 金属”组合中都成立（4.4 节），比绝对值可靠。"
                 if n_order == n_combo else
                 f"三步的相对顺序（正向 1 &lt; 2 &lt; 3，逆向 1 &gt; 2 &gt; 3）只在 {n_order}/{n_combo} 个“设定 × 金属”组合中成立（4.4 节）。")
    sol_len = max(num(r["t_end_ps"]) for r in d.screen if r["system"] == "solution")
    c1 = d.kmeta["c_single_M"]
    oh_other = [num(r["REE_O_hydroxide_mean"]) for r in d.cn
                if r["model"] != "DPA4C-Plus" and r["variant"] in ("ALL", "box_pbc") and r["collapsed"] != "1"]
    unclear = ([f"石墨 {MODEL_LABEL[r['model']]} {r['metal']}" for r in d.gcut if r["OH_X_shell_min_clear"] != "True"]
               + [f"溶液 {MODEL_LABEL[r['model']]} {r['metal']}" for r in d.scut if r["OH_X_shell_min_clear"] != "True"])
    items = [
        f"<b>表观量：</b>committed 速率依赖过滤窗口与截断，例如石墨 DPA4C-Plus Dy 3→2 的 ΔG<sup>‡</sup> 在 W = 0.2 / 2.0 ps 时为 {fmt(w02, 1)} / {fmt(w20, 1)} kJ/mol。{order_txt}",
        f"<b>溶液统计量少：</b>每个温度只有 1 条 {fmt(sol_len, 0)} ps 轨迹，DPA4C-Plus 在 300 K 每个方向只有 {min(sol_n)}–{max(sol_n)} 个事件，E<sub>a</sub> 的 95% CI 很宽（见 4.3 节）。",
        f"<b>pKa 等效值的假设：</b>pK<sub>w</sub> = 14、理想溶液、标准态换算用 c<sub>1</sub> = 1 个分子 / box = {fmt(c1['solution'], 2)} M（溶液）/ {fmt(c1['graphite'], 2)} M（石墨），游离 OH<sup>−</sup> 浓度 (4 − n)c<sub>1</sub> 达到约 1 M，没有活度与离子对校正；300 K 与 298.15 K 的差别忽略不计。",
        "<b>起始构型记忆：</b>石墨 DPA4C-Plus 的占据数到 ~90 ps 才收敛（5.4 节）。",
        f"<b>UMA 坍塌：</b>UMA Dy 溶液全部排除；石墨 UMA {uma_txt}。UMA 的结果只作参考。",
        f"<b>势函数之间定性不同：</b>DPA4-Pro、MACE、UMA 中 3 个 OH 基本始终在 Ln 上（300 K REE–O(OH) = {rng(oh_other)}），DPA4C-Plus 中 OH 会频繁离开；哪个更接近真实，需要 AIMD 或实验结构来独立验证。",
        "<b>论文的可比性：</b>论文的 pKa<sub>1</sub> 是小 box（64 个水、每点 15 ps）约束 AIMD 的 PMF 终点，pKa<sub>2</sub>、pKa<sub>3</sub> 是固定配位数（Nd 9、Dy 8）团簇 + CPCM 的相对自由能；只有 298 K 热力学（且第 2、3 步实验值的读法有歧义，6.4 节），没有能垒。我们的 E<sub>a</sub> 无法与论文直接对比，只能和论文的 ΔG 比量级。",
        (f"<b>OH 溶剂化壳层：</b>{'、'.join(unclear)} 的 O–X(总) 壳层极小不清晰（3.2 节），其壳层 O–X(总) 不宜下结论，建议改看 O–REE(壳层) 与 REE–O(OH)/3。"
         if unclear else "<b>OH 溶剂化壳层：</b>壳层极小都清晰，但 O–X(总) 依赖截断定义（2.3 节），建议同时看 O–REE(壳层) 与 REE–O(OH)/3。"),
    ]
    return card("limits", "7 局限性", "解读数字前请先看这里", '<ul class="tight">' + "".join(f"<li>{x}</li>" for x in items) + "</ul>")


# ------------------------------------------------------------------ 8 downloads
FILE_INFO = {
    "README.md": "说明文档：定义、流程与输出文件",
    "analyze.py": "任务 (a)(b)：g(r)、截断、CN 与 g(r) 图（石墨 + 溶液）",
    "kinetics.py": "轨迹提取、坍塌筛查、committed 跃迁、速率、ΔG‡、Arrhenius、占据数、van 't Hoff",
    "literature.py": "Shiery 2021 表格转录、SI 结构与 PMF 提取（自动下载 SI PDF）",
    "compare.py": "与 Shiery 2021 的逐步/累积 pKa、占据数、结构、PMF 对比",
    "build_page.py": "生成本网页与 zip（只读取上面脚本的输出）",
    "cn_summary_300K.csv": "300 K CN 汇总（任务 a、b）：石墨 4 种起始构型 + ALL，溶液 box_pbc；含截断与标准误",
    "graphite_cn_summary_box_pbc.csv": "石墨 CN（势函数 × 金属 × 起始构型，300 K），含截断与两种标准误",
    "graphite_cn_by_replica_box_pbc.csv": "石墨 CN 逐条轨迹（筛查状态、可用时长、采样数）",
    "graphite_cutoffs_box_pbc.csv": "石墨截断：REE–O 极小，OH 中心严格/壳层截断与峰→谷降幅",
    "graphite_gr_box_pbc.csv": "石墨 g(r)：REE–O 与 OH 中心各原子对（Δr = 0.02 Å）",
    "graphite_rdf_accumulated.npz": "石墨 RDF 直方图累加器（analyze.py --reuse-rdf 可直接复用）",
    "graphite_rdf_accumulated.json": "上面 npz 的键与说明",
    "solution_cn_box_pbc.csv": "溶液 CN：全部势函数 × 300–500 K（逐帧计数，5 段 block 标准误）",
    "solution_gr_box_pbc_300K.csv": "溶液 300 K g(r)：REE–O 与 OH 中心各原子对",
    "solution_cutoffs_box_pbc_300K.csv": "溶液 300 K 截断与清晰度（峰→谷降幅）",
    "collapse_screen.csv": "每条轨迹的坍塌筛查：势能/几何起点、可用时长、状态",
    "graphite_REE_O_total_box_pbc_300K.svg": "图：石墨 REE–O g(r) 与截断",
    "graphite_OH_X_total_box_pbc_300K.svg": "图：石墨 OH 中心 O–X g(r) 与截断",
    "graphite_CN_box_pbc_300K.svg": "图：石墨 CN",
    "solution_REE_O_total_box_pbc_300K.svg": "图：溶液 REE–O g(r) 与截断",
    "solution_OH_X_total_box_pbc_300K.svg": "图：溶液 OH 中心 O–X g(r) 与截断",
    "solution_CN_box_pbc_300K.svg": "图：溶液 CN",
    "kinetics/kinetics_rates.csv": "全部设定 × 体系 × 势函数 × 金属 × 温度 × 步 × 方向：事件数、停留时间、k、ΔG‡ 与 95% CI",
    "kinetics/kinetics_arrhenius.csv": "Arrhenius / Eyring 拟合：E<sub>a</sub>、ln A、ΔH‡、ΔS‡（含 95% CI）与拟合质量",
    "kinetics/kinetics_occupancy.csv": "平衡占据数 P<sub>0</sub>–P<sub>3</sub>（95% CI）与 ⟨n⟩（全部设定）",
    "kinetics/kinetics_thermo.csv": "逐步 ΔG<sub>occ</sub>、ΔG°、pK、pKa 等效值（全部设定）",
    "kinetics/kinetics_vanthoff.csv": "van 't Hoff：ΔH°、ΔS°（失 OH 方向）",
    "kinetics/kinetics_cutoffs.csv": "动力学用 Ln–O(OH) 截断 r<sub>c</sub>（径向密度极小）",
    "kinetics/kinetics_events_by_traj.csv": "逐条轨迹：各类跃迁次数 N<sub>i→j</sub>、各状态停留时间、原始 P<sub>n</sub>、OH 多于 3 个 / 壳层内有 O<sup>2−</sup> 的帧比例",
    "kinetics/kinetics_occupancy_vs_time.csv": "P<sub>n</sub> 随时间变化（10 ps 分块）",
    "kinetics/kinetics_oh_density.csv": "Ln–O(OH) 径向数密度（用于定 r<sub>c</sub>）",
    "kinetics/kinetics_traces_example.csv": "示例轨迹 n(t)：原始值与 committed 值",
    "kinetics/kinetics_metadata.json": "分析设定与常数（pK<sub>w</sub>、c<sub>1</sub>、采样间隔、bootstrap 次数等）",
    "kinetics/fig_kinetics_dGact_300K.svg": "图：300 K 三步 ΔG‡",
    "kinetics/fig_kinetics_arrhenius_DPA4C_solution.svg": "图：溶液 DPA4C-Plus Arrhenius",
    "kinetics/fig_kinetics_sensitivity.svg": "图：分析设定敏感性",
    "kinetics/fig_kinetics_occupancy.svg": "图：占据数 P<sub>n</sub>",
    "kinetics/fig_kinetics_traces.svg": "图：示例轨迹 n(t)",
    "kinetics/fig_kinetics_OH_density_cutoff.svg": "图：Ln–O(OH) 径向密度与 r<sub>c</sub>",
    "comparison/compare_activation.csv": "三步 ΔG‡（300 K）与 E<sub>a</sub> 总表（体系 × 势函数 × 金属 × 方向）",
    "comparison/compare_pKa_stepwise.csv": "逐步 pKa、ΔG<sub>acid</sub>、ΔΔG：论文计算、Table 2 实验、MD 等效值",
    "comparison/compare_pKa_cumulative.csv": "累积常数读法：论文求和、实验、MD",
    "comparison/compare_box_occupancy.csv": "论文 pKa 在本 box 中预测的占据数 vs MD 占据数",
    "comparison/compare_structure_by_nOH.csv": "按 n<sub>OH</sub> 分开的 CN 与 Ln–O 距离：SI 团簇 vs MD",
    "comparison/compare_PMF_step1.csv": "第一步：论文 PMF 终点 / ΔG<sub>1</sub> vs MD ΔG<sub>1</sub> 等效与 ΔG‡",
    "comparison/fig_compare_pKa_stepwise.svg": "图：逐步 pKa 对比",
    "comparison/fig_compare_pKa_cumulative.svg": "图：累积常数对比",
    "comparison/fig_compare_box_occupancy.svg": "图：占据数对比",
    "comparison/fig_compare_structure.svg": "图：结构对比",
    "comparison/fig_compare_PMF_step1.svg": "图：第一步 PMF 与 MD",
    "literature/lit_table1_calc_pKa.csv": "论文 Table 1：15 个 Ln 的计算 pKa<sub>1–3</sub>、MAD、ΔG（298.15 K）",
    "literature/lit_table2_pKa_all.csv": "论文 Table 2：全部实验 pKa 与引用",
    "literature/lit_table2_mad_check.csv": "用 Table 2 重算 Table 1 的 MAD（核对转录）",
    "literature/lit_tableS1_multiplicity.csv": "SI Table S1：电子组态与自旋多重度",
    "literature/lit_tableS2_basis_set_MAD.csv": "SI Table S2：ma-def2-TZVPP 优化结构的 MAD",
    "literature/lit_SI_structures_metrics.csv": "SI 60 个优化结构的几何指标",
    "literature/lit_SI_structures_Nd.xyz": "SI 中 Nd 的 4 个 M06 优化结构（n<sub>OH</sub> = 0–3）",
    "literature/lit_SI_structures_Dy.xyz": "SI 中 Dy 的 4 个 M06 优化结构（n<sub>OH</sub> = 0–3）",
    "literature/lit_SI_PMF_first_hydrolysis_digitized.csv": "SI 图 S1–S15 第一步 PMF 数字化（矢量路径读取）",
    "literature/lit_SI_PMF_endpoint_check.csv": "PMF 终点换算的 pKa 与 Table 1 pKa<sub>1</sub> 对照",
    "literature/lit_references.csv": "Table 2 引用文献编号与条目",
    "literature/lit_metadata.json": "转录来源、核对方式与单位换算",
    "literature/source/table2_references.json": "Table 2 引用文献列表（literature.py 的输入，内容同 lit_references.csv）",
}
DOWNLOAD_GROUPS = (
    ("CN 与 g(r)（任务 a、b）", lambda rel: "/" not in rel and not rel.endswith((".py", ".md"))),
    ("动力学 kinetics/", lambda rel: rel.startswith("kinetics/")),
    ("与 Shiery 2021 的对比 comparison/", lambda rel: rel.startswith("comparison/")),
    ("文献转录 literature/", lambda rel: rel.startswith("literature/")),
    ("脚本与说明", lambda rel: rel.endswith((".py", ".md"))),
)


def bundle_files() -> list[str]:
    out = []
    for root, dirs, files in os.walk(HERE):
        dirs[:] = sorted(x for x in dirs if x not in SKIP_DIR_NAMES and not x.startswith("."))
        for name in sorted(files):
            rel = (Path(root) / name).relative_to(HERE).as_posix()
            if (name.startswith(".") or name.endswith((".log", ".tmp")) or name in SKIP_FILES
                    or rel in SKIP_RELS or rel == OUT_HTML.name):
                continue
            out.append(rel)
    return out


def file_sort_key(rel: str):
    ext = rel.rsplit(".", 1)[-1].lower()
    return ({"csv": 0, "json": 1, "npz": 1, "xyz": 1, "svg": 2}.get(ext, 3), rel.lower())


def describe(rel: str) -> str:
    if rel in FILE_INFO:
        return FILE_INFO[rel]
    ext = rel.rsplit(".", 1)[-1].lower()
    return {"csv": "CSV 数据", "svg": "SVG 图", "json": "JSON 元数据", "py": "Python 脚本"}.get(ext, "文件")


def dir_size(path: Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) if path.is_dir() else 0


def cache_size_label() -> str:
    n = dir_size(HERE / "kinetics_cache")
    return size_label(n) if n else "约 370 MB"


def sec_downloads(files: list[str], zip_label: str) -> str:
    blocks, used = [], set()
    for title, pred in DOWNLOAD_GROUPS:
        sel = sorted((f for f in files if f not in used and pred(f)), key=file_sort_key)
        used.update(sel)
        if sel:
            blocks.append(f"<h3>{title}（{len(sel)} 个文件）</h3>"
                          + table(["文件", "内容"], [[dl(f, f), describe(f)] for f in sel], cls="compact"))
    rest = sorted(f for f in files if f not in used)
    if rest:
        blocks.append("<h3>其它</h3>" + table(["文件", "内容"], [[dl(f, f), describe(f)] for f in rest], cls="compact"))
    total = sum((HERE / f).stat().st_size for f in files)
    body = f"""
<p><a class="dl" href="{OUT_ZIP.name}" download>{OUT_ZIP.name}</a><span class="size">{zip_label}</span>
：下面全部 {len(files)} 个文件（未压缩 {size_label(total)}）+ 本网页；解压后直接用浏览器打开 <code>{ZIP_ROOT}/index.html</code> 即可离线浏览，所有链接都是相对路径。</p>
{''.join(blocks)}
<p class="muted">没有打包：<code>kinetics_cache/</code>（{cache_size_label()}，<code>python kinetics.py extract 8</code> 可从原始轨迹重新生成）；
出版社的 SI 原文件 <code>literature/source/ic1c00662_si_001.pdf</code>（请从 <a href="{SI_URL}" target="_blank" rel="noopener">figshare</a> 下载，<code>literature.py</code> 也会自动下载）；以及日志和 Python 缓存。</p>"""
    return card("downloads", "8 数据下载", "每个文件都可单独点击下载；CSV 为 UTF-8、逗号分隔，第一行是列名", body)


# ------------------------------------------------------------------ 9 sources
REPRO = """module load python/3.12-26.1.0
python kinetics.py extract 8    # 原始轨迹 -> kinetics_cache/（{cache}）
python kinetics.py screen       # 坍塌筛查 -> collapse_screen.csv
python analyze.py               # 任务 (a)(b)：g(r)、截断、CN、图（--reuse-rdf 复用石墨 RDF）
python kinetics.py analyze 8    # 速率、ΔG‡、Arrhenius、占据数、pKa 等效值
python kinetics.py figures
python literature.py            # 论文表格、SI 结构与 PMF（自动下载 SI PDF）
python compare.py 8             # 与论文对比
python build_page.py            # 本网页 + zip"""


SOURCE_ZH = {
    "Table 1/2": ("Table 1/2", "转录自博士论文第 5 章 Table 5-1/5-2（与文章 Table 1/2 相同）；用 Table 2 重算 Table 1 的 MAD 核对（lit_table2_mad_check.csv）"),
    "Table S1/S2": ("Table S1/S2", "SI 第 S2、S9 页"),
    "structures": ("SI 优化结构", "SI 第 S10–S39 页（pypdf 读取文字层）"),
    "PMF": ("第一步 PMF", "SI 图 S1–S15：从 PDF 矢量路径读取数据点中心；坐标轴 0.95–1.45 Å、0–60 kJ/mol"),
}


def sec_sources(d) -> str:
    rows = []
    for k, v in d.lmeta.get("table_sources", {}).items():
        if k not in SOURCE_ZH:
            UNTRANSLATED.add(f"{k}: {v}")
        key, text = SOURCE_ZH.get(k, (k, v))
        rows.append([esc(key), esc(text)])
    conv = d.lmeta.get("conversion", "")
    m = re.search(r"=\s*([\d.]+)\s*kJ/mol at ([\d.]+) K", conv)
    if m:
        conv = f"ΔG = RT ln10 · pKa；{m.group(2)} K 时 RT ln10 = {m.group(1)} kJ/mol"
    else:
        UNTRANSLATED.add(conv)
    rows.append(["单位换算", esc(conv)])
    body = f"""
<ul class="tight">
<li>文章：R. C. Shiery, K. A. Cooper, D. C. Cantu, “Computational Prediction of All Lanthanide Aqua Ion Acidity Constants”, <i>Inorg. Chem.</i> 2021, 60, 10257–10266，<a href="{ARTICLE_URL}" target="_blank" rel="noopener">doi:10.1021/acs.inorgchem.1c00662</a>。</li>
<li>SI：<a href="{SI_URL}" target="_blank" rel="noopener">ic1c00662_si_001.pdf</a>（Table S1、S2，图 S1–S15 第一步 PMF，M06 优化结构坐标）。</li>
<li>R. C. Shiery, <i>Ab initio molecular dynamics simulations of lanthanide coordination structures in water and in faujasite</i>，博士学位论文，University of Nevada, Reno, 2022，<a href="{THESIS_URL}" target="_blank" rel="noopener">PDF</a>。第 5 章即这篇文章（Table 5-1/5-2 = 文章 Table 1/2），方法细节取自该章。</li>
</ul>
<h3>论文数据如何得到</h3>
{table(["内容", "来源与核对"], rows, cls="compact")}
<h3>复现</h3>
<p>脚本按相对路径读取 <code>genesis/</code> 下的原始数据：石墨 <code>Graphite-dy-nd-50-denser-dumping/</code>；溶液动力学 <code>md/ln-150ps-denser-dumping/</code>、<code>md/ln-150ps-every-step-dumping/</code>；溶液 CN/g(r) <code>md/dense-gr-cn/</code>。在本目录依次运行：</p>
<pre>{esc(REPRO.replace("{cache}", cache_size_label()))}</pre>"""
    return card("sources", "9 来源与复现", "文献链接、转录与核对方式、运行顺序", body)


# ------------------------------------------------------------------ assemble
def render(d, files: list[str], zip_label: str, stamp: str) -> str:
    sections = [toc(), metrics(d), notice(d), sec_summary(d), sec_methods(d), sec_cn(d), sec_kinetics(d),
                sec_occupancy(d), sec_paper(d), sec_limits(d), sec_downloads(files, zip_label), sec_sources(d)]
    body = "\n".join(sections)
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Dy、Nd 水解：配位数、三步表观活化能与 Shiery 2021 对比</title>
<style>{CSS}</style>
</head>
<body>
{page_header(zip_label)}
<main>
{body}
<div class="foot">生成于 {stamp}（build_page.py）。网页中的数字全部来自本目录的 CSV/JSON；页面只做取范围、作差之类的简单汇总。</div>
</main>
</body>
</html>
"""


def write_zip(files: list[str]) -> None:
    tmp = OUT_ZIP.with_suffix(".zip.tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        zf.write(OUT_HTML, f"{ZIP_ROOT}/{OUT_HTML.name}")
        for rel in files:
            zf.write(HERE / rel, f"{ZIP_ROOT}/{rel}")
    tmp.replace(OUT_ZIP)


# ------------------------------------------------------------------ validation
VOID = {"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "source", "track", "wbr"}
NOT_IN_P = {"div", "ul", "ol", "table", "h1", "h2", "h3", "h4", "figure", "details", "pre", "section", "p",
            "header", "main", "nav", "blockquote", "form"}
SUSPICIOUS = re.compile(r"(?<![A-Za-z])(nan|NaN|None|inf|null)(?![A-Za-z])|[{}]|−0\.0+(?![0-9])")


class PageChecker(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack, self.errors, self.links, self.text = [], [], [], []
        self.ids: Counter = Counter()
        self.skip = 0

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag in NOT_IN_P and "p" in self.stack:
            self.errors.append(f"<{tag}> inside <p> (line {self.getpos()[0]})")
        if a.get("id"):
            self.ids[a["id"]] += 1
        self.links += [a[k] for k in ("href", "src") if a.get(k)]
        if tag in ("style", "script"):
            self.skip += 1
        if tag not in VOID:
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if not self.stack or self.stack[-1] != tag:
            self.errors.append(f"</{tag}> but open <{self.stack[-1] if self.stack else '-'}> (line {self.getpos()[0]})")
            if tag not in self.stack:
                return
            while self.stack[-1] != tag:
                self.stack.pop()
        self.stack.pop()
        if tag in ("style", "script"):
            self.skip -= 1

    def handle_data(self, data):
        if not self.skip:
            self.text.append(data)


def validate(page: str, files: list[str]) -> list[str]:
    chk = PageChecker()
    chk.feed(page)
    chk.close()
    problems = list(chk.errors)
    problems += [f"untranslated note: {n}" for n in sorted(UNTRANSLATED)]
    if chk.stack:
        problems.append(f"unclosed tags: {chk.stack}")
    problems += [f"duplicate id {k}" for k, v in chk.ids.items() if v > 1]
    for link in chk.links:
        if link.startswith("#"):
            if link[1:] not in chk.ids:
                problems.append(f"broken anchor {link}")
        elif link.startswith(("http://", "https://")):
            if link not in EXTERNAL_URLS:
                problems.append(f"unexpected external link {link}")
        elif not (HERE / unquote(link.split("#")[0])).is_file():
            problems.append(f"missing file {link}")
    if set(files) - LINKED:
        problems.append(f"bundle files without a link: {sorted(set(files) - LINKED)}")
    if LINKED - set(files):
        problems.append(f"linked files missing from the bundle: {sorted(LINKED - set(files))}")
    text = " ".join(" ".join(chk.text).split())
    for m in SUSPICIOUS.finditer(text):
        problems.append(f"suspicious text {m.group(0)!r}: ...{text[max(0, m.start() - 60):m.end() + 30]}...")
    with zipfile.ZipFile(OUT_ZIP) as zf:
        bad, names = zf.testzip(), set(zf.namelist())
    if bad:
        problems.append(f"corrupt zip member {bad}")
    want = {f"{ZIP_ROOT}/{f}" for f in files} | {f"{ZIP_ROOT}/{OUT_HTML.name}"}
    if names != want:
        problems.append(f"zip content mismatch: {sorted(names ^ want)[:6]}")
    return problems


def main() -> None:
    d = load()
    files = bundle_files()
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    label, page = "…", ""
    for _ in range(5):
        LINKED.clear()
        page = render(d, files, label, stamp)
        OUT_HTML.write_text(page, encoding="utf-8")
        write_zip(files)
        new = size_label(OUT_ZIP.stat().st_size)
        if new == label:
            break
        label = new
    problems = validate(page, files)
    print(f"index.html {size_label(OUT_HTML.stat().st_size)} | {OUT_ZIP.name} {label} "
          f"({len(files) + 1} files) | linked files {len(LINKED)}")
    for p in problems:
        print("PROBLEM:", p)
    if problems:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
