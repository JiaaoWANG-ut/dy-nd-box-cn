#!/usr/bin/env python3
"""Machine-readable tables from Shiery, Cooper & Cantu, Inorg. Chem. 2021, 60, 10257.

"Computational Prediction of All Lanthanide Aqua Ion Acidity Constants",
DOI 10.1021/acs.inorgchem.1c00662.

* Tables 1 and 2 are transcribed from R. C. Shiery's open-access dissertation (University of
  Nevada, Reno, 2022), Chapter 5, Tables 5-1 and 5-2, which reproduce the article tables.
  The transcription is checked by recomputing every MAD of Table 1 from Table 2.
* Tables S1/S2, the 60 optimized structures (pp. S10-S39) and the 15 first-hydrolysis PMFs
  (Figures S1-S15, digitized from the vector drawing) come from the Supporting Information
  (figshare file 28696985).
"""

from __future__ import annotations

import csv
import json
import math
import re
import urllib.request
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "literature"
SOURCE = OUT / "source"
SI_PDF = SOURCE / "ic1c00662_si_001.pdf"
SI_URL = "https://ndownloader.figshare.com/files/28696985"
ARTICLE_URL = "https://pubs.acs.org/doi/10.1021/acs.inorgchem.1c00662"
DISSERTATION_URL = "https://scholarwolf.unr.edu/bitstreams/ae32a3d1-8c2e-40ab-a18f-9499bf7f21e1/download"
LN = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
RT_LN10_298 = 8.314462618e-3 * 298.15 * math.log(10.0)

# Table 1: Ln -> (pKa1, error, MAD1, pKa2, MAD2, pKa3, MAD3); None = "-" in the table.
TABLE1 = {
    "La": (8.4, 0.9, 0.7, 13.4, 4.5, 19.6, 7.1),
    "Ce": (8.8, 0.3, 0.5, 17.1, 0.8, 20.0, 5.3),
    "Pr": (9.0, 0.5, 0.6, 14.0, None, 20.8, 4.9),
    "Nd": (8.3, 0.5, 0.5, 12.9, 3.8, 17.8, 6.6),
    "Pm": (8.7, 0.5, 2.2, 12.8, None, 19.9, None),
    "Sm": (7.0, 0.1, 1.1, 10.9, 4.5, 16.7, 6.9),
    "Eu": (10.1, 0.2, 2.4, 19.8, None, 23.9, 0.4),
    "Gd": (8.4, 0.4, 0.6, 15.8, 0.3, 17.2, 5.5),
    "Tb": (9.1, 0.5, 1.6, 15.4, None, 20.4, 3.0),
    "Dy": (9.1, 0.1, 1.6, 16.1, None, 21.2, 1.7),
    "Ho": (8.4, 0.5, 0.9, 15.2, 1.1, 18.8, 3.9),
    "Er": (7.2, 0.9, 0.9, 11.8, 4.2, 15.9, 6.8),
    "Tm": (8.9, 0.1, 1.7, 12.1, 3.9, 15.8, 6.1),
    "Yb": (8.2, 0.3, 0.9, 15.1, 0.2, 18.3, 3.6),
    "Lu": (7.9, 0.9, 0.9, 14.0, 1.2, 17.1, 4.3),
}

# Table 2 in table order: "value(ref)"; TW = "This Work" (calculated in the article).
TABLE2 = {
    "La": (
        "7.4(73) 7.77(74) 8.14(75) 8.4±0.9(TW) 8.4(76) 8.5(2) 8.5(77) 8.52(78) 8.53(79) 8.6(81) 8.81(82) "
        "8.87(83) 8.89(84) 9.1(85) 9.2(86) 9.33(87) 9.94(88) 10(89) 10.0(90) 10.1(91) 10.3(92)",
        "13.4(TW) 17.9(81)",
        "19.6(TW) 26.12(80) 26.56(79) 26.84(78) 27.3(81)",
    ),
    "Ce": (
        "8.02(75) 8.1(93) 8.31(84) 8.34(82) 8.34(77) 8.43(79) 8.8±0.3(TW) 9.1(95) 9.3(96)",
        "16.3(93) 17.1(TW)",
        "20.0(TW) 23.65(94) 26.0(93) 26.25(79)",
    ),
    "Pr": (
        "7.1(73) 7.86(75) 8.1(77) 8.3(84) 8.32(82) 8.32(79) 8.5(96) 8.5(97) 8.54(83) 8.82(87) 8.94(98) "
        "9±0.5(TW) 9.0(90) 9.5(92)",
        "14.0(TW)",
        "20.8(TW) 25.66(79)",
    ),
    "Nd": (
        "7(73) 7.76(75) 8.0(2) 8.0(77) 8.1(99) 8.13(84) 8.18(82) 8.24(79) 8.3±0.5(TW) 8.5(96) 8.5(97) "
        "8.7(87) 9.0(90) 9.4(104) 9.46(88)",
        "9.2(100) 12.9(TW) 16.2(99) 17.35(103)",
        "17.8(TW) 23.54(80) 23.8(101) 23.9(102) 24.0(100) 24.3(99) 25.16(79) 26.2(103)",
    ),
    "Pm": ("6.5(73) 8.7±0.5(TW)", "12.8(TW)", "19.9(TW)"),
    "Sm": (
        "4.4(73) 7.0±0.1(TW) 7.11(78) 7.5(105) 7.6(75) 7.84(82) 7.84(84) 7.9(77) 8.02(79) 8.61(87) 8.9(90)",
        "10.9(TW) 15.0(105) 15.84(78)",
        "16.7(TW) 22.7(105) 23.44(78) 24.62(79)",
    ),
    "Eu": (
        "4.8(73) 6.7(106) 7.3(107) 7.49(75) 7.49(108) 7.66(84) 7.76(82) 7.8(77) 7.91(79) 8(85) 8.08(109) "
        "8.1(110) 8.2(111) 8.31(112) 8.34(83) 8.58(87) 8.7(113) 10.1±0.2(TW)",
        "19.8(TW)",
        "23.9(TW) 24.34(79)",
    ),
    "Gd": (
        "7.1(73) 7.3(114) 7.36(75) 7.83(82) 7.87(84) 7.87(79) 7.87(115) 8(77) 8.2(89) 8.3(85) 8.44(116) "
        "8.4±0.4(TW) 8.62(87) 9.2(117)",
        "15.16(115) 15.8(TW) 15.82(116)",
        "17.2(TW) 21.9(114) 22.16(80) 23.92(79)",
    ),
    "Tb": (
        "5.2(73) 7.18(75) 7.55(115) 7.6(84) 7.64(82) 7.74(79) 7.9(77) 8.43(87) 8.44(116) 9.1±0.5(TW)",
        "15.4(TW)",
        "20.4(TW) 23.43(79)",
    ),
    "Dy": (
        "5.6(73) 7.5(118) 7.53(84) 7.59(82) 7.72(79) 8.0(77) 8.37(87) 9.1±0.1(TW)",
        "16.1(TW)",
        "21.2(TW) 22.92(79)",
    ),
    "Ho": (
        "5.7(73) 6.8(75) 7.43(84) 7.56(82) 7.66(79) 7.85(74) 8(77) 8.14(119) 8.31(87) 8.4±0.5(TW)",
        "15.2(TW) 16.25(79)",
        "18.8(TW) 22.71(79)",
    ),
    "Er": (
        "5.5(73) 6.3(120) 6.56(75) 7.2±0.9(TW) 7.46(84) 7.52(82) 7.63(79) 7.9(77) 8.16(83) 8.26(87) "
        "8.4(88) 9(121)",
        "11.8(TW) 14.5(2) 16.15(79) 17.4(89)",
        "15.9(TW) 22.23(79) 23.1(2)",
    ),
    "Tm": (
        "4.4(73) 6.36(75) 7.34(84) 7.39(82) 7.51(79) 7.7(77) 8.22(87) 8.58(122) 8.9±0.0(TW)",
        "12.1(TW) 15.98(79)",
        "15.8(TW) 21.91(79)",
    ),
    "Yb": (
        "4.3(73) 6.24(75) 7.24(82) 7.31(84) 7.45(79) 7.7(123) 7.7(77) 8.0(85) 8.19(87) 8.4(124) "
        "8.2±0.3(TW) 8.6(92)",
        "15.1(TW) 15.5(123) 15.59(79)",
        "18.3(TW) 20.68(80) 21.72(79) 23.2(123)",
    ),
    "Lu": (
        "3.5(73) 6.04(75) 6.6(96) 7.27(82) 7.33(84) 7.41(79) 7.6(77) 7.7(85) 7.9±0.9(TW) 8.0(125) "
        "8.11(83) 8.17(87)",
        "14.0(TW) 15.21(79)",
        "17.1(TW) 21.43(79)",
    ),
}

TABLE_S1 = {
    "La": ("5s2 5p6", 1), "Ce": ("4f1", 2), "Pr": ("4f2", 3), "Nd": ("4f3", 4), "Pm": ("4f4", 5),
    "Sm": ("4f5", 6), "Eu": ("4f6", 7), "Gd": ("4f7", 8), "Tb": ("4f8", 7), "Dy": ("4f9", 6),
    "Ho": ("4f10", 5), "Er": ("4f11", 4), "Tm": ("4f12", 3), "Yb": ("4f13", 2), "Lu": ("4f14", 1),
}

# Table S2: pKa2/pKa3 MAD in the manuscript vs. with ma-def2-TZVPP in the optimizations.
TABLE_S2 = [
    ("La", 4.5, 5.5, 7.1, 8.6),
    ("Ce", 0.8, 5.1, 5.3, 5.9),
    ("Pr", None, None, 4.9, 9.6),
    ("Nd", 3.8, 4.0, 6.6, 5.7),
    ("Ho", 1.1, 2.1, 3.9, 2.7),
    ("Er", 4.2, 8.7, 6.8, 5.8),
    ("Lu", 1.2, 2.3, 4.3, 1.4),
]


def write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def fmt(value, digits: int = 4):
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return value


def parse_table2() -> list[dict]:
    rows = []
    for ln, steps in TABLE2.items():
        for step, text in enumerate(steps, 1):
            for order, token in enumerate(text.split()):
                match = re.fullmatch(r"([\d.]+)(?:±([\d.]+))?\((\w+)\)", token)
                if not match:
                    raise ValueError(f"cannot parse Table 2 token {token!r}")
                value, error, ref = match.groups()
                rows.append({
                    "Ln": ln, "step": step, "order": order, "pKa": float(value),
                    "error": float(error) if error else None,
                    "source": "this_work_calc" if ref == "TW" else "experiment",
                    "ref": "" if ref == "TW" else int(ref),
                })
    return rows


def mad_check(table2: list[dict]) -> list[list]:
    out = []
    for ln in LN:
        calc = TABLE1[ln]
        for step, (value, mad_paper) in enumerate(((calc[0], calc[2]), (calc[3], calc[4]), (calc[5], calc[6])), 1):
            exp = [row["pKa"] for row in table2 if row["Ln"] == ln and row["step"] == step and row["source"] == "experiment"]
            mad = float(np.mean(np.abs(np.asarray(exp) - value))) if exp else None
            exp_mean = float(np.mean(exp)) if exp else None
            exp_median = float(np.median(exp)) if exp else None
            ok = "" if mad is None or mad_paper is None else ("yes" if abs(mad - mad_paper) <= 0.051 else "NO")
            out.append([ln, step, value, len(exp), fmt(exp_mean, 3), fmt(exp_median, 3),
                        fmt(min(exp) if exp else None, 2), fmt(max(exp) if exp else None, 2),
                        fmt(mad, 3), fmt(mad_paper, 1), ok])
    return out


def ensure_si() -> Path:
    if not SI_PDF.is_file():
        SI_PDF.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(SI_URL, SI_PDF)
    return SI_PDF


def si_reader():
    import logging

    import pypdf

    logging.getLogger("pypdf").setLevel(logging.ERROR)
    return pypdf.PdfReader(str(ensure_si()))


def parse_structures(reader) -> list[dict]:
    text = "\n".join(page.extract_text() or "" for page in reader.pages)
    structures = []
    current = None
    for line in text.splitlines():
        header = re.match(r"^\s*\[(La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu)\((.*)$", line)
        if header:
            label = line.strip()
            n_oh = 0
            if label.startswith(f"[{header.group(1)}(OH"):
                count = re.match(r"\[\w+\(OH\)?(\d?)", label).group(1)
                n_oh = int(count) if count else 1
            current = {"Ln": header.group(1), "label": label, "n_oh_label": n_oh, "atoms": []}
            structures.append(current)
            continue
        atom = re.match(r"^\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", line)
        if atom and current is not None:
            current["atoms"].append((atom.group(1), *map(float, atom.groups()[1:])))
    if len(structures) != 60:
        raise ValueError(f"expected 60 SI structures, found {len(structures)}")
    return structures


def structure_metrics(structure: dict) -> dict:
    symbols = np.array([atom[0] for atom in structure["atoms"]])
    xyz = np.array([atom[1:] for atom in structure["atoms"]], dtype=float)
    ln_index = int(np.flatnonzero(symbols == structure["Ln"])[0])
    o_idx = np.flatnonzero(symbols == "O")
    h_idx = np.flatnonzero(symbols == "H")
    d_oh = np.linalg.norm(xyz[o_idx][:, None, :] - xyz[h_idx][None, :, :], axis=2)
    host = np.argmin(d_oh, axis=0)
    n_h = np.bincount(host, minlength=o_idx.size)
    d_ln_o = np.linalg.norm(xyz[o_idx] - xyz[ln_index], axis=1)
    water = n_h == 2
    hydroxide = n_h == 1
    covalent_oh = d_oh[np.arange(o_idx.size)[host], np.arange(h_idx.size)]
    hb_counts = []
    for k in np.flatnonzero(hydroxide):
        foreign = host != k
        hb_counts.append(int(np.count_nonzero(d_oh[k, foreign] < 2.5)))
    return {
        "Ln": structure["Ln"], "label": structure["label"], "n_atoms": int(symbols.size),
        "n_oh_label": structure["n_oh_label"], "n_oh_geom": int(np.count_nonzero(hydroxide)),
        "n_water": int(np.count_nonzero(water)), "n_other_O": int(np.count_nonzero(~(water | hydroxide))),
        "CN_O_3A": int(np.count_nonzero(d_ln_o < 3.0)),
        "Ln_O_water_mean": float(d_ln_o[water].mean()) if water.any() else None,
        "Ln_O_water_min": float(d_ln_o[water].min()) if water.any() else None,
        "Ln_O_water_max": float(d_ln_o[water].max()) if water.any() else None,
        "Ln_O_OH_mean": float(d_ln_o[hydroxide].mean()) if hydroxide.any() else None,
        "Ln_O_OH_min": float(d_ln_o[hydroxide].min()) if hydroxide.any() else None,
        "Ln_O_OH_max": float(d_ln_o[hydroxide].max()) if hydroxide.any() else None,
        "Ln_O_all_mean": float(d_ln_o.mean()),
        "O_H_covalent_max": float(covalent_oh.max()),
        "OH_O_Hwater_lt2.5_mean": float(np.mean(hb_counts)) if hb_counts else None,
        "symbols": symbols, "xyz": xyz,
    }


def write_structures(metrics: list[dict]) -> None:
    keys = ["Ln", "label", "n_atoms", "n_oh_label", "n_oh_geom", "n_water", "n_other_O", "CN_O_3A",
            "Ln_O_water_mean", "Ln_O_water_min", "Ln_O_water_max", "Ln_O_OH_mean", "Ln_O_OH_min",
            "Ln_O_OH_max", "Ln_O_all_mean", "O_H_covalent_max", "OH_O_Hwater_lt2.5_mean"]
    write_csv(OUT / "lit_SI_structures_metrics.csv", keys, [[fmt(m[k]) for k in keys] for m in metrics])
    for ln in ("Nd", "Dy"):
        lines = []
        for m in (m for m in metrics if m["Ln"] == ln):
            lines.append(str(m["symbols"].size))
            lines.append(f"{m['label']}  (Shiery et al. Inorg. Chem. 2021 SI, M06 optimized)")
            for symbol, (x, y, z) in zip(m["symbols"], m["xyz"]):
                lines.append(f"{symbol:2s} {x:16.10f} {y:16.10f} {z:16.10f}")
        (OUT / f"lit_SI_structures_{ln}.xyz").write_text("\n".join(lines) + "\n")


def digitize_pmfs(reader) -> list[list]:
    marker = re.compile(r"(-?[\d.]+) (-?[\d.]+) m((?: -?[\d.]+){6} c){4} h f")
    axis_x = re.compile(r"q ([\d.]+) 0 0 (-?[\d.]+) ([\d.]+) ([\d.]+) cm 0 0 m (\d+) 1 l S")
    axis_y = re.compile(r"q ([\d.]+) 0 0 (-?[\d.]+) ([\d.]+) ([\d.]+) cm 0 0 m 1 (\d+) l S")
    rows = []
    for page_index in range(2, 7):
        page = reader.pages[page_index]
        content = re.sub(r"\s+", " ", page.get_contents().get_data().decode("latin1"))
        text = page.extract_text()
        titles = re.findall(r"(\w\w) First Hydrolysis Potential of Mean Force", text)
        x_ticks = [list(map(float, m.split())) for m in re.findall(r"((?:\d\.\d\d ?){5,})", text)]
        x_axes = [(float(m[2]), float(m[3]), float(m[0]) * float(m[4])) for m in axis_x.findall(content)]
        y_axes = [(float(m[2]), float(m[3]), float(m[0]) * float(m[4])) for m in axis_y.findall(content)]
        points = []
        for m in marker.finditer(content):
            numbers = list(map(float, re.findall(r"-?[\d.]+", m.group(0))))
            xs, ys = numbers[0::2], numbers[1::2]
            points.append(((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2))
        if not (len(titles) == len(x_axes) == len(y_axes) == len(x_ticks) == 3):
            raise ValueError(f"page {page_index + 1}: unexpected chart layout")
        charts = sorted(zip(x_axes, y_axes), key=lambda pair: -pair[0][1])
        for title, ticks, ((x0, y0, length_x), (_, _, length_y)) in zip(titles, x_ticks, charts):
            mine = sorted(p for p in points if x0 - 1 <= p[0] <= x0 + length_x + 1 and y0 - 3 <= p[1] <= y0 + length_y + 3)
            r_min, r_max = ticks[0], ticks[-1]
            g_min, g_max = 0.0, 60.0
            for k, (px, py) in enumerate(mine):
                r = r_min + (px - x0) / length_x * (r_max - r_min)
                g = g_min + (py - y0) / length_y * (g_max - g_min)
                rows.append([title, f"S{LN.index(title) + 1}", k, round(r, 4), round(g, 3)])
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    references = {int(k): v for k, v in json.loads((SOURCE / "table2_references.json").read_text()).items()}

    write_csv(
        OUT / "lit_table1_calc_pKa.csv",
        ["Ln", "pKa1_calc", "pKa1_error", "pKa1_MAD", "pKa2_calc", "pKa2_MAD", "pKa3_calc", "pKa3_MAD",
         "dG1_kJmol_298K", "dG2_kJmol_298K", "dG3_kJmol_298K", "assumed_CN_aqua_ion"],
        [[ln, *[fmt(v, 1) for v in TABLE1[ln]],
          *[fmt(RT_LN10_298 * TABLE1[ln][i], 2) for i in (0, 3, 5)],
          9 if LN.index(ln) <= LN.index("Sm") else 8] for ln in LN],
    )

    table2 = parse_table2()
    write_csv(
        OUT / "lit_table2_pKa_all.csv",
        ["Ln", "step", "table_order", "pKa", "pKa_error", "source", "ref_number", "citation", "dG_kJmol_298K"],
        [[r["Ln"], r["step"], r["order"], r["pKa"], fmt(r["error"], 1), r["source"], r["ref"],
          references.get(r["ref"], "Shiery, Cooper, Cantu, Inorg. Chem. 2021, 60, 10257 (calculated)") if r["ref"] != "" else
          "Shiery, Cooper, Cantu, Inorg. Chem. 2021, 60, 10257 (calculated)",
          fmt(RT_LN10_298 * r["pKa"], 2)] for r in table2],
    )
    write_csv(
        OUT / "lit_table2_mad_check.csv",
        ["Ln", "step", "pKa_calc", "n_exp", "exp_mean", "exp_median", "exp_min", "exp_max",
         "MAD_recomputed", "MAD_table1", "agrees_within_0.05"],
        mad_check(table2),
    )
    used = sorted({r["ref"] for r in table2 if r["ref"] != ""})
    write_csv(OUT / "lit_references.csv", ["ref_number", "citation"], [[k, references[k]] for k in used])
    write_csv(OUT / "lit_tableS1_multiplicity.csv", ["Ln", "electronic_configuration", "multiplicity"],
              [[ln, *TABLE_S1[ln]] for ln in LN])
    write_csv(OUT / "lit_tableS2_basis_set_MAD.csv",
              ["Ln", "pKa2_MAD_manuscript", "pKa2_MAD_ma_def2_TZVPP_opt", "pKa3_MAD_manuscript", "pKa3_MAD_ma_def2_TZVPP_opt"],
              [[row[0], *[fmt(v, 1) for v in row[1:]]] for row in TABLE_S2])

    reader = si_reader()
    metrics = [structure_metrics(s) for s in parse_structures(reader)]
    for m in metrics:
        if m["n_oh_geom"] != m["n_oh_label"] or m["n_other_O"] or m["O_H_covalent_max"] > 1.1:
            raise ValueError(f"inconsistent SI structure {m['label']}: {m['n_oh_geom']} OH, max O-H {m['O_H_covalent_max']:.3f}")
    write_structures(metrics)

    pmf = digitize_pmfs(reader)
    first = {row[0]: row[4] for row in pmf if row[2] == 0}
    write_csv(OUT / "lit_SI_PMF_first_hydrolysis_digitized.csv",
              ["Ln", "SI_figure", "point", "r_OH_A", "dG_kJmol_marker", "dG_kJmol_rel_first_point"],
              [[*row, round(row[4] - first[row[0]], 3)] for row in pmf])
    check = []
    for ln in LN:
        pts = [row for row in pmf if row[0] == ln]
        end = pts[-1][4] - pts[0][4]
        check.append([ln, len(pts), pts[0][3], pts[-1][3], fmt(pts[0][4], 2), fmt(pts[-1][4], 2), fmt(end, 2),
                      fmt(end / RT_LN10_298, 2), TABLE1[ln][0]])
    write_csv(OUT / "lit_SI_PMF_endpoint_check.csv",
              ["Ln", "n_points", "r_first_A", "r_last_A", "dG_first", "dG_last", "dG_last_minus_first",
               "pKa_from_PMF_endpoint", "pKa1_table1"], check)

    meta = {
        "article": "Shiery, R. C.; Cooper, K. A.; Cantu, D. C. Computational Prediction of All Lanthanide Aqua Ion "
                   "Acidity Constants. Inorg. Chem. 2021, 60 (14), 10257-10266. DOI 10.1021/acs.inorgchem.1c00662",
        "article_url": ARTICLE_URL, "si_url": SI_URL, "dissertation_url": DISSERTATION_URL,
        "table_sources": {
            "Table 1/2": "transcribed from the dissertation Chapter 5 Tables 5-1/5-2 (reproduce the article tables); "
                         "verified by recomputing the Table 1 MADs (lit_table2_mad_check.csv)",
            "Table S1/S2": "SI p. S2 and S9",
            "structures": "SI pp. S10-S39 (pypdf text layer)",
            "PMF": "SI Figures S1-S15, marker centres read from the PDF vector paths; axes 0.95-1.45 A and 0-60 kJ/mol",
        },
        "conversion": f"dG = RT ln(10) pKa with RT ln(10) = {RT_LN10_298:.4f} kJ/mol at 298.15 K",
    }
    (OUT / "lit_metadata.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False))
    print(f"wrote literature tables to {OUT}")


if __name__ == "__main__":
    main()
