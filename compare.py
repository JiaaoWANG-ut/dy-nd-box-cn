#!/usr/bin/env python3
"""Compare the box-PBC MD Ln-OH results with Shiery, Cooper & Cantu, Inorg. Chem. 2021, 60, 10257.

Literature inputs are the transcribed tables in literature/ (Table 1, Table 2, Tables S1/S2, SI structures and
SI PMFs).  MD inputs are kinetics/*.csv (run `kinetics.py analyze` first), the total first-shell cutoffs in
cn_summary_300K.csv / solution_cn_box_pbc.csv, and kinetics_cache/ for the n_OH-resolved first-shell structure.

MD pKa-equivalents describe OH- binding in the box, pKa_eq = pKw - log10 K_OH with pKw = 14; the paper values are
298.15 K acidity constants from constrained (Bluemoon) AIMD in explicit water for pKa1 and from M06/CPCM
cluster-continuum relative free energies anchored on pKa1 for pKa2 and pKa3.  Writes comparison/*.csv and
comparison/*.svg.

Usage: python compare.py [workers]
"""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import kinetics as K

HERE = Path(__file__).resolve().parent
LIT = HERE / "literature"
KIN = HERE / "kinetics"
CMP = HERE / "comparison"
T_PAPER = 298.15
T_MD = 300.0
METALS = ("Nd", "Dy")
STEPS = (1, 2, 3)
MD_SERIES = (
    ("solution", "DPA4C-Plus"), ("graphite", "DPA4C-Plus"), ("graphite", "uma-s-1.2.1"),
    ("solution", "DPA4-Pro"), ("solution", "MACE-POLAR-1-M"), ("solution", "uma-s-1.2.1"),
)
SOURCE_CALC = "Shiery2021 Table 1 (DFT, 298.15 K)"
SOURCE_EXP = "Shiery2021 Table 2 (experiment, as listed)"
KRAGTEN_ND = (8.1, 16.2, 24.3)  # Table 2 ref 99 (Kragten & Decnop-Weever 1984), steps 1-3 as listed
STRUCT_START_PS = 10.0
REACTION = {
    1: "pKa1: Ln(3+) -> Ln(OH)(2+) + H+  (MD Ln-OH loss 1->0)",
    2: "pKa2: Ln(OH)(2+) -> Ln(OH)2(+) + H+  (MD 2->1)",
    3: "pKa3: Ln(OH)2(+) -> Ln(OH)3 + H+  (MD 3->2)",
}
STATUS_NOTE = {
    "measured": "both states visited with committed dwell >= 1 ps",
    "transient": "one state only visited transiently (committed dwell < 1 ps); low confidence",
    "upper_bound": "n-1 never visited; bound assumes true P(n-1) <= P(n)",
    "lower_bound": "n never visited; bound assumes true P(n) <= P(n-1)",
    "undetermined": "neither state visited",
}


def read(path: Path) -> list[dict]:
    return K.read_rows(path)


def num(x) -> float:
    return K.num(x)


def rnd(x, nd: int = 2):
    return "" if x is None or not math.isfinite(x) else round(float(x), nd)


def rtln10(temperature: float) -> float:
    return K.R_KJ * temperature * K.LN10


def label(system: str, model: str) -> str:
    return f"{K.MODEL_LABEL[model]} ({system})"


# --------------------------------------------------------------------------- inputs

def load_literature() -> dict:
    return {
        "t1": {r["Ln"]: r for r in read(LIT / "lit_table1_calc_pKa.csv")},
        "t2": read(LIT / "lit_table2_pKa_all.csv"),
        "mad": {(r["Ln"], int(r["step"])): r for r in read(LIT / "lit_table2_mad_check.csv")},
        "pmf": read(LIT / "lit_SI_PMF_first_hydrolysis_digitized.csv"),
        "pmf_end": {r["Ln"]: r for r in read(LIT / "lit_SI_PMF_endpoint_check.csv")},
        "si": read(LIT / "lit_SI_structures_metrics.csv"),
    }


def paper_pka(lit: dict, metal: str, step: int) -> float:
    return num(lit["t1"][metal][f"pKa{step}_calc"])


def exp_values(lit: dict, metal: str, step: int) -> list[float]:
    return [num(r["pKa"]) for r in lit["t2"]
            if r["Ln"] == metal and int(r["step"]) == step and r["source"] == "experiment"]


def md_rows(name: str, scenario: str = "primary") -> dict:
    out = {}
    for r in read(KIN / name):
        if r["scenario"] == scenario and r["variant"] in ("ALL", "box_pbc") and num(r["T_K"]) == T_MD:
            key = (r["system"], r["model"], r["metal"])
            if "hydrolysis_step" in r:
                key += (int(r["hydrolysis_step"]),)
            if "direction" in r:
                key += (r["direction"],)
            out[key] = r
    return out


def classify(row: dict) -> tuple[float, str, float]:
    """(pKa_eq, status, bound) of one kinetics_thermo row; bound is pKa_eq at P(n-1) = P(n)."""
    step = int(row["hydrolysis_step"])
    p_n, p_m = num(row["P_n"]), num(row["P_n_minus_1"])
    edge = K.PKW + math.log10((4 - step) * num(row["c_single_M"]))
    if p_n > 0 and p_m > 0:
        status = "measured" if row["both_states_committed"] == "True" else "transient"
        return num(row["pKa_equivalent"]), status, math.nan
    if p_n > 0:
        return math.nan, "upper_bound", edge
    if p_m > 0:
        return math.nan, "lower_bound", edge
    return math.nan, "undetermined", math.nan


# --------------------------------------------------------------------------- tables

def stepwise_rows(lit: dict, thermo: dict, thermo90: dict) -> list[dict]:
    rows = []
    for metal in METALS:
        for step in STEPS:
            calc = paper_pka(lit, metal, step)
            exps = exp_values(lit, metal, step)
            median = float(np.median(exps)) if exps else math.nan
            base = {"metal": metal, "hydrolysis_step": step, "reaction": REACTION[step]}

            def derived(pka: float, temperature: float) -> dict:
                return {
                    "log10_K_OH_binding": rnd(K.PKW - pka),
                    "dG_acid_kJmol": rnd(rtln10(temperature) * pka),
                    "dG0_OH_binding_kJmol": rnd(-rtln10(temperature) * (K.PKW - pka)),
                }

            err = num(lit["t1"][metal].get("pKa1_error")) if step == 1 else math.nan
            mad_t1 = num(lit["t1"][metal].get(f"pKa{step}_MAD"))
            paper_system = ("constrained AIMD, explicit water" if step == 1
                            else "cluster + implicit solvent (CPCM), relative to pKa1")
            rows.append({**base, "source": SOURCE_CALC, "system": paper_system,
                         "model": "AIMD-PBE" if step == 1 else "DFT-M06",
                         "status": "calc", "T_K": T_PAPER, "pKa": calc,
                         "pKa_CI95_lo": rnd(calc - err), "pKa_CI95_hi": rnd(calc + err), "pKa_bound": "",
                         **derived(calc, T_PAPER), "n_values": 1,
                         "dpKa_vs_calc": 0.0, "ddG_vs_calc_kJmol_298K": 0.0,
                         "dpKa_vs_exp_median": rnd(calc - median),
                         "MAD_vs_table2_as_listed": rnd(mad_t1),
                         "note": "pKa1 error as given in Table 1; MAD as given in Table 1"})
            if exps:
                rows.append({**base, "source": SOURCE_EXP, "system": "experiment", "model": "",
                             "status": "experiment", "T_K": T_PAPER, "pKa": rnd(median),
                             "pKa_CI95_lo": rnd(min(exps)), "pKa_CI95_hi": rnd(max(exps)), "pKa_bound": "",
                             **derived(median, T_PAPER), "n_values": len(exps),
                             "dpKa_vs_calc": rnd(median - calc),
                             "ddG_vs_calc_kJmol_298K": rnd(rtln10(T_PAPER) * (median - calc)),
                             "dpKa_vs_exp_median": 0.0, "MAD_vs_table2_as_listed": "",
                             "note": "median (lo/hi = min/max) of Table 2 entries for this step; steps 2/3 are "
                                     "numerically consistent with cumulative -log*beta_n"})
            for system, model in MD_SERIES:
                row = thermo.get((system, model, metal, step))
                if row is None:
                    continue
                pka, status, bound = classify(row)
                alt = thermo90.get((system, model, metal, step)) if system == "graphite" else None
                alt_pka = classify(alt)[0] if alt else math.nan
                ref = pka if math.isfinite(pka) else bound
                rows.append({
                    **base, "source": "this work (MD, box PBC)", "system": system, "model": model,
                    "status": status, "T_K": T_MD, "pKa": rnd(pka),
                    "pKa_CI95_lo": rnd(num(row["pKa_equivalent_CI95_lo"])) if math.isfinite(pka) else "",
                    "pKa_CI95_hi": rnd(num(row["pKa_equivalent_CI95_hi"])) if math.isfinite(pka) else "",
                    "pKa_bound": rnd(bound), "pKa_t_ge_90ps": rnd(alt_pka),
                    **(derived(pka, T_MD) if math.isfinite(pka) else {}),
                    "dG0_dissociation_kJmol": rnd(num(row["dG0_dissociation_kJmol"])) if math.isfinite(pka) else "",
                    "P_n": row["P_n"], "P_n_minus_1": row["P_n_minus_1"],
                    "n_values": row["n_trajectories"],
                    "dpKa_vs_calc": rnd(ref - calc),
                    "ddG_vs_calc_kJmol_298K": rnd(rtln10(T_PAPER) * (ref - calc)),
                    "dpKa_vs_exp_median": rnd(ref - median),
                    "MAD_vs_table2_as_listed": rnd(float(np.mean([abs(pka - e) for e in exps])))
                    if exps and math.isfinite(pka) else "",
                    "note": STATUS_NOTE[status] + ("; differences use the bound" if not math.isfinite(pka)
                                                   and math.isfinite(bound) else ""),
                })
    return rows


def cumulative_rows(lit: dict, thermo: dict) -> list[dict]:
    rows = []
    for metal in METALS:
        for n in STEPS:
            calc = sum(paper_pka(lit, metal, s) for s in range(1, n + 1))
            base = {"metal": metal, "n_OH": n, "quantity": f"-log*beta_{n} (cumulative, Ln3+ + n H2O)"}
            rows.append({**base, "source": SOURCE_CALC + ", summed", "value": rnd(calc), "lo": "", "hi": "",
                         "diff_vs_calc_cumulative": 0.0, "note": "sum of stepwise pKa1..pKa_n"})
            exps = exp_values(lit, metal, n)
            if exps:
                med = float(np.median(exps))
                rows.append({**base, "source": SOURCE_EXP, "value": rnd(med), "lo": rnd(min(exps)),
                             "hi": rnd(max(exps)), "diff_vs_calc_cumulative": rnd(med - calc),
                             "note": f"median of {len(exps)} Table 2 step-{n} entries, read as cumulative"})
            if metal == "Nd":
                rows.append({**base, "source": "Table 2 ref 99 (Kragten 1984) as listed",
                             "value": KRAGTEN_ND[n - 1], "lo": "", "hi": "",
                             "diff_vs_calc_cumulative": rnd(KRAGTEN_ND[n - 1] - calc),
                             "note": "8.1/16.2/24.3 = 8.1 x n, i.e. cumulative constants"})
            for system, model in MD_SERIES:
                parts = [thermo.get((system, model, metal, s)) for s in range(1, n + 1)]
                if any(p is None for p in parts):
                    continue
                cls = [classify(p) for p in parts]
                if not all(math.isfinite(c[0]) for c in cls):
                    continue
                value = sum(c[0] for c in cls)
                lo = value - math.sqrt(sum((c[0] - num(p["pKa_equivalent_CI95_lo"])) ** 2 for c, p in zip(cls, parts)))
                hi = value + math.sqrt(sum((num(p["pKa_equivalent_CI95_hi"]) - c[0]) ** 2 for c, p in zip(cls, parts)))
                weakest = "transient" if any(c[1] == "transient" for c in cls) else "measured"
                rows.append({**base, "source": f"this work: {label(system, model)}", "value": rnd(value),
                             "lo": rnd(lo), "hi": rnd(hi), "diff_vs_calc_cumulative": rnd(value - calc),
                             "note": f"sum of stepwise pKa_eq ({weakest}); CI combined in quadrature (approximate)"})
    return rows


def predicted_occupancy(pkas: tuple[float, float, float], c1: float) -> list[float]:
    ratio = {k: 10.0 ** (pkas[k - 1] - K.PKW) / ((4 - k) * c1) for k in STEPS}  # P(k-1)/P(k)
    weights = [ratio[3] * ratio[2] * ratio[1], ratio[3] * ratio[2], ratio[3], 1.0]
    total = sum(weights)
    return [w / total for w in weights]


def occupancy_rows(lit: dict, thermo: dict) -> list[dict]:
    c1 = {}
    for (system, *_), row in thermo.items():
        c1[system] = num(row["c_single_M"])
    occ = md_rows("kinetics_occupancy.csv")
    rows = []
    for metal in METALS:
        for system in ("solution", "graphite"):
            sources = [(SOURCE_CALC, tuple(paper_pka(lit, metal, s) for s in STEPS))]
            if metal == "Nd":
                sources += [("Kragten 1984 (Table 2 ref 99) as listed", KRAGTEN_ND),
                            ("Kragten 1984 read as cumulative (stepwise 8.1/8.1/8.1)", (8.1, 8.1, 8.1))]
            for source, pkas in sources:
                probs = predicted_occupancy(pkas, c1[system])
                rows.append({"metal": metal, "system": system, "c_single_M": rnd(c1[system], 4),
                             "source": source, "kind": "predicted from pKa",
                             "pKa1": pkas[0], "pKa2": pkas[1], "pKa3": pkas[2],
                             **{f"P{n}": f"{p:.3g}" for n, p in enumerate(probs)},
                             "mean_n_OH": rnd(sum(n * p for n, p in enumerate(probs)), 3)})
            for sys_md, model in MD_SERIES:
                row = occ.get((sys_md, model, metal))
                if sys_md != system or row is None:
                    continue
                rows.append({"metal": metal, "system": system, "c_single_M": rnd(c1[system], 4),
                             "source": f"this work: {label(system, model)}", "kind": "MD occupancy (300 K)",
                             **{f"P{n}": row[f"P{n}"] for n in range(4)},
                             "mean_n_OH": row["mean_n_OH_bound"]})
    return rows


# --------------------------------------------------------------------------- n_OH-resolved structure

def structure_job(args: tuple) -> tuple[tuple, dict]:
    key, path, start, end, rc_oh, r_shell = args
    with np.load(path) as data:
        t = data["t_ps"]
        sel = (t >= start - 1e-9) & (t <= end + 1e-9)
        near_d = data["near_d"][sel].astype(np.float64)
        near_type = data["near_type"][sel]
        oh_d = data["oh_d"][sel]
    n_oh = (oh_d < rc_oh).sum(axis=1)
    shell = near_d < r_shell
    water = shell & (near_type == 0)
    hydroxide = (near_d < rc_oh) & (near_type == 1)
    other = shell & (near_type == 2)
    acc = {}
    for n in np.unique(n_oh):
        m = n_oh == n
        d_w, d_h, d_a = near_d[m][water[m]], near_d[m][hydroxide[m]], near_d[m][shell[m]]
        acc[int(n)] = np.array([m.sum(), water[m].sum(), hydroxide[m].sum(), other[m].sum(), shell[m].sum(),
                                d_w.sum(), (d_w ** 2).sum(), d_w.size, d_h.sum(), (d_h ** 2).sum(), d_h.size,
                                d_a.sum(), d_a.size], dtype=np.float64)
    return key, acc


def structure_rows(lit: dict, workers: int) -> list[dict]:
    rc_oh = {(r["system"], r["model"], r["metal"]): num(r["rc_A"]) for r in read(KIN / "kinetics_cutoffs.csv")}
    shell = {}
    for r in read(HERE / "cn_summary_300K.csv"):
        if r["variant"] == "ALL":
            shell[("graphite", r["model"], r["metal"], 300)] = num(r["rc_REE_O_A"])
    for r in read(HERE / "solution_cn_box_pbc.csv"):
        if r["collapsed"] == "0":
            shell[("solution", r["model"], r["metal"], int(num(r["T_K"])))] = num(r["rc_Ln_O_total_A"])
    screen = K.load_screen()
    tasks = []
    for job in K.all_jobs():
        key = (job["system"], job["model"], job["metal"], int(job["T_K"]))
        srow = screen[(job["system"], job["model"], job["variant"], job["metal"], int(job["replica"]), int(job["T_K"]))]
        if srow["status"] == "excluded" or key not in shell or not math.isfinite(shell[key]):
            continue
        tasks.append((key, str(K.cache_path(job)), STRUCT_START_PS, float(srow["usable_end_ps"]),
                      rc_oh[key[:3]], shell[key]))
    totals: dict = defaultdict(lambda: defaultdict(lambda: np.zeros(13)))
    n_traj: dict = defaultdict(int)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for key, acc in pool.map(structure_job, tasks, chunksize=8):
            n_traj[key] += 1
            for n, vec in acc.items():
                totals[key][n] += vec
    rows = []
    for key in sorted(totals):
        system, model, metal, temperature = key
        frames_all = sum(v[0] for v in totals[key].values())
        for n in sorted(totals[key]):
            v = totals[key][n]
            frames = v[0]

            def mean_sd(s: float, s2: float, cnt: float) -> tuple:
                if cnt <= 0:
                    return "", ""
                mu = s / cnt
                return rnd(mu, 3), rnd(math.sqrt(max(s2 / cnt - mu * mu, 0.0)), 3)

            w_mu, w_sd = mean_sd(v[5], v[6], v[7])
            h_mu, h_sd = mean_sd(v[8], v[9], v[10])
            rows.append({
                "source": "this work (MD)", "system": system, "model": model, "metal": metal, "T_K": temperature,
                "label": f"Ln(OH){n} state", "n_OH": n, "n_trajectories": n_traj[key], "frames": int(frames),
                "fraction": rnd(frames / frames_all, 4), "CN_water": rnd(v[1] / frames, 3),
                "CN_OH": rnd(v[2] / frames, 3), "CN_other": rnd(v[3] / frames, 3), "CN_total": rnd(v[4] / frames, 3),
                "Ln_O_water_mean_A": w_mu, "Ln_O_water_sd_A": w_sd, "Ln_O_OH_mean_A": h_mu, "Ln_O_OH_sd_A": h_sd,
                "Ln_O_all_mean_A": rnd(v[11] / v[12], 3) if v[12] > 0 else "",
                "rc_OH_A": rc_oh[key[:3]], "r_shell_A": rnd(shell[key], 3),
            })
    for r in lit["si"]:
        if r["Ln"] not in METALS:
            continue
        n = int(r["n_oh_geom"])
        rows.append({
            "source": "Shiery2021 SI optimized cluster", "system": "cluster", "model": "DFT", "metal": r["Ln"],
            "T_K": "", "label": r["label"], "n_OH": n, "n_trajectories": "", "frames": "", "fraction": "",
            "CN_water": r["n_water"], "CN_OH": n, "CN_other": r["n_other_O"], "CN_total": r["CN_O_3A"],
            "Ln_O_water_mean_A": rnd(num(r["Ln_O_water_mean"]), 3), "Ln_O_water_sd_A": "",
            "Ln_O_OH_mean_A": rnd(num(r["Ln_O_OH_mean"]), 3), "Ln_O_OH_sd_A": "",
            "Ln_O_all_mean_A": rnd(num(r["Ln_O_all_mean"]), 3), "rc_OH_A": "", "r_shell_A": 3.0,
        })
    return rows


# --------------------------------------------------------------------------- PMF and activation

def pmf_rows(lit: dict, thermo: dict, rates: dict) -> list[dict]:
    rows = []
    for metal in METALS:
        pts = sorted((num(r["r_OH_A"]), num(r["dG_kJmol_rel_first_point"])) for r in lit["pmf"] if r["Ln"] == metal)
        drops = [b[1] - a[1] for a, b in zip(pts, pts[1:])]
        end = lit["pmf_end"][metal]
        rows.append({"metal": metal, "source": "Shiery2021 SI PMF (step 1, r(O-H) 0.99->1.40 A)",
                     "quantity": "dG(end) - dG(start) = deprotonation free energy", "value_kJmol": end["dG_last_minus_first"],
                     "lo": "", "hi": "", "T_K": T_PAPER,
                     "note": f"monotonic ({len(pts)} points, smallest increment {min(drops):.2f} kJ/mol): "
                             "no barrier, so no activation energy is defined; pKa from endpoint "
                             f"{end['pKa_from_PMF_endpoint']}"})
        rows.append({"metal": metal, "source": SOURCE_CALC, "quantity": "dG1 (acid, step 1)",
                     "value_kJmol": lit["t1"][metal]["dG1_kJmol_298K"], "lo": "", "hi": "", "T_K": T_PAPER,
                     "note": "Table 1"})
        for system, model in MD_SERIES:
            row = thermo.get((system, model, metal, 1))
            if row is None:
                continue
            pka, status, bound = classify(row)
            if math.isfinite(pka):
                rows.append({"metal": metal, "source": f"this work: {label(system, model)}",
                             "quantity": "dG1 acid-equivalent = RT ln10 pKa1_eq", "value_kJmol": rnd(rtln10(T_MD) * pka),
                             "lo": rnd(rtln10(T_MD) * num(row["pKa_equivalent_CI95_lo"])),
                             "hi": rnd(rtln10(T_MD) * num(row["pKa_equivalent_CI95_hi"])), "T_K": T_MD,
                             "note": STATUS_NOTE[status]})
            for direction, what in (("dissociation", "dG_act Ln-OH loss 1->0 (Ln-OH + H2O -> Ln-OH2 + OH-)"),
                                    ("association", "dG_act OH uptake 0->1")):
                rate = rates.get((system, model, metal, 1, direction))
                if rate is None or rate["dG_act_kJmol"] in ("", "inf"):
                    continue
                rows.append({"metal": metal, "source": f"this work: {label(system, model)}", "quantity": what,
                             "value_kJmol": rate["dG_act_kJmol"], "lo": rate["dG_act_CI95_lo"],
                             "hi": rate["dG_act_CI95_hi"], "T_K": T_MD,
                             "note": f"apparent Eyring barrier, N = {rate['N_events']} events; no paper counterpart"})
    return rows


def activation_rows(lit: dict, rates: dict, arr: dict) -> list[dict]:
    rows = []
    for (system, model, metal, step, direction), r in sorted(rates.items()):
        if metal not in METALS:
            continue
        a = arr.get((system, model, metal, step, direction), {})
        has_ea = a.get("Ea_kJmol", "") not in ("", None)
        rows.append({
            "system": system, "model": model, "metal": metal, "hydrolysis_step": step,
            "md_transition": K.STEP_LABEL[step].replace("→", "->"), "direction": direction,
            "process": "Ln-OH loss (n -> n-1)" if direction == "dissociation" else "OH- uptake (n-1 -> n)",
            "N_events_300K": r["N_events"], "dwell_origin_state_ps": r["T_origin_state_ps"],
            "k_300K_ps-1": r["k_ps-1"], "dG_act_300K_kJmol": r["dG_act_kJmol"],
            "dG_act_CI95_lo": r["dG_act_CI95_lo"], "dG_act_CI95_hi": r["dG_act_CI95_hi"],
            "dG_act_lower_bound_kJmol": r["dG_act_lower_bound_if_no_events"],
            "frac_proton_transfer_like": r["frac_proton_transfer_like"],
            "Ea_kJmol": a.get("Ea_kJmol", "") if has_ea else "", "Ea_CI95_lo": a.get("Ea_CI95_lo", "") if has_ea else "",
            "Ea_CI95_hi": a.get("Ea_CI95_hi", "") if has_ea else "", "Ea_quality": a.get("quality", ""),
            "Ea_T_used": a.get("T_used", ""), "N_per_T": a.get("N_per_T", ""),
            "dH_act_kJmol": a.get("dH_act_kJmol", "") if has_ea else "",
            "dS_act_JmolK": a.get("dS_act_JmolK", "") if has_ea else "",
            "paper_dG_acid_step_kJmol_298K": lit["t1"][metal][f"dG{step}_kJmol_298K"],
            "paper_barrier": "none reported; step-1 PMF is monotonic (thermodynamics only)",
        })
    return rows


# --------------------------------------------------------------------------- figures

def series_style(system: str, model: str) -> dict:
    return {"color": K.SERIES_COLOR[(system, model)], "marker": "o" if system == "solution" else "s"}


def fig_pka(lit: dict, step_rows: list[dict]) -> None:
    plt = K.pyplot()
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6), sharey=True)
    md = {(r["metal"], r["hydrolysis_step"], r["system"], r["model"]): r for r in step_rows
          if r["source"].startswith("this work")}
    for ax, metal in zip(axes, METALS):
        for step in STEPS:
            exps = exp_values(lit, metal, step)
            if exps:
                xs = step - 0.34 + np.linspace(-0.05, 0.05, len(exps))
                ax.scatter(xs, exps, s=14, color="#9a9a9a", zorder=3)
        calc = [paper_pka(lit, metal, s) for s in STEPS]
        err = num(lit["t1"][metal]["pKa1_error"])
        ax.plot([s - 0.24 for s in STEPS], calc, "-D", color="black", ms=6, zorder=4)
        ax.errorbar(1 - 0.24, calc[0], yerr=err, color="black", capsize=3, lw=1)
        for i, (system, model) in enumerate(MD_SERIES):
            st = series_style(system, model)
            for step in STEPS:
                r = md.get((metal, step, system, model))
                if r is None:
                    continue
                x = step - 0.12 + i * 0.09
                if r["status"] in ("measured", "transient"):
                    y, lo, hi = num(r["pKa"]), num(r["pKa_CI95_lo"]), num(r["pKa_CI95_hi"])
                    lo = lo if math.isfinite(lo) else y
                    hi = hi if math.isfinite(hi) else y
                    ax.errorbar(x, y, yerr=[[y - lo], [hi - y]], fmt=st["marker"], color=st["color"], ms=6,
                                mfc=st["color"] if r["status"] == "measured" else "white", capsize=2, lw=1, zorder=5)
                elif r["status"] in ("upper_bound", "lower_bound"):
                    y = num(r["pKa_bound"])
                    sign = -1 if r["status"] == "upper_bound" else 1
                    ax.annotate("", xy=(x, y + sign * 2.2), xytext=(x, y),
                                arrowprops={"arrowstyle": "-|>", "color": st["color"], "lw": 1.2})
                    ax.plot(x, y, "_", color=st["color"], ms=10, mew=2)
        ax.axhline(K.PKW, color="#888", ls="--", lw=0.8)
        ax.set_xticks(STEPS, ["pKa1\n(MD 1→0)", "pKa2\n(MD 2→1)", "pKa3\n(MD 3→2)"])
        ax.set_xlim(0.45, 3.55)
        ax.set_ylim(4, 28)
        ax.set_title(f"{metal}: stepwise pKa (paper) vs pKa-equivalent from box occupancy (MD, 300 K)")
        ax.grid(axis="y", color="#e8ecf2", lw=0.7)
    axes[0].set_ylabel("pKa  (MD: 14 + log10 K_d, OH- dissociation)")
    sec = axes[1].secondary_yaxis("right", functions=(lambda p: p * rtln10(T_PAPER), lambda g: g / rtln10(T_PAPER)))
    sec.set_ylabel("ΔG_acid = RT ln10 pKa (kJ/mol, 298 K)")
    handles = [Line2D([], [], color="#9a9a9a", marker="o", ls="", label="Table 2 experiment (each entry, as listed)"),
               Line2D([], [], color="black", marker="D", label="Table 1 DFT (this paper)")]
    handles += [Line2D([], [], color=series_style(*s)["color"], marker=series_style(*s)["marker"], ls="",
                       label=label(*s)) for s in MD_SERIES]
    handles += [Line2D([], [], color="#555", marker="o", mfc="white", ls="", label="MD transient (low confidence)"),
                Line2D([], [], color="#555", marker="_", ms=10, mew=2, ls="", label="MD bound (arrow: direction)"),
                Line2D([], [], color="#888", ls="--", lw=0.8, label="pKa = pKw = 14 (K_d = 1 M)")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(CMP / "fig_compare_pKa_stepwise.svg")
    plt.close(fig)


def fig_cumulative(lit: dict, cum_rows: list[dict]) -> None:
    plt = K.pyplot()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, metal in zip(axes, METALS):
        calc = np.cumsum([paper_pka(lit, metal, s) for s in STEPS])
        ax.plot(STEPS, calc, "-D", color="black", label="Table 1 DFT, cumulative sum")
        ax.plot(STEPS, [paper_pka(lit, metal, s) for s in STEPS], ":D", color="black", mfc="white",
                label="Table 1 DFT, stepwise")
        for n in STEPS:
            exps = exp_values(lit, metal, n)
            if exps:
                ax.scatter(n + np.linspace(-0.05, 0.05, len(exps)), exps, s=14, color="#9a9a9a", zorder=3,
                           label="Table 2 experiment (as listed)" if n == 1 else None)
        for system, model in MD_SERIES:
            pts = [(r["n_OH"], num(r["value"]), num(r["lo"]), num(r["hi"])) for r in cum_rows
                   if r["metal"] == metal and r["source"] == f"this work: {label(system, model)}"]
            if not pts:
                continue
            st = series_style(system, model)
            xs = [p[0] + (0.08 if system == "graphite" else -0.08) for p in pts]
            ax.errorbar(xs, [p[1] for p in pts], yerr=[[p[1] - p[2] for p in pts], [p[3] - p[1] for p in pts]],
                        fmt=st["marker"] + "-", color=st["color"], ms=5, capsize=2, lw=1, label=label(system, model))
        ax.set_xticks(STEPS, ["n = 1", "n = 2", "n = 3"])
        ax.set_title(f"{metal}: cumulative -log*β_n")
        ax.grid(axis="y", color="#e8ecf2", lw=0.7)
    axes[0].set_ylabel("-log*β_n  (sum of stepwise pKa)")
    axes[0].legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(CMP / "fig_compare_pKa_cumulative.svg")
    plt.close(fig)


def fig_occupancy(occ_rows: list[dict]) -> None:
    plt = K.pyplot()
    colors = ["#fde0c5", "#facba6", "#f59e72", "#d9534f"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    for ax, metal in zip(axes, METALS):
        rows = [r for r in occ_rows if r["metal"] == metal]
        names = []
        for i, r in enumerate(rows):
            bottom = 0.0
            for n in range(4):
                p = num(r[f"P{n}"])
                ax.bar(i, p, 0.75, bottom=bottom, color=colors[n], edgecolor="white", lw=0.4,
                       label=f"n = {n}" if i == 0 else None)
                bottom += p
            src = r["source"].replace("this work: ", "MD ").replace(SOURCE_CALC, "paper Table 1")
            src = src.replace("Kragten 1984 (Table 2 ref 99) as listed", "Kragten as listed")
            src = src.replace("Kragten 1984 read as cumulative (stepwise 8.1/8.1/8.1)", "Kragten cumulative")
            names.append(src if src.startswith("MD ") else f"{src}\n[{r['system']} box]")
        ax.set_xticks(range(len(rows)), names, rotation=90, fontsize=6.5)
        ax.set_title(f"{metal}: P(n bound OH) in our box, predicted from pKa vs MD")
    axes[0].set_ylabel("P(n)")
    axes[0].set_ylim(0, 1)
    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_, fontsize=8, ncol=4, frameon=False, loc="upper center")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(CMP / "fig_compare_box_occupancy.svg")
    plt.close(fig)


def fig_structure(struct_rows: list[dict]) -> None:
    plt = K.pyplot()
    from matplotlib.lines import Line2D
    metrics = (("Ln_O_OH_mean_A", "Ln-O(OH) mean (Å)"), ("Ln_O_water_mean_A", "Ln-O(water) mean (Å)"),
               ("CN_total", "first-shell CN (all O)"))
    fig, axes = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
    for j, metal in enumerate(METALS):
        si = sorted((r for r in struct_rows if r["source"].startswith("Shiery") and r["metal"] == metal),
                    key=lambda r: r["n_OH"])
        for i, (col, ylabel) in enumerate(metrics):
            ax = axes[i][j]
            pts = [(r["n_OH"], num(r[col])) for r in si if math.isfinite(num(r[col]))]
            ax.plot(*zip(*pts), "-D", color="black", ms=6, zorder=5)
            for k, (system, model) in enumerate(MD_SERIES):
                st = series_style(system, model)
                rows = [r for r in struct_rows if r["source"] == "this work (MD)" and r["system"] == system
                        and r["model"] == model and r["metal"] == metal and r["T_K"] == 300
                        and num(r["fraction"]) >= 0.005 and r["n_OH"] <= 3]
                xs = [r["n_OH"] + (k - 2.5) * 0.06 for r in rows if math.isfinite(num(r[col]))]
                ys = [num(r[col]) for r in rows if math.isfinite(num(r[col]))]
                sd_col = col.replace("_mean_A", "_sd_A")
                sds = [num(r[sd_col]) if sd_col != col and sd_col in r else 0.0
                       for r in rows if math.isfinite(num(r[col]))]
                if xs:
                    ax.errorbar(xs, ys, yerr=[s if math.isfinite(s) else 0 for s in sds], fmt=st["marker"],
                                color=st["color"], ms=5, capsize=1.5, lw=0.7, alpha=0.9)
            ax.set_ylabel(ylabel)
            ax.grid(color="#eef1f5", lw=0.6)
            if i == 0:
                ax.set_title(f"{metal}: SI cluster (black) vs MD state-resolved, 300 K (bars = thermal SD)",
                             fontsize=10)
            if i == 2:
                ax.set_xlabel("n_OH bound")
                ax.set_xticks(range(4))
    handles = [Line2D([], [], color="black", marker="D", label="Shiery2021 SI optimized cluster")]
    handles += [Line2D([], [], color=series_style(*s)["color"], marker=series_style(*s)["marker"], ls="",
                       label=label(*s)) for s in MD_SERIES]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(CMP / "fig_compare_structure.svg")
    plt.close(fig)


def fig_pmf(lit: dict, pmf: list[dict]) -> None:
    plt = K.pyplot()
    fig, ax = plt.subplots(figsize=(8.5, 5))
    curves = defaultdict(list)
    for r in lit["pmf"]:
        curves[r["Ln"]].append((num(r["r_OH_A"]), num(r["dG_kJmol_rel_first_point"])))
    for ln, pts in curves.items():
        if ln not in METALS:
            ax.plot(*zip(*sorted(pts)), color="#cccccc", lw=0.8, zorder=1)
    for metal, color in (("Nd", "#d95f02"), ("Dy", "#0072b2")):
        ax.plot(*zip(*sorted(curves[metal])), "-o", color=color, ms=4, label=f"{metal} PMF (paper SI)", zorder=3)
    md = [r for r in pmf if r["quantity"].startswith("dG1 acid-equivalent")]
    label_y = {}
    for k in sorted(range(len(md)), key=lambda k: num(md[k]["value_kJmol"])):
        y = num(md[k]["value_kJmol"])
        prev = max(label_y.values(), default=-1e9)
        label_y[k] = max(y, prev + 3.2)
    for k, r in enumerate(md):
        color = "#d95f02" if r["metal"] == "Nd" else "#0072b2"
        system = "solution" if "(solution)" in r["source"] else "graphite"
        x = 1.43 + 0.012 * k
        y, lo, hi = num(r["value_kJmol"]), num(r["lo"]), num(r["hi"])
        ax.errorbar(x, y, yerr=[[y - lo], [hi - y]], fmt="o" if system == "solution" else "s", color=color,
                    mfc="white" if "transient" in r["note"] else color, ms=6, capsize=2)
        text = f"{r['metal']} " + r["source"].replace("this work: ", "")
        if "transient" in r["note"]:
            text += " [transient]"
        ax.annotate(text, (x, y), xytext=(1.47 + 0.012 * len(md), label_y[k]), textcoords="data",
                    fontsize=6.5, color=color, va="center",
                    arrowprops={"arrowstyle": "-", "color": color, "lw": 0.5, "alpha": 0.6})
    ax.axvline(1.415, color="#bbb", lw=0.6, ls=":")
    ax.set_xlim(0.95, 1.74)
    ax.set_ylim(-2, 90)
    ax.set_xlabel("r(O-H) of the deprotonating Ln-bound water (Å)")
    ax.set_ylabel("ΔG (kJ/mol)")
    ax.set_title("Step-1 PMF (paper SI; grey = other Ln) vs MD acid-equivalent ΔG1 (right, 300 K)", fontsize=9)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.grid(color="#eef1f5", lw=0.6)
    fig.tight_layout()
    fig.savefig(CMP / "fig_compare_PMF_step1.svg")
    plt.close(fig)


def main(workers: int) -> None:
    CMP.mkdir(exist_ok=True)
    lit = load_literature()
    thermo = md_rows("kinetics_thermo.csv")
    thermo90 = md_rows("kinetics_thermo.csv", "start90")
    rates = md_rows("kinetics_rates.csv")
    arr = {}
    for r in read(KIN / "kinetics_arrhenius.csv"):
        if r["scenario"] == "primary":
            arr[(r["system"], r["model"], r["metal"], int(r["hydrolysis_step"]), r["direction"])] = r

    step = stepwise_rows(lit, thermo, thermo90)
    cum = cumulative_rows(lit, thermo)
    occ = occupancy_rows(lit, thermo)
    pmf = pmf_rows(lit, thermo, rates)
    act = activation_rows(lit, rates, arr)
    print("computing n_OH-resolved structure ...", flush=True)
    struct = structure_rows(lit, workers)
    for name, rows in (("compare_pKa_stepwise.csv", step), ("compare_pKa_cumulative.csv", cum),
                       ("compare_box_occupancy.csv", occ), ("compare_structure_by_nOH.csv", struct),
                       ("compare_PMF_step1.csv", pmf), ("compare_activation.csv", act)):
        K.write_rows(CMP / name, rows)
    fig_pka(lit, step)
    fig_cumulative(lit, cum)
    fig_occupancy(occ)
    fig_structure(struct)
    fig_pmf(lit, pmf)
    print("comparison written to", CMP, flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
