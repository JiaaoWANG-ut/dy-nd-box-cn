#!/usr/bin/env python3
"""Stepwise Ln-OH (de)coordination kinetics for Dy/Nd (solution box_pbc and graphite).

Stages (run in order):
  extract  dense 20 fs per-trajectory time series -> kinetics_cache/<group>/<run>.npz
  screen   per-trajectory collapse onset            -> collapse_screen.csv
  analyze  committed transitions, rates, Eyring/Arrhenius and thermodynamics -> kinetics/*.csv
  figures  SVG figures (+ plotted data) from the analyze tables and the cache -> kinetics/*.svg

Step n->n-1 is the loss of one Ln-bound hydroxide:
  Ln(OH)n(3-n)+  ->  Ln(OH)(n-1)(4-n)+  +  OH-
so 1->0, 2->1 and 3->2 are the reverse of the first, second and third hydrolysis steps.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
CACHE = OUT / "kinetics_cache"
GRAPHITE = ROOT / "Graphite-dy-nd-50-denser-dumping"
SOLUTION_RUNS = {
    "uma-s-1.2.1": ROOT / "md" / "ln-150ps-denser-dumping" / "runs",
    "DPA4C-Plus": ROOT / "md" / "ln-150ps-denser-dumping" / "runs",
    "DPA4-Pro": ROOT / "md" / "ln-150ps-every-step-dumping" / "runs",
    "MACE-POLAR-1-M": ROOT / "md" / "ln-150ps-every-step-dumping" / "runs",
}
SAMPLE_FS = 20.0
CHUNK = 500
N_NEAR_O = 12
N_OH_KEEP = 4


def solution_jobs() -> list[dict]:
    jobs = []
    for model, runs in SOLUTION_RUNS.items():
        for metal in ("Dy", "Nd"):
            for temperature in (300, 350, 400, 450, 500):
                base = runs / model / "box_pbc" / metal / f"T{temperature}"
                prefix = next(base.glob("traj-*-meta.json")).name[: -len("-meta.json")]
                jobs.append({
                    "system": "solution", "model": model, "metal": metal, "variant": "box_pbc",
                    "replica": 0, "T_K": temperature, "dir": str(base), "prefix": prefix,
                })
    return jobs


def graphite_jobs() -> list[dict]:
    jobs = []
    for job in json.loads((GRAPHITE / "jobs.json").read_text())["jobs"]:
        base = GRAPHITE / "runs" / job["model"] / job["variant"] / job["metal"] / f"r{job['replica']:02d}"
        jobs.append({
            "system": "graphite", "model": job["model"], "metal": job["metal"],
            "variant": job["variant"], "replica": int(job["replica"]), "T_K": 300,
            "dir": str(base), "prefix": "traj-denser-dumping",
        })
    return jobs


def all_jobs() -> list[dict]:
    return solution_jobs() + graphite_jobs()


def cache_path(job: dict) -> Path:
    run = f"T{job['T_K']}" if job["system"] == "solution" else f"r{job['replica']:02d}"
    return CACHE / job["system"] / job["model"] / job["variant"] / job["metal"] / f"{run}.npz"


def extract(job: dict) -> str:
    target = cache_path(job)
    if target.is_file():
        return f"skip {target}"
    base = Path(job["dir"])
    meta = json.loads((base / f"{job['prefix']}-meta.json").read_text())
    positions = np.load(base / f"{job['prefix']}-positions.npy", mmap_mode="r")
    energies = np.load(base / f"{job['prefix']}-energies.npy", mmap_mode="r")
    symbols = np.asarray(meta["symbols"])
    oxygen = np.flatnonzero(symbols == "O")
    hydrogen = np.flatnonzero(symbols == "H")
    metal = int(np.flatnonzero(symbols == job["metal"])[0])
    if not all(meta["pbc"]):
        raise ValueError(f"{base}: expected 3D PBC")
    cell = np.asarray(meta["cell_A"], dtype=np.float64)
    if np.count_nonzero(cell - np.diag(np.diag(cell))) != 0:
        raise ValueError(f"{base}: expected an orthorhombic cell")
    box = np.diag(cell).astype(np.float64)
    dump_fs = float(meta["dump_interval_fs"])
    stride = int(round(SAMPLE_FS / dump_fs))
    n_written = int(meta["frames_written"])
    frames = np.arange(0, n_written, stride, dtype=np.int64)
    n_t = frames.size

    near_d = np.empty((n_t, N_NEAR_O), dtype=np.float32)
    near_type = np.empty((n_t, N_NEAR_O), dtype=np.int8)
    oh_d = np.full((n_t, N_OH_KEEP), np.inf, dtype=np.float32)
    n_oh = np.empty(n_t, dtype=np.int8)
    n_h3o = np.empty(n_t, dtype=np.int8)
    n_o_bare = np.empty(n_t, dtype=np.int8)
    min_ln_h = np.empty(n_t, dtype=np.float32)
    epot = np.asarray(energies[frames, 0], dtype=np.float64)

    for start in range(0, n_t, CHUNK):
        sel = frames[start : start + CHUNK]
        pos = np.asarray(positions[sel], dtype=np.float64)
        for k in range(pos.shape[0]):
            frame = pos[k]
            o_xyz = np.mod(frame[oxygen], box)
            h_xyz = np.mod(frame[hydrogen], box)
            o_xyz = np.where(o_xyz >= box, 0.0, o_xyz)
            h_xyz = np.where(h_xyz >= box, 0.0, h_xyz)
            tree = cKDTree(o_xyz, boxsize=box)
            _, host = tree.query(h_xyz, k=1)
            n_h = np.bincount(host, minlength=oxygen.size)
            o_type = np.where(n_h == 2, 0, np.where(n_h == 1, 1, 2)).astype(np.int8)

            delta = frame[oxygen] - frame[metal]
            delta -= box * np.rint(delta / box)
            d_o = np.sqrt(np.sum(delta * delta, axis=1))
            order = np.argsort(d_o)[:N_NEAR_O]
            row = start + k
            near_d[row] = d_o[order]
            near_type[row] = o_type[order]
            d_oh = np.sort(d_o[n_h == 1])[:N_OH_KEEP]
            oh_d[row, : d_oh.size] = d_oh
            n_oh[row] = int(np.count_nonzero(n_h == 1))
            n_h3o[row] = int(np.count_nonzero(n_h == 3))
            n_o_bare[row] = int(np.count_nonzero(n_h == 0))

            delta_h = frame[hydrogen] - frame[metal]
            delta_h -= box * np.rint(delta_h / box)
            min_ln_h[row] = float(np.sqrt(np.min(np.sum(delta_h * delta_h, axis=1))))

    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(".tmp.npz")
    np.savez_compressed(
        tmp,
        t_ps=(frames * dump_fs / 1000.0).astype(np.float64),
        near_d=near_d, near_type=near_type, oh_d=oh_d, n_oh=n_oh, n_h3o=n_h3o,
        n_o_bare=n_o_bare, min_ln_h=min_ln_h, epot_eV=epot,
        box_A=box, n_atoms=np.int32(symbols.size), n_oxygen=np.int32(oxygen.size),
        n_hydrogen=np.int32(hydrogen.size), n_carbon=np.int32(np.count_nonzero(symbols == "C")),
        dump_interval_fs=np.float64(dump_fs), sample_fs=np.float64(SAMPLE_FS),
    )
    os.replace(tmp, target)
    return f"done {target}"


BLOCK_SAMPLES = 50
DROP_THRESHOLD_EV = {"graphite": -25.0, "solution": -15.0}
REFINE_THRESHOLD_EV = -10.0
GEOM_LN_O_A = 1.85
GEOM_LN_H_A = 1.90
GEOM_SUSTAIN = 10
ONSET_MARGIN_PS = 1.0
SCREEN_CSV = OUT / "collapse_screen.csv"


def first_sustained(mask: np.ndarray, length: int) -> int:
    padded = np.concatenate([[0], mask.astype(np.int8), [0]])
    edges = np.diff(padded)
    starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
    long_runs = np.flatnonzero(ends - starts >= length)
    return int(starts[long_runs[0]]) if long_runs.size else -1


def screen_one(job: dict) -> dict:
    """Collapse onset from abrupt potential-energy drops, cross-checked with short Ln-O/Ln-H contacts."""
    data = np.load(cache_path(job))
    t, epot = data["t_ps"], data["epot_eV"]
    ln_o, ln_h = data["near_d"][:, 0], data["min_ln_h"]
    n_blocks = t.size // BLOCK_SAMPLES
    blocks = epot[: n_blocks * BLOCK_SAMPLES].reshape(n_blocks, BLOCK_SAMPLES).mean(axis=1)
    drops = np.array([blocks[b] - blocks[max(0, b - 5) : b].min() for b in range(1, n_blocks)])
    hit = np.flatnonzero(drops < DROP_THRESHOLD_EV[job["system"]])
    energy_onset = math.inf
    if hit.size:
        block = int(hit[0]) + 1
        reference = blocks[max(0, block - 5) : block].min()
        lo = (block - 1) * BLOCK_SAMPLES
        below = np.flatnonzero(epot[lo:] < reference + REFINE_THRESHOLD_EV)
        energy_onset = float(t[lo + below[0]]) if below.size else float(t[block * BLOCK_SAMPLES])
    i_o = first_sustained(ln_o < GEOM_LN_O_A, GEOM_SUSTAIN)
    i_h = first_sustained(ln_h < GEOM_LN_H_A, GEOM_SUSTAIN)
    geom_onset = min(float(t[i_o]) if i_o >= 0 else math.inf, float(t[i_h]) if i_h >= 0 else math.inf)
    onset = min(energy_onset, geom_onset)
    t_end = float(t[-1])
    usable_end = t_end if not math.isfinite(onset) else max(0.0, onset - ONSET_MARGIN_PS)
    keep = t <= usable_end
    status = "ok" if not math.isfinite(onset) else ("truncated" if usable_end >= 20.0 else "excluded")
    return {
        "system": job["system"], "model": job["model"], "variant": job["variant"], "metal": job["metal"],
        "replica": job["replica"], "T_K": job["T_K"], "t_end_ps": round(t_end, 3),
        "energy_onset_ps": "" if not math.isfinite(energy_onset) else round(energy_onset, 3),
        "geom_onset_ps": "" if not math.isfinite(geom_onset) else round(geom_onset, 3),
        "usable_end_ps": round(usable_end, 3),
        "max_block_drop_eV": round(float(drops.min()), 2),
        "min_LnO_usable_A": round(float(ln_o[keep].min()), 3) if keep.any() else "",
        "min_LnH_usable_A": round(float(ln_h[keep].min()), 3) if keep.any() else "",
        "status": status,
    }


def run_screen() -> None:
    rows = [screen_one(job) for job in all_jobs()]
    with SCREEN_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = defaultdict(lambda: defaultdict(int))
    for row in rows:
        summary[(row["system"], row["model"], row["metal"])][row["status"]] += 1
    for key, counts in sorted(summary.items()):
        print(key, dict(counts))


def load_screen() -> dict[tuple, dict]:
    with SCREEN_CSV.open() as handle:
        return {
            (r["system"], r["model"], r["variant"], r["metal"], int(r["replica"]), int(r["T_K"])): r
            for r in csv.DictReader(handle)
        }


def run_extract(workers: int) -> None:
    jobs = all_jobs()
    todo = [job for job in jobs if not cache_path(job).is_file()]
    print(f"extract: {len(jobs)} trajectories, {len(todo)} to process with {workers} workers", flush=True)
    if workers <= 1:
        for i, job in enumerate(todo, 1):
            extract(job)
            if i % 20 == 0:
                print(f"  {i}/{len(todo)}", flush=True)
        return
    with ProcessPoolExecutor(workers) as pool:
        for i, _ in enumerate(pool.map(extract, todo, chunksize=1), 1):
            if i % 20 == 0:
                print(f"  {i}/{len(todo)}", flush=True)


# --------------------------------------------------------------------------- analyze
#
# State n(t) = number of hydroxide O (exactly one nearest H) within r_c of Ln, from the 20 fs series.
# Raw n(t) flickers on 20-100 fs because protons shuttle between Ln-bound OH/H2O and neighbouring
# waters, and implied timescales of a lag-time Markov model keep growing with lag (no Markovian
# plateau).  Committed transitions are therefore defined at an explicit time resolution: a centred
# sliding-window majority (mode) filter of width W followed by absorption of filtered runs shorter
# than W/2 (cf. the t* = 0.5 ps convention for ligand-exchange counting).  Occupancies (thermodynamics)
# come from the unfiltered frames and do not depend on W.

KIN_OUT = OUT / "kinetics"
R_KJ = 8.314462618e-3
KB_OVER_H = 2.083661912e10
N_AVOGADRO = 6.02214076e23
PKW = 14.0
LN10 = math.log(10.0)
V_SOLUTION_A3 = 15.781866235591567**3
N_O_SOLUTION = 131
DT_PS = SAMPLE_FS / 1000.0
N_STATE = 4
MIN_WINDOW_PS = 5.0
BLOCK_PS = 10.0
N_BOOT = 2000
THERMO_START_PS = {"solution": 10.0, "graphite": 50.0}
MECH_HALF_PS = 0.5
VH_MIN_DWELL_PS = 1.0
SCENARIOS = (
    # name, filter window W (ps), r_c override (A) or None for the group value, analysis start (ps)
    ("primary", 0.5, None, 10.0),
    ("W0.2", 0.2, None, 10.0),
    ("W1.0", 1.0, None, 10.0),
    ("W2.0", 2.0, None, 10.0),
    ("rc2.90", 0.5, 2.90, 10.0),
    ("rc3.30", 0.5, 3.30, 10.0),
    ("start1", 0.5, None, 1.0),
    ("start50", 0.5, None, 50.0),
    ("start90", 0.5, None, 90.0),
)
MODEL_LABEL = {
    "DPA4C-Plus": "DPA4C-Plus", "uma-s-1.2.1": "UMA-s-1.2.1",
    "DPA4-Pro": "DPA4-Pro", "MACE-POLAR-1-M": "MACE-POLAR-1-M",
}
STEP_LABEL = {3: "3→2", 2: "2→1", 1: "1→0"}


def group_key(job: dict) -> tuple[str, str, str]:
    return (job["system"], job["model"], job["metal"])


def odd_samples(width_ps: float) -> int:
    return 2 * int(round(width_ps / DT_PS / 2.0)) + 1


def mode_filter(s: np.ndarray, width: int) -> np.ndarray:
    """Centred sliding-window majority; ties keep the raw centre state when it is among the maxima."""
    onehot = np.zeros((s.size, N_STATE + 1), dtype=np.int32)
    onehot[np.arange(s.size), np.clip(s, 0, N_STATE)] = 1
    csum = np.vstack([np.zeros((1, N_STATE + 1), dtype=np.int32), np.cumsum(onehot, axis=0)])
    idx = np.arange(s.size)
    half = width // 2
    counts = csum[np.clip(idx + half + 1, 0, s.size)] - csum[np.clip(idx - half, 0, s.size)]
    best = counts.argmax(axis=1)
    keep_raw = counts[idx, s] == counts[idx, best]
    return np.where(keep_raw, s, best)


def runs_of(s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    change = np.flatnonzero(np.diff(s) != 0) + 1
    return np.r_[0, change], np.r_[change, s.size]


def absorb_short(s: np.ndarray, min_len: int) -> np.ndarray:
    """Runs shorter than min_len samples take the last committed state (leading short runs take the first)."""
    starts, ends = runs_of(s)
    long_runs = np.flatnonzero(ends - starts >= min_len)
    if long_runs.size == 0:
        return np.full_like(s, np.bincount(s, minlength=N_STATE).argmax())
    out = s.copy()
    current = s[starts[long_runs[0]]]
    for a, b in zip(starts, ends):
        if b - a >= min_len:
            current = s[a]
        out[a:b] = current
    return out


def committed_states(raw: np.ndarray, width_ps: float) -> np.ndarray:
    width = odd_samples(width_ps)
    return absorb_short(mode_filter(raw, width), max(1, width // 2))


def oh_density(files_windows: list[tuple[Path, float, float]]) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(1.5, 6.0 + 1e-9, 0.02)
    hist = np.zeros(edges.size - 1)
    for path, start, end in files_windows:
        with np.load(path) as data:
            t = data["t_ps"]
            sel = (t >= start) & (t <= end)
            d = data["oh_d"][sel].ravel()
        hist += np.histogram(d[np.isfinite(d)], edges)[0]
    r = 0.5 * (edges[1:] + edges[:-1])
    return r, hist / r**2


def choose_cutoff(r: np.ndarray, density: np.ndarray) -> dict:
    smooth = np.convolve(density, np.ones(5) / 5.0, mode="same")
    first = np.flatnonzero((r >= 1.9) & (r <= 2.8))
    peak = first[smooth[first].argmax()]
    gap = np.flatnonzero((r >= 2.7) & (r <= 3.6))
    low = smooth[gap].min()
    at_min = gap[smooth[gap] <= low + 1e-12 * max(smooth[peak], 1.0)]
    # Centre of the longest run of minimal bins (the gap is often exactly empty).
    splits = np.split(at_min, np.flatnonzero(np.diff(at_min) > 1) + 1)
    longest = max(splits, key=len)
    rc = round(float(r[longest].mean()), 2)
    second = np.flatnonzero((r >= 3.6) & (r <= 5.5))
    return {
        "rc_A": rc,
        "first_peak_A": float(r[peak]),
        "second_shell_peak_A": float(r[second[smooth[second].argmax()]]) if smooth[second].max() > 0 else math.nan,
        "density_min_over_peak": float(low / smooth[peak]) if smooth[peak] > 0 else math.nan,
    }


def per_trajectory(args: tuple[dict, dict, dict]) -> dict:
    job, screen_row, rc_by_group = args
    usable_end = float(screen_row["usable_end_ps"])
    with np.load(cache_path(job)) as data:
        t = data["t_ps"]
        oh_d = data["oh_d"]
        near_d = data["near_d"]
        near_type = data["near_type"]
    out = {"job": {k: job[k] for k in ("system", "model", "metal", "variant", "replica", "T_K")},
           "status": screen_row["status"], "usable_end_ps": usable_end, "scenarios": {}}
    thermo_start = THERMO_START_PS[job["system"]]
    for name, width_ps, rc_override, start_ps in SCENARIOS:
        rc = rc_override if rc_override is not None else rc_by_group[group_key(job)]
        sel = np.flatnonzero((t >= start_ps - 1e-9) & (t <= usable_end + 1e-9))
        if sel.size * DT_PS < MIN_WINDOW_PS:
            out["scenarios"][name] = None
            continue
        tt = t[sel]
        raw_all = (oh_d[sel] < rc).sum(axis=1)
        raw = np.minimum(raw_all, N_STATE - 1)
        cn_total = (near_d[sel] < rc).sum(axis=1)
        committed = committed_states(raw, width_ps)
        block = np.minimum(((tt - tt[0]) // BLOCK_PS).astype(int), max(0, int((tt[-1] - tt[0]) // BLOCK_PS) - 1)) \
            if (tt[-1] - tt[0]) >= BLOCK_PS else np.zeros(tt.size, dtype=int)
        n_blocks = int(block.max()) + 1
        n_mat = np.zeros((n_blocks, N_STATE, N_STATE))
        mech = np.zeros((N_STATE, N_STATE, 3))
        change = np.flatnonzero(np.diff(committed) != 0)
        half = max(1, int(round(MECH_HALF_PS / DT_PS)))
        for k in change:
            a, b = int(committed[k]), int(committed[k + 1])
            n_mat[block[k + 1], a, b] += 1
            before = cn_total[max(0, k + 1 - half) : k + 1].mean()
            after = cn_total[k + 1 : k + 1 + half].mean()
            d_cn = int(round(after - before))
            expected = 1 if b > a else -1
            mech[a, b, 0 if d_cn == 0 else (1 if d_cn == expected * abs(b - a) else 2)] += 1
        dwell = np.zeros((n_blocks, N_STATE))
        occ = np.zeros((n_blocks, N_STATE))
        np.add.at(dwell, (block, committed), DT_PS)
        np.add.at(occ, (block, raw), 1)
        thermo = tt >= thermo_start - 1e-9
        occ_thermo = np.bincount(raw[thermo], minlength=N_STATE)[:N_STATE].astype(float)
        dwell_thermo = np.bincount(committed[thermo], minlength=N_STATE)[:N_STATE] * DT_PS
        out["scenarios"][name] = {
            "rc_A": rc, "t_start_ps": float(tt[0]), "t_end_ps": float(tt[-1]),
            "N_blocks": n_mat, "T_blocks": dwell, "occ_blocks": occ, "occ_thermo": occ_thermo,
            "T_thermo": dwell_thermo,
            "mech": mech,
            "n_frames": int(sel.size), "frac_raw_gt3": float(np.mean(raw_all > N_STATE - 1)),
            "frac_oxo_in_shell": float(np.mean(((near_type[sel] == 2) & (near_d[sel] < rc)).any(axis=1))),
            "n_raw_changes": int(np.count_nonzero(np.diff(raw))),
        }
    return out


def garwood(n: float, conf: float = 0.95) -> tuple[float, float]:
    from scipy.stats import chi2
    alpha = 1.0 - conf
    lo = 0.0 if n == 0 else chi2.ppf(alpha / 2, 2 * n) / 2
    hi = chi2.ppf(1 - alpha / 2, 2 * n + 2) / 2
    return float(lo), float(hi)


def eyring_dg(k_ps: float, temperature: float) -> float:
    if not k_ps > 0:
        return math.inf
    return R_KJ * temperature * math.log(KB_OVER_H * temperature / (k_ps * 1e12))


def step_rates(n_mat: np.ndarray, dwell: np.ndarray, step: int) -> dict:
    """n_mat/dwell may carry leading bootstrap axes; step n: n -> lower (dissociation), n-1 -> higher (association)."""
    lower = n_mat[..., step, :step].sum(axis=-1)
    higher = n_mat[..., step - 1, step:].sum(axis=-1)
    t_hi = dwell[..., step]
    t_lo = dwell[..., step - 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        k_d = np.where(t_hi > 0, lower / np.where(t_hi > 0, t_hi, 1.0), np.nan)
        k_a = np.where(t_lo > 0, higher / np.where(t_lo > 0, t_lo, 1.0), np.nan)
    return {"N_d": lower, "N_a": higher, "T_hi": t_hi, "T_lo": t_lo, "k_d": k_d, "k_a": k_a}


def concentration_single(system: str, n_oxygen: int) -> float:
    """Concentration (M) of one solute particle in the liquid volume (graphite: solution molar volume per O)."""
    volume_a3 = V_SOLUTION_A3 if system == "solution" else V_SOLUTION_A3 * n_oxygen / N_O_SOLUTION
    return 1.0 / (N_AVOGADRO * volume_a3 * 1e-27)


def thermo_from_occ(occ: np.ndarray, step: int, temperature: float, c1: float) -> dict:
    p = occ / occ.sum(axis=-1, keepdims=True)
    hi, lo = p[..., step], p[..., step - 1]
    rt = R_KJ * temperature
    with np.errstate(divide="ignore", invalid="ignore"):
        dg = -rt * np.log(lo / hi)
    dg0 = dg - rt * math.log((N_STATE - step) * c1)
    pkd = dg0 / (rt * LN10)
    return {"P_n": hi, "P_n_minus_1": lo, "dG": dg, "dG0_d": dg0, "pK_d": pkd, "pKa_eq": PKW - pkd}


def percentile_ci(samples: np.ndarray) -> tuple[float, float]:
    s = np.asarray(samples, dtype=float)
    s = s[~np.isnan(s)]
    if s.size == 0:
        return (math.nan, math.nan)
    s = np.sort(s)
    lo = s[int(math.floor(0.025 * (s.size - 1)))]
    hi = s[int(math.ceil(0.975 * (s.size - 1)))]
    return (float(lo), float(hi))


def wls_line(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    sw = w.sum()
    xm, ym = (w * x).sum() / sw, (w * y).sum() / sw
    sxx = (w * (x - xm) ** 2).sum()
    slope = (w * (x - xm) * (y - ym)).sum() / sxx
    return float(slope), float(ym - slope * xm)


def fmt(x: float, nd: int = 3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    if isinstance(x, float) and math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return round(float(x), nd)


def bootstrap_units(units: list[dict], strata: list[str], rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Resample units (blocks or trajectories) with replacement inside each stratum; returns summed arrays."""
    keys = ("N", "T", "occ", "occ_thermo")
    stacked = {k: np.stack([u[k] for u in units]) for k in keys}
    totals = {k: np.zeros((N_BOOT,) + stacked[k].shape[1:]) for k in keys}
    strata_arr = np.asarray(strata)
    for stratum in np.unique(strata_arr):
        members = np.flatnonzero(strata_arr == stratum)
        draws = members[rng.integers(0, members.size, size=(N_BOOT, members.size))]
        for k in keys:
            totals[k] += stacked[k][draws].sum(axis=1)
    return totals


def write_rows(path: Path, rows: list[dict]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def unit_list(results: list[dict], scenario: str, per_block: bool) -> tuple[list[dict], list[str]]:
    units, strata = [], []
    for res in results:
        sc = res["scenarios"].get(scenario)
        if sc is None:
            continue
        if per_block:
            for b in range(sc["N_blocks"].shape[0]):
                units.append({"N": sc["N_blocks"][b], "T": sc["T_blocks"][b],
                              "occ": sc["occ_blocks"][b], "occ_thermo": sc["occ_blocks"][b],
                              "T_thermo": sc["T_blocks"][b],
                              "mech": sc["mech"] if b == 0 else np.zeros_like(sc["mech"])})
                strata.append("all")
        else:
            units.append({"N": sc["N_blocks"].sum(0), "T": sc["T_blocks"].sum(0),
                          "occ": sc["occ_blocks"].sum(0), "occ_thermo": sc["occ_thermo"],
                          "T_thermo": sc["T_thermo"], "mech": sc["mech"]})
            strata.append(res["job"]["variant"])
    return units, strata


def summarize_group(meta: dict, units: list[dict], strata: list[str], temperature: float,
                    c1: float, ci_method: str, rng: np.random.Generator) -> tuple[list[dict], list[dict], dict]:
    total = {k: np.sum([u[k] for u in units], axis=0) for k in ("N", "T", "occ", "occ_thermo", "T_thermo", "mech")}
    boot = bootstrap_units(units, strata, rng)
    rate_rows, thermo_rows, boot_store = [], [], {}
    rt = R_KJ * temperature
    for step in (3, 2, 1):
        point = step_rates(total["N"], total["T"], step)
        bs = step_rates(boot["N"], boot["T"], step)
        boot_store[step] = bs
        mech_d = total["mech"][step, :step].sum(axis=0)
        mech_a = total["mech"][step - 1, step:].sum(axis=0)
        k_kin = {}
        for direction, n_key, k_key, t_key, mech, multi, per in (
            ("dissociation", "N_d", "k_d", "T_hi", mech_d, total["N"][step, : step - 1].sum(), step),
            ("association", "N_a", "k_a", "T_lo", mech_a, total["N"][step - 1, step + 1 :].sum(), N_STATE - step),
        ):
            n_ev = float(point[n_key])
            t_origin = float(point[t_key])
            k = float(point[k_key]) if t_origin > 0 else math.nan
            k_kin[direction] = k
            b_lo, b_hi = percentile_ci(bs[k_key])
            g_lo, g_hi = garwood(n_ev)
            g_lo, g_hi = (g_lo / t_origin, g_hi / t_origin) if t_origin > 0 else (math.nan, math.nan)
            if t_origin > 0:
                lo, hi = float(np.nanmin([b_lo, g_lo])), float(np.nanmax([b_hi, g_hi]))
            else:
                lo = hi = math.nan
            dg = eyring_dg(k, temperature) if t_origin > 0 else math.nan
            n_mech = mech.sum()
            rate_rows.append({
                **meta, "T_K": temperature, "step": STEP_LABEL[step], "hydrolysis_step": step,
                "direction": direction,
                "process": (f"Ln(OH){step} -> Ln(OH){step - 1} + OH-" if direction == "dissociation"
                            else f"Ln(OH){step - 1} + OH- -> Ln(OH){step}"),
                "N_events": int(n_ev), "N_multistep": int(multi), "T_origin_state_ps": fmt(t_origin, 2),
                "k_ps-1": fmt(k, 5), "k_CI95_lo": fmt(lo, 5), "k_CI95_hi": fmt(hi, 5),
                "k_boot_lo": fmt(b_lo, 5), "k_boot_hi": fmt(b_hi, 5),
                "k_garwood_lo": fmt(g_lo, 5), "k_garwood_hi": fmt(g_hi, 5),
                "dG_act_kJmol": fmt(dg, 2),
                "dG_act_CI95_lo": fmt(eyring_dg(hi, temperature), 2) if t_origin > 0 else "",
                "dG_act_CI95_hi": fmt(eyring_dg(lo, temperature), 2) if t_origin > 0 else "",
                "dG_act_lower_bound_if_no_events": fmt(eyring_dg(g_hi, temperature), 2) if n_ev == 0 and t_origin > 0 else "",
                "dG_act_per_OH_kJmol": fmt(dg + rt * math.log(per), 2) if math.isfinite(dg) else "",
                "frac_proton_transfer_like": fmt(mech[0] / n_mech, 3) if n_mech else "",
                "frac_ligand_exchange_like": fmt(mech[1] / n_mech, 3) if n_mech else "",
                "ci_method": ci_method,
            })
        occ = total["occ_thermo"]
        th = thermo_from_occ(occ, step, temperature, c1)
        th_b = thermo_from_occ(boot["occ_thermo"], step, temperature, c1)
        cis = {key: percentile_ci(np.where(np.isfinite(th_b[key]), th_b[key], np.nan)) for key in ("dG", "dG0_d", "pKa_eq")}
        undefined = float(np.mean(~np.isfinite(th_b["dG"])))
        k_d, k_a = k_kin["dissociation"], k_kin["association"]
        dg_kin = -rt * math.log(k_d / k_a) if (k_d and k_a and k_d > 0 and k_a > 0) else math.nan
        thermo_rows.append({
            **meta, "T_K": temperature, "step": STEP_LABEL[step], "hydrolysis_step": step,
            "P_n": fmt(float(th["P_n"]), 4), "P_n_minus_1": fmt(float(th["P_n_minus_1"]), 4),
            "frames": int(occ.sum()),
            "T_committed_n_ps": fmt(float(total["T_thermo"][step]), 2),
            "T_committed_n_minus_1_ps": fmt(float(total["T_thermo"][step - 1]), 2),
            "both_states_committed": bool(min(total["T_thermo"][step], total["T_thermo"][step - 1]) >= VH_MIN_DWELL_PS),
            "dG_occ_kJmol": fmt(float(th["dG"]), 2), "dG_occ_CI95_lo": fmt(cis["dG"][0], 2),
            "dG_occ_CI95_hi": fmt(cis["dG"][1], 2), "bootstrap_undefined_frac": fmt(undefined, 3),
            "dG_kinetic_kJmol": fmt(dg_kin, 2),
            "c_single_M": fmt(c1, 4), "free_OH_conc_M": fmt((N_STATE - step) * c1, 4),
            "dG0_dissociation_kJmol": fmt(float(th["dG0_d"]), 2),
            "dG0_dissociation_CI95_lo": fmt(cis["dG0_d"][0], 2), "dG0_dissociation_CI95_hi": fmt(cis["dG0_d"][1], 2),
            "pK_dissociation": fmt(float(th["pK_d"]), 2),
            "pKa_equivalent": fmt(float(th["pKa_eq"]), 2),
            "pKa_equivalent_CI95_lo": fmt(cis["pKa_eq"][0], 2), "pKa_equivalent_CI95_hi": fmt(cis["pKa_eq"][1], 2),
            "dG_acid_equivalent_kJmol": fmt(float(th["pKa_eq"]) * rt * LN10, 2),
        })
    p_all = total["occ_thermo"] / max(total["occ_thermo"].sum(), 1)
    p_boot = boot["occ_thermo"] / np.maximum(boot["occ_thermo"].sum(axis=1, keepdims=True), 1)
    occ_row = {**meta, "T_K": temperature, "frames": int(total["occ_thermo"].sum())}
    for n in range(N_STATE):
        lo, hi = percentile_ci(p_boot[:, n])
        occ_row.update({f"P{n}": fmt(float(p_all[n]), 4), f"P{n}_CI95_lo": fmt(lo, 4), f"P{n}_CI95_hi": fmt(hi, 4)})
    occ_row["mean_n_OH_bound"] = fmt(float((p_all * np.arange(N_STATE)).sum()), 3)
    window_p = total["occ"] / max(total["occ"].sum(), 1)
    occ_row.update({f"P{n}_from_analysis_start": fmt(float(window_p[n]), 4) for n in range(N_STATE)})
    filtered_p = total["T"] / max(total["T"].sum(), 1e-12)
    occ_row.update({f"P{n}_committed": fmt(float(filtered_p[n]), 4) for n in range(N_STATE)})
    return rate_rows, thermo_rows, {"boot": boot_store, "boot_occ": boot["occ_thermo"], "total": total, "occ_row": occ_row}


def arrhenius(points: list[dict]) -> dict:
    """points: per-T dicts with T, k, N, k_boot, N_boot.  Fits of ln k and ln(k/T) vs 1/T weighted by event
    counts (Poisson: var ln k ~ 1/N).  CI = envelope of the block/trajectory bootstrap and the Poisson-parametric
    95% interval of the weighted fit."""
    use = [p for p in points if p["N"] > 0 and p["k"] > 0]
    n_total = int(sum(p["N"] for p in use))
    out = {"T_used": ";".join(str(int(p["T"])) for p in use), "n_T": len(use), "N_total": n_total,
           "N_per_T": ";".join(str(int(p["N"])) for p in use)}
    if len(use) < 2:
        out.update({"quality": "none", "note": "fewer than 2 temperatures with events"})
        return out
    x = np.array([1.0 / p["T"] for p in use])
    w = np.array([p["N"] for p in use], dtype=float)
    y = np.log([p["k"] for p in use])
    y_e = y - np.log(1.0 / x)
    slope, icpt = wls_line(x, y, w)
    slope_e, icpt_e = wls_line(x, y_e, w)
    xm = (w * x).sum() / w.sum()
    sxx = (w * (x - xm) ** 2).sum()
    se_slope = math.sqrt(1.0 / sxx)
    se_icpt = math.sqrt(1.0 / w.sum() + xm**2 / sxx)
    ln_kbh = math.log(KB_OVER_H * 1e-12)
    ea_b, dh_b, ds_b = [], [], []
    n_boot = len(use[0]["k_boot"])
    for b in range(n_boot):
        kb = np.array([p["k_boot"][b] for p in use])
        nb = np.array([p["N_boot"][b] for p in use], dtype=float)
        ok = (nb > 0) & (kb > 0) & np.isfinite(kb)
        if ok.sum() < 2:
            continue
        sb, _ = wls_line(x[ok], np.log(kb[ok]), nb[ok])
        se, ie = wls_line(x[ok], np.log(kb[ok]) - np.log(1.0 / x[ok]), nb[ok])
        ea_b.append(-sb * R_KJ)
        dh_b.append(-se * R_KJ)
        ds_b.append((ie - ln_kbh) * R_KJ * 1000.0)
    ea, dh, ds = -slope * R_KJ, -slope_e * R_KJ, (icpt_e - ln_kbh) * R_KJ * 1000.0
    z = 1.959963984540054
    param = {
        "Ea": (ea - z * se_slope * R_KJ, ea + z * se_slope * R_KJ),
        "dH": (dh - z * se_slope * R_KJ, dh + z * se_slope * R_KJ),
        "dS": (ds - z * se_icpt * R_KJ * 1000.0, ds + z * se_icpt * R_KJ * 1000.0),
    }
    boot = {"Ea": percentile_ci(ea_b), "dH": percentile_ci(dh_b), "dS": percentile_ci(ds_b)}
    env = {key: (float(np.nanmin([boot[key][0], param[key][0]])), float(np.nanmax([boot[key][1], param[key][1]])))
           for key in param}
    min_n = min(p["N"] for p in use)
    quality = "good" if (len(use) >= 4 and n_total >= 20 and min_n >= 3) else (
        "fair" if (len(use) >= 3 and n_total >= 8) else "poor")
    out.update({
        "quality": quality,
        "Ea_kJmol": fmt(ea, 2), "Ea_CI95_lo": fmt(env["Ea"][0], 2), "Ea_CI95_hi": fmt(env["Ea"][1], 2),
        "Ea_boot_lo": fmt(boot["Ea"][0], 2), "Ea_boot_hi": fmt(boot["Ea"][1], 2),
        "Ea_poisson_lo": fmt(param["Ea"][0], 2), "Ea_poisson_hi": fmt(param["Ea"][1], 2),
        "lnA_ps-1": fmt(icpt, 3), "A_ps-1": fmt(math.exp(icpt), 4),
        "dH_act_kJmol": fmt(dh, 2), "dH_act_CI95_lo": fmt(env["dH"][0], 2), "dH_act_CI95_hi": fmt(env["dH"][1], 2),
        "dS_act_JmolK": fmt(ds, 1), "dS_act_CI95_lo": fmt(env["dS"][0], 1), "dS_act_CI95_hi": fmt(env["dS"][1], 1),
        "dG_act_300K_from_fit_kJmol": fmt(dh - 300.0 * ds / 1000.0, 2),
        "bootstrap_valid_frac": fmt(len(ea_b) / n_boot, 3),
        "note": "2-point slope" if len(use) == 2 else "",
    })
    return out


def vant_hoff(points: list[dict]) -> dict:
    use = [p for p in points if np.isfinite(p["lnK"])]
    out = {"T_used": ";".join(str(int(p["T"])) for p in use), "n_T": len(use)}
    if len(use) < 2:
        out["note"] = "fewer than 2 temperatures with both states committed (dwell >= 1 ps)"
        return out
    x = np.array([1.0 / p["T"] for p in use])
    y = np.array([p["lnK"] for p in use])
    slope, icpt = wls_line(x, y, np.ones_like(x))
    dh_b, ds_b = [], []
    for b in range(len(use[0]["lnK_boot"])):
        yb = np.array([p["lnK_boot"][b] for p in use])
        ok = np.isfinite(yb)
        if ok.sum() < 2:
            continue
        sb, ib = wls_line(x[ok], yb[ok], np.ones(ok.sum()))
        dh_b.append(-sb * R_KJ)
        ds_b.append(ib * R_KJ * 1000.0)
    dh_ci, ds_ci = percentile_ci(dh_b), percentile_ci(ds_b)
    dh, ds = -slope * R_KJ, icpt * R_KJ * 1000.0
    out.update({
        "dH0_dissociation_kJmol": fmt(dh, 2), "dH0_CI95_lo": fmt(dh_ci[0], 2), "dH0_CI95_hi": fmt(dh_ci[1], 2),
        "dS0_dissociation_JmolK": fmt(ds, 1), "dS0_CI95_lo": fmt(ds_ci[0], 1), "dS0_CI95_hi": fmt(ds_ci[1], 1),
        "dG0_dissociation_300K_fit_kJmol": fmt(dh - 300.0 * ds / 1000.0, 2),
        "note": "" if len(use) > 2 else "2-point line",
    })
    return out


def run_analyze(workers: int) -> None:
    KIN_OUT.mkdir(parents=True, exist_ok=True)
    screen = load_screen()
    selected, skipped = [], []
    for job in all_jobs():
        row = screen[(job["system"], job["model"], job["variant"], job["metal"], int(job["replica"]), int(job["T_K"]))]
        (skipped if row["status"] == "excluded" else selected).append((job, row))
    print(f"analyze: {len(selected)} trajectories ({len(skipped)} excluded by collapse screen)", flush=True)

    windows = defaultdict(list)
    for job, row in selected:
        windows[group_key(job)].append((cache_path(job), 10.0, float(row["usable_end_ps"])))
    rc_by_group, cutoff_rows = {}, []
    for key, items in sorted(windows.items()):
        r, density = oh_density(items)
        info = choose_cutoff(r, density)
        rc_by_group[key] = info["rc_A"]
        cutoff_rows.append({"system": key[0], "model": key[1], "metal": key[2], "n_trajectories": len(items),
                            **{k: fmt(v, 4) for k, v in info.items()}})
        print("  cutoff", key, info, flush=True)
    write_rows(KIN_OUT / "kinetics_cutoffs.csv", cutoff_rows)

    tasks = [(job, row, rc_by_group) for job, row in selected]
    with ProcessPoolExecutor(workers) as pool:
        results = list(pool.map(per_trajectory, tasks, chunksize=4))
    print("  per-trajectory pass done", flush=True)

    traj_rows = []
    for res in results:
        sc = res["scenarios"]["primary"]
        row = {**res["job"], "status": res["status"], "usable_end_ps": res["usable_end_ps"]}
        if sc is None:
            row["note"] = "window shorter than 5 ps after 10 ps"
            traj_rows.append(row)
            continue
        n_mat, dwell, occ = sc["N_blocks"].sum(0), sc["T_blocks"].sum(0), sc["occ_blocks"].sum(0)
        row.update({"rc_A": sc["rc_A"], "t_start_ps": sc["t_start_ps"], "t_end_ps": sc["t_end_ps"],
                    "n_raw_state_changes": sc["n_raw_changes"], "n_committed_transitions": int(n_mat.sum())})
        for a in range(N_STATE):
            for b in range(N_STATE):
                if a != b:
                    row[f"N_{a}to{b}"] = int(n_mat[a, b])
        row.update({f"T{n}_ps": round(float(dwell[n]), 2) for n in range(N_STATE)})
        row.update({f"P{n}_raw": round(float(occ[n] / occ.sum()), 4) for n in range(N_STATE)})
        row.update({"frac_frames_more_than_3_OH": round(sc["frac_raw_gt3"], 5),
                    "frac_frames_oxo_in_shell": round(sc["frac_oxo_in_shell"], 5)})
        traj_rows.append(row)
    write_rows(KIN_OUT / "kinetics_events_by_traj.csv", traj_rows)

    rng = np.random.default_rng(20240607)
    n_oxygen = {}
    for job, _ in selected:
        if job["system"] not in n_oxygen:
            with np.load(cache_path(job)) as data:
                n_oxygen[job["system"]] = int(data["n_oxygen"])
    rate_rows, thermo_rows, occ_rows, arr_rows, vh_rows = [], [], [], [], []
    for scenario, width_ps, rc_override, start_ps in SCENARIOS:
        sol_store = defaultdict(dict)
        by_group = defaultdict(list)
        for res in results:
            j = res["job"]
            if j["system"] == "solution":
                by_group[("solution", j["model"], j["metal"], "box_pbc", j["T_K"])].append(res)
            else:
                by_group[("graphite", j["model"], j["metal"], "ALL", 300)].append(res)
                if scenario == "primary":
                    by_group[("graphite", j["model"], j["metal"], j["variant"], 300)].append(res)
        for (system, model, metal, variant, temperature), members in sorted(by_group.items(), key=lambda kv: str(kv[0])):
            per_block = system == "solution"
            units, strata = unit_list(members, scenario, per_block)
            if not units:
                continue
            sc0 = next(m["scenarios"][scenario] for m in members if m["scenarios"].get(scenario) is not None)
            thermo_start = max(start_ps, THERMO_START_PS[system]) if system == "graphite" else start_ps
            meta = {"scenario": scenario, "filter_window_ps": width_ps, "rc_A": sc0["rc_A"],
                    "analysis_start_ps": start_ps, "thermo_start_ps": thermo_start,
                    "system": system, "model": model, "metal": metal, "variant": variant,
                    "n_trajectories": len([m for m in members if m["scenarios"].get(scenario) is not None])}
            c1 = concentration_single(system, n_oxygen[system])
            rr, tr, extra = summarize_group(meta, units, strata, float(temperature), c1,
                                            "block bootstrap (10 ps blocks)" if per_block
                                            else "trajectory bootstrap stratified by start configuration", rng)
            rate_rows += rr
            thermo_rows += tr
            occ_rows.append(extra["occ_row"])
            if system == "solution":
                sol_store[(model, metal)][int(temperature)] = (extra, c1)
        for (model, metal), per_t in sorted(sol_store.items()):
            for step in (3, 2, 1):
                for direction, k_key, n_key in (("dissociation", "k_d", "N_d"), ("association", "k_a", "N_a")):
                    points = []
                    for temperature, (extra, _) in sorted(per_t.items()):
                        point = step_rates(extra["total"]["N"], extra["total"]["T"], step)
                        bs = extra["boot"][step]
                        points.append({"T": temperature, "k": float(point[k_key]) if np.isfinite(point[k_key]) else 0.0,
                                       "N": float(point[n_key]), "k_boot": np.nan_to_num(bs[k_key]), "N_boot": bs[n_key]})
                    arr_rows.append({"scenario": scenario, "filter_window_ps": width_ps, "analysis_start_ps": start_ps,
                                     "system": "solution", "model": model, "metal": metal,
                                     "step": STEP_LABEL[step], "hydrolysis_step": step, "direction": direction,
                                     **arrhenius(points)})
                if scenario == "primary":
                    points = []
                    for temperature, (extra, c1) in sorted(per_t.items()):
                        dwell = extra["total"]["T_thermo"]
                        if min(dwell[step], dwell[step - 1]) < VH_MIN_DWELL_PS:
                            continue
                        th = thermo_from_occ(extra["total"]["occ_thermo"], step, temperature, c1)
                        th_b = thermo_from_occ(extra["boot_occ"], step, temperature, c1)
                        rt = R_KJ * temperature
                        points.append({"T": temperature, "lnK": float(-th["dG0_d"] / rt),
                                       "lnK_boot": -th_b["dG0_d"] / rt})
                    vh_rows.append({"system": "solution", "model": model, "metal": metal,
                                    "step": STEP_LABEL[step], "hydrolysis_step": step, **vant_hoff(points)})
        print(f"  scenario {scenario} done", flush=True)
    write_rows(KIN_OUT / "kinetics_rates.csv", rate_rows)
    write_rows(KIN_OUT / "kinetics_thermo.csv", thermo_rows)
    write_rows(KIN_OUT / "kinetics_occupancy.csv", occ_rows)
    write_rows(KIN_OUT / "kinetics_arrhenius.csv", arr_rows)
    write_rows(KIN_OUT / "kinetics_vanthoff.csv", vh_rows)

    # Occupancy versus time (10 ps blocks of the primary window) to document equilibration.
    time_rows = []
    pooled = defaultdict(lambda: defaultdict(lambda: np.zeros(N_STATE)))
    for res in results:
        sc = res["scenarios"]["primary"]
        if sc is None:
            continue
        j = res["job"]
        for b in range(sc["occ_blocks"].shape[0]):
            mid = sc["t_start_ps"] + (b + 0.5) * BLOCK_PS
            for variant in ((j["variant"], "ALL") if j["system"] == "graphite" else (j["variant"],)):
                pooled[(j["system"], j["model"], j["metal"], variant, j["T_K"])][round(mid, 1)] += sc["occ_blocks"][b]
    for key, blocks in sorted(pooled.items(), key=lambda kv: str(kv[0])):
        for mid, occ in sorted(blocks.items()):
            time_rows.append({"system": key[0], "model": key[1], "metal": key[2], "variant": key[3], "T_K": key[4],
                              "block_mid_ps": mid, "frames": int(occ.sum()),
                              **{f"P{n}": round(float(occ[n] / occ.sum()), 4) for n in range(N_STATE)},
                              "mean_n": round(float((occ * np.arange(N_STATE)).sum() / occ.sum()), 3)})
    write_rows(KIN_OUT / "kinetics_occupancy_vs_time.csv", time_rows)
    meta = {
        "description": "Stepwise Ln-OH (de)coordination kinetics and occupancy thermodynamics",
        "state": "n = number of hydroxide O (one nearest H) within r_c of Ln (clamped to 3)",
        "sampling_fs": SAMPLE_FS,
        "committed_transitions": "centred mode filter of width W then runs shorter than W/2 absorbed",
        "scenarios": [dict(zip(("name", "filter_window_ps", "rc_override_A", "analysis_start_ps"), s)) for s in SCENARIOS],
        "thermo_start_ps": THERMO_START_PS, "block_ps": BLOCK_PS, "n_bootstrap": N_BOOT,
        "eyring": "dG_act = RT ln(kB T / (h k))", "pKw": PKW,
        "standard_state": "dG0_d = dG_occ - RT ln((4-n) c1 / 1 M); pKa_eq = pKw - dG0_d/(RT ln10)",
        "c_single_M": {s: concentration_single(s, n) for s, n in n_oxygen.items()},
        "excluded_trajectories": len(skipped),
    }
    (KIN_OUT / "kinetics_metadata.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    print(f"analyze: wrote tables to {KIN_OUT}", flush=True)


# --------------------------------------------------------------------------- figures

SERIES = (
    ("graphite", "DPA4C-Plus"), ("solution", "DPA4C-Plus"), ("graphite", "uma-s-1.2.1"),
    ("solution", "DPA4-Pro"), ("solution", "MACE-POLAR-1-M"), ("solution", "uma-s-1.2.1"),
)
SERIES_COLOR = {
    ("graphite", "DPA4C-Plus"): "#1b9e77", ("solution", "DPA4C-Plus"): "#d95f02",
    ("graphite", "uma-s-1.2.1"): "#7570b3", ("solution", "DPA4-Pro"): "#e7298a",
    ("solution", "MACE-POLAR-1-M"): "#66a61e", ("solution", "uma-s-1.2.1"): "#a6761d",
}
SCENARIO_ORDER = ("W0.2", "primary", "W1.0", "W2.0", "start1", "start50", "start90", "rc2.90", "rc3.30")
SCENARIO_LABEL = {"W0.2": "W 0.2 ps", "primary": "W 0.5 ps\n(primary)", "W1.0": "W 1 ps", "W2.0": "W 2 ps",
                  "start1": "start 1 ps", "start50": "start 50 ps", "start90": "start 90 ps", "rc2.90": "r_c 2.90", "rc3.30": "r_c 3.30"}
STEP_COLOR = {3: "#1f78b4", 2: "#33a02c", 1: "#e31a1c"}


def series_label(system: str, model: str) -> str:
    return f"{MODEL_LABEL[model]} ({system})"


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def num(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"svg.fonttype": "none", "font.size": 9, "axes.titlesize": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    return plt


def fig_dgact(rates: list[dict]) -> None:
    plt = pyplot()
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    sel = {(r["system"], r["model"], r["metal"], int(r["hydrolysis_step"]), r["direction"]): r for r in rates
           if r["scenario"] == "primary" and num(r["T_K"]) == 300 and r["variant"] in ("ALL", "box_pbc")}
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), sharey=True)
    width = 0.13
    top = 30.0
    for i, metal in enumerate(("Dy", "Nd")):
        for j, direction in enumerate(("dissociation", "association")):
            ax = axes[i][j]
            for s_idx, key in enumerate(SERIES):
                color = SERIES_COLOR[key]
                for st_idx, step in enumerate((3, 2, 1)):
                    row = sel.get((key[0], key[1], metal, step, direction))
                    if row is None:
                        continue
                    x = st_idx + (s_idx - (len(SERIES) - 1) / 2) * width
                    dg, lo, hi = num(row["dG_act_kJmol"]), num(row["dG_act_CI95_lo"]), num(row["dG_act_CI95_hi"])
                    bound = num(row["dG_act_lower_bound_if_no_events"])
                    if math.isfinite(dg):
                        ax.bar(x, dg, width * 0.92, color=color, alpha=0.85, edgecolor="none")
                        hi_draw = hi if math.isfinite(hi) else top
                        ax.errorbar(x, dg, yerr=[[dg - lo], [hi_draw - dg]], fmt="none", ecolor="#333", lw=0.9, capsize=2)
                        if not math.isfinite(hi):
                            ax.plot(x, top - 0.6, marker="^", color="#333", ms=4)
                        ax.text(x, 0.4, row["N_events"], ha="center", va="bottom", fontsize=6, color="white", rotation=90)
                    elif math.isfinite(bound):
                        ax.plot(x, bound, marker="^", color=color, ms=7, mec="#333", mew=0.6)
                        ax.vlines(x, bound, top, color=color, lw=1.0, linestyles=":")
            ax.set_xticks(range(3), [f"{STEP_LABEL[s]}" for s in (3, 2, 1)])
            ax.set_ylim(0, top)
            ax.set_title(f"{metal}: {'Ln-OH loss (n -> n-1)' if direction == 'dissociation' else 'OH- uptake (n-1 -> n)'}")
            if j == 0:
                ax.set_ylabel("apparent ΔG‡ (kJ/mol), 300 K")
            ax.grid(axis="y", color="#e5e9f0", lw=0.7)
    handles = [Patch(color=SERIES_COLOR[k], label=series_label(*k)) for k in SERIES]
    handles.append(Line2D([], [], marker="^", ls="", color="#888", mec="#333", label="lower bound (no committed event)"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8)
    fig.suptitle("Stepwise apparent activation free energies (Eyring, W = 0.5 ps committed transitions; numbers in bars = events)", fontsize=10)
    fig.tight_layout(rect=(0, 0.08, 1, 0.96))
    fig.savefig(KIN_OUT / "fig_kinetics_dGact_300K.svg")
    plt.close(fig)


def fig_arrhenius(rates: list[dict], arr: list[dict]) -> None:
    plt = pyplot()
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True)
    for i, metal in enumerate(("Dy", "Nd")):
        for j, step in enumerate((3, 2, 1)):
            ax = axes[i][j]
            for direction, marker, face in (("dissociation", "o", None), ("association", "s", "white")):
                color = "#c0392b" if direction == "dissociation" else "#2471a3"
                rows = sorted((r for r in rates if r["scenario"] == "primary" and r["system"] == "solution"
                               and r["model"] == "DPA4C-Plus" and r["metal"] == metal
                               and int(r["hydrolysis_step"]) == step and r["direction"] == direction),
                              key=lambda r: num(r["T_K"]))
                for r in rows:
                    temperature, k = num(r["T_K"]), num(r["k_ps-1"])
                    lo, hi, n_ev = num(r["k_CI95_lo"]), num(r["k_CI95_hi"]), int(r["N_events"])
                    x = 1000.0 / temperature
                    if n_ev > 0 and k > 0:
                        lo_draw = lo if lo > 0 else k / 20.0
                        ax.errorbar(x, k, yerr=[[k - lo_draw], [hi - k]], fmt=marker, color=color, mfc=face or color,
                                    ms=5, capsize=2, lw=0.9)
                        ax.annotate(str(n_ev), (x, k), textcoords="offset points", xytext=(4, 3), fontsize=6, color=color)
                    elif math.isfinite(hi) and hi > 0:
                        ax.plot(x, hi, marker="v", color=color, ms=5, mfc="white")
                fit = next((a for a in arr if a["scenario"] == "primary" and a["model"] == "DPA4C-Plus" and a["metal"] == metal
                            and int(a["hydrolysis_step"]) == step and a["direction"] == direction), None)
                if fit and fit.get("Ea_kJmol") not in (None, ""):
                    temps = [float(t) for t in fit["T_used"].split(";")]
                    xs = np.linspace(1000.0 / max(temps), 1000.0 / min(temps), 20)
                    ys = np.exp(num(fit["lnA_ps-1"]) - num(fit["Ea_kJmol"]) / (R_KJ * 1000.0 / xs))
                    ax.plot(xs, ys, color=color, lw=1.2, ls="-" if direction == "dissociation" else "--")
                    label = "loss" if direction == "dissociation" else "uptake"
                    y_text = 0.95 if direction == "dissociation" else 0.86
                    ax.text(0.03, y_text, f"{label}: Ea = {fit['Ea_kJmol']} [{fit['Ea_CI95_lo']}, {fit['Ea_CI95_hi']}] ({fit['quality']})",
                            transform=ax.transAxes, fontsize=7, color=color, va="top")
            ax.set_yscale("log")
            ax.set_ylim(5e-3, 20)
            ax.set_title(f"{metal}  {STEP_LABEL[step]}")
            ax.grid(color="#e5e9f0", lw=0.6)
            if i == 1:
                ax.set_xlabel("1000 / T (1/K)")
            if j == 0:
                ax.set_ylabel("k (1/ps)")
            twin = ax.secondary_xaxis("top", functions=(lambda x: 1000.0 / np.maximum(x, 1e-9), lambda t: 1000.0 / np.maximum(t, 1e-9)))
            twin.set_xticks([300, 350, 400, 450, 500])
            twin.tick_params(labelsize=7)
    fig.suptitle("DPA4C-Plus solution: Arrhenius plots (filled circles = OH loss n->n-1, open squares = OH uptake; "
                 "v = Garwood upper bound when no event; labels = event counts)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(KIN_OUT / "fig_kinetics_arrhenius_DPA4C_solution.svg")
    plt.close(fig)


def fig_sensitivity(rates: list[dict], arr: list[dict]) -> None:
    plt = pyplot()
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))
    xs = np.arange(len(SCENARIO_ORDER))
    combos = [(s, d) for s in (3, 2, 1) for d in ("dissociation", "association")]
    for i, metal in enumerate(("Dy", "Nd")):
        ax = axes[0][i]
        for c_idx, (step, direction) in enumerate(combos):
            off = (c_idx - 2.5) * 0.1
            for s_idx, scenario in enumerate(SCENARIO_ORDER):
                fit = next((a for a in arr if a["scenario"] == scenario and a["model"] == "DPA4C-Plus" and a["metal"] == metal
                            and int(a["hydrolysis_step"]) == step and a["direction"] == direction), None)
                if not fit or fit.get("Ea_kJmol") in (None, ""):
                    continue
                ea, lo, hi = num(fit["Ea_kJmol"]), num(fit["Ea_CI95_lo"]), num(fit["Ea_CI95_hi"])
                filled = fit["quality"] in ("good", "fair")
                ax.errorbar(s_idx + off, ea, yerr=[[ea - max(lo, -45)], [min(hi, 45) - ea]], fmt="o" if direction == "dissociation" else "s",
                            color=STEP_COLOR[step], mfc=STEP_COLOR[step] if filled else "white", ms=4, lw=0.8, capsize=1.5,
                            label=f"{STEP_LABEL[step]} {'loss' if direction == 'dissociation' else 'uptake'}" if s_idx == 1 else None)
        ax.axhline(0, color="#999", lw=0.8)
        ax.set_ylim(-45, 45)
        ax.set_xticks(xs, [SCENARIO_LABEL[s] for s in SCENARIO_ORDER], fontsize=7)
        ax.set_ylabel("Arrhenius Ea (kJ/mol)")
        ax.set_title(f"DPA4C-Plus solution {metal}: Ea vs analysis choice (open = poor fit)")
        ax.grid(axis="y", color="#e5e9f0", lw=0.6)
        if i == 0:
            ax.legend(fontsize=7, ncol=3, frameon=False, loc="lower left")
        ax = axes[1][i]
        for c_idx, (step, direction) in enumerate(combos):
            off = (c_idx - 2.5) * 0.1
            for s_idx, scenario in enumerate(SCENARIO_ORDER):
                row = next((r for r in rates if r["scenario"] == scenario and r["system"] == "graphite" and r["model"] == "DPA4C-Plus"
                            and r["metal"] == metal and r["variant"] == "ALL" and int(r["hydrolysis_step"]) == step
                            and r["direction"] == direction), None)
                if row is None or row["dG_act_kJmol"] in ("", "inf"):
                    continue
                dg, lo, hi = num(row["dG_act_kJmol"]), num(row["dG_act_CI95_lo"]), num(row["dG_act_CI95_hi"])
                ax.errorbar(s_idx + off, dg, yerr=[[dg - lo], [hi - dg]], fmt="o" if direction == "dissociation" else "s",
                            color=STEP_COLOR[step], mfc=STEP_COLOR[step] if direction == "dissociation" else "white", ms=4, lw=0.8, capsize=1.5)
        ax.set_xticks(xs, [SCENARIO_LABEL[s] for s in SCENARIO_ORDER], fontsize=7)
        ax.set_ylabel("ΔG‡ at 300 K (kJ/mol)")
        ax.set_title(f"DPA4C-Plus graphite {metal}: ΔG‡ vs analysis choice")
        ax.grid(axis="y", color="#e5e9f0", lw=0.6)
    fig.tight_layout()
    fig.savefig(KIN_OUT / "fig_kinetics_sensitivity.svg")
    plt.close(fig)


def fig_occupancy(occ_rows: list[dict], time_rows: list[dict]) -> None:
    plt = pyplot()
    colors = ["#fde0c5", "#facba6", "#f59e72", "#d9534f"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    bars = []
    for metal in ("Dy", "Nd"):
        for temperature in (300, 350, 400, 450, 500):
            bars.append((f"{metal} {temperature} K\nDPA4C sol", "solution", "DPA4C-Plus", metal, temperature))
        bars.append((f"{metal} 300 K\nDPA4C graph", "graphite", "DPA4C-Plus", metal, 300))
        bars.append((f"{metal} 300 K\nUMA graph", "graphite", "uma-s-1.2.1", metal, 300))
    for b_idx, (label, system, model, metal, temperature) in enumerate(bars):
        row = next((r for r in occ_rows if r["scenario"] == "primary" and r["system"] == system and r["model"] == model
                    and r["metal"] == metal and num(r["T_K"]) == temperature and r["variant"] in ("ALL", "box_pbc")), None)
        if row is None:
            continue
        bottom = 0.0
        for n in range(N_STATE):
            p = num(row[f"P{n}"])
            ax.bar(b_idx, p, 0.8, bottom=bottom, color=colors[n], edgecolor="white", lw=0.4, label=f"n = {n}" if b_idx == 0 else None)
            bottom += p
    ax.set_xticks(range(len(bars)), [b[0] for b in bars], fontsize=6.5, rotation=90)
    ax.set_ylabel("P(n)  (n = bound OH)")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=7, ncol=4, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), borderaxespad=0.2)
    ax.set_title("Occupancy of Ln(OH)n states (solution: t >= 10 ps; graphite: t >= 50 ps)", pad=22)
    ax = axes[1]
    styles = {"3near": "-", "3far": "--", "1near2far": "-.", "2near1far": ":"}
    for metal, color in (("Dy", "#0072b2"), ("Nd", "#d95f02")):
        for variant, ls in styles.items():
            pts = sorted((num(r["block_mid_ps"]), num(r["mean_n"])) for r in time_rows if r["system"] == "graphite"
                         and r["model"] == "DPA4C-Plus" and r["metal"] == metal and r["variant"] == variant)
            if pts:
                ax.plot(*zip(*pts), ls=ls, color=color, lw=1.1, label=f"{metal} {variant}")
    ax.axvline(THERMO_START_PS["graphite"], color="#888", lw=0.8, ls=":")
    ax.text(THERMO_START_PS["graphite"] + 2, 2.85, "primary thermo window", fontsize=7, color="#666")
    ax.axvline(90.0, color="#bbb", lw=0.8, ls="--")
    ax.text(92, 2.65, "start90 sensitivity", fontsize=7, color="#888")
    ax.set_xlabel("t (ps, 10 ps blocks)")
    ax.set_ylabel("<n> bound OH")
    ax.set_ylim(0, 3)
    ax.set_title("DPA4C-Plus graphite: start-configuration memory")
    ax.legend(fontsize=6.5, ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(KIN_OUT / "fig_kinetics_occupancy.svg")
    plt.close(fig)


def fig_traces() -> None:
    plt = pyplot()
    screen = load_screen()
    cutoffs = {(r["system"], r["model"], r["metal"]): float(r["rc_A"]) for r in read_rows(KIN_OUT / "kinetics_cutoffs.csv")}
    picks = [
        ("DPA4C-Plus solution Nd 300 K", "solution", "DPA4C-Plus", "box_pbc", "Nd", 0, 300),
        ("DPA4C-Plus solution Dy 300 K", "solution", "DPA4C-Plus", "box_pbc", "Dy", 0, 300),
        ("DPA4C-Plus graphite Dy 3near r00", "graphite", "DPA4C-Plus", "3near", "Dy", 0, 300),
        ("UMA-s-1.2.1 graphite Nd 3far r00", "graphite", "uma-s-1.2.1", "3far", "Nd", 0, 300),
    ]
    jobs = {(j["system"], j["model"], j["variant"], j["metal"], j["replica"], j["T_K"]): j for j in all_jobs()}
    fig, axes = plt.subplots(len(picks) + 1, 1, figsize=(12, 9.5))
    rows_out = []
    zoom = None
    for ax, (label, system, model, variant, metal, replica, temperature) in zip(axes, picks):
        job = jobs[(system, model, variant, metal, replica, temperature)]
        srow = screen[(system, model, variant, metal, replica, temperature)]
        rc = cutoffs[(system, model, metal)]
        with np.load(cache_path(job)) as data:
            t, oh_d = data["t_ps"], data["oh_d"]
        keep = t <= float(srow["usable_end_ps"]) + 1e-9
        t, oh_d = t[keep], oh_d[keep]
        raw = np.minimum((oh_d < rc).sum(axis=1), N_STATE - 1)
        committed = committed_states(raw, 0.5)
        ax.plot(t, raw, color="#bbbbbb", lw=0.4, drawstyle="steps-post", label="raw (20 fs)")
        ax.plot(t, committed, color="#c0392b", lw=1.2, drawstyle="steps-post", label="committed (W = 0.5 ps)")
        ax.axvspan(0, 10, color="#eeeeee", zorder=0)
        ax.set_ylim(-0.3, 3.3)
        ax.set_yticks(range(4))
        ax.set_ylabel("n")
        ax.set_title(f"{label} (r_c = {rc:.2f} Å, {srow['status']})", fontsize=9, loc="left")
        for ti, ri, ci in zip(t, raw, committed):
            rows_out.append({"trajectory": label, "t_ps": round(float(ti), 3), "n_raw": int(ri), "n_committed_W0.5ps": int(ci)})
        if zoom is None:
            zoom = (label, t, raw, committed)
    axes[0].legend(fontsize=7, frameon=False, loc="upper right", ncol=2)
    label, t, raw, committed = zoom
    ax = axes[-1]
    window = (t >= 40) & (t <= 50)
    ax.plot(t[window], raw[window], color="#999999", lw=0.7, drawstyle="steps-post", label="raw")
    ax.plot(t[window], committed[window], color="#c0392b", lw=1.4, drawstyle="steps-post", label="committed")
    ax.set_ylim(-0.3, 3.3)
    ax.set_yticks(range(4))
    ax.set_xlabel("t (ps)")
    ax.set_ylabel("n")
    ax.set_title(f"zoom 40-50 ps: {label} - fast proton-shuttling flicker is filtered", fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig(KIN_OUT / "fig_kinetics_traces.svg")
    plt.close(fig)
    write_rows(KIN_OUT / "kinetics_traces_example.csv", rows_out)


def fig_density() -> None:
    plt = pyplot()
    screen = load_screen()
    cutoffs = {(r["system"], r["model"], r["metal"]): r for r in read_rows(KIN_OUT / "kinetics_cutoffs.csv")}
    windows = defaultdict(list)
    for job in all_jobs():
        row = screen[(job["system"], job["model"], job["variant"], job["metal"], int(job["replica"]), int(job["T_K"]))]
        if row["status"] != "excluded":
            windows[group_key(job)].append((cache_path(job), 10.0, float(row["usable_end_ps"])))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    rows_out = []
    for key, items in sorted(windows.items()):
        r, density = oh_density(items)
        rel = density / density.max()
        ax = axes[0] if key[0] == "graphite" else axes[1]
        color = SERIES_COLOR[(key[0], key[1])]
        ls = "-" if key[2] == "Dy" else "--"
        ax.semilogy(r, np.maximum(rel, 1e-7), color=color, ls=ls, lw=1.0, label=f"{MODEL_LABEL[key[1]]} {key[2]}")
        rc = float(cutoffs[key]["rc_A"])
        ax.axvline(rc, color=color, ls=ls, lw=0.6, alpha=0.7)
        for ri, di, gi in zip(r, density, rel):
            rows_out.append({"system": key[0], "model": key[1], "metal": key[2], "r_A": round(float(ri), 3),
                             "density_counts_over_r2": float(di), "relative_to_max": float(gi), "rc_A": rc})
    for ax, title in zip(axes, ("graphite (300 K)", "solution (all T)")):
        ax.set_xlim(1.5, 6.0)
        ax.set_ylim(1e-7, 2)
        ax.set_xlabel("Ln-O(hydroxide) distance (Å)")
        ax.set_title(f"OH- oxygen radial density around Ln, {title}; vertical lines = r_c")
        ax.legend(fontsize=7, frameon=False, loc="upper right")
        ax.grid(color="#eef1f5", lw=0.6)
    axes[0].set_ylabel("relative density (counts / r², floor 1e-7)")
    fig.tight_layout()
    fig.savefig(KIN_OUT / "fig_kinetics_OH_density_cutoff.svg")
    plt.close(fig)
    write_rows(KIN_OUT / "kinetics_oh_density.csv", rows_out)


def run_figures() -> None:
    rates = read_rows(KIN_OUT / "kinetics_rates.csv")
    arr = read_rows(KIN_OUT / "kinetics_arrhenius.csv")
    fig_dgact(rates)
    fig_arrhenius(rates, arr)
    fig_sensitivity(rates, arr)
    fig_occupancy(read_rows(KIN_OUT / "kinetics_occupancy.csv"), read_rows(KIN_OUT / "kinetics_occupancy_vs_time.csv"))
    fig_traces()
    fig_density()
    print("figures written to", KIN_OUT, flush=True)


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "extract"
    if stage == "extract":
        workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
        run_extract(workers)
    elif stage == "extract-one":
        index = int(sys.argv[2])
        print(extract(all_jobs()[index]))
    elif stage == "screen":
        run_screen()
    elif stage == "analyze":
        run_analyze(int(sys.argv[2]) if len(sys.argv) > 2 else 8)
    elif stage == "figures":
        run_figures()
    else:
        raise SystemExit(f"unknown stage {stage}")
