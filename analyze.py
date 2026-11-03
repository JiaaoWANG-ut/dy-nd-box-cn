#!/usr/bin/env python3
"""Box-PBC RDF/CN comparison for Dy and Nd in solution and on graphite."""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
GRAPHITE = ROOT / "Graphite-dy-nd-50-denser-dumping"
SOLUTION = ROOT / "md" / "dense-gr-cn"
SCREEN_CSV = OUT / "collapse_screen.csv"

SAMPLE_FS = 1000.0
EQUIL_FS = 10_000.0
DR = 0.02
R_MAX = 8.50
EDGES = np.arange(0.0, R_MAX + DR / 2, DR)
R = (EDGES[:-1] + EDGES[1:]) / 2
SHELL_VOLUME = 4.0 * np.pi * R**2 * DR
CHUNK_FRAMES = 8
O_W, O_OH, O_OTHER = 0, 1, 2
REE, CARBON, O_W_LABEL, O_OH_LABEL, O_OTHER_LABEL, H_W_LABEL, H_OH_LABEL, H_OTHER_LABEL = range(8)
LABEL_NAMES = (
    "REE", "C", "O_water", "O_hydroxide", "O_other",
    "H_water", "H_hydroxide", "H_other",
)
CN_NAMES = (
    "REE_O_water", "REE_O_hydroxide", "REE_O_other", "REE_O_total",
    "OH_REE_first", "OH_Hwater_first", "OH_Owater_first", "OH_Xtotal_first",
    "OH_REE_shell", "OH_Hwater_shell", "OH_Owater_shell", "OH_Xtotal_shell",
)


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def smooth(y: np.ndarray, width: int = 5) -> np.ndarray:
    return np.convolve(y, np.ones(width) / width, mode="same")


def load_collapse_screen() -> dict[tuple[str, str, str, int], dict]:
    """Per-trajectory usable window written by `kinetics.py screen`."""
    if not SCREEN_CSV.is_file():
        raise FileNotFoundError(f"{SCREEN_CSV} missing; run `kinetics.py extract` then `kinetics.py screen`")
    screen = {}
    with SCREEN_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["system"] == "graphite":
                screen[(row["model"], row["variant"], row["metal"], int(row["replica"]))] = row
    return screen


def screen_row(job: dict, screen: dict) -> dict:
    return screen[(job["model"], job["variant"], job["metal"], int(job["replica"]))]


def analysis_frames(meta: dict, usable_end_ps: float) -> np.ndarray:
    """1 ps samples from EQUIL_FS up to the last pre-collapse time (frame time = index * dump interval)."""
    dump_fs = float(meta["dump_interval_fs"])
    stride = int(round(SAMPLE_FS / dump_fs))
    start = int(round(EQUIL_FS / dump_fs))
    frame_ids = np.arange(start, meta["frames_written"], stride, dtype=np.int64)
    return frame_ids[frame_ids * dump_fs <= usable_end_ps * 1000.0 + 1e-6]


def first_min_after(y: np.ndarray, peak: int, max_r: float) -> int | None:
    ys = smooth(y)
    end = min(int(np.searchsorted(R, max_r)), len(R) - 1)
    for i in range(peak + 1, end):
        if ys[i] <= ys[i - 1] and ys[i] < ys[i + 1]:
            return i
    return None


def deep_minimum(y: np.ndarray, lo: float, hi: float, window: float) -> tuple[float, float, float]:
    ys = smooth(y)
    candidates = np.flatnonzero((R >= lo) & (R <= hi))
    if candidates.size == 0 or np.max(ys[candidates]) <= 0:
        return (math.nan, math.nan, math.nan)
    peak = int(candidates[np.argmax(ys[candidates])])
    end = min(int(np.searchsorted(R, R[peak] + window)), len(R) - 1)
    minimum = peak + int(np.argmin(ys[peak : end + 1]))
    return (float(R[peak]), float(ys[peak]), float(R[minimum]))


def real_first_peak(y: np.ndarray, lo: float, hi: float) -> int | None:
    """First peak with g > 0.5 that drops by at least 10% before rising again (as in md/dense-gr-cn/gr_cn.py)."""
    ys = smooth(y)
    for i in np.flatnonzero((R >= lo) & (R <= hi))[1:-1]:
        if ys[i] > 0.5 and ys[i] >= ys[i - 1] and ys[i] >= ys[i + 1]:
            j = i
            while j + 1 < len(ys) and ys[j + 1] <= ys[j]:
                j += 1
            if ys[j] < 0.9 * ys[i]:
                return int(i)
    return None


def minimum_after(y: np.ndarray, r0: float, r_hi: float) -> tuple[int, int]:
    """Climb to the local maximum at or after r0, then descend to the next minimum."""
    ys = smooth(y)
    i = int(np.searchsorted(R, r0))
    last = int(np.searchsorted(R, r_hi))
    while i + 1 < last and ys[i + 1] >= ys[i]:
        i += 1
    peak = i
    while i + 1 < last and ys[i + 1] <= ys[i]:
        i += 1
    return peak, i


def cutoff_data(rdf: dict[str, np.ndarray]) -> dict[str, float]:
    rp, rg, rc_ree = deep_minimum(rdf["REE_O_total"], 1.6, 3.4, 1.6)
    if not math.isfinite(rc_ree):
        raise ValueError("Could not locate the first-shell minimum of REE-O total g(r)")

    p_oh_total = real_first_peak(rdf["OH_X_total"], 1.2, 4.5)
    m_oh_first = first_min_after(rdf["OH_X_total"], p_oh_total, 5.5) if p_oh_total is not None else None
    r_oh_ree, g_oh_ree, _ = deep_minimum(rdf["OH_REE"], 1.9, 2.8, 1.5)
    r_ree_peak = r_oh_ree if math.isfinite(g_oh_ree) and g_oh_ree > 0.5 else 2.25
    p_oh_ree, m_oh_shell = minimum_after(rdf["OH_X_total"], r_ree_peak - 0.1, 3.0)
    g_oh_total = smooth(rdf["OH_X_total"])

    def drop_fraction(peak: int | None, minimum: int | None) -> float:
        if peak is None or minimum is None or g_oh_total[peak] <= 0:
            return math.nan
        return float(1.0 - g_oh_total[minimum] / g_oh_total[peak])

    first_drop = drop_fraction(p_oh_total, m_oh_first)
    shell_drop = drop_fraction(p_oh_ree, m_oh_shell)
    return {
        "REE_O_peak_A": rp,
        "REE_O_peak_g": rg,
        "REE_O_total_min_A": rc_ree,
        "OH_X_first_peak_A": float(R[p_oh_total]) if p_oh_total is not None else math.nan,
        "OH_X_first_min_A": float(R[m_oh_first]) if m_oh_first is not None else math.nan,
        "OH_X_first_drop_fraction": first_drop,
        "OH_X_first_min_clear": bool(first_drop >= 0.10) if math.isfinite(first_drop) else False,
        "OH_REE_peak_A": r_oh_ree,
        "OH_REE_peak_g": g_oh_ree,
        "OH_REE_peak_used_A": r_ree_peak,
        "OH_X_shell_peak_A": float(R[p_oh_ree]),
        "OH_X_shell_min_A": float(R[m_oh_shell]) if m_oh_shell is not None else math.nan,
        "OH_X_shell_drop_fraction": shell_drop,
        "OH_X_shell_min_clear": bool(shell_drop >= 0.10) if math.isfinite(shell_drop) else False,
    }


def min_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box * np.rint(delta / box)


def classify_chunk(
    pos: np.ndarray, oxygen: np.ndarray, hydrogen: np.ndarray, metal_idx: int, carbon: np.ndarray, box: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_frames = pos.shape[0]
    o_pos = pos[:, oxygen, :]
    h_pos = pos[:, hydrogen, :]
    delta = min_image(h_pos[:, :, None, :] - o_pos[:, None, :, :], box)
    host = np.sum(delta * delta, axis=-1).argmin(axis=2)
    n_h = np.zeros((n_frames, oxygen.size), dtype=np.int16)
    frame_ids = np.arange(n_frames)
    for frame in range(n_frames):
        np.add.at(n_h[frame], host[frame], 1)
    o_type = np.where(n_h == 2, O_W, np.where(n_h == 1, O_OH, O_OTHER)).astype(np.int8)

    labels = np.full((n_frames, pos.shape[1]), -1, dtype=np.int8)
    labels[:, metal_idx] = REE
    labels[:, carbon] = CARBON
    labels[:, oxygen] = np.where(
        o_type == O_W, O_W_LABEL, np.where(o_type == O_OH, O_OH_LABEL, O_OTHER_LABEL)
    )
    h_type = o_type[frame_ids[:, None], host]
    labels[:, hydrogen] = np.where(
        h_type == O_W, H_W_LABEL, np.where(h_type == O_OH, H_OH_LABEL, H_OTHER_LABEL)
    )
    return labels, host, o_type


def bins_for(dist: np.ndarray) -> np.ndarray:
    return np.floor(dist / DR).astype(np.int32)


def add_hist(destination: np.ndarray, bins: np.ndarray, mask: np.ndarray | None = None) -> None:
    if mask is not None:
        bins = bins[mask]
    bins = bins[(bins >= 0) & (bins < len(R))]
    destination += np.bincount(bins, minlength=len(R))


def accumulate_rdf_for_job(job: dict, stats: dict, usable_end_ps: float) -> None:
    base = GRAPHITE / "runs" / job["model"] / job["variant"] / job["metal"] / f"r{job['replica']:02d}"
    meta = json.loads((base / "traj-denser-dumping-meta.json").read_text())
    positions = np.load(base / "traj-denser-dumping-positions.npy", mmap_mode="r")
    symbols = np.asarray(meta["symbols"])
    oxygen = np.flatnonzero(symbols == "O")
    hydrogen = np.flatnonzero(symbols == "H")
    carbon = np.flatnonzero(symbols == "C")
    metal_idx = int(np.flatnonzero(symbols == job["metal"])[0])
    box = np.asarray([meta["cell_A"][i][i] for i in range(3)], dtype=np.float64)
    if not np.all(meta["pbc"]):
        raise ValueError(f"Expected box PBC for {base}")
    frame_ids = analysis_frames(meta, usable_end_ps)
    key = (job["model"], job["metal"])
    item = stats.setdefault(
        key,
        {
            "volume": float(np.prod(box)),
            "ree_hist": np.zeros((3, len(R)), dtype=np.float64),
            "ree_exposure": np.zeros(3, dtype=np.float64),
            "oh_hist": np.zeros((9, len(R)), dtype=np.float64),
            "oh_exposure": np.zeros(9, dtype=np.float64),
            "n_oh_centres": 0,
            "n_frames": 0,
        },
    )
    if not np.allclose(item["volume"], np.prod(box), atol=1e-5):
        raise ValueError(f"Inconsistent box volumes in {key}")

    for start_chunk in range(0, frame_ids.size, CHUNK_FRAMES):
        chosen = frame_ids[start_chunk : start_chunk + CHUNK_FRAMES]
        pos = np.asarray(positions[chosen], dtype=np.float64)
        labels, host, o_type = classify_chunk(pos, oxygen, hydrogen, metal_idx, carbon, box)
        o_xyz = pos[:, oxygen, :]
        delta_ree_o = min_image(o_xyz - pos[:, metal_idx : metal_idx + 1, :], box)
        dist_ree_o = np.sqrt(np.sum(delta_ree_o * delta_ree_o, axis=-1))
        b_ree_o = bins_for(dist_ree_o)
        for typ in (O_W, O_OH, O_OTHER):
            mask = o_type == typ
            add_hist(item["ree_hist"][typ], b_ree_o, mask)
            item["ree_exposure"][typ] += int(mask.sum())

        f_idx, o_local = np.nonzero(o_type == O_OH)
        if f_idx.size:
            o_global = oxygen[o_local]
            centers = pos[f_idx, o_global, :]
            delta = min_image(pos[f_idx] - centers[:, None, :], box)
            dist = np.sqrt(np.sum(delta * delta, axis=-1))
            bins = bins_for(dist)
            partner_labels = labels[f_idx]
            valid = np.ones(partner_labels.shape, dtype=bool)
            valid[np.arange(f_idx.size), o_global] = False
            valid[:, hydrogen] = host[f_idx] != o_local[:, None]
            partner_labels = partner_labels[valid]
            partner_bins = bins[valid]
            in_range = (partner_bins >= 0) & (partner_bins < len(R))
            combo = partner_labels[in_range].astype(np.int32) * len(R) + partner_bins[in_range]
            item["oh_hist"][:8] += np.bincount(
                combo, minlength=8 * len(R)
            ).reshape(8, len(R))
            counts = np.bincount(partner_labels, minlength=8)
            item["oh_exposure"][:8] += counts
            item["n_oh_centres"] += f_idx.size
        item["n_frames"] += len(chosen)


def make_rdf(stats: dict) -> dict[str, np.ndarray]:
    volume = stats["volume"]
    ree_h, ree_n = stats["ree_hist"], stats["ree_exposure"]
    oh_h, oh_n = stats["oh_hist"], stats["oh_exposure"]
    ree_total_h, ree_total_n = ree_h.sum(axis=0), ree_n.sum()
    oh_total_h, oh_total_n = oh_h[:8].sum(axis=0), oh_n[:8].sum()

    def normalized(hist: np.ndarray, exposure: float) -> np.ndarray:
        if exposure <= 0:
            return np.zeros_like(R)
        return hist * volume / (exposure * SHELL_VOLUME)

    out = {
        "REE_O_water": normalized(ree_h[O_W], ree_n[O_W]),
        "REE_O_hydroxide": normalized(ree_h[O_OH], ree_n[O_OH]),
        "REE_O_other": normalized(ree_h[O_OTHER], ree_n[O_OTHER]),
        "REE_O_total": normalized(ree_total_h, ree_total_n),
        "OH_X_total": normalized(oh_total_h, oh_total_n),
        "OH_REE": normalized(oh_h[REE], oh_n[REE]),
        "OH_H_water": normalized(oh_h[H_W_LABEL], oh_n[H_W_LABEL]),
        "OH_O_water": normalized(oh_h[O_W_LABEL], oh_n[O_W_LABEL]),
    }
    return out


def count_cn_for_job(job: dict, cutoffs: dict[str, float], screen: dict) -> dict:
    base = GRAPHITE / "runs" / job["model"] / job["variant"] / job["metal"] / f"r{job['replica']:02d}"
    meta = json.loads((base / "traj-denser-dumping-meta.json").read_text())
    positions = np.load(base / "traj-denser-dumping-positions.npy", mmap_mode="r")
    symbols = np.asarray(meta["symbols"])
    oxygen = np.flatnonzero(symbols == "O")
    hydrogen = np.flatnonzero(symbols == "H")
    carbon = np.flatnonzero(symbols == "C")
    metal_idx = int(np.flatnonzero(symbols == job["metal"])[0])
    box = np.asarray([meta["cell_A"][i][i] for i in range(3)], dtype=np.float64)
    frame_ids = analysis_frames(meta, float(screen["usable_end_ps"]))
    values = {name: np.zeros(frame_ids.size, dtype=np.float64) for name in CN_NAMES}
    n_oh_frames = np.zeros(frame_ids.size, dtype=np.int64)
    rc_ree = cutoffs["REE_O_total_min_A"]
    rc_first = cutoffs["OH_X_first_min_A"]
    rc_shell = cutoffs["OH_X_shell_min_A"]

    for start_chunk in range(0, frame_ids.size, CHUNK_FRAMES):
        sl = slice(start_chunk, min(start_chunk + CHUNK_FRAMES, frame_ids.size))
        chosen = frame_ids[sl]
        pos = np.asarray(positions[chosen], dtype=np.float64)
        labels, host, o_type = classify_chunk(pos, oxygen, hydrogen, metal_idx, carbon, box)
        d = min_image(pos[:, oxygen, :] - pos[:, metal_idx : metal_idx + 1, :], box)
        d = np.sqrt(np.sum(d * d, axis=-1))
        near = d < rc_ree
        values["REE_O_water"][sl] = (near & (o_type == O_W)).sum(axis=1)
        values["REE_O_hydroxide"][sl] = (near & (o_type == O_OH)).sum(axis=1)
        values["REE_O_other"][sl] = (near & (o_type == O_OTHER)).sum(axis=1)
        values["REE_O_total"][sl] = near.sum(axis=1)

        f_idx, o_local = np.nonzero(o_type == O_OH)
        n_oh = (o_type == O_OH).sum(axis=1)
        if f_idx.size:
            o_global = oxygen[o_local]
            centers = pos[f_idx, o_global, :]
            delta = min_image(pos[f_idx] - centers[:, None, :], box)
            dist = np.sqrt(np.sum(delta * delta, axis=-1))
            neighbor_labels = labels[f_idx]
            valid = np.ones(neighbor_labels.shape, dtype=bool)
            valid[np.arange(f_idx.size), o_global] = False
            valid[:, hydrogen] = host[f_idx] != o_local[:, None]
            within_first = (dist < rc_first) & valid
            within_shell = (dist < rc_shell) & valid
            per_center = {
                "OH_REE_first": (within_first & (neighbor_labels == REE)).sum(axis=1),
                "OH_Hwater_first": (within_first & (neighbor_labels == H_W_LABEL)).sum(axis=1),
                "OH_Owater_first": (within_first & (neighbor_labels == O_W_LABEL)).sum(axis=1),
                "OH_Xtotal_first": within_first.sum(axis=1),
                "OH_REE_shell": (within_shell & (neighbor_labels == REE)).sum(axis=1),
                "OH_Hwater_shell": (within_shell & (neighbor_labels == H_W_LABEL)).sum(axis=1),
                "OH_Owater_shell": (within_shell & (neighbor_labels == O_W_LABEL)).sum(axis=1),
                "OH_Xtotal_shell": within_shell.sum(axis=1),
            }
            for name, count in per_center.items():
                np.add.at(values[name], f_idx + start_chunk, count)
        n_oh_frames[sl] = n_oh
    n_oh_total = int(n_oh_frames.sum())
    return {
        "model": job["model"],
        "metal": job["metal"],
        "variant": job["variant"],
        "replica": job["replica"],
        "screen_status": screen["status"],
        "usable_end_ps": float(screen["usable_end_ps"]),
        "n_samples": len(frame_ids),
        "n_oh_centres": n_oh_total,
        "rc_REE_O_A": rc_ree,
        "rc_OX_first_A": rc_first,
        "rc_OX_shell_A": rc_shell,
        **{name: float(values[name].mean()) for name in CN_NAMES[:4]},
        # OH-centred CN: average over all OH centres (not over frames).
        **{name: float(values[name].sum() / n_oh_total) if n_oh_total else math.nan for name in CN_NAMES[4:]},
    }


def weighted_mean_sem(values: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    """Weighted mean with a cluster (per-trajectory) ratio-estimator standard error."""
    ok = np.isfinite(values) & (weights > 0)
    x, w = values[ok], weights[ok].astype(np.float64)
    if x.size == 0:
        return math.nan, math.nan
    mean = float(np.sum(w * x) / np.sum(w))
    if x.size < 2:
        return mean, 0.0
    sem = math.sqrt(x.size / (x.size - 1) * np.sum((w * (x - mean)) ** 2)) / float(np.sum(w))
    return mean, sem


def summarize_replicas(replica_rows: list[dict], excluded: dict[tuple[str, str, str], int]) -> list[dict]:
    groups: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in replica_rows:
        groups[(row["model"], row["metal"], row["variant"])].append(row)
        groups[(row["model"], row["metal"], "ALL")].append(row)
    summaries = []
    for (model, metal, variant), rows in sorted(groups.items()):
        summary = {
            "model": model,
            "metal": metal,
            "variant": variant,
            "collapsed": 0,
            "n_excluded": excluded.get((model, metal, variant), 0),
            "n_replicas": len(rows),
            "n_truncated": sum(1 for row in rows if row["screen_status"] == "truncated"),
            "n_samples": sum(int(row["n_samples"]) for row in rows),
            "rc_REE_O_A": rows[0]["rc_REE_O_A"],
            "rc_OX_first_A": rows[0]["rc_OX_first_A"],
            "rc_OX_shell_A": rows[0]["rc_OX_shell_A"],
        }
        for metric in CN_NAMES:
            data = np.asarray([row[metric] for row in rows], dtype=np.float64)
            weight_col = "n_samples" if metric.startswith("REE_") else "n_oh_centres"
            weights = np.asarray([row[weight_col] for row in rows], dtype=np.float64)
            mean, sem = weighted_mean_sem(data, weights)
            summary[f"{metric}_mean"] = mean
            summary[f"{metric}_sem_replica"] = sem
            summary[f"{metric}_sem"] = sem
        summaries.append(summary)
    return summaries


def load_solution() -> tuple[list[dict], dict[tuple[str, str], dict[str, np.ndarray]]]:
    with (SOLUTION / "cn_table.csv").open(newline="") as handle:
        all_rows = list(csv.DictReader(handle))
    box_rows = [row for row in all_rows if row["variant"] == "box_pbc"]
    derived = (
        "OH_REE_first_from_gr", "OH_Hwater_first_from_gr",
        "OH_Xtotal_first_from_gr", "OH_REE_shell_from_gr",
        "OH_Hwater_shell_from_gr", "OH_Xtotal_shell_from_gr",
    )
    for row in box_rows:
        row.update({name: "" for name in derived})
        if row["collapsed"] != "0":
            continue
        model, metal, temperature = row["model"], row["metal"], row["T_K"]
        stem = f"{model}__box_pbc__{metal}__T{temperature}"
        gr_path = SOLUTION / "results" / f"{stem}_gr.npz"
        result_path = SOLUTION / "results" / f"{stem}.json"
        result = json.loads(result_path.read_text())
        with np.load(gr_path) as data:
            radii = data["r"]
            dr = float(radii[1] - radii[0])
            shell_volume = 4.0 * np.pi * radii**2 * dr
            means = result["mean_counts"]
            neighbor_counts = {
                "OH_REE": means["REE"],
                "OH_Hwater": means["H_w"],
                "OH_Xtotal": sum(means.values()) - 2.0,
            }
            curves = {
                "OH_REE": data["OH__REE"],
                "OH_Hwater": data["OH__H_w"],
                "OH_Xtotal": data["OH__X_total"],
            }
            for suffix, cutoff_col in (
                ("first", "rc_OH_X_first_A"),
                ("shell", "rc_OH_X_shell_A"),
            ):
                within = radii < float(row[cutoff_col])
                for name, curve in curves.items():
                    row[f"{name}_{suffix}_from_gr"] = float(
                        np.sum(curve[within] * shell_volume[within])
                        * neighbor_counts[name] / (15.781866235591567**3)
                    )
    write_csv(OUT / "solution_cn_box_pbc.csv", box_rows)

    gr_data = {}
    gr_rows = []
    for row in box_rows:
        if row["T_K"] != "300" or row["collapsed"] != "0":
            continue
        model, metal = row["model"], row["metal"]
        path = SOLUTION / "results" / f"{model}__box_pbc__{metal}__T300_gr.npz"
        with np.load(path) as data:
            gr_data[(model, metal)] = {
                "r": data["r"],
                "REE_O_water": data["REE__O_w"],
                "REE_O_hydroxide": data["REE__O_OH"],
                "REE_O_total": data["REE__O_total"],
                "OH_X_total": data["OH__X_total"],
                "OH_REE": data["OH__REE"],
                "OH_H_water": data["OH__H_w"],
            }
            for pair, curve in gr_data[(model, metal)].items():
                if pair == "r":
                    continue
                gr_rows.extend(
                    {
                        "model": model, "metal": metal, "T_K": 300,
                        "pair": pair, "r_A": float(radius), "g_r": float(value),
                    }
                    for radius, value in zip(data["r"], curve)
                )
    write_csv(OUT / "solution_gr_box_pbc_300K.csv", gr_rows)
    return box_rows, gr_data


def write_solution_cutoff_clarity(box_rows: list[dict], gr_data: dict) -> None:
    """Shell clarity of the 300 K solution O(OH)-X total g(r) at the source cutoffs, using the graphite 10% rule.

    first: highest smoothed g(r) between 0.9 A and rc_first versus g(rc_first).
    shell: local maximum reached by climbing from the O(OH)-REE peak - 0.1 A versus g(rc_shell).
    """
    rows = []
    for row in box_rows:
        key = (row["model"], row["metal"])
        if row["T_K"] != "300" or row["collapsed"] != "0" or key not in gr_data:
            continue
        radii = gr_data[key]["r"]
        g_total = smooth(gr_data[key]["OH_X_total"])
        g_ree = smooth(gr_data[key]["OH_REE"])
        rc_first, rc_shell = float(row["rc_OH_X_first_A"]), float(row["rc_OH_X_shell_A"])
        at = lambda r: int(np.argmin(np.abs(radii - r)))  # noqa: E731
        window = np.flatnonzero((radii >= 0.9) & (radii <= rc_first))
        first_peak = int(window[np.argmax(g_total[window])])
        ree_window = np.flatnonzero((radii >= 1.9) & (radii <= 2.8))
        ree_peak = int(ree_window[np.argmax(g_ree[ree_window])])
        r_ree_peak = float(radii[ree_peak]) if g_ree[ree_peak] > 0.5 else 2.25
        shell_peak = at(r_ree_peak - 0.1)
        while shell_peak + 1 < at(rc_shell) and g_total[shell_peak + 1] >= g_total[shell_peak]:
            shell_peak += 1
        out = {"model": row["model"], "metal": row["metal"], "T_K": 300,
               "rc_REE_O_A": float(row["rc_Ln_O_total_A"]), "OH_REE_peak_A": float(radii[ree_peak]),
               "OH_REE_peak_g": float(g_ree[ree_peak])}
        for suffix, peak, rc in (("first", first_peak, rc_first), ("shell", shell_peak, rc_shell)):
            drop = float(1.0 - g_total[at(rc)] / g_total[peak]) if g_total[peak] > 0 else math.nan
            out.update({
                f"OH_X_{suffix}_peak_A": float(radii[peak]),
                f"OH_X_{suffix}_min_A": rc,
                f"OH_X_{suffix}_drop_fraction": drop,
                f"OH_X_{suffix}_min_clear": bool(drop >= 0.10) if math.isfinite(drop) else False,
            })
        rows.append(out)
    write_csv(OUT / "solution_cutoffs_box_pbc_300K.csv", rows)


def write_rdf_csv(rdf_by_group: dict, cutoff_by_group: dict) -> None:
    rows = []
    for (model, metal), rdf in sorted(rdf_by_group.items()):
        cuts = cutoff_by_group[(model, metal)]
        for pair, curve in rdf.items():
            for radius, value in zip(R, curve):
                rows.append({
                    "system": "graphite_box_pbc",
                    "model": model,
                    "metal": metal,
                    "pair": pair,
                    "r_A": float(radius),
                    "g_r": float(value),
                    "rc_REE_O_A": cuts["REE_O_total_min_A"],
                    "rc_OX_first_A": cuts["OH_X_first_min_A"],
                    "rc_OX_shell_A": cuts["OH_X_shell_min_A"],
                    "sample_fs": SAMPLE_FS,
                    "equilibration_ps": EQUIL_FS / 1000.0,
                })
    write_csv(OUT / "graphite_gr_box_pbc.csv", rows)


def write_cn_figure(
    path: Path,
    summaries: list[dict],
    models: list[str],
    metrics: list[tuple[str, str, str]],
    variant: str,
) -> None:
    width, height = 1540, 1250
    left, top = 72, 48
    panel_w, panel_h = 330, 350
    gap_x, gap_y = 20, 45
    colors = {"Dy": "#2878b5", "Nd": "#d04a35"}
    rows = [r for r in summaries if r["variant"] == variant]
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:21px;font-weight:bold}.label{font-size:12px}.small{font-size:10px}.paneltitle{font-size:14px;font-weight:bold}</style>',
        f'<text x="{width/2}" y="25" text-anchor="middle" class="title">{path.stem.replace("_", " ")} — {variant}</text>',
    ]
    for p, (metric, label, err_col) in enumerate(metrics):
        row_i, col_i = divmod(p, 4)
        x0 = left + col_i * (panel_w + gap_x)
        y0 = top + row_i * (panel_h + gap_y)
        plot_x, plot_y = x0 + 48, y0 + 35
        plot_w, plot_h = panel_w - 62, panel_h - 94
        subset = [
            r for r in rows
            if r["model"] in models and str(r.get("collapsed", "0")) != "1"
        ]
        vals = [r[f"{metric}_mean"] for r in subset if math.isfinite(r[f"{metric}_mean"])]
        max_y = max(vals, default=1.0)
        max_y = max(max_y * 1.20, 0.5)
        svg.append(f'<text x="{x0+panel_w/2}" y="{y0+16}" text-anchor="middle" class="paneltitle">{label}</text>')
        svg.append(f'<line x1="{plot_x}" y1="{plot_y+plot_h}" x2="{plot_x+plot_w}" y2="{plot_y+plot_h}" stroke="#444"/>')
        svg.append(f'<line x1="{plot_x}" y1="{plot_y}" x2="{plot_x}" y2="{plot_y+plot_h}" stroke="#444"/>')
        for tick in range(1, 5):
            yy = plot_y + plot_h * (1 - tick / 4)
            svg.append(f'<line x1="{plot_x}" y1="{yy:.1f}" x2="{plot_x+plot_w}" y2="{yy:.1f}" stroke="#ddd"/>')
            svg.append(f'<text x="{plot_x-6}" y="{yy+4:.1f}" text-anchor="end" class="small">{max_y*tick/4:.2g}</text>')
        slot = plot_w / max(len(models), 1)
        bar_w = min(24, slot / 4)
        for i, model in enumerate(models):
            center = plot_x + slot * (i + 0.5)
            svg.append(f'<text x="{center:.1f}" y="{plot_y+plot_h+16}" text-anchor="middle" class="small">{model.replace("uma-s-1.2.1","UMA").replace("MACE-POLAR-1-M","MACE-POLAR").replace("DPA4C-Plus","DPA4C").replace("DPA4-Pro","DPA4-Pro")}</text>')
            for j, metal in enumerate(("Dy", "Nd")):
                match = next((r for r in subset if r["model"] == model and r["metal"] == metal), None)
                if match is None:
                    continue
                value = match[f"{metric}_mean"]
                error = match.get(err_col, 0.0) if err_col else 0.0
                if not math.isfinite(value):
                    continue
                bx = center + (j - 0.5) * (bar_w + 4) - bar_w / 2
                bh = value / max_y * plot_h
                by = plot_y + plot_h - bh
                svg.append(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" fill="{colors[metal]}"/>')
                if error > 0:
                    ey = error / max_y * plot_h
                    cx = bx + bar_w / 2
                    svg.append(f'<line x1="{cx:.1f}" y1="{by-ey:.1f}" x2="{cx:.1f}" y2="{by+ey:.1f}" stroke="#222"/>')
                    svg.append(f'<line x1="{cx-4:.1f}" y1="{by-ey:.1f}" x2="{cx+4:.1f}" y2="{by-ey:.1f}" stroke="#222"/>')
                    svg.append(f'<line x1="{cx-4:.1f}" y1="{by+ey:.1f}" x2="{cx+4:.1f}" y2="{by+ey:.1f}" stroke="#222"/>')
        svg.append(f'<text x="{plot_x+plot_w/2}" y="{plot_y+plot_h+35}" text-anchor="middle" class="small">potential</text>')
    legend_y = height - 23
    for i, metal in enumerate(("Dy", "Nd")):
        x = width / 2 - 60 + i * 120
        svg.append(f'<rect x="{x}" y="{legend_y-10}" width="12" height="12" fill="{colors[metal]}"/>')
        svg.append(f'<text x="{x+18}" y="{legend_y}" class="label">{metal}</text>')
    svg.append("</svg>")
    path.write_text("\n".join(svg))


def write_gr_figure(
    path: Path, datasets: dict[tuple[str, str], dict[str, np.ndarray]], pair: str, cuts: dict | None = None
) -> None:
    width, height = 1120, 470
    margin_l, margin_t, panel_w, panel_h = 72, 56, 450, 335
    gap = 70
    colors = ("#2878b5", "#d04a35", "#2a9d55", "#8a55a2")
    models = sorted({key[0] for key in datasets})
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:20px;font-weight:bold}.label{font-size:12px}.small{font-size:10px}.paneltitle{font-size:14px;font-weight:bold}</style>',
        f'<text x="{width/2}" y="27" text-anchor="middle" class="title">{path.stem.replace("_", " ")}</text>',
    ]
    for pi, metal in enumerate(("Dy", "Nd")):
        x0 = margin_l + pi * (panel_w + gap)
        y0 = margin_t
        plot_x, plot_y = x0 + 45, y0 + 24
        plot_w, plot_h = panel_w - 60, panel_h - 65
        curves = [(model, datasets[(model, metal)]) for model in models if (model, metal) in datasets]
        max_y = max(
            (
                float(np.max(curve[pair][np.asarray(curve.get("r", R)) < 5.5]))
                for _, curve in curves
            ),
            default=1.0,
        )
        max_y = max(max_y * 1.12, 1.0)
        svg.append(f'<text x="{x0+panel_w/2}" y="{y0+10}" text-anchor="middle" class="paneltitle">{metal}</text>')
        svg.append(f'<line x1="{plot_x}" y1="{plot_y+plot_h}" x2="{plot_x+plot_w}" y2="{plot_y+plot_h}" stroke="#444"/>')
        svg.append(f'<line x1="{plot_x}" y1="{plot_y}" x2="{plot_x}" y2="{plot_y+plot_h}" stroke="#444"/>')
        for tick in range(1, 5):
            yy = plot_y + plot_h * (1 - tick / 4)
            svg.append(f'<line x1="{plot_x}" y1="{yy:.1f}" x2="{plot_x+plot_w}" y2="{yy:.1f}" stroke="#ddd"/>')
            svg.append(f'<text x="{plot_x-5}" y="{yy+4:.1f}" text-anchor="end" class="small">{max_y*tick/4:.1f}</text>')
        for j, (model, curve) in enumerate(curves):
            pts = []
            radii = np.asarray(curve.get("r", R))
            for radius, value in zip(radii, curve[pair]):
                if radius > 5.5:
                    break
                x = plot_x + radius / 5.5 * plot_w
                y = plot_y + plot_h - min(float(value), max_y) / max_y * plot_h
                pts.append(f"{x:.1f},{y:.1f}")
            svg.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{colors[j % len(colors)]}" stroke-width="1.8"/>')
            svg.append(f'<line x1="{plot_x+8}" y1="{plot_y+18+j*16}" x2="{plot_x+28}" y2="{plot_y+18+j*16}" stroke="{colors[j % len(colors)]}" stroke-width="2"/>')
            svg.append(f'<text x="{plot_x+33}" y="{plot_y+22+j*16}" class="small">{model}</text>')
        if cuts:
            for j, model in enumerate(models):
                if (model, metal) not in cuts:
                    continue
                rc = cuts[(model, metal)].get("REE_O_total_min_A") if pair == "REE_O_total" else cuts[(model, metal)].get("OH_X_shell_min_A")
                if rc and math.isfinite(rc):
                    x = plot_x + rc / 5.5 * plot_w
                    svg.append(f'<line x1="{x:.1f}" y1="{plot_y}" x2="{x:.1f}" y2="{plot_y+plot_h}" stroke="{colors[j % len(colors)]}" stroke-dasharray="4,3" opacity="0.55"/>')
        svg.append(f'<text x="{plot_x+plot_w/2}" y="{plot_y+plot_h+35}" text-anchor="middle" class="label">r (Å)</text>')
    svg.append("</svg>")
    path.write_text("\n".join(svg))


def solution_summary_rows(solution_rows: list[dict]) -> list[dict]:
    result = []
    fields = {
        "REE_O_water": ("REE_Ow_mean", "REE_Ow_sem"),
        "REE_O_hydroxide": ("REE_OOH_mean", "REE_OOH_sem"),
        "REE_O_total": ("REE_Ototal_mean", "REE_Ototal_sem"),
        "OH_REE_first": ("OH_REE_first_from_gr", ""),
        "OH_Hwater_first": ("OH_Hwater_first_from_gr", ""),
        "OH_REE_shell": ("OH_REE_shell_mean", "OH_REE_shell_sem"),
        "OH_Hwater_shell": ("OH_Hw_shell_mean", "OH_Hw_shell_sem"),
        "OH_Xtotal_first": ("OH_Xtotal_mean", "OH_Xtotal_sem"),
        "OH_Xtotal_shell": ("OH_Xtotal_shell_mean", "OH_Xtotal_shell_sem"),
    }
    for row in solution_rows:
        if row["T_K"] != "300":
            continue
        item = {
            "model": row["model"], "metal": row["metal"], "variant": "box_pbc",
            "collapsed": row["collapsed"],
            "n_replicas": 1 if row["collapsed"] == "0" else 0,
            "rc_REE_O_A": row["rc_Ln_O_total_A"],
            "rc_OX_first_A": row["rc_OH_X_first_A"],
            "rc_OX_shell_A": row["rc_OH_X_shell_A"],
        }
        usable = row["collapsed"] == "0"
        for metric, (source, sem) in fields.items():
            value = row[source] if usable else None
            item[f"{metric}_mean"] = float(value) if value is not None and value != "" else math.nan
            error = row[sem] if (sem and usable) else None
            item[f"{metric}_sem"] = float(error) if error is not None and error != "" else math.nan
        result.append(item)
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    solution_rows, solution_gr = load_solution()
    write_solution_cutoff_clarity(solution_rows, solution_gr)

    jobs_doc = json.loads((GRAPHITE / "jobs.json").read_text())
    jobs = jobs_doc["jobs"]
    valid_jobs = []
    for job in jobs:
        base = GRAPHITE / "runs" / job["model"] / job["variant"] / job["metal"] / f"r{job['replica']:02d}"
        if (base / "done.json").is_file() and (base / "traj-denser-dumping-positions.npy").is_file():
            valid_jobs.append(job)
    if len(valid_jobs) != 800:
        raise RuntimeError(f"Expected 800 complete graphite trajectories; found {len(valid_jobs)}")
    screen = load_collapse_screen()
    excluded: dict[tuple[str, str, str], int] = defaultdict(int)
    analysed_jobs = []
    for job in valid_jobs:
        if screen_row(job, screen)["status"] == "excluded":
            excluded[(job["model"], job["metal"], job["variant"])] += 1
            excluded[(job["model"], job["metal"], "ALL")] += 1
        else:
            analysed_jobs.append(job)
    print(f"Collapse screen: analysing {len(analysed_jobs)}/{len(valid_jobs)} trajectories "
          f"(excluded {len(valid_jobs) - len(analysed_jobs)}; others truncated at usable_end_ps)", flush=True)

    cache_path = OUT / "graphite_rdf_accumulated.npz"
    cache_meta = OUT / "graphite_rdf_accumulated.json"
    rdf_stats = {}
    if "--reuse-rdf" in sys.argv:
        cache_doc = json.loads(cache_meta.read_text())
        if not isinstance(cache_doc, dict) or cache_doc.get("truncated_by") != SCREEN_CSV.name:
            raise RuntimeError("RDF cache predates the collapse screen; rerun without --reuse-rdf")
        cache_keys = cache_doc["keys"]
        with np.load(cache_path) as cache:
            for i, key in enumerate(cache_keys):
                rdf_stats[tuple(key)] = {
                    "volume": float(cache[f"volume_{i}"]),
                    "ree_hist": cache[f"ree_hist_{i}"].copy(),
                    "ree_exposure": cache[f"ree_exposure_{i}"].copy(),
                    "oh_hist": cache[f"oh_hist_{i}"].copy(),
                    "oh_exposure": cache[f"oh_exposure_{i}"].copy(),
                    "n_oh_centres": int(cache[f"n_oh_centres_{i}"]),
                    "n_frames": int(cache[f"n_frames_{i}"]),
                }
    else:
        print(f"Accumulating graphite RDFs from {len(analysed_jobs)} screened trajectories", flush=True)
        for i, job in enumerate(analysed_jobs, start=1):
            accumulate_rdf_for_job(job, rdf_stats, float(screen_row(job, screen)["usable_end_ps"]))
            if i % 40 == 0:
                print(f"  RDF trajectories: {i}/{len(analysed_jobs)}", flush=True)
        cache_keys = [list(key) for key in sorted(rdf_stats)]
        cache_data = {}
        for i, key in enumerate(sorted(rdf_stats)):
            item = rdf_stats[key]
            for name in ("volume", "ree_hist", "ree_exposure", "oh_hist", "oh_exposure",
                         "n_oh_centres", "n_frames"):
                cache_data[f"{name}_{i}"] = item[name]
        np.savez_compressed(cache_path, **cache_data)
        cache_meta.write_text(json.dumps({"keys": cache_keys, "truncated_by": SCREEN_CSV.name}, indent=2) + "\n")

    rdf_by_group = {key: make_rdf(value) for key, value in rdf_stats.items()}
    cutoff_by_group = {key: cutoff_data(rdf) for key, rdf in rdf_by_group.items()}
    cutoff_rows = []
    for (model, metal), cuts in sorted(cutoff_by_group.items()):
        cutoff_rows.append({
            "model": model, "metal": metal, "variant": "ALL",
            "n_trajectories": sum(1 for job in analysed_jobs if job["model"] == model and job["metal"] == metal),
            "n_excluded": excluded.get((model, metal, "ALL"), 0),
            "n_samples_total": rdf_stats[(model, metal)]["n_frames"],
            **cuts,
        })
    write_csv(OUT / "graphite_cutoffs_box_pbc.csv", cutoff_rows)
    write_rdf_csv(rdf_by_group, cutoff_by_group)
    if "--rdf-only" in sys.argv:
        print(f"RDF and cutoff tables written to {OUT}; skipped CN pass by request", flush=True)
        return
    for key, cuts in cutoff_by_group.items():
        if not all(math.isfinite(cuts[name]) for name in (
            "REE_O_total_min_A", "OH_X_first_min_A", "OH_X_shell_min_A",
        )):
            raise ValueError(f"No required RDF minimum for {key}; refusing to report CN.")

    replica_path = OUT / "graphite_cn_by_replica_box_pbc.csv"
    if "--reuse-cn" in sys.argv and replica_path.is_file():
        print("Reusing per-trajectory CN table", flush=True)
        with replica_path.open(newline="") as handle:
            replica_rows = list(csv.DictReader(handle))
        integer_columns = ("replica", "n_samples", "n_oh_centres")
        float_columns = ("usable_end_ps", "rc_REE_O_A", "rc_OX_first_A", "rc_OX_shell_A") + CN_NAMES
        for row in replica_rows:
            for name in integer_columns:
                row[name] = int(row[name])
            for name in float_columns:
                row[name] = float(row[name])
    else:
        print("Counting per-trajectory CN using pooled cutoffs", flush=True)
        replica_rows = []
        for i, job in enumerate(analysed_jobs, start=1):
            replica_rows.append(count_cn_for_job(
                job, cutoff_by_group[(job["model"], job["metal"])], screen_row(job, screen),
            ))
            if i % 40 == 0:
                print(f"  CN trajectories: {i}/{len(analysed_jobs)}", flush=True)
        write_csv(replica_path, replica_rows)
    graphite_summary = summarize_replicas(replica_rows, excluded)
    write_csv(OUT / "graphite_cn_summary_box_pbc.csv", graphite_summary)

    solution_summary = solution_summary_rows(solution_rows)
    combined_summary = graphite_summary + solution_summary
    write_csv(OUT / "cn_summary_300K.csv", combined_summary)

    # Export an aligned set of 300 K RDF curves for the two systems.
    for system, datasets in (("solution", solution_gr), ("graphite", rdf_by_group)):
        for pair in ("REE_O_total", "OH_X_total"):
            path = OUT / f"{system}_{pair}_box_pbc_300K.svg"
            write_gr_figure(
                path, datasets, pair,
                cutoff_by_group if system == "graphite" else None,
            )

    plot_metrics = [
        ("REE_O_water", "REE–O water", "REE_O_water_sem"),
        ("REE_O_hydroxide", "REE–O hydroxide", "REE_O_hydroxide_sem"),
        ("REE_O_total", "REE–O total", "REE_O_total_sem"),
        ("OH_REE_first", "OH–REE (first minimum)", "OH_REE_first_sem"),
        ("OH_Hwater_first", "OH–H water (first minimum)", "OH_Hwater_first_sem"),
        ("OH_REE_shell", "OH–REE (shell cutoff)", "OH_REE_shell_sem"),
        ("OH_Hwater_shell", "OH–H water (shell cutoff)", "OH_Hwater_shell_sem"),
        ("OH_Xtotal_first", "OH–X total (first minimum)", "OH_Xtotal_first_sem"),
        ("OH_Xtotal_shell", "OH–X total (shell minimum)", "OH_Xtotal_shell_sem"),
    ]
    model_order = ["uma-s-1.2.1", "DPA4C-Plus", "DPA4-Pro", "MACE-POLAR-1-M"]
    write_cn_figure(
        OUT / "graphite_CN_box_pbc_300K.svg", graphite_summary,
        ["uma-s-1.2.1", "DPA4C-Plus"], plot_metrics, "ALL",
    )
    write_cn_figure(
        OUT / "solution_CN_box_pbc_300K.svg", solution_summary,
        model_order, plot_metrics, "box_pbc",
    )
    print(f"Analysis outputs written to {OUT}", flush=True)


if __name__ == "__main__":
    main()
