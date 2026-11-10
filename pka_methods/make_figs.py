"""Step-by-step figures of the pKa1-3 protocol (Chinese labels).  Output: figs/*.png"""
import csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Circle

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent / "shiery_dft_repro"
ALL = ROOT / "all_ree"
OUT = HERE / "figs"
FONTS = HERE.parent / "genesis_hydration_figure" / "fonts"
for f in FONTS.glob("*.otf"):
    font_manager.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Noto Sans SC", "DejaVu Sans"], "axes.unicode_minus": False,
                     "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 200, "savefig.bbox": "tight"})
CONV = 2625.499639 / 0.529177210903  # Ha/bohr -> kJ/mol/A
RTLN10 = 8.314462618e-3 * 298.15 * np.log(10.0)
ORDER = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Sc", "Y"]
SUM = {r["metal"]: r for r in csv.DictReader(open(ALL / "results/ree_pKa_summary.csv"))}
SPEC = {r["species"]: r for r in csv.DictReader(open(ALL / "results/ree_cluster_species.csv"))}
SENS = list(csv.DictReader(open(ALL / "results/ree_pKa_sensitivity.csv")))
C1, C2, C3 = "#1f5fa8", "#d1631d", "#2a9d6f"
CPAP, CEXP = "#7a4fb0", "#555555"
EL_COL = {"H": "#f2f2f2", "O": "#d62828", "Cl": "#3fae49"}
EL_R = {"H": 0.28, "O": 0.5, "Cl": 0.7}


def fnum(x):
    try:
        v = float(x)
        return None if np.isnan(v) else v
    except (TypeError, ValueError):
        return None


def read_xyz(path, frame=-1):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0]); nf = len(lines) // (n + 2)
    blocks = [lines[i * (n + 2):(i + 1) * (n + 2)] for i in range(nf)]
    blocks = [b for b in blocks if all(len(l.split()) == 4 for l in b[2:])]
    b = blocks[frame]
    sym = [l.split()[0] for l in b[2:]]
    xyz = np.array([[float(v) for v in l.split()[1:]] for l in b[2:]])
    return sym, xyz


def draw_mol(ax, sym, xyz, metal, hl=(), scale=1.0, bond_cut=1.25, alpha_far=None, lw=1.6):
    """Ball-and-stick, projected on x-y, painter's order by z."""
    order = np.argsort(xyz[:, 2])
    z = xyz[:, 2]; zn = (z - z.min()) / (np.ptp(z) + 1e-9)
    for i in range(len(sym)):
        for j in range(i + 1, len(sym)):
            if {sym[i], sym[j]} == {"O", "H"} and np.linalg.norm(xyz[i] - xyz[j]) < bond_cut:
                a = 1.0 if alpha_far is None else 0.25 + 0.75 * min(zn[i], zn[j])
                ax.plot(*xyz[[i, j], :2].T, color="#888", lw=lw, zorder=2 + 10 * min(zn[i], zn[j]), alpha=a)
    for i in order:
        s = sym[i]
        if s == metal:
            col, r = "#1f5fa8", 0.8
        else:
            col, r = EL_COL.get(s, "#999"), EL_R.get(s, 0.5)
        a = 1.0 if alpha_far is None else 0.25 + 0.75 * zn[i]
        ec = "k"
        lwi = 0.6
        if i in hl:
            ec, lwi = "#f5b700", 2.4
        ax.add_patch(Circle(xyz[i, :2], r * scale, fc=col, ec=ec, lw=lwi, alpha=a, zorder=3 + 10 * zn[i]))
    ax.set_aspect("equal"); ax.axis("off")


# ---------------- Fig 0: workflow ----------------
def fig0():
    fig, ax = plt.subplots(figsize=(13, 5.6)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 5.6)

    def box(x, y, w, h, txt, fc, fs=10.5, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06,rounding_size=0.12", fc=fc, ec="#333", lw=1))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, weight="bold" if bold else None,
                linespacing=1.35)

    def arr(x0, y0, x1, y1):
        ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", lw=1.4, color="#333"))

    ax.text(0.1, 5.35, "pKa1：显式水 AIMD + Blue-Moon 约束 PMF（CP2K, PBE-D3）", fontsize=13, weight="bold", color=C1)
    y = 3.45; h = 1.5
    box(0.1, y, 2.3, h, "① 建盒子\nLn³⁺ + 3Cl⁻ + 64 H₂O\n12.42 Å 立方 PBC\nCN 9 (La–Sm) / 8 (Eu–Lu)", "#e3eefa")
    box(2.75, y, 2.3, h, "② 平衡\nGEO_OPT → NVT 298 K\n≥ 2 ps，检查 CN\n(Eu/Dy/Ho 推出一个水)", "#e3eefa")
    box(5.4, y, 2.45, h, "③ 9 个约束窗口\n第一层水 O–H 固定为\nr = 1.00, 1.05 … 1.40 Å\n各 8 ps，弃前 2 ps", "#e3eefa")
    box(8.2, y, 2.3, h, "④ 平均约束力\ndA/dr = −⟨λ⟩\n（SHAKE 乘子）", "#e3eefa")
    box(10.85, y, 2.05, h, "⑤ 积分 → pKa1\nA(1.40 Å)\npKa1 = ΔA / (RT ln10)", "#cfe0f7", bold=True)
    for x0, x1 in [(2.4, 2.75), (5.05, 5.4), (7.85, 8.2), (10.5, 10.85)]:
        arr(x0, y + h / 2, x1, y + h / 2)
    ax.text(0.1, 2.85, "pKa2、pKa3：第一层团簇 + C-PCM，以 pKa1 为锚逐级相加（PySCF）", fontsize=13, weight="bold", color=C2)
    y = 0.95
    box(0.1, y, 2.3, h, "⑥ 4 个团簇\nM(H₂O)ₓ³⁺, …OH²⁺,\n…(OH)₂⁺, …(OH)₃\n起始：Shiery SI 结构", "#fbe7d9")
    box(2.75, y, 2.3, h, "⑦ 结构 + 频率\nM06 / 4f-in-core ECP\n/ cc-pVTZ / C-PCM\n→ G_corr (qRRHO)", "#fbe7d9")
    box(5.4, y, 2.45, h, "⑧ 相对论单点\nM06 / sfX2C-1e / C-PCM\nSARC-DKH2(Ln)\nma-def2-TZVPP(O,H)\n4f 平均占据 → E", "#fbe7d9")
    box(8.2, y, 2.3, h, "⑨ 二阶差分\nΔΔG = G[n+1] − 2G[n]\n + G[n−1]\n（H⁺ 自由能抵消）", "#fbe7d9")
    box(10.85, y, 2.05, h, "⑩ 逐级相加\npKa(n+1) = pKa(n)\n + ΔΔG / (RT ln10)", "#f7d3ba", bold=True)
    for x0, x1 in [(2.4, 2.75), (5.05, 5.4), (7.85, 8.2), (10.5, 10.85)]:
        arr(x0, y + h / 2, x1, y + h / 2)
    ax.annotate("", (11.9, y + h + 0.05), (11.9, 3.45 - 0.05),
                arrowprops=dict(arrowstyle="-|>", lw=2, color=C1, ls="--"))
    ax.text(12.0, 2.65, "pKa1 作锚", color=C1, fontsize=10, va="center")
    ax.text(0.1, 0.3, "流程同 Shiery, Cooper & Cantu, Inorg. Chem. 2021, 60, 10257；本工作对 La–Lu 全部重算，Sc、Y 为扩展。"
            "RT ln10 = 5.708 kJ/mol（298.15 K），即 1 个 pKa 单位。", fontsize=9.5, color="#444")
    fig.savefig(OUT / "fig0_workflow.png"); plt.close(fig)


# ---------------- Fig 1: AIMD box + constrained O-H ----------------
def fig1():
    d = ROOT / "aimd/Nd/02_windows"
    don = json.loads((d / "donor.json").read_text())
    io, ih, ia = don["donor_O"] - 1, don["H"] - 1, don["acceptor_O"] - 1
    fig, axs = plt.subplots(1, 3, figsize=(14, 5.2), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
    sym, xyz = read_xyz(d / "r1.00/init.xyz", 0)
    L = 12.42
    m = sym.index("Nd")
    # wrap molecules around the box (keep waters whole)
    xyz = xyz.copy()
    for i, s in enumerate(sym):
        if s == "O":
            for j, t in enumerate(sym):
                if t == "H":
                    dv = xyz[j] - xyz[i]
                    if np.linalg.norm(dv - L * np.round(dv / L)) < 1.25:
                        xyz[j] = xyz[i] + dv - L * np.round(dv / L)
    ax = axs[0]
    draw_mol(ax, sym, xyz, "Nd", hl={io, ih, ia}, scale=0.75, alpha_far=True, lw=1.2)
    ax.add_patch(plt.Rectangle((0, 0), L, L, fill=False, ls="--", lw=1.2, ec="#333"))
    rO = np.linalg.norm(xyz[[i for i, s in enumerate(sym) if s == "O"]] - xyz[m], axis=1)
    ax.set_xlim(-1.5, L + 1.5); ax.set_ylim(-1.5, L + 1.5)
    ax.set_title("(a) Nd³⁺ + 3Cl⁻ + 64 H₂O，12.42 Å 周期盒\n（xy 投影，浅色 = 远离观察者）", fontsize=11)
    for k, (fr, title) in enumerate([(0, "(b) r(O–H) = 1.00 Å 窗口：起点"), (-1, "(c) r(O–H) = 1.40 Å 窗口：末帧")]):
        win = "r1.00" if k == 0 else "r1.40"
        s2, x2 = read_xyz(d / win / "window-pos-1.xyz", fr)
        mm = s2.index("Nd")
        dv = x2 - x2[mm]; dv -= L * np.round(dv / L); x2 = x2[mm] + dv
        # Nd, first shell O + their H, acceptor water
        Oidx = [i for i, s in enumerate(s2) if s == "O" and np.linalg.norm(x2[i] - x2[mm]) < 3.3]
        keep = set([mm] + Oidx + [ia])
        for i in list(keep):
            if s2[i] == "O":
                keep |= {j for j, t in enumerate(s2) if t == "H" and np.linalg.norm(x2[j] - x2[i]) < 1.25}
        keep |= {ih}
        keep = sorted(keep)
        sub = x2[keep]
        # rotate so that Nd -> donor O lies along +x
        v = x2[io] - x2[mm]; v /= np.linalg.norm(v)
        zt = np.cross(v, [0, 0, 1.0]); zt /= np.linalg.norm(zt); yt = np.cross(zt, v)
        R = np.array([v, yt, zt])
        sub = (sub - x2[mm]) @ R.T
        idx = {g: n for n, g in enumerate(keep)}
        draw_mol(axs[k + 1], [s2[i] for i in keep], sub, "Nd", hl={idx[io], idx[ih], idx[ia]}, scale=0.85)
        rOH = np.linalg.norm(x2[ih] - x2[io]); rHA = np.linalg.norm(x2[ih] - x2[ia]); rLnO = np.linalg.norm(x2[io] - x2[mm])
        p = sub[idx[ih]]
        axs[k + 1].text(0.5, 0.0, f"O–H = {rOH:.2f} Å　　H 到受体 O = {rHA:.2f} Å\nNd–O(供体) = {rLnO:.2f} Å",
                        transform=axs[k + 1].transAxes, ha="center", va="bottom", fontsize=9.5,
                        bbox=dict(fc="white", ec="#f5b700"))
        axs[k + 1].set_title(title, fontsize=11)
        axs[k + 1].autoscale_view()
        axs[k + 1].set_xlim(sub[:, 0].min() - 1.2, sub[:, 0].max() + 1.2); axs[k + 1].set_ylim(sub[:, 1].min() - 3.2, sub[:, 1].max() + 1.0)
    fig.text(0.5, -0.02, "黄框：被约束的第一层水 O（供体）、它的 H 和第二层受体 O（氢键最短者）。蓝 = Nd，红 = O，白 = H，绿 = Cl。"
             "(b)(c) 只画 Nd、第一层水和受体水，Nd→供体 O 方向朝右。", ha="center", fontsize=9.5, color="#444")
    fig.savefig(OUT / "fig1_aimd_box.png"); plt.close(fig)


# ---------------- Fig 2: per-window constraint force vs time ----------------
def read_lambda(path):
    out = []
    for l in open(path):
        w = l.split()
        if l.startswith("Shake") and len(w) == 4:
            try:
                out.append(float(w[-1]))
            except ValueError:
                pass
    return np.array(out)


def fig2():
    d = ROOT / "aimd/Nd/02_windows"
    rs = [1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40]
    pmf = {round(float(r["r"]), 2): r for r in csv.DictReader(open(ALL / "results/pmf/pmf_Nd.csv"))}
    fig, axs = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw=dict(width_ratios=[1.6, 1]))
    cmap = plt.get_cmap("viridis")
    for k, r in enumerate(rs):
        f = -read_lambda(d / f"r{r:.2f}/window-1.LagrangeMultLog") * CONV
        t = np.arange(len(f)) * 1e-3
        nb = 250
        nblk = len(f) // nb
        fb = f[:nblk * nb].reshape(nblk, nb).mean(1); tb = (np.arange(nblk) + 0.5) * nb * 1e-3
        axs[0].plot(tb, fb, color=cmap(k / 8), lw=1.1, label=f"{r:.2f} Å")
        axs[0].hlines(float(pmf[r]["force"]), 2, t[-1], color=cmap(k / 8), lw=2.2, ls="--")
    axs[0].axvspan(0, 2, color="#ddd", alpha=0.6)
    axs[0].text(1.0, 0.97, "弃去\n(平衡)", transform=axs[0].get_xaxis_transform(), ha="center", va="top", fontsize=9.5)
    axs[0].set_xlabel("模拟时间 t (ps)"); axs[0].set_ylabel("约束力 −λ (kJ mol⁻¹ Å⁻¹)")
    axs[0].set_title("(a) Nd：9 个窗口的约束力（0.25 ps 块平均；虚线 = 2 ps 后平均值）", fontsize=11)
    axs[0].legend(title="r(O–H)", fontsize=8.5, ncol=3, loc="lower right", frameon=False)
    # quarters per window
    ax = axs[1]
    for k, r in enumerate(rs):
        q = [float(pmf[r][f"force_q{i}"]) for i in range(1, 5)]
        ax.plot([r + (i - 1.5) * 0.007 for i in range(4)], q, "o", ms=4, color=cmap(k / 8))
        ax.errorbar(r, float(pmf[r]["force"]), yerr=float(pmf[r]["force_sem"]), fmt="s", ms=7, color="k", mfc="white", capsize=3, zorder=5)
    ax.set_xlabel("r(O–H) (Å)"); ax.set_ylabel("平均约束力 (kJ mol⁻¹ Å⁻¹)")
    ax.set_title("(b) 每窗口的平均力：方块 ± SEM；\n彩色点 = 生产段四等分各自的平均（看收敛）", fontsize=11)
    n_unb = int(pmf[1.0]["n_donor_unbound"]); n_prod = int(pmf[1.0]["n_prod"])
    ax.text(0.98, 0.70, f"注：Nd 1.00 Å 窗口中被约束的水\n有 {n_unb}/{n_prod} 步离开第一层；\n平均力只用仍在第一层的步。",
            transform=ax.transAxes, fontsize=8.8, color="#a33", ha="right", va="top")
    fig.savefig(OUT / "fig2_windows_force.png"); plt.close(fig)


# ---------------- Fig 3: PMF integration ----------------
def fig3():
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    for el, col in [("Nd", C1), ("Dy", C2)]:
        rows = list(csv.DictReader(open(ALL / f"results/pmf/pmf_{el}.csv")))
        r = np.array([float(x["r"]) for x in rows]); F = np.array([float(x["force"]) for x in rows])
        se = np.array([float(x["force_sem"]) for x in rows]); A = np.array([float(x["A"]) for x in rows])
        lab = f"{el}（盒子 CN {SUM[el]['box_cn']}）"
        axs[0].errorbar(r, F, yerr=se, fmt="o-", color=col, capsize=3, label=lab)
        axs[0].fill_between(r, 0, F, color=col, alpha=0.08)
        axs[1].plot(r, A, "o-", color=col, label=lab)
        pk = A[-1] / RTLN10
        axs[1].annotate(f"ΔA = {A[-1]:.1f} kJ/mol\npKa1 = {A[-1]:.1f} / 5.708 = {pk:.2f}",
                        (r[-1], A[-1]), (r[-1] - 0.2, A[-1] + 6) if el == "Dy" else (r[-1] - 0.12, A[-1] - 30), fontsize=10, color=col,
                        arrowprops=dict(arrowstyle="->", color=col))
    axs[0].axhline(0, color="#999", lw=0.8)
    axs[0].set_xlabel("r(O–H) (Å)"); axs[0].set_ylabel("dA/dr = −⟨λ⟩ (kJ mol⁻¹ Å⁻¹)")
    axs[0].set_title("(a) 平均力 ± SEM（阴影面积 = 做功）", fontsize=11); axs[0].legend(frameon=False)
    axs[1].set_xlabel("r(O–H) (Å)"); axs[1].set_ylabel("A(r) − A(1.00 Å) (kJ/mol)")
    axs[1].set_title("(b) 梯形积分得到 PMF；端点 1.40 Å 的 ΔA 换算为 pKa1", fontsize=11); axs[1].legend(frameon=False, loc="upper left")
    axs[1].set_ylim(-3, 62)
    fig.text(0.5, -0.04, "A(r) = ∫₁.₀₀ʳ ⟨dA/dr'⟩ dr'；pKa1 = A(1.40 Å) / (RT ln10)，RT ln10 = 5.708 kJ/mol。"
             "误差（论文同法）：生产段前半与后半分别积分得到的 pKa1 之差。", ha="center", fontsize=9.5, color="#444")
    fig.savefig(OUT / "fig3_pmf.png"); plt.close(fig)


def exp_range(s):
    if not s:
        return None
    a, b = s.rsplit("-", 1) if s.count("-") == 1 else s.split("-", 1)
    return float(a), float(b)


# ---------------- Fig 4: pKa1 all elements ----------------
def fig4():
    fig, ax = plt.subplots(figsize=(13, 4.8))
    x = np.arange(len(ORDER))
    for i, el in enumerate(ORDER):
        s = SUM[el]
        er = exp_range(s["exp_range1"])
        if er:
            ax.plot([i, i], er, color="#bbb", lw=7, solid_capstyle="round", zorder=1)
            ax.plot(i, float(s["exp_median1"]), "_", color=CEXP, ms=16, mew=2, zorder=2)
        if fnum(s["paper_pKa1"]):
            ax.plot(i + 0.18, float(s["paper_pKa1"]), "D", color=CPAP, ms=6, zorder=3)
        ax.errorbar(i - 0.18, float(s["pKa1"]), yerr=float(s["pKa1_err"]), fmt="o", color=C1, ms=7, capsize=3, zorder=4)
    ax.plot([], [], color="#bbb", lw=7, label="实验 / 汇编值范围（Shiery Table 2）"); ax.plot([], [], "_", color=CEXP, ms=16, mew=2, label="实验中位数")
    ax.plot([], [], "D", color=CPAP, label="Shiery 2021 计算"); ax.errorbar([], [], yerr=[], fmt="o", color=C1, label="本工作 AIMD-PMF ± 半段差")
    ax.axvspan(14.5, 16.5, color="#f4f4f4", zorder=0); ax.text(15.5, 3.4, "扩展\n(论文无)", ha="center", fontsize=9, color="#777")
    ax.set_xticks(x, ORDER); ax.set_ylabel("pKa1"); ax.set_ylim(3, 11)
    ax.legend(ncol=4, frameon=False, fontsize=9.5, loc="upper center", bbox_to_anchor=(0.5, 1.12))
    m = [abs(float(SUM[e]["pKa1"]) - float(SUM[e]["paper_pKa1"])) for e in ORDER[:15]]
    ax.set_title(f"pKa1 全元素结果（与 Shiery 平均绝对差 {np.mean(m):.2f}；元素 MAD vs 实验 1.05，论文 1.15）", fontsize=11, pad=40)
    fig.savefig(OUT / "fig4_pka1_all.png"); plt.close(fig)


# ---------------- Fig 5: Nd clusters ----------------
def fig5():
    names = [("aq", "Nd(H₂O)₉³⁺"), ("OH1", "Nd(H₂O)₈(OH)²⁺"), ("OH2", "Nd(H₂O)₇(OH)₂⁺"), ("OH3", "Nd(H₂O)₆(OH)₃")]
    fig, axs = plt.subplots(1, 4, figsize=(15, 4.0))
    for ax, (k, lab) in zip(axs, names):
        sym, xyz = read_xyz(ROOT / f"cluster/lc/Nd_{k}/opt_final.xyz", 0)
        m = sym.index("Nd"); xyz = xyz - xyz[m]
        Hs = [i for i, s in enumerate(sym) if s == "H"]
        ohO = [i for i, s in enumerate(sym) if s == "O" and sum(np.linalg.norm(xyz[j] - xyz[i]) < 1.25 for j in Hs) == 1]
        # PCA view
        xyz = best_view(xyz, sym)
        draw_mol(ax, sym, xyz, "Nd", hl=set(ohO), scale=0.9)
        sp = SPEC[f"Nd_{k}"]
        ax.set_xlim(-4.3, 4.3); ax.set_ylim(-4.3, 4.3)
        ax.set_title(f"{lab}\n电荷 {sp['charge']}，CN {sp['cn_final']}，虚频 {sp['n_imag']} 个", fontsize=11)
    fig.text(0.5, 0.02, "M06 / C-PCM(水) 优化结构（Nd，来自 Shiery SI 起始结构）。黄框 = OH⁻ 的 O。每一级：一个配位水失去一个 H⁺，变成 OH⁻，总电荷 −1。",
             ha="center", fontsize=9.5, color="#444")
    fig.savefig(OUT / "fig5_clusters_nd.png"); plt.close(fig)


# ---------------- Fig 6: free-energy staircase and decomposition ----------------
def best_view(xyz, sym):
    """Rotation whose xy projection keeps heavy atoms apart (least visual overlap)."""
    rng = np.random.default_rng(0)
    heavy = [i for i, t in enumerate(sym) if t != "H"]
    best, bs = None, -1
    for _ in range(3000):
        q = rng.normal(size=4); q /= np.linalg.norm(q); a, b, c, d = q
        R = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])
        p = (xyz @ R.T)[heavy, :2]
        dd = np.linalg.norm(p[:, None] - p[None], axis=-1)[np.triu_indices(len(heavy), 1)]
        sc = np.sort(dd)[:3].sum()
        if sc > bs:
            bs, best = sc, R
    return xyz @ best.T


def species_set(el):
    return SUM[el]["cluster_set"]


def ddg_parts(el):
    st = species_set(el)
    g = [float(SPEC[f"{st}_{k}"]["Gcorr_qrrho_inv_kJ"]) for k in ("aq", "OH1", "OH2", "OH3")]
    dd21 = float(SUM[el]["ddG21_kJ"]); dd32 = float(SUM[el]["ddG32_kJ"])
    gc21 = g[2] - 2 * g[1] + g[0]; gc32 = g[3] - 2 * g[2] + g[1]
    return dict(dd21=dd21, dd32=dd32, gc21=gc21, gc32=gc32, de21=dd21 - gc21, de32=dd32 - gc32, g=g)


def fig6():
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.2), gridspec_kw=dict(width_ratios=[1, 1.35]))
    ax = axs[0]
    for el, col, off in [("Nd", C1, -0.17), ("Dy", C2, 0.17)]:
        s = SUM[el]; dg = [0, float(s["dG1_kJ"]), float(s["dG2_kJ"]), float(s["dG3_kJ"])]
        for n in range(1, 4):
            ax.bar(n + off, dg[n], width=0.3, color=col, alpha=0.85 if n == 1 else 0.45, edgecolor=col,
                   label=el if n == 1 else None)
            pk = float(s[f"pKa{n}"])
            ax.text(n + off, dg[n] + (1.5 if el == "Nd" else 7.5), f"{dg[n]:.1f}\npKa {pk:.2f}", ha="center", fontsize=8.6, color=col)
        p = ddg_parts(el)
        ax.annotate("", (2 + off, dg[2]), (2 + off, dg[1]), arrowprops=dict(arrowstyle="->", color="k", lw=1))
        ax.annotate("", (3 + off, dg[3]), (3 + off, dg[2]), arrowprops=dict(arrowstyle="->", color="k", lw=1))
    ax.set_xticks([1, 2, 3], ["第 1 步\nAIMD-PMF", "第 2 步\n= 第 1 步 + ΔΔG₂₁", "第 3 步\n= 第 2 步 + ΔΔG₃₂"])
    ax.set_ylabel("逐步脱质子自由能 ΔGₙ (kJ/mol)"); ax.set_ylim(0, 118)
    ax.set_title("(a) 爬楼梯：深色 = AIMD 直接算；浅色 = 团簇 ΔΔG 加上去", fontsize=11); ax.legend(frameon=False, loc="upper left")
    ax = axs[1]
    els = ORDER
    x = np.arange(len(els)); w = 0.38
    P = [ddg_parts(e) for e in els]
    for j, (key, lab, col, dx) in enumerate([("21", "ΔΔG₂₁ → pKa2", C2, -w / 2), ("32", "ΔΔG₃₂ → pKa3", C3, w / 2)]):
        de = [p["de" + key] for p in P]; gc = [p["gc" + key] for p in P]
        ax.bar(x + dx, de, w, color=col, alpha=0.75, label=f"{lab}：电子能 ΔΔE (X2C)")
        ax.bar(x + dx, gc, w, bottom=de, color="none", edgecolor=col, hatch="///", lw=0.8, label=f"{lab}：热校正 ΔΔG_corr")
        ax.plot(x + dx, [p["dd" + key] for p in P], "k_", ms=10, mew=2)
    ax.plot([], [], "k_", ms=10, mew=2, label="合计 ΔΔG")
    ax.axhline(0, color="#999", lw=0.8)
    ax.set_xticks(x, els); ax.set_ylabel("ΔΔG (kJ/mol)；5.708 kJ/mol = 1 pKa")
    ax.set_title("(b) 全元素 ΔΔG 拆分：G = E(X2C 单点) + G_corr(M06 频率)", fontsize=11)
    ax.legend(fontsize=8.3, ncol=2, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.0))
    ax.set_ylim(min(0, min(min(p["de21"], p["de32"], p["dd21"], p["dd32"]) for p in P)) - 2, 72)
    fig.text(0.5, -0.03, "ΔΔG(n+1, n) = G[n+1] − 2G[n] + G[n−1]，[n] = 含 n 个 OH⁻ 的团簇。H⁺ 与溶剂化 H⁺ 的自由能在相邻两步相减时抵消。"
             "Ce/Sm/Eu/Tm/Yb 用 4f 冻结模型，Sc 用 CN 7 团簇。", ha="center", fontsize=9.5, color="#444")
    fig.savefig(OUT / "fig6_ddg.png"); plt.close(fig)


# ---------------- Fig 7: sensitivity ----------------
def fig7():
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), sharey=False)
    for ax, n, col in [(axs[0], 2, C2), (axs[1], 3, C3)]:
        for i, el in enumerate(ORDER):
            rows = [r for r in SENS if r["metal"] == el and r["pKa1_anchor"] == "PMF"]
            vals = {(r["set_kind"], r["G_corr"]): float(r[f"pKa{n}"]) for r in rows}
            allv = list(vals.values())
            ax.plot([i, i], [min(allv), max(allv)], color=col, lw=6, alpha=0.25, solid_capstyle="round")
            mk = {"main": "o", "dropimag": "^", "rrho": "s"}
            for (kind, gc), v in vals.items():
                ax.plot(i + (0.12 if kind == "frozen" else -0.12), v, mk[gc], ms=4.5,
                        mfc="white" if kind == "frozen" else col, mec=col)
            ax.plot(i, float(SUM[el][f"pKa{n}"]), "k*", ms=11, zorder=5)
        ax.set_xticks(range(len(ORDER)), ORDER); ax.set_ylabel(f"pKa{n}")
        ax.set_title(f"({'ab'[n - 2]}) pKa{n} 对处理方式的敏感性", fontsize=11)
    h = [plt.Line2D([], [], marker="*", color="k", ls="", ms=11, label="主结果"),
         plt.Line2D([], [], marker="o", color="#666", ls="", label="虚频按 |ν| 计（主）"),
         plt.Line2D([], [], marker="^", color="#666", ls="", label="去掉虚频"),
         plt.Line2D([], [], marker="s", color="#666", ls="", label="纯谐振 RRHO"),
         plt.Line2D([], [], marker="o", mfc="white", mec="#666", ls="", label="空心 = 4f 冻结模型"),
         plt.Line2D([], [], color="#999", lw=6, alpha=0.4, label="所有变体的范围")]
    fig.legend(handles=h, ncol=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.06), fontsize=9.5)
    fig.text(0.5, -0.03, "同一 pKa1 锚（AIMD-PMF），只改变团簇热校正（虚频处理、qRRHO/RRHO）和 4f 电子模型。范围常达 2–3 个 pKa 单位：pKa2/3 是半定量结果。",
             ha="center", fontsize=9.5, color="#444")
    fig.savefig(OUT / "fig7_sensitivity.png"); plt.close(fig)


# ---------------- Fig 8: final pKa1-3 ----------------
def fig8():
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, n, col in zip(axs, (1, 2, 3), (C1, C2, C3)):
        for i, el in enumerate(ORDER):
            s = SUM[el]
            er = exp_range(s[f"exp_range{n}"])
            if er:
                ax.plot([i, i], [er[0] - (0.08 if er[0] == er[1] else 0), er[1] + (0.08 if er[0] == er[1] else 0)],
                        color="#bbb", lw=6, solid_capstyle="round", zorder=1)
            if fnum(s[f"paper_pKa{n}"]):
                ax.plot(i + 0.15, float(s[f"paper_pKa{n}"]), "D", color=CPAP, ms=5, zorder=3)
            ax.plot(i - 0.15, float(s[f"pKa{n}"]), "o", color=col, ms=6.5, zorder=4)
        ax.set_xticks(range(len(ORDER)), ORDER, fontsize=8.5)
        ax.set_title(f"pKa{n}：" + ("AIMD-PMF" if n == 1 else f"pKa{n - 1} + ΔΔG（团簇）"), fontsize=11)
    axs[0].set_ylabel("pKa")
    h = [plt.Line2D([], [], marker="o", color="#333", ls="", label="本工作"),
         plt.Line2D([], [], marker="D", color=CPAP, ls="", label="Shiery 2021 计算"),
         plt.Line2D([], [], color="#bbb", lw=6, label="文献汇编值范围（Shiery Table 2）")]
    fig.legend(handles=h, ncol=3, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.05))
    fig.text(0.5, -0.04, "注意：Table 2 中 pKa2/pKa3 的许多数值在数量级上像累积常数 −log β*ₙ（如 pKa3 ≈ 21–27），而计算值是逐步常数，"
             "两者不能直接比较；pKa1 才是可比的基准。", ha="center", fontsize=9.5, color="#a33")
    fig.savefig(OUT / "fig8_final.png"); plt.close(fig)


if __name__ == "__main__":
    for f in (fig0, fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8):
        f(); print("done", f.__name__)
