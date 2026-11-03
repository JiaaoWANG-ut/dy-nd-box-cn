# Dy and Nd hydrolysis: box-PBC coordination numbers, three-step apparent activation energies, and comparison with Shiery 2021

**English** | [中文](README_zh.md)

- Web pages: https://jiaaowang-ut.github.io/dy-nd-box-cn/ (English, homepage) and https://jiaaowang-ut.github.io/dy-nd-box-cn/index_zh.html (Chinese)
- Bundle download (GitHub release): https://github.com/JiaaoWANG-ut/dy-nd-box-cn/releases/latest
- Offline: open `index.html` (English) or `index_zh.html` (Chinese), or unzip `dy_nd_box_cn_data.zip` and open `dy_nd_box_cn/index.html`. The English / 中文 buttons at the top right of each page switch the language. All links are relative, so the pages also work offline.

Every table and figure on the page has download links to its raw data underneath, and section 8 lists all files. Every number on the page comes from the CSV/JSON files in this directory; the page itself only takes ranges, differences and similar simple summaries.

## Contents

1. Task (a), REE-centred: REE–O(water), REE–O(OH), REE–O(total). The cutoff is the minimum after the first shell of the total REE–O g(r).
2. Task (b), OH-centred: O–REE, O–H(water), O–X(total). Two treatments, a strict cutoff and a shell cutoff; see below.
3. Three-step apparent activation energies: committed transition rates first, then the Eyring ΔG‡ (300 K) and the Arrhenius Ea (solution, 300–500 K).
4. Stepwise pKa equivalents and ΔG from the occupancies, compared with Shiery, Cooper, Cantu, *Inorg. Chem.* 2021, 60, 10257: Table 1, Table 2, SI Tables S1/S2, the SI optimized structures and the SI first-step PMFs.

Only three-dimensional periodic box (box_pbc) trajectories are used. The system is Ln(OH)3 + water; the graphite system also contains one graphite layer.

- Graphite: DPA4C-Plus, UMA-s-1.2.1; 300 K; 4 OH starting configurations (3near / 2near1far / 1near2far / 3far) × 50 trajectories × 150 ps.
- Solution: DPA4-Pro, DPA4C-Plus, MACE-POLAR-1-M, UMA-s-1.2.1; 300/350/400/450/500 K; 1 trajectory × 150 ps per temperature.

## Run order

Run in this directory, in this order. The numeric argument 8 is the number of parallel processes.

```
module load python/3.12-26.1.0
python kinetics.py extract 8    # raw trajectories -> kinetics_cache/ (20 fs time series, ~372 MB)
python kinetics.py screen       # collapse screening -> collapse_screen.csv
python analyze.py               # tasks (a)(b): g(r), cutoffs, CN, figures
python kinetics.py analyze 8    # rates, ΔG‡, Arrhenius, occupancies, pKa equivalents -> kinetics/
python kinetics.py figures      # kinetics figures -> kinetics/*.svg
python literature.py            # paper tables, SI structures and PMFs -> literature/ (downloads the SI PDF automatically)
python compare.py 8             # comparison with the paper -> comparison/
python build_page.py            # pages index.html (English) and index_zh.html (Chinese) + dy_nd_box_cn_data.zip; checks links and untranslated notes
```

Options of `analyze.py`:

- `--reuse-rdf`: reuse `graphite_rdf_accumulated.npz` and skip the graphite RDF accumulation.
- `--rdf-only`: write only g(r) and the cutoffs.
- `--reuse-cn`: reuse the existing per-trajectory CN tables.

Raw data (paths relative to `genesis/`):

- Graphite: `Graphite-dy-nd-50-denser-dumping/`
- Solution trajectories:
  - `md/ln-150ps-denser-dumping/`: DPA4C-Plus and UMA, one frame every 2 fs;
  - `md/ln-150ps-every-step-dumping/`: DPA4-Pro and MACE, one frame every 1 fs.
- Solution CN and g(r): already computed from the solution trajectories above by `md/dense-gr-cn/gr_cn.py` (`cn_table.csv`, `results/*_gr.npz`, `results/*.json`); `analyze.py` reads them directly.

## Methods

### Collapse screening (`kinetics.py screen` → `collapse_screen.csv`)

Some UMA trajectories collapse: the potential energy drops abruptly and H/O are pulled onto Ln. A trajectory is flagged if either criterion holds:

- Energy: a 1 ps block average falls 25 eV (graphite) or 15 eV (solution) below the lowest of the previous 5 blocks; the onset is refined at −10 eV.
- Geometry: Ln–O < 1.85 Å or Ln–H < 1.90 Å for 10 consecutive samples.

Trajectories are cut 1 ps before the collapse onset. If ≥ 20 ps remain, the trajectory is marked "truncated" and still used; otherwise it is "excluded".

- 91 trajectories are excluded in total.
- Graphite UMA: Dy 79 excluded, 116 truncated; Nd 7 excluded, 116 truncated.
- Solution: all 5 UMA Dy trajectories are excluded; UMA Nd 500 K is truncated.
- DPA4-Pro, DPA4C-Plus and MACE never collapse.

Graphite CN/g(r) and all kinetics use only the usable time windows. Solution CN/g(r) are taken from `md/dense-gr-cn`, which gives no CN for collapsed trajectories, so UMA Nd 500 K in solution has kinetics results only (from the first 66 ps).

### Coordination numbers (`analyze.py`)

- **O species:** in every frame each H is assigned to its nearest O (minimum-image distance). An O with 2 H is water, with 1 H it is OH⁻, and with 0 or ≥ 3 H it is "other" (O²⁻/H₃O⁺). After a proton transfer, the OH and water identities swap accordingly.
- **REE–O(other):** 0 on graphite; |total − water − OH| ≤ 0.002 in the non-collapsed solution trajectories.
- **(a) REE-centred:** r_c is the minimum after the first shell of the total REE–O g(r). g(r) grid 0.02 Å, 5-point moving average.
  - Procedure: first find the main peak (the highest point within 1.6–3.4 Å on graphite, 1.8–3.0 Å in solution), then take the lowest point within 1.6 Å (solution: 1.5 Å) after the main peak.
  - This skips the small dip between the OH (≈ 2.25 Å) and water (≈ 2.5 Å) sub-peaks.
  - Water, OH and total are all counted with this one r_c.
- **(b) OH-centred:** the centres are the O atoms of all OH in each frame. X is every atom except the central O and its own H: REE, water O/H, O/H of the other OH, "other" O/H, and C in the graphite system.
- **Sampling and errors:**
  - Graphite: one frame every 1 ps from t ≥ 10 ps. Errors are standard errors across trajectories: REE-centred quantities are weighted by frame count, OH-centred quantities by the number of OH centres. The `ALL` row combines the 4 starting configurations.
  - Solution: per-frame counts from `md/dense-gr-cn` (t ≥ 10 ps, every 2 fs); errors are 5-block standard errors. The exception is O–REE and O–H(water) with the strict cutoff: they are integrated from partial g(r), have no error bars, and differ from direct counts by at most 0.05.

### When the OH solvation shell is unclear

The O–X(total) g(r) is made of several overlapping peaks: hydrogen-bonded water H (1.6–1.8 Å; DPA4C-Plus also has a shared proton at 1.23–1.25 Å), REE (2.2–2.3 Å) and water O (2.4–2.7 Å). The "first minimum of the total g(r)" is therefore not unique. Treatment:

1. **Strict cutoff:** the first minimum after the first real peak of O–X(total) (g > 0.5, followed by a drop of ≥ 10%). It always lies before the O–REE peak, so O–REE ≈ 0.
2. **Shell cutoff:** start 0.1 Å below the O–REE peak position (2.25 Å if the O–REE peak has g ≤ 0.5), climb along O–X(total) to the local maximum, then descend to the next minimum (at most 3.0 Å), so that the REE is included.
3. **Clarity:** a peak-to-minimum drop 1 − g_min/g_peak ≥ 10% counts as clear (`*_min_clear`). Only the shell minimum of graphite DPA4C-Plus Nd is unclear.
4. **Auxiliary checks:**
   - Remainder: inside the shell, X − REE − H(water) should be ≈ 0.
   - Reciprocal reference: REE–O(OH)/3 does not depend on the OH-centred cutoff and can be compared with the shell O–REE.

The strict and shell cutoffs for solution come from `md/dense-gr-cn`, which uses the same rules; `analyze.py` adds the clarity flag with the same 10% criterion.

### Kinetics (`kinetics.py analyze`)

- **State:** n(t) = number of OH⁻ O atoms within r_c of Ln, capped at 3, one frame every 20 fs. r_c is the minimum after the first peak of the Ln–O(OH) radial number density (3.03–3.47 Å, see `kinetics_cutoffs.csv`).
- **Mapping to the paper's steps:**
  - MD forward n−1 → n (Ln gains an OH⁻) corresponds to the paper's step-n acid dissociation Ln(OH)ₙ₋₁ + H₂O → Ln(OH)ₙ + H⁺.
  - MD reverse n → n−1 (Ln loses an OH⁻) is its reverse.
- **Committed transitions:** a centred mode filter of width W = 0.5 ps, after which segments shorter than W/2 are absorbed. Rates are counted from t ≥ 10 ps.
- **Rates and barriers:**
  - k = N / T_origin, ΔG‡ = RT ln(k_B T / h k).
  - The 95% CI is the envelope of a bootstrap (2000 resamples) and the Garwood Poisson interval. Graphite is resampled by trajectory; solution, with only 1 trajectory per temperature, is resampled in 10 ps blocks.
  - For N = 0, the Poisson upper limit 3.689 / T_origin gives a lower bound on ΔG‡.
- **Arrhenius/Eyring:** fits weighted by N. Fit quality:
  - good: ≥ 4 temperatures, ≥ 20 events, and ≥ 3 events at every temperature;
  - fair: ≥ 3 temperatures and ≥ 8 events;
  - poor: everything else;
  - none: fewer than 2 temperatures with events.
- **Sensitivity:** W = 0.2 / 1.0 / 2.0 ps; r_c = 2.90 / 3.30 Å; start at t ≥ 1 / 50 / 90 ps. See the `scenario` column of `kinetics_rates.csv`.

### Occupancies and pKa equivalents

- Occupancies Pₙ are counted from t ≥ 10 ps (solution) or t ≥ 50 ps (graphite).
- Dissociation constant: K_d = (Pₙ₋₁/Pₙ)·(4 − n)·c₁, where c₁ is the concentration of one molecule in the box: 0.4224 M in solution, 0.3594 M on graphite.
- pKaₙ(equiv.) = pK_w + log₁₀ K_d with pK_w = 14, treated as an ideal solution without activity corrections.
- Status classes (`status` column of `compare_pKa_stepwise.csv`; the display on the page is given in parentheses):
  - `measured` (bold value): the committed dwells in both states are ≥ 1 ps;
  - `transient` ("transient"): both states occur, but one of them has committed dwells < 1 ps, so the value is an order of magnitude only;
  - `upper_bound` / `lower_bound` ("≤ / ≥ … bound"): only one state occurs. Assuming the unseen state is not more populated than the seen one, the pKa at Pₙ₋₁ = Pₙ is given as a bound;
  - `undetermined` ("undetermined"): neither state occurs.
- The graphite DPA4C-Plus occupancies remember the starting configuration and converge only after about 90 ps; see the `start90` setting.

### Comparison with the paper (`literature.py`, `compare.py`)

- **Methods of the paper** (from chapter 5 of the first author's PhD dissertation, which is this paper):
  - pKa₁: constrained AIMD gives a PMF, and its end point is used.
    - Setup: Bluemoon, CP2K/PBE, 64 waters + 3 Cl⁻.
    - 9 r(O–H) constraint points, 15 ps NVT each; the constraint force is averaged over 5–15 ps.
  - pKa₂, pKa₃: reaction free energies of M06 clusters (explicit first-shell water and OH) + CPCM implicit solvent (ORCA), taken relative to pKa₁.
    - The cluster coordination number is fixed: 9 for La–Sm and 8 for Eu–Lu, i.e. 9 for Nd and 8 for Dy.
  - The paper gives only 298.15 K thermodynamics, with no barriers or activation energies. The first-step PMF rises monotonically and has no transition state.
- **Tables 1/2:** transcribed from dissertation Tables 5-1/5-2, which are the same as Tables 1/2 of the article. Check: recomputing the Table 1 MAD from Table 2 gives 38/39 agreeing entries (±0.05); only Yb step 2 disagrees, and another 6 entries have no experimental value.
- **SI:** Tables S1/S2, the 60 M06-optimized structures, and the PMFs of Figures S1–S15 (read from the PDF vector paths).
- **Unit conversion:** ΔG = RT ln10 · pKa.
- **How to read the step-2 and step-3 experimental values:**
  - For Nd, Kragten 1984 (Table 2 ref 99) lists 8.1/16.2/24.3 = 8.1 × n, which looks like cumulative constants −log *βₙ rather than stepwise constants.
  - The other experimental values for steps 2 and 3 are also numerically closer to cumulative values.
  - The reading decides the conclusion (predicted ⟨n⟩ ≈ 1 or 3). `compare_pKa_cumulative.csv` gives the comparison read as cumulative constants.

## Output files

- This directory:
  - Pages: `index.html` (English, homepage), `index_zh.html` (Chinese); documentation: `README.md` (English), `README_zh.md` (Chinese)
  - CN summary: `cn_summary_300K.csv`
  - Graphite: `graphite_cn_summary_box_pbc.csv`, `graphite_cn_by_replica_box_pbc.csv`, `graphite_cutoffs_box_pbc.csv`, `graphite_gr_box_pbc.csv`, `graphite_rdf_accumulated.npz/.json`
  - Solution: `solution_cn_box_pbc.csv`, `solution_cutoffs_box_pbc_300K.csv`, `solution_gr_box_pbc_300K.csv`
  - Collapse screening: `collapse_screen.csv`
  - Figures: `*_box_pbc_300K.svg`
- `kinetics/`:
  - Rates and barriers: `kinetics_rates.csv`, `kinetics_arrhenius.csv`
  - Thermodynamics: `kinetics_occupancy.csv`, `kinetics_thermo.csv`, `kinetics_vanthoff.csv`
  - Cutoffs and diagnostics: `kinetics_cutoffs.csv`, `kinetics_oh_density.csv`, `kinetics_events_by_traj.csv`, `kinetics_occupancy_vs_time.csv`, `kinetics_traces_example.csv`
  - Settings: `kinetics_metadata.json`
  - Figures: `fig_*.svg`
- `comparison/`:
  - `compare_activation.csv` (master table of the three-step ΔG‡ and Ea)
  - `compare_pKa_stepwise.csv` (stepwise pKa, ΔG, ΔΔG)
  - `compare_pKa_cumulative.csv`, `compare_box_occupancy.csv`, `compare_structure_by_nOH.csv`, `compare_PMF_step1.csv`
  - Figures: `fig_compare_*.svg`
- `literature/`: transcriptions of the paper's Tables 1/2 and S1/S2, SI structures and geometry metrics, PMF digitization and end-point check, `lit_metadata.json`.

Section 8 of the page gives a one-line description of every file.

## Not included in the bundle

- `kinetics_cache/` (~372 MB): regenerate it from the raw trajectories with `python kinetics.py extract 8`.
- The publisher's SI file `literature/source/ic1c00662_si_001.pdf`: `literature.py` downloads it automatically.
- Logs and Python caches.

## Main limitations

See section 7 of the page for details.

- **Apparent quantities:** ΔG‡ and Ea depend on the filter window and the cutoff.
  - For example, ΔG‡ of graphite DPA4C-Plus Dy 3→2 is 6.2 / 10.1 kJ/mol at W = 0.2 / 2.0 ps.
  - The relative order of the three steps (forward step 1 < 2 < 3, reverse the opposite) holds in all 9 settings × 2 metals and is more reliable than the absolute values.
- **Few statistics in solution:** only one 150 ps trajectory per temperature, so the 95% CIs of Ea are wide. Some Ea values are negative, which shows that these are composite apparent rate constants rather than single elementary steps.
- **Assumptions behind the pKa equivalents:** pK_w = 14, ideal solution, and a standard state of c₁ = 1 molecule per box. The free OH⁻ concentration (4 − n)c₁ reaches about 1 M, with no activity or ion-pair corrections.
- **Starting-configuration memory:** the graphite DPA4C-Plus occupancies converge only after about 90 ps.
- **UMA collapse:** UMA trajectories often collapse, so the UMA results are for reference only.
- **Qualitative differences between potentials:** with DPA4-Pro, MACE and UMA the 3 OH stay on Ln almost all the time, whereas with DPA4C-Plus OH leaves frequently. Which is closer to reality needs independent validation by AIMD or experimental structures.
- **Comparability with the paper:** the paper has only 298 K thermodynamics and no barriers, so our Ea can be compared with the paper's ΔG only in order of magnitude.
- **OH solvation shell:** the O–X(total) shell minimum of graphite DPA4C-Plus Nd is unclear, so its shell O–X(total) should not be used for conclusions; use O–REE(shell) and REE–O(OH)/3 instead.
