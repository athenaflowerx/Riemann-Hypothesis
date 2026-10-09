# RH Proof Package — Code Instruments (README)

Authors: Guang Yang (ORCID: 0000-0003-0599-2881) · Yueting Xiao (ORCID: 0009-0002-0268-0655)
Date: 2026-10-09 · Package v3.6 + addendum v3.7

Nine local numerical instruments, a terminal-form suite (26 checks, sections A–I),
the v3.7 terminal-block-II suite (22 checks), a nine-section instrument suite, an
unreachability-transformation closed-form suite, plus a one-click regression driver
(14 jobs, legacy regression included). Everything runs offline: no network access,
no data download, no external services. Each instrument is self-contained, prints
its own PASS/FAIL verdict, and exits nonzero on failure.

> 中文要点见文末（Chinese summary at the end）。

---

## 1. Requirements

- Python ≥ 3.10
- `numpy`, `scipy`, `mpmath` (`pip install numpy scipy mpmath`)

No other dependencies. No GPU, no network.

## 2. Quick start

```bash
cd code
python3 run_all.py                        # 14 jobs: 9 instruments + 2 legacy + 3 v3.7 suites
python3 final_instruments_2026-09-30.py   # terminal-form suite A–I (26 checks)
python3 v37_verification.py               # terminal block II suite (22 checks)
```

Expected final lines: `== all regression passed ==` (≈ 60 s),
`== RESULT: 26/26 PASS ==` (≈ 30 s), and `22/22 PASS` (≈ 1 s).
Any instrument can also be run standalone, e.g. `python3 weld_detector.py`.

## 3. Instrument manifest

| File | Role | What it verifies / produces |
|---|---|---|
| `theta_RS_guarded.py` | Base module | Guarded Riemann–Siegel θ(t) with three guards (branch, precision, range); self-test PASS. Imported by `instrument_discipline.py`. |
| `weld_detector.py` | Weld detector | Unwrapped-phase (honest-record) instrument: pin gap π (m = 1) vs weld gap mπ (m = 2, 3, 5) via the argument principle with Richardson extrapolation; welds invisible under the mod-π/2 reading = projection loss at the weld (W2). Data source of Seam Theorems B.3/B.4. |
| `extrapolation_crossing.py` | Quantum crossing-zone extrapolation | Crossing of the height-quantum ladder and the discrimination cone at γₓ = 2π·e^{4π} ≈ 1.80 × 10⁶; N(γ) quantum cells, seam-theorem validity window. |
| `blindness_limits.py` | Cross-layer blindness limits | What each instrument layer cannot see; data source of Appendix A. |
| `rotation_scan.py` | Rotation-scan discriminant | Lock-line reading under rotation; on-line/off-line discrimination. |
| `rh_discriminator.py` | Five-layer unified discriminator | The counterexample-exclusion chain: synthetic implants (β₀ = 0.6 → slope 1.1) are caught; genuine zeros pass all five layers. |
| `verify_0i_program.py` | 0i-program numerical verification | Phase of ζ′(ρ₁) on the lock lattice; ±0i seam readings at machine precision. |
| `verify_V1_V6.py` | Checklist V1–V6 | M5 pin differences, skeleton-vertical heights, near-skeleton-vertical zeros (γ = 711.900289 …). |
| `instrument_discipline.py` | Discipline automation | Instrument-layer discipline checks (guards active, precision policy, no silent fallbacks). |
| `final_instruments_2026-09-30.py` | Terminal-form suite A–I | Reproduces every numerical claim of the v3.6 final increment: layer ladder (half-turn law / Gram deviation); configuration issuing + data-wall law; weld parity readout; correct/incorrect discrimination (m-integrality); convolution readout c₂; FE closure H-condition (vacuous in strip); FE residue + seam faces; Euclid increment + i-vortex (1+i)-power tower; face book (seam two archives in the Arg ξ′ readout). **26/26 PASS.** |
| `v37_verification.py` | Terminal block II suite | Every numerical anchor of the v3.7 addendum (TB-II): the issuing formula and the p̄ readout; the i^i discriminator ĉ_k(t) = exp(t·e^{−ikπ/2}) with archive readout at t = π/2; the canonical recursion (fixed point π/2 exact; rate λ = π/(4 ln 2) − 1/2 = 0.633090… closed form); the pair-coordinate map and the settled-representative convention; the zero anchor ξ′(ρ) = E(ρ)ζ′(ρ); the FE residue at dps = 40. **22/22 PASS.** |
| `RH_instruments_v6_2026-10-05.py` | Nine-section instrument suite (A–I) | The five-archive-operations instrument layer: layer ladder; configuration issuing + data-wall law; weld-parity readout; m-discrimination battery; convolution readout; H-condition vacuity; FE + seam faces; tower + Euclid increment; face book. Exits `ALL PASS`. |
| `unreachability_closedforms.py` | Unreachability transformation, closed forms | The four closed forms of the five archive operations: the generating series G_n = Σ (iπ/2)^k/k! vs its closed form; the remainder R_n = i·γ(n+1, iπ/2)/Γ(n+1) (exact |R_32| = 3.41 × 10⁻³¹); the recursion rate λ = π/(4 ln 2) − 1/2; the basin rates. All assertions PASS. |
| `run_all.py` | Regression driver | Runs the 9 instruments, 2 legacy regressions (`legacy_v1.9/verify_merged_v5.py`, `legacy_v1.9/rc_T126_T130_experiments.py`), and the 3 v3.7 suites; reports per-job PASS/FAIL and timings. |

## 4. Expected regression output

```
[PASS] theta_RS three-guard self-test            (~0 s)
[PASS] cross-layer blindness-limit experiment    (~0 s)
[PASS] rotation-scan discriminant                (~3 s)
[PASS] weld detector                             (~1 s)
[PASS] quantum crossing-zone extrapolation       (~1 s)
[PASS] five-layer unified discriminator          (~21 s)
[PASS] instrument-layer discipline automation    (~1 s)
[PASS] 0i program numerical verification         (~1 s)
[PASS] checklist V1-V6 verification              (~1 s)
[PASS] [legacy] merged v5                        (~1 s)
[PASS] [legacy] T126-T130 experiments            (~12 s)
[PASS] terminal-form suite A-I (26 checks)       (~30 s)
[PASS] terminal block II suite (22 checks)       (~1 s)
[PASS] unreachability closed forms               (~1 s)
== all regression passed ==
```

Terminal-form suite: `== RESULT: 26/26 PASS ==` (160-zero table; FE relative
residual 6.55 × 10⁻⁵¹ at dps = 50). Terminal block II suite: `22/22 PASS`.

## 5. Key numerical results (reproduced by the suites)

- Weld gaps: ζ³@ρ₁ gap 2π (m = 2), ζ²⁵@ρ₁ gap 3π (m = 3), synthetic (s − ρ₁)⁵ gap 5π;
  pin gap π (m = 1) — order exact, by the argument principle.
- Quantum crossing: γₓ = 2π·e^{4π} ≈ 1.80 × 10⁶.
- Blind-region crossing min |ζ| = 0.12203.
- ±0i pair sums: 7.9 × 10⁻⁶ / 7.9 × 10⁻⁸.
- Implanted off-line zero β₀ = 0.6 reads walk slope ≈ 1.1 (theory: ½ + 0.6); real zeros read 1.0.
- Terminal-form suite: 26/26 PASS, incl. the face book showing only the two
  [+0i, −0i] archive rows, and the FE closure H-condition vacuous inside the strip.
- v3.7 suite: 22/22 PASS, incl. the fixed-point identity 2^{α/2}·e^{−π/4} = 1
  (α = π/(2 ln 2)) exact, the recursion rate λ = 0.633090… (closed form), and the
  period–quantum discrimination (200 certified zeros: deviation 0.00; 60 random
  in-strip points: none on-lattice).

## 6. Notes

- `legacy_v1.9/` (archival originals) ships with the package and is exercised by
  `run_all.py` as two guarded regression jobs; the 9 instruments plus the suites
  above are the current maintained versions.
- All computations are local; floating-point policy is stated inside each instrument
  (mpmath at ≥ 30 dps where branch decisions are made; FE-residue checks at dps = 50
  in the terminal-form suite and dps = 40 in the v3.7 suite, as labeled).
- Numerics-compatible statement: all numerical claims unchanged in substance since
  v3.6; the terminal-form suite count is 26 (the v3.6 paper's "29 executed checks"
  line is superseded by the merge list in record §191); the v3.7 suite is 22.
- License/citation: cite the package as *RH Proof Package v3.6 + v3.7 addendum
  (2026-10-09)*.

---

## Citation

- v3.7 addendum (2026-10-09, terminal block II): https://doi.org/10.5281/zenodo.23261294
- v3.6 base (2026-09-30): https://doi.org/10.5281/zenodo.23068324
- v2.1 (2026-09-29): https://doi.org/10.5281/zenodo.23027848
- v2.0 (2026-09-28): https://doi.org/10.5281/zenodo.23019558

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23261294.svg)](https://doi.org/10.5281/zenodo.23261294)

---

## 中文摘要

9 个本地数值仪器 + 终形核验套件（A–I 共 **26** 项）+ v3.7 终块 II 套件（**22** 项）+
九节仪器套件 + 不可达转化闭式套件 + 一键回归（`run_all.py`，**14** 个任务，
含 2 项 legacy 回归），全程离线、零联网。

- **依赖**：Python ≥ 3.10，`numpy / scipy / mpmath`。
- **运行**：`cd code && python3 run_all.py`，末行 `== all regression passed ==`（约 60 秒）；
  `final_instruments_2026-09-30.py` 末行 `== RESULT: 26/26 PASS ==`（约 30 秒）；
  `v37_verification.py` 末行 `22/22 PASS`（约 1 秒）；每个仪器可单独运行。
- **v3.7 新增**：终块 II 套件（签发公式与 p̄ 读出；i^i 判别式 ĉ_k(t) = exp(t·e^{−ikπ/2})，
  t = π/2 时偶通道读出双档案；典范递归不动点 π/2 恒等、速率闭式
  λ = π/(4ln2) − 1/2 = 0.633090…；对坐标映射与定居代表约定；零锚定
  ξ′(ρ) = E(ρ)ζ′(ρ)；FE 残差 dps=40）；九节仪器套件；不可达转化四闭式
  （G_n 级数对闭式、余项 R_n = i·γ(n+1, iπ/2)/Γ(n+1)，|R_32| = 3.41×10⁻³¹）。
- **关键复现数**：焊 gap 2π/3π/5π（阶精确）；盲区穿越 min|ζ| = 0.12203；
  ±0i 对和 7.9×10⁻⁶ / 7.9×10⁻⁸；植入零点读 slope ≈ 1.1，真零点读 1.0；
  FE 相对残差 6.55×10⁻⁵¹（dps=50）；不动点恒等 2^{α/2}·e^{−π/4} = 1 精确；
  周期-量子判别：200 认证零点偏差 0.00，带内 60 随机点无一落格。
- `legacy_v1.9/` 随包附存并纳入回归守护；分支判定处 mpmath ≥ 30 位，
  FE 残差检查按套件标注 dps（50 / 40）。
