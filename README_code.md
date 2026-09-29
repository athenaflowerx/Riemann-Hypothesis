# RH Proof Package — Code Instruments (README)

Authors: Guang Yang (ORCID: 0000-0003-0599-2881) · Yueting Xiao (ORCID: 0009-0002-0268-0655)
Date: 2026-09-28 · Package v3-13.1

Nine local numerical instruments plus a one-click regression driver. Everything runs
offline: no network access, no data download, no external services. Each instrument is
self-contained, prints its own PASS/FAIL verdict, and exits nonzero on failure.

> 中文要点见文末（Chinese summary at the end）。

---

## 1. Requirements

- Python ≥ 3.10
- `numpy`, `scipy`, `mpmath` (`pip install numpy scipy mpmath`)

No other dependencies. No GPU, no network.

## 2. Quick start

```bash
cd code
python3 run_all.py        # runs all 9 instruments sequentially
```

Expected final line: `== all regression passed ==` (total runtime ≈ 30 s).
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
| `run_all.py` | Regression driver | Runs the 9 instruments above, reports per-job PASS/FAIL and timings. |

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
== all regression passed ==
```

## 5. Key numerical results (reproduced by the suite)

- Weld gaps: ζ³@ρ₁ gap 2π (m = 2), ζ²⁵@ρ₁ gap 3π (m = 3), synthetic (s − ρ₁)⁵ gap 5π;
  pin gap π (m = 1) — order exact, by the argument principle.
- Quantum crossing: γₓ = 2π·e^{4π} ≈ 1.80 × 10⁶.
- Blind-region crossing min |ζ| = 0.12203.
- ±0i pair sums: 7.9 × 10⁻⁶ / 7.9 × 10⁻⁸.
- Implanted off-line zero β₀ = 0.6 reads walk slope ≈ 1.1 (theory: ½ + 0.6); real zeros read 1.0.

## 6. Notes

- `legacy_v1.9/` (archival originals) is intentionally **not** part of this release;
  the 9 instruments above are the current maintained versions.
- All computations are local; floating-point policy is stated inside each instrument
  (mpmath at ≥ 30 dps where branch decisions are made).
- License/citation: cite the package as *RH Proof Package v3-13.1 (2026-09-28)*.

---

## Citation

- Latest version (always points to newest): https://doi.org/10.5281/zenodo.23027848
- v2.0 (2026-09-28): https://doi.org/10.5281/zenodo.23019558

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23027848.svg)](https://doi.org/10.5281/zenodo.23027848)

---

## 中文摘要

9 个本地数值仪器 + 一键回归（`run_all.py`），全程离线、零联网。

- **依赖**：Python ≥ 3.10，`numpy / scipy / mpmath`。
- **运行**：`cd code && python3 run_all.py`，末行 `== all regression passed ==` 即全过（约 30 秒）；每个仪器也可单独运行。
- **构成**：θ 三守卫基础模块、焊接探测器（B.3/B.4 数据来源：钉 gap π / 焊 gap mπ，mod-π/2 下焊不可见 = 焊接投影损失 W2）、量子交叉外推（γₓ = 2π·e⁴ᵖⁱ ≈ 1.80×10⁶）、盲区极限（附录 A 数据源）、旋转扫描、五层统一判别器（植入 β₀=0.6 被捕获、真零点全过）、0i 纲领验证、V1–V6 清单核验、仪器层纪律自动化。
- **关键复现数**：焊 gap 2π/3π/5π（阶精确）；盲区穿越 min|ζ| = 0.12203；±0i 对和 7.9×10⁻⁶ / 7.9×10⁻⁸；植入零点读 slope ≈ 1.1，真零点读 1.0。
- **不含** `legacy_v1.9/`（归档旧版）；所有计算本地完成，分支判定处使用 mpmath ≥ 30 位精度。
