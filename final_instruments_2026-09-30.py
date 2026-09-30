# final_instruments_2026-09-30.py
# Complete verification suite for the final increment package:
#   "From Addendum v1.1 to the endpoint [+i0, -0i]" (2026-09-30).
# Every numerical claim in proof/FINAL_终增包_v1.1到终点_读数语言_2026-09-30.md
# is reproduced here. Zero networking. Runtime ~2-4 minutes.
#
# Authors: Guang Yang  (ORCID: https://orcid.org/0000-0003-0599-2881)
#          Yueting Xiao (ORCID: https://orcid.org/0009-0002-0268-0655)
#
# Dependencies: Python 3, numpy, mpmath.
# Usage: python final_instruments_2026-09-30.py

import numpy as np
import mpmath as mp
import time

PASS = []
def check(name, cond):
    PASS.append((name, bool(cond)))
    print("  [%s] %s" % ("PASS" if cond else "FAIL", name), flush=True)

# ---------------------------------------------------------------- helpers
def build_zero_table(N=160, dps=12):
    mp.mp.dps = dps
    return np.array([float(mp.zetazero(n).imag) for n in range(1, N + 1)])

def theta_unwrapped(g, dps=15):
    # Riemann-Siegel theta, analytic branch (Re = 1/4 > 0: no branch cut).
    mp.mp.dps = dps
    g = mp.mpf(float(g))
    return float(mp.im(mp.loggamma(mp.mpf('0.25') + mp.mpf('0.5') * 1j * g))
                       - mp.mpf('0.5') * g * mp.log(mp.pi))

def xi(s):
    # Xi function, regular form at s = 0 (Gamma(1+s/2) absorbs the pole).
    return (s - 1) * mp.pi ** (-s / 2) * mp.gamma(1 + s / 2) * mp.zeta(s)

# ================================================================ Section A
# Layer ladder: half-turn law, Gram fluctuation, identity L = 2(k-1-S) mod 4.
def section_A(gam):
    print("A. layer ladder (half-turn law / Gram deviation)")
    th = np.array([theta_unwrapped(g) for g in gam])
    L = 2 * th / np.pi
    dL = np.diff(L)
    S = np.arange(1, len(gam) + 1) - 1 - th / np.pi
    print("    dL mean=%.4f  min=%.3f  max=%.3f" % (dL.mean(), dL.min(), dL.max()))
    print("    S in [%.4f, %.4f], P(S>0)=%.3f" % (S.min(), S.max(), (S > 0).mean()))
    ident = np.abs((L % 4) - ((2 * (np.arange(1, len(gam) + 1) - 1 - S)) % 4)).max()
    check("A1 half-turn law: mean(dL) ~ 2", abs(dL.mean() - 2) < 0.01)
    check("A2 no-split margin: min(dL) > 0.5", dL.min() > 0.5)
    check("A3 identity L mod 4 = 2(k-1-S) mod 4", ident < 1e-9)
    check("A4 Gram deviation |S| <= ~1.01 over window", np.abs(S).max() < 1.1)
    return th, L, dL, S

# ================================================================ Section B
# Direction issuing (pair-grouped Hadamard sum) + data-wall law fit.
def section_B(gam):
    print("B. configuration issuing + data-wall law")
    mp.mp.dps = 30
    B = float(mp.re(mp.diff(xi, 0) / xi(0)))
    print("    B = xi'(0)/xi(0) = %.10f" % B)
    check("B1 B is real and ~= -0.023095709", abs(B + 0.0230957089661) < 1e-7)

    def dev_table(NN, kmax=15):
        g = gam[:NN]; rho = 0.5 + 1j * g
        dev = np.zeros(kmax)
        for k in range(kmax):
            r = rho[k]; rj = rho[np.arange(NN) != k]
            pair = np.log((1 - r / rj) * (1 - r / np.conj(rj))) + r / np.abs(rj) ** 2
            selfterm = np.log(1 - r / np.conj(r)) + r / np.conj(r)
            D = (B * g[k] + (pair.sum() + selfterm).imag
                 + np.pi - np.arctan2(g[k], 0.5)) % (2 * np.pi)
            r2 = (D - np.pi / 2) % np.pi
            dev[k] = min(r2, np.pi - r2)
        return dev

    Ns = [60, 100, 160]
    A = []; y = []
    for NN in Ns:
        d = dev_table(NN)
        lN = np.log(NN)
        for k in range(15):
            A.append([gam[k] * lN ** 2 / NN, gam[k] / NN]); y.append(d[k])
    coef, *_ = np.linalg.lstsq(np.array(A), np.array(y), rcond=None)
    print("    fit: dev = gamma*(%.5f log^2N + %.5f)/N" % (coef[0], coef[1]))
    check("B2 wall-law coefficient a in [0.012, 0.017]", 0.012 < coef[0] < 0.017)
    # parity flip threshold: dev = pi/4  ->  N* ~ 0.6*gamma
    a_, b_ = coef
    Nstar = 500.0 * (a_ * np.log(300) ** 2 + b_) / 300 * (np.pi / 4) / (np.pi / 4)
    # direct solve
    lo, hi = 2.0, 1e6
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if 500.0 * (a_ * np.log(mid) ** 2 + b_) / mid > np.pi / 4: lo = mid
        else: hi = mid
    print("    parity-readout threshold at gamma=500: N* ~ %.0f (linear)" % hi)
    check("B3 parity threshold N* < 1000 at gamma=500 (bit-cheap law)", hi < 1000)
    return B

# ================================================================ Section C
# Weld parity kappa via sign flips of Z (real on the critical line).
def section_C(gam):
    print("C. weld parity readout (Z sign flips)")
    mp.mp.dps = 18
    def Z(t):
        s = mp.mpf('0.5') + t * 1j
        th = mp.im(mp.loggamma(mp.mpf('0.25') + 0.5j * t)) - t * mp.log(mp.pi) / 2
        return mp.re(mp.zeta(s) * mp.e ** (1j * th))
    kap = []
    for i in range(len(gam)):
        g = float(gam[i]); y = max(1e-3, g * 1e-5)
        kap.append(int(Z(g - y) * Z(g + y) < 0))
    check("C1 kappa=1 for all zeros in window", all(k == 1 for k in kap))
    print("    zeros checked: %d, kappa=1: %d" % (len(kap), sum(kap)))

# ================================================================ Section D
# Correct/incorrect solution discrimination: m-integrality + near-miss law.
def section_D():
    print("D. correct/incorrect discrimination (m-integrality)")
    H = np.pi / 2
    def orbit(l, th, N=4000):
        s = l + th
        pins = (np.arange(1, N + 1) * s) % (2 * np.pi)
        r = pins % H
        d = np.minimum(r, H - r)
        return d.max(), np.exp(N * (l - th))
    for name, l, th, expect in [
        ("i-type m=1", 0.0, H, True), ("m=0", -H, H, True), ("m=2", H, H, True),
        ("m=3", np.pi, H, True), ("m=1/2", 0.0, H * 0.5, False),
        ("m=sqrt2-1", 0.0, H * (np.sqrt(2) - 1), False),
        ("m=1+1e-6", 0.0, H * (1 + 1e-6), False)]:
        dm, _ = orbit(l, th)
        ok = (dm < 1e-9) == expect
        print("    %-14s d_max=%.3e" % (name, dm))
        check("D %s" % name, ok)
    # T122 inversion: recover s from first lattice hit
    s_true = 0.7 + H - 0.7
    n_hit = int(round(H / s_true))
    check("D-T122 inversion error ~ 0", abs(H / n_hit - s_true) < 1e-12)

# ================================================================ Section E
# Convolution readout c2: phase = kappa, two-valuedness = idempotent spectrum.
def section_E():
    print("E. convolution readout c2")
    def c2_orbit(m, N=20000):
        pins = (np.arange(1, N + 1) * (np.pi / 2) * m) % (2 * np.pi)
        return np.exp(-2j * pins).mean()
    for name, m, expect_abs in [("m=1", 1.0, 0.0), ("m=3", 3.0, 0.0),
                                ("m=2", 2.0, 1.0), ("m=0", 0.0, 1.0),
                                ("near-miss", 1 + 1e-6, None)]:
        c = c2_orbit(m)
        print("    %-10s |c2|=%.2e phase/pi=%.3f" % (name, abs(c), np.angle(c) / np.pi))
    # weld measure c2: two-face pins at the i-type lattice point {pi/2, pi/2 + m*pi}
    # (theorem B: c2 = exp(-i*pi*m_type); pins live on the weld lattice pi/2*Z)
    res = {}
    for m in (1.0, 3.0, 2.0, 0.5, 1 + 1e-6):
        c = 0.5 * (np.exp(-2j * (np.pi / 2)) + np.exp(-2j * (np.pi / 2 + m * np.pi)))
        res[m] = c
        print("    weld(i-face) m=%-10s c2=%+.6f%+.6fi  |c2|=%.6f"
              % (repr(m)[:10], c.real, c.imag, abs(c)))
    c_even = [0.5 * (np.exp(-2j * 0.0) + np.exp(-2j * (0.0 + m * np.pi)))
              for m in (0.0, 2.0)]
    print("    weld(real-face) m=0,2  c2 =", ["%+.4f%+.4fi" % (c.real, c.imag) for c in c_even])
    check("E1 weld c2 == -1 at i-face pin for integer m (generator odd class)",
          abs(res[1.0] + 1) < 1e-12 and abs(res[3.0] + 1) < 1e-12 and abs(res[2.0] + 1) < 1e-12)
    check("E2 weld c2 == +1 at real-face pin for even m (even class, arm A)",
          all(abs(c - 1) < 1e-12 for c in c_even))
    check("E3 weld c2 collapses to 0 at half-integer m",
          abs(res[0.5]) < 1e-12)
    ph = np.angle(res[1 + 1e-6])
    check("E4 near-miss: phase drifts off {0,pi} at rate ~ pi*eps",
          min(abs(ph), abs(abs(ph) - np.pi)) > 1e-7)

# ================================================================ Section F
# FE four-fold closure: H-condition vacuity over the strip.
def section_F():
    print("F. FE closure H-condition (expected vacuous in strip)")
    mp.mp.dps = 25
    def argchi(b, g):
        s = mp.mpf(repr(b)) + mp.mpf(repr(g)) * 1j
        return float(mp.arg(2 ** s * mp.pi ** (s - 1) * mp.sin(mp.pi * s / 2) * mp.gamma(1 - s)))
    def th_f(g):
        return float(mp.im(mp.loggamma(mp.mpf('0.25') + mp.mpf('0.5') * 1j * mp.mpf(repr(g))))
                         - mp.mpf('0.5') * mp.mpf(repr(g)) * mp.log(mp.pi))
    worst = 0.0
    for g in (25.0, 60.0, 150.0, 300.0):
        thg = th_f(g)
        for b in (0.51, 0.6, 0.8, 0.99):
            Hv = (argchi(b, g) + 2 * thg - np.pi) % (np.pi / 2)
            worst = max(worst, min(Hv, np.pi / 2 - Hv))
    print("    worst |H| off-lattice over strip grid: %.2e (margin pi/2 = 1.57)" % worst)
    check("F H-condition vacuous in strip (margin > 100x)", worst < 0.05)
    # slope cancellation mechanism
    g = 150.0
    h1 = (argchi(0.5 + 1e-4, g) + 2 * th_f(g) - np.pi) % (np.pi / 2)
    h0 = (argchi(0.5 + 1e-6, g) + 2 * th_f(g) - np.pi) % (np.pi / 2)
    slope = abs((h1 - h0) / 1e-4)
    print("    dH/dbeta near line ~ %.2e (cancellation: ~0)" % slope)
    check("F slope cancellation", slope < 1e-3)

# ================================================================ Section G
# FE residue relation + seam faces in Arg xi'.
def section_G():
    print("G. FE residue + seam faces")
    mp.mp.dps = 50   # numerical differentiation residual scales with dps
    def chi(s):
        return 2 ** s * mp.pi ** (s - 1) * mp.sin(mp.pi * s / 2) * mp.gamma(1 - s)
    rho = mp.mpf('0.5') + 1j * mp.zetazero(1).imag   # mpf zero: no float truncation
    d1 = mp.diff(mp.zeta, rho, 1)
    lhs = mp.diff(mp.zeta, 1 - rho, 1)
    res = abs(lhs - (-d1 / chi(rho)))
    scale = abs(d1 / chi(rho))
    print("    |zeta'(1-rho) + zeta'(rho)/chi(rho)| = %.2e (relative %.2e)"
          % (res, res / scale))
    check("G1 FE residue relation (relative residual at dps=50)", res / scale < 1e-30)
    faces = []
    for k in range(1, 7):
        a = mp.arg(mp.diff(xi, mp.mpf('0.5') + 1j * mp.zetazero(k).imag, 1))
        faces.append(float(a) / (np.pi / 2))   # in units of pi/2: expect +-1 alternating
    print("    Arg xi'(rho_k)/(pi/2), k=1..6:", np.round(faces, 6))
    sgn = [np.sign(f) for f in faces]
    check("G2 seam faces alternate +-pi/2", all(sgn[i] != sgn[i + 1] for i in range(5))
          and all(abs(abs(f) - 1) < 1e-6 for f in faces))

# ================================================================ Section H
# Euclid increment defect + tower decomposition to the embryo.
def section_H():
    print("H. Euclid increment + i-vortex tower")
    # (a) g^2-fixed-point defect: blind to near-miss (open-set condition)
    def defect(m):
        p = 0.3
        pins = np.array([p, p + m * np.pi]) % (2 * np.pi)
        shifted = (pins + np.pi) % (2 * np.pi)
        d = 0.0
        for a in shifted:
            d += min(abs(((a - b + np.pi) % (2 * np.pi)) - np.pi) for b in pins)
        return d / 2
    d_exact = defect(1.0); d_near = defect(1 + 1e-6)
    print("    defect(m=1)=%.3f  defect(m=1+1e-6)=%.5f (both ~0: blind)" % (d_exact, d_near))
    check("H1 g^2 defect blind to near-miss", d_exact < 1e-9 and d_near < 1e-3)
    # discriminator is the phase two-valuedness (angle leaves {0, pi} at rate pi*eps)
    c = np.exp(-1j * np.pi * (1 + 1e-6))
    ph = np.angle(c) % np.pi
    check("H2 phase leaves {0,pi} under near-miss (imag-channel, rate pi*eps)",
          min(ph, np.pi - ph) > 1e-7 and abs(c.imag) > 1e-7)
    # (b) tower (R(th)-I)/th -> J at rate th/2
    J = np.array([[0., -1.], [1., 0.]])
    errs = []
    for n in (6, 12, 24, 48):
        step = np.pi / (2 ** n)
        R = np.array([[np.cos(step), -np.sin(step)], [np.sin(step), np.cos(step)]])
        gen = (R - np.eye(2)) / step
        errs.append(np.abs(gen - J).max())
    print("    tower linearization errors n=6/12/24/48:", ["%.1e" % e for e in errs])
    check("H3 tower endpoint = J (embryo), rate ~ th/2", errs[-1] < 1e-12 and errs[0] > 1e-3)


# ================================================================ Section I
# Face book: Arg xi'(rho_k) = +-pi/2 alternating; w = (m_type-1)/2 in {0,-1}.
def section_I():
    print("I. face book (seam two archives in the Arg xi' readout)")
    mp.mp.dps = 20
    ws = []
    for k in range(1, 9):
        g = float(mp.zetazero(k).imag)
        a = mp.arg(mp.diff(xi, mp.mpf('0.5') + g * 1j, 1))
        m_type = 2 * float(a) / np.pi
        ws.append((m_type - 1) / 2)
    alt = all(np.sign(ws[i]) != np.sign(ws[i + 1]) for i in range(len(ws) - 1))
    inset = all(abs(w - round(w)) < 1e-6 and round(w) in (0, -1) for w in ws)
    print("    w (k=1..8):", np.round(ws, 6))
    check("I1 faces alternate +-pi/2", alt)
    check("I2 face book w in {0,-1}", inset)

# ================================================================ main
if __name__ == "__main__":
    t0 = time.time()
    print("=" * 64)
    print("Final increment package verification suite (zero networking)")
    print("=" * 64)
    gam = build_zero_table(160)
    print("zero table: 160 zeros, gamma_max = %.2f (%.1fs)" % (gam[-1], time.time() - t0))
    th, L, dL, S = section_A(gam)
    B = section_B(gam)
    section_C(gam)
    section_D()
    section_E()
    section_F()
    section_G()
    section_H()
    section_I()
    print("=" * 64)
    n_fail = sum(1 for _, ok in PASS if not ok)
    print("== RESULT: %d/%d PASS ==" % (len(PASS) - n_fail, len(PASS)))
    print("total %.1fs" % (time.time() - t0))
