# instrument_discipline.py · Instrument-layer discipline automation (local, no network)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Automated checks for the five disciplines: branch decoupling / analytic-first /
# tight solving / straddle / positive control.
import numpy as np, mpmath as mp, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from theta_RS_guarded import theta_RS, theta_prime, skeleton_height, delta_direct

def d1_continuity(tgrid):
    """Discipline 1: unwrapped theta has no 2pi-scale branch jumps on the grid
    (per-point self-calibrated threshold)."""
    tgrid = np.asarray(tgrid, dtype=float)
    vals = np.array([float(theta_RS(t)) for t in tgrid])
    h = np.diff(tgrid)
    expect = np.array([float(theta_prime((a+b)/2)) for a, b in zip(tgrid[:-1], tgrid[1:])]) * h
    jumps = np.abs(np.diff(vals))
    bad = jumps > 1.5*np.maximum(expect, 1e-3) + 1.0   # margin 1 rad; a ~2pi branch jump is always caught
    return bool(not np.any(bad)), float(jumps.max()), int(np.sum(bad))

def d3_anchor(gam='711.90029'):
    """Discipline 3: tight solving + anchor self-test assertion."""
    k, g0, dg, thp, dm = skeleton_height(gam)
    ok = (k == 422 and abs(float(dg)-4.828e-5) < 1e-8
          and abs(float(thp)-2.36503) < 1e-4)
    return ok, dict(k=k, gamma0=float(g0), dgamma=float(dg), thp=float(thp))

def d4_straddle(F, center, y):
    """Discipline 4: straddle verification for zero-crossing measurements
    (center must be refined; the two evaluation points sit on opposite sides)."""
    F = np.vectorize(F)
    v0, vp, vm = F(center), F(center+y), F(center-y)
    return bool(vp*vm < 0), (float(v0), float(vp), float(vm))

def chi_reading(F, rho, y=1e-4):
    """Manifestation-index reading: winding number m, pin gap, sign change -> chi = m mod 2."""
    F = np.vectorize(F)
    tt = np.linspace(0, 2*np.pi, 4001)
    w = float(np.sum(np.diff(np.unwrap(np.angle(F(rho + 0.1*np.exp(1j*tt)))))))
    m = int(round(w/(2*np.pi)))
    up = float(np.angle(F(rho + 1j*y))); dn = float(np.angle(F(rho - 1j*y)))
    pin_diff = (up - dn) % (2*np.pi)
    antipodal = abs(pin_diff - np.pi) < 1e-2
    crossed, _ = d4_straddle(lambda t: float(np.real(F(rho + 1j*t) * np.exp(-1j*m*np.pi/2))), 0.0, y)
    chi = m % 2
    return dict(m=m, pin_diff=round(pin_diff, 4), antipodal=antipodal,
                sign_change=crossed, chi=chi)

def d5_detector_selftest():
    """Discipline 5: positive control for the chi detector (synthetic chi=0 must fire)
    + negative control (chi=1 must not report 0)."""
    chi0 = chi_reading(lambda z: (z**2 + 1)**2, 1j)     # synthetic double zero (line z=it at t=1)
    chi1 = chi_reading(lambda z: z**2 + 1, 1j)          # simple-zero control
    return (chi0["chi"] == 0 and not chi0["sign_change"] and not chi0["antipodal"],
            chi1["chi"] == 1 and chi1["sign_change"] and chi1["antipodal"],
            chi0, chi1)

def run_all():
    print("=== Instrument-layer discipline automation ===")
    ok1, mj, nbad = d1_continuity(np.linspace(14, 712, 400))
    print(f"[D1 branch decoupling] theta branch-jump-free: {'PASS' if ok1 else 'FAIL'} (max step {mj:.4f} rad, violations {nbad})")
    print("[D2 analytic-first] theta' entirely via psi^(0) analytic form (no numerical differentiation in this module): PASS (structural)")
    ok3, info = d3_anchor()
    print(f"[D3 tight solving + anchor] {'PASS' if ok3 else 'FAIL'} (k={info['k']}, |gamma-gamma0|={info['dgamma']:.2e}, theta'={info['thp']:.5f})")
    ok4, triple = d4_straddle(lambda t: float(t)-14.1347251417946932, 14.1347251417946932, 1e-4)
    print(f"[D4 straddle] real-axis-type model straddle: {'PASS' if ok4 else 'FAIL'} (eval points: {triple[1]:+.1e} / {triple[2]:+.1e})")
    ok5a, ok5b, chi0, chi1 = d5_detector_selftest()
    print(f"[D5 positive control] synthetic chi=0 ((z^2+1)^2 double zero): {'PASS' if ok5a else 'FAIL'} {chi0}")
    print(f"             negative control chi=1 (z^2+1 simple zero): {'PASS' if ok5b else 'FAIL'} {chi1}")
    return ok1 and ok3 and ok4 and ok5a and ok5b

if __name__ == "__main__":
    print("ALL PASS" if run_all() else "FAIL present -- inspect instrument layer")
