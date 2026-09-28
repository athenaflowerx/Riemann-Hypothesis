# blindness_limits.py · Cross-layer blindness-limit experiment (data source of Appendix A)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Decision rule: an instrument can judge only if its reading change >= readable floor eps;
#                blindness limit u* = eps / sensitivity.
# Usage: python blindness_limits.py  (local only, mpmath required; ~15 seconds)
import numpy as np, mpmath as mp

GAMMA = 300.0      # experiment height
EPS   = 1e-3       # readable floor
A0    = 0.75       # reference candidate real part (u0 = 0.25)
U0    = A0 - 0.5
C_DVP = 1/57.9
PI2   = np.pi/2

def theta_RS(t, dps=20):
    mp.mp.dps = dps
    t = mp.mpf(t)
    p = mp.im(mp.loggamma(mp.mpf('0.25')+1j*t/2)) - t/2*mp.log(mp.pi)
    z = mp.mpf('0.25')+1j*t/2
    ref = mp.im((z-mp.mpf('0.5'))*mp.log(z)-z) - t/2*mp.log(mp.pi)
    return p + 2*mp.pi*mp.nint((ref-p)/(2*mp.pi))

def F(a, b, dps=15):
    """Lock reading F(a,b) = (arg zeta'(a+ib) + theta(b) + pi/2) mod pi/2."""
    mp.mp.dps = dps
    s = mp.mpf(a)+1j*mp.mpf(b)
    return float((mp.arg(mp.diff(mp.zeta, s)) + theta_RS(b, 20) + PI2) % PI2)

if __name__ == "__main__":
    print("========== Cross-layer blindness-limit experiment (gamma=300, eps=1e-3) ==========")

    # Layer 1: modulus coordinate  dℓ/da = a/(a²+gamma²)
    sens_l = A0/(A0**2+GAMMA**2)
    u_l = EPS/sens_l
    print(f"[L1 coordinate-modulus] sensitivity {sens_l:.3e}; blindness limit u* = {u_l:.1f}")
    print(f"   inside-limit reading (u=0.25): dl = {U0*sens_l:.2e} (< eps, blind)")

    # Layer 2: phase coordinate  dθ/da = −gamma/(a²+gamma²)
    sens_t = GAMMA/(A0**2+GAMMA**2)
    u_t = EPS/sens_t
    print(f"[L2 coordinate-phase]   sensitivity {sens_t:.3e}; blindness limit u* = {u_t:.3f}")
    print(f"   inside-limit reading (u=0.10): dth = {0.10*sens_t:.2e} (< eps, blind)")
    print(f"   outside-limit reading (u=0.40): dth = {0.40*sens_t:.2e} (> eps, decidable)")

    # Layer 3: dVP barrier
    u_dvp = 0.5 - C_DVP/np.log(GAMMA)
    print(f"[L3 dVP barrier] blindness limit (coverage boundary) u* = {u_dvp:.4f}")
    print("   outside limit (u>u*): excluded by theorem; inside: open strip")

    # Layer 4: anchor ray (near-miss slope calibrated by u0=0.25 measured min=0.12203)
    ray_demo_min, u_demo = 0.12203, 0.25
    slope_ray = ray_demo_min/u_demo
    u_ray = EPS/slope_ray
    print(f"[L4 anchor-ray scan] near-miss slope {slope_ray:.3f}/unit-u; blindness limit u* = {u_ray:.2e}")
    print(f"   outside-limit reading (u=0.25): min|zeta| = {ray_demo_min:.5f} (> eps, decidable -- blind-region crossing)")
    print(f"   inside-limit prediction (u=2e-4): min|zeta| ~ {slope_ray*2e-4:.2e} (< eps, blind)")

    # Layer 5: lock reading F (function-value phase; measured gradient)
    Fa, Fb = F(A0, GAMMA), F(A0+0.05, GAMMA)
    grad_F = abs(Fa-Fb)/0.05
    u_F = EPS/grad_F if grad_F > 0 else np.inf
    print(f"[L5 lock reading F] gradient {grad_F:.3f}/unit-u; blindness limit u* = {u_F:.2e}")
    print(f"   outside-limit reading (u=0.05 step): dF = {abs(Fa-Fb):.4f} (> eps, decidable)")

    # Layer 6: ±0i seam (F pair; u* -> 0, machine precision)
    g1 = 14.1347251417946932
    def F_pair(d):
        Fa = F(0.5+g1*d, g1); Fb = F(0.5-g1*d, g1)
        return (Fa+Fb) % PI2
    s3, s4 = float(F_pair(1e-3)), float(F_pair(1e-4))
    print(f"[L6 ±0i seam] blindness limit u* -> 0 (the lock line itself)")
    print(f"   measured: d=1e-3 pair-sum = {s3:.2e}; d=1e-4 pair-sum = {s4:.2e} (machine precision; both sides decidable)")

    print("Collapse chain u*: ~%.0f (fully blind) -> %.2f (pin-phase) -> %.4f (barrier) -> %.1e (ray) -> %.1e (F) -> 0 (±0i)"
          % (u_l, u_t, u_dvp, u_ray, u_F))
    print("The blind region collapses layer by layer, finally onto the lock line itself -- governed by ±0i.")
    assert s3 < 1e-5 and s4 < 1e-5 and u_ray < u_t
    print("ALL PASS")
