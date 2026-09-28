# rotation_scan.py · Rotation-scan discriminant (best expression form of T-LOCK)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Reproducible verification. Usage: python rotation_scan.py
# (local only, mpmath required; ~40 seconds)
import numpy as np, mpmath as mp

C_DVP = 1/57.9   # classical de la Vallee Poussin constant
G1 = 14.1347251417946932

def theta_RS(t, dps=25):
    mp.mp.dps = dps
    t = mp.mpf(t)
    p = mp.im(mp.loggamma(mp.mpf('0.25')+1j*t/2)) - t/2*mp.log(mp.pi)
    z = mp.mpf('0.25')+1j*t/2
    ref = mp.im((z-mp.mpf('0.5'))*mp.log(z)-z) - t/2*mp.log(mp.pi)
    return p + 2*mp.pi*mp.nint((ref-p)/(2*mp.pi))

def gamma_star(delta):
    """Finite-segment theorem: uncertified length of a non-vertical anchor ray ~ 1/(2*delta)."""
    return 0.5/delta

def ray_scan(delta, g0=0.5, g1=700, n=350, dps=12):
    """Anchor-ray scan: s = 1/2 + gamma*delta + i*gamma; returns (min|zeta|, argmin gamma)."""
    mp.mp.dps = dps
    gs = np.linspace(g0, g1, n)
    vals = [float(abs(mp.zeta(mp.mpf(0.5+g*delta)+1j*mp.mpf(g)))) for g in gs]
    i = int(np.argmin(vals))
    return min(vals), gs[i]

def near_miss_pred(gam_k, delta, dps=15):
    """Near-miss formula: min ~ |zeta'(rho_k)| * gamma_k * delta."""
    mp.mp.dps = dps
    return float(abs(mp.diff(mp.zeta, mp.mpf('0.5')+1j*mp.mpf(gam_k)))) * gam_k * delta

def F_pair(gam_k, delta, dps=25):
    """±0i seam reading: F(±delta) and their (mod pi/2) sum."""
    mp.mp.dps = dps
    Fa = float((mp.arg(mp.diff(mp.zeta, mp.mpf(0.5+gam_k*delta)+1j*mp.mpf(gam_k)))
                + theta_RS(gam_k) + mp.pi/2) % (mp.pi/2))
    Fb = float((mp.arg(mp.diff(mp.zeta, mp.mpf(0.5-gam_k*delta)+1j*mp.mpf(gam_k)))
                + theta_RS(gam_k) + mp.pi/2) % (mp.pi/2))
    return Fa, Fb, (Fa+Fb) % (mp.pi/2)

if __name__ == "__main__":
    print("=== Showcase A: blind-region crossing (point 0.75+300i) ===")
    delta = float(np.arctan(0.25/300))
    m, g = ray_scan(delta)
    print(f" delta={delta:.3e}, gamma*={gamma_star(delta):.0f}; ray min|zeta|={m:.5f} (at gamma={g:.1f})")
    print(f" near-miss formula (zero #1 local): {near_miss_pred(G1, delta):.5f}")
    print("=== Showcase B: ±0i seam (zero #1, delta -> 0) ===")
    for d in [1e-3, 1e-4]:
        Fa, Fb, s = F_pair(G1, d)
        print(f" delta={d:g}: F(+d)={Fa:.6f}, F(-d)={Fb:.6f}, sum mod pi/2 = {float(s):.2e}")
    print("Assertion: sum < 1e-5 (machine precision of the ±0i pair)",
          all(F_pair(G1, d)[2] < 1e-5 for d in [1e-3, 1e-4]))
