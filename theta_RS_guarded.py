# theta_RS_guarded.py · Guarded Riemann-Siegel theta module (three guards)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Local only, no network.
# Guard 1: Stirling continuous reference for branch unwrapping;
# Guard 2: analytic psi derivative (no numerical differentiation);
# Guard 3: tight findroot.
# Self-test: gamma = 711.90029 -> k=422, |gamma-gamma0|=4.828e-5, theta'=2.36503
import mpmath as mp

def theta_RS(t, dps=25):
    """Continuous (unwrapped) Riemann-Siegel theta(t)."""
    mp.mp.dps = dps
    t = mp.mpf(t)
    p   = mp.im(mp.loggamma(mp.mpf('0.25')+1j*t/2)) - t/2*mp.log(mp.pi)
    z   = mp.mpf('0.25') + 1j*t/2
    ref = mp.im((z-mp.mpf('0.5'))*mp.log(z) - z) - t/2*mp.log(mp.pi)
    return p + 2*mp.pi*mp.nint((ref - p)/(2*mp.pi))

def theta_prime(t, dps=25):
    """Analytic derivative theta'(t) = 0.5*Re psi(1/4+it/2) - 0.5*ln pi."""
    mp.mp.dps = dps
    t = mp.mpf(t)
    return 0.5*mp.re(mp.psi(0, mp.mpf('0.25')+1j*t/2)) - 0.5*mp.log(mp.pi)

def skeleton_height(gam, dps=25):
    """Nearest skeleton height: gamma0 with theta(gamma0) = pi/2 + k*pi.
    Returns (k, gamma0, |gamma-gamma0|, theta'(gamma0), delta_model)."""
    mp.mp.dps = dps
    gam = mp.mpf(gam); thv = theta_RS(gam)
    k = int(round((thv - mp.pi/2)/mp.pi))
    g0 = mp.findroot(lambda t: theta_RS(t)-(mp.pi/2+k*mp.pi), (gam-1, gam+1),
                     tol=mp.mpf('1e-12'), solver='anderson')
    thp = theta_prime(g0)
    return k, g0, abs(gam-g0), thp, abs(gam-g0)*thp

def delta_direct(gam, dps=25):
    """delta_direct = dist(theta(gamma) mod pi, pi/2)."""
    mp.mp.dps = dps
    return abs(float(theta_RS(gam)) % float(mp.pi) - float(mp.pi)/2)

if __name__ == "__main__":
    k, g0, dg, thp, dm = skeleton_height('711.90029')
    dd = delta_direct('711.90029')
    assert k == 422, f"k={k} != 422"
    assert abs(float(dg) - 4.828e-5) < 1e-8
    assert abs(float(thp) - 2.36503) < 1e-4
    assert abs(dm - dd) < 1e-9
    print("Self-test passed: k=422, |gamma-gamma0|=4.828e-5, theta'=2.36503, delta_model=delta_direct")
