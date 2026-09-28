# weld_detector.py · Weld detector (unwrapped-phase instrument / 4D honest record)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# [0,0i] language: a simple zero is an antipodal pin pair (pin gap pi); a multiple
# zero is a weld (pin gap m*pi). The F-reading (mod pi/2, a 3D projection) cannot
# distinguish them (W2 projection loss). This instrument reads unwrapped phase (mod 2pi)
# and additionally gives the exact order via the argument principle.
#
# Two readings:
#   pin gap  Δ = continuous arg variation along the upper semicircle from (rho-u) to
#              (rho+u)  ->  m*pi  (Richardson-extrapolated, u -> 0)
#   order    m = (1/2pi) * contour integral of d arg g around |s-rho| = r
# Self-test: zeta@rho1 (1) · zeta^2@rho1 (2) · zeta^3@rho1 (3) · (s-rho1)^5 (5) ·
#            toy f_eps on-line zero (1).  Dependency: mpmath. Local, no network.
import mpmath as mp

mp.mp.dps = 40
PI  = mp.pi
GAMMA1 = mp.mpf('14.134725141734693790457251983562470270784257115699')
RHO1   = mp.mpf('0.5') + 1j*GAMMA1

def pin_gap_raw(g, rho, M=32):
    """Pin gap: continuous arg variation along the upper semicircle -> m*pi.
    Richardson extrapolation removes the O(u) drift."""
    def at(u):
        total = mp.mpf('0')
        prev = mp.arg(g(rho + u*mp.e**(1j*PI)))
        for j in range(1, M + 1):
            th = PI*(1 - j/M)
            cur = mp.arg(g(rho + u*mp.e**(1j*th)))
            d = cur - prev
            while d >  PI: d -= 2*PI
            while d < -PI: d += 2*PI
            total += d
            prev = cur
        return abs(total)
    g1, g2 = at(mp.mpf('1e-3')), at(mp.mpf('1e-4'))
    return (10*g2 - g1)/9        # linear extrapolation u->0, remainder O(u^2)

def winding_order(g, rho, r=mp.mpf('1e-3'), M=64):
    """Exact order by the argument principle: total arg variation / 2pi on a small circle."""
    total = mp.mpf('0')
    prev  = mp.arg(g(rho + r))
    for j in range(1, M + 1):
        th = 2*PI*j/M
        cur = mp.arg(g(rho + r*mp.e**(1j*th)))
        d = cur - prev
        while d >  PI: d -= 2*PI
        while d < -PI: d += 2*PI
        total += d
        prev = cur
    return int(round(float(total / (2*PI))))

def detect(g, rho, label, expect):
    m  = winding_order(g, rho)
    gap = pin_gap_raw(g, rho)
    ok = (m == expect) and (abs(gap - expect*PI) < mp.mpf('1e-7'))
    parity = 'antipodal(odd)' if (abs((gap/PI) % 2 - 1) < mp.mpf('1e-6')) else 'weld(even)'
    print(f"  {label:<30} order m={m} (expect {expect})  gap/pi={mp.nstr(gap/PI,10)}  [{parity}]  {'OK' if ok else 'FAIL'}")
    return ok

def main():
    print("=== Weld detector self-check (unwrapped phase + argument principle) ===")
    ok = True
    ok &= detect(mp.zeta,                 RHO1, "zeta @ rho1 (simple)",       1)
    ok &= detect(lambda s: mp.zeta(s)**2, RHO1, "zeta^2 @ rho1 (double weld)", 2)
    ok &= detect(lambda s: mp.zeta(s)**3, RHO1, "zeta^3 @ rho1 (triple weld)", 3)
    ok &= detect(lambda s: (s - RHO1)**5, RHO1, "(s-rho1)^5 (synthetic 5)",    5)

    eps = mp.mpf('1e-4')
    def f_eps(s):
        z = s - mp.mpf('0.5')
        return mp.cosh(PI*z) - eps*mp.cos(PI*z)
    z0 = mp.findroot(lambda z: mp.cosh(PI*z) - eps*mp.cos(PI*z), 0.5j)
    ok &= detect(f_eps, mp.mpf('0.5') + z0, "f_eps on-line zero (toy world)", 1)

    print("\n=== Contrast: weld blindness of the F-reading (mod pi/2) ===")
    th1 = mp.im(mp.loggamma(mp.mpf('0.25') + 0.5j*GAMMA1)) - 0.5*GAMMA1*mp.log(PI)
    PI2 = PI/2
    def F(v, phi): return (mp.arg(v) + phi + PI2) % PI2
    u = mp.mpf('1e-3')
    def seam(g, phi, m, rho):
        return F(mp.diff(g, rho + u, m), phi) + F(mp.diff(g, rho - u, m), phi)
    s1 = seam(mp.zeta, th1, 1, RHO1)
    s2 = seam(lambda s: mp.zeta(s)**2, 2*th1, 2, RHO1)
    print(f"  zeta  simple zero: F pair-sum = {mp.nstr(s1,12)} (|.-pi/2|={mp.nstr(abs(s1-PI2),8)})")
    print(f"  zeta^2 double zero: F pair-sum = {mp.nstr(s2,12)} (|.-pi/2|={mp.nstr(abs(s2-PI2),8)})")
    print("  -> identical structure mod pi/2 (weld invisible); mod-2pi unwrapped reading distinguishes pi/2pi/3pi/5pi")
    print("\n" + ("ALL PASS -- weld detector self-check complete" if ok else "FAIL"))
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
