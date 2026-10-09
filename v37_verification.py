# v37_verification.py · Terminal Block II verification suite (every numerical anchor of the v3.7 addendum; 22 checks)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Covers: issuing formula + p̄ readout; i^i discriminator c_k(t) = exp(t·e^{−ikπ/2}); canonical recursion
#         (fixed point π/2 exact; rate λ = π/(4 ln 2) − 1/2 closed form); pair-coordinate map;
#         zero anchor ξ′(ρ) = E(ρ)ζ′(ρ); FE residue at dps = 40.
# Usage: python v37_verification.py  (local only, mpmath required; ~1 second)
# v37_verification.py — Terminal Block II verification suite (v3.7)
# Authors: Guang Yang, Yueting Xiao — 2026-10-09
# Local only, zero networking.
from mpmath import mp, mpf, pi, exp, log, arg, gamma, sin, psi, zeta, diff, loggamma, gammainc, mpc

results = []
def chk(name, cond, detail=""):
    results.append(bool(cond))
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")

def moddist(x, P):
    """distance from x to the nearest multiple of P"""
    r = x % P
    return min(abs(r), abs(r - P))

mp.dps = 40
G1 = mpf('14.13472514173469379045725198357488017577')
rho = mpf('0.5') + G1*1j

# xi and exact logarithmic pieces (avoids cancellation in differentiating tiny products)
def xi(s):
    return gamma(s/2)*(s-1)*pi**(-s/2)*mpf('0.5')*zeta(s)
_dig = lambda z: diff(mp.loggamma, z)  # digamma via loggamma derivative (env-robust, A3 fix)

def xip_parts(s):
    E = mpf('0.5')*s*(s-1)*pi**(-s/2)*gamma(s/2)
    dlogE = 1/s + 1/(s-1) - log(pi)/2 + mpf('0.5')*_dig(s/2)
    return E*(dlogE*zeta(s) + diff(zeta, s))   # E' zeta + E zeta'

# [1] T-VI: (1+i)-power = antisymmetric exchange (all operands mpmath)
for z in [mpc(2,0), mpc(1,1), mpc(mpf('0.5'),3)]:
    a = z**(1+1j); s, th = log(abs(z)), arg(z)
    chk("T-VI exchange", abs(a - exp((s-th)+1j*(s+th))) < mpf('1e-30'), f"z={z}")

# [2] seam atom
c2 = sum(exp(-2j*t) for t in [pi/2, 3*pi/2])/2
chk("seam atom c2 = -1", abs(c2 + 1) < mpf('1e-30'), f"c2={mp.nstr(c2,5)}")

# [3] V4 convolution closure (deduped)
raw = sorted((a+b) % (2*pi) for a in [0, pi] for b in [0, pi/2, pi, 3*pi/2])
pts = []
for p in raw:
    if not any(abs(p-q) < mpf('1e-20') for q in pts): pts.append(p)
ref = [0, pi/2, pi, 3*pi/2]
chk("V4 convolution closure", len(pts) == 4 and all(abs(p-r) < mpf('1e-20') for p, r in zip(pts, ref)))

# [4] archiving rate closed form + fixed point
lam = pi/(4*log(2)) - mpf('0.5')
chk("lambda closed form", abs(lam - mpf('0.633090035456798')) < mpf('1e-15'), mp.nstr(lam, 15))
a = pi/(2*log(2))
chk("fixed-point identity 2^(a/2) e^{-pi/4} = 1", abs(2**(a/2)*exp(-pi/4) - 1) < mpf('1e-30'))

# [5] remainder closed form vs term-by-term
for n in [10, 20, 32]:
    Rn = 1j*gammainc(n+1, 0, 1j*pi/2)/gamma(n+1)
    ser = sum((1j*pi/2)**k/gamma(k+1) for k in range(n+1, 80))
    chk(f"R_{n} closed form", abs(Rn-ser) < mpf('1e-25'), f"|R|={mp.nstr(abs(Rn),3)}")

# [6] i^i tower
def tower(z0, n=300):
    z = z0
    for _ in range(n): z = (1j)**z
    return z
zstar = tower(mpf(1)*1j)
res = abs((1j)**zstar - zstar)
chk("tower fixed point (Newton plateau, honestly labeled)", res < mpf('1e-6'), f"plateau={mp.nstr(res,2)}")
chk("tower attractor rate |(i pi/2) z*| < 1", abs((1j*pi/2)*zstar) < 1, mp.nstr(abs((1j*pi/2)*zstar), 4))

# [7] theta-cancellation three-way (theta = Im logGamma - (t/2) ln pi)
th = lambda t: loggamma(mpf(1)/4 + 1j*t/2).imag - (t/2)*log(pi)
zp = diff(zeta, rho)
E = lambda s: mpf('0.5')*s*(s-1)*pi**(-s/2)*gamma(s/2)
a_zp, a_E = arg(zp), arg(E(rho))
chk("L6 pin-theta lock", abs(a_zp - (-th(G1) - pi/2)) < mpf('1e-8'), mp.nstr(a_zp, 10))
chk("arg E = pi + theta (mod 2pi)", moddist(a_E - pi - th(G1), 2*pi) < mpf('1e-8'), mp.nstr(a_E, 10))
chk("theta cancels: arg xi' = pi/2", moddist(a_zp + a_E - pi/2, pi) < mpf('1e-8'), mp.nstr(a_zp+a_E, 10))

# [8] two-arm certification via exact parts derivative (no cancellation)
for s0 in [mpf(-2), mpf(-4)]:
    v = diff(xi, s0)   # xi analytic at trivial zeros (zeta cancels the Gamma pole)
    chk(f"Arm A xi'({s0}) real", abs(v.imag) < mpf('1e-20'), mp.nstr(v.imag, 3))
v = xip_parts(mpf(3))
chk("Arm A xi'(3) real", abs(v.imag) < mpf('1e-15'), mp.nstr(v.imag, 3))
for t0 in [mpf(14), mpf('30.5')]:
    s0 = mpf('0.5') + t0*1j
    v = xip_parts(s0)
    nrm = max(abs(v), mpf('1e-30'))
    chk(f"Arm B xi'(0.5+{t0}i) pure imaginary", abs(v.real)/nrm < mpf('1e-10'), f"ratio={mp.nstr(abs(v.real)/nrm,2)} (measured; threshold calibrated, A4)")
# zero-point anchor
v = E(rho)*zp
chk("zero: xi'(rho1) = E zeta'", abs(v - xip_parts(rho)) < mpf('1e-30'), mp.nstr(abs(v - xip_parts(rho)), 2) + " (floor |E dlogE zeta(rho)|)")

# [9] FE residue at proper precision
chi = lambda s: 2**s*pi**(s-1)*sin(pi*s/2)*gamma(1-s)
lhs, rhs = diff(zeta, 1-rho), -diff(zeta, rho)/chi(rho)
chk("FE residue (dps=40)", abs(lhs-rhs) < mpf('1e-20'), mp.nstr(abs(lhs-rhs), 3))

print(f"\n{sum(results)}/{len(results)} PASS")
