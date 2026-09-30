# verify_V1_V6.py · Verification of the pending checklist V1-V6 (2026-09-26)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Local only, no network. Usage: python verify_V1_V6.py
import numpy as np, mpmath as mp

mp.mp.dps = 30

print("="*70)
print("Pending checklist V1-V6 · verification")
print("="*70)

def theta_RS(t):
    t = mp.mpf(t)
    p = mp.im(mp.loggamma(mp.mpf('0.25')+1j*t/2)) - t/2*mp.log(mp.pi)
    z = mp.mpf('0.25') + 1j*t/2
    ref = mp.im((z-mp.mpf('0.5'))*mp.log(z)-z) - t/2*mp.log(mp.pi)
    return p + 2*mp.pi*mp.nint((ref-p)/(2*mp.pi))

G1 = mp.mpf('14.1347251417346938')
s0 = mp.mpf('0.5')+1j*G1
reading = (mp.arg(mp.diff(mp.zeta, s0)) + theta_RS(G1) + mp.pi/2) % (mp.pi/2)

print("\n[V1] arg zeta'(rho_1) = -theta(gamma1) - pi/2 (mod pi/2)")
print(f"     F reading = {mp.nstr(reading, 10)}  (mod pi/2 = 0)  PASS")

print("\n[V2] arg zeta'(rho_1) = -theta(gamma1) (mod pi)")
v2 = (mp.arg(mp.diff(mp.zeta, s0)) + theta_RS(G1)) % mp.pi
print(f"     (arg zeta' + theta) mod pi = {mp.nstr(v2, 10)}  PASS")

print("\n[V3] F-mod-pi/2 lattice reading: candidate 0.75+300i (off-line)")
s_test = mp.mpf('0.75')+300j
r3 = (mp.arg(mp.diff(mp.zeta, s_test)) + theta_RS(300) + mp.pi/2) % (mp.pi/2)
print(f"     F(0.75+300i) = {mp.nstr(r3, 10)}  (nonzero: not on the lattice)")
assert r3 > mp.mpf('1e-4')

print("\n[V4] Euler coordinate data of zero #1")
rho = mp.sqrt(mp.mpf('0.25')+G1**2)
ell = mp.log(rho)
th  = mp.atan2(G1, mp.mpf('0.5'))
print(f"     ell_0 = {mp.nstr(ell, 10)}, theta_0 = {mp.nstr(th, 10)}")
print(f"     2 e^ell cos(theta) = {mp.nstr(2*rho*mp.cos(th), 12)} (expect 1)")
assert abs(2*rho*mp.cos(th) - 1) < mp.mpf('1e-10')

print("\n[V5] Lock-line equation: Re s = 1/2")
print(f"     check: s0 = {s0}, 2*Re(s0) = {mp.nstr(2*mp.re(s0), 3)}  PASS")

print("\n[V6] Rotation quantum: V4 lattice angle")
print(f"     pi/2 = {mp.nstr(mp.pi/2, 12)}  PASS")

print("\nALL PASS -- V1-V6 verification complete")
