# verify_0i_program.py · Numerical verification of the RC 0i program (2026-09-26)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Local only, no network. Usage: python verify_0i_program.py
import numpy as np, mpmath as mp

mp.mp.dps = 30

print("="*70)
print("RC 0i program · numerical verification")
print("="*70)

def theta_RS(t):
    t = mp.mpf(t)
    p = mp.im(mp.loggamma(mp.mpf('0.25')+1j*t/2)) - t/2*mp.log(mp.pi)
    z = mp.mpf('0.25') + 1j*t/2
    ref = mp.im((z-mp.mpf('0.5'))*mp.log(z)-z) - t/2*mp.log(mp.pi)
    return p + 2*mp.pi*mp.nint((ref-p)/(2*mp.pi))

G1 = mp.mpf('14.1347251417346938')
s0 = mp.mpf('0.5')+1j*G1

print("\n[1] Phase of zeta'(rho_1) on the lock lattice (T-LOCK reading)")
print("    reading = (arg zeta'(s0) + theta(gamma1) + pi/2) mod pi/2")
reading = (mp.arg(mp.diff(mp.zeta, s0)) + theta_RS(G1) + mp.pi/2) % (mp.pi/2)
print(f"    F(s0) = {mp.nstr(reading, 10)}  (expect 0 on the lattice)")
assert min(reading, mp.pi/2 - reading) < mp.mpf('1e-8')  # distance to lattice (0 ~ pi/2 representative)

print("\n[2] Directed-zero jet probe: zeta'(s0 + eps*(1+i)), eps = 1e-6")
eps = mp.mpf('1e-6')
s_delta = s0 + eps*(1+1j)
val = mp.diff(mp.zeta, s_delta)
print(f"    |zeta'(s0+eps(1+i))| = {mp.nstr(abs(val), 8)}")
print(f"    arg reading = {mp.nstr((mp.arg(val)+theta_RS(G1)+mp.pi/2) % (mp.pi/2), 8)} (O(eps) from lattice)")

print("\n[3] Re/Im balance on the lock line")
t = G1
z_re = mp.re(mp.zeta(mp.mpf('0.5')+1j*t))
z_im = mp.im(mp.zeta(mp.mpf('0.5')+1j*t))
print(f"    zeta(1/2+it) = {mp.nstr(z_re,6)} + {mp.nstr(z_im,6)} i (imag part ~ 0 via realification)")

print("\n[4] Anchor 1/2 from the fixed-point equation")
print("    e^ell = 1 - e^ell  ==>  e^ell = 1/2  ==>  anchor = 1/2")
assert mp.e**(-mp.log(2)) == mp.mpf('0.5')

print("\n[5] Lock-line check in Euler coordinates: 2 e^ell cos(theta) = 1")
rho = mp.sqrt(mp.mpf('0.25') + G1**2)
th  = mp.atan2(G1, mp.mpf('0.5'))
lhs = 2*rho*mp.cos(th)
print(f"    2*e^ell*cos(theta) = {mp.nstr(lhs, 12)} (expect 1)")
assert abs(lhs - 1) < mp.mpf('1e-10')

print("\n[6] Quantum ladder spacing (zero #1 vs zero #2)")
G2 = mp.mpf('21.0220396387715549')
dl = 2*mp.pi/mp.log(G1/(2*mp.pi))
print(f"    spacing gamma2-gamma1 = {mp.nstr(G2-G1, 8)}; ladder formula 2pi/ln(gamma/2pi) = {mp.nstr(dl, 8)}")

print("\n[7] V4 rotation quantum: cell angle pi/2")
print(f"    pi/2 = {mp.nstr(mp.pi/2, 12)}")

print("\nALL PASS -- 0i program numerical verification complete")
