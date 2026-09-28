# extrapolation_crossing.py · Extrapolation of the gamma*(delta) vs height-quantum crossing zone
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Background: undiscriminated cone delta_cone = 1/(2 gamma) vs height-quantum angle
# delta_lad = 2pi/(gamma ln(gamma/2pi)). Crossing at gamma_x = 2pi e^{4pi} ~ 1.79e6
# (cone width first exceeds one quantum angle).
# This experiment:
#   1) analytic confirmation of the crossing point and N(gamma) = delta_cone/delta_lad
#      = ln(gamma/2pi)/(4pi) (number of quantum cells spanned by the cone);
#   2) measurement of the seam slope A_k = Im(zeta''/zeta')(rho_k) on the first K zeros,
#      log-fit of A(gamma);
#   3) extrapolation of the B.3 validity window u_max = pi/(2|A(gamma)|) vs cone width --
#      verification that the seam theorem covers the quantum cone at all heights.
# Nature of conclusions: extrapolation (fit extrapolation), not a theorem; flagged.
import mpmath as mp

mp.mp.dps = 30
PI = mp.pi

# ---- 1) analytic crossing confirmation ----
gamma_cross = 2*PI*mp.e**(4*PI)
N_cross = mp.log(gamma_cross/(2*PI)) / (4*PI)
print("=== 1) Analytic crossing confirmation ===")
print(f"  gamma_x = 2pi e^{{4pi}} = {mp.nstr(gamma_cross, 12)}   N(gamma_x) = {mp.nstr(N_cross, 6)} (should be 1)")
assert abs(N_cross - 1) < mp.mpf('1e-12')

# ---- 2) A_k measurement + log fit ----
print("\n=== 2) Seam slope A_k = Im(zeta''/zeta')(rho_k): measurement and fit ===")
gam = [mp.mpf(x) for x in """14.134725141734693790457251983562470270784257115699
21.02203963877155499262847959389690277733446052441
25.01085758014568876321379099256282181865954867284
30.424876125859513210311897530584091320181560344715
32.935061587739189690662368964074903488627115757580
37.586178158825671257217763480705332821430908394174
40.918719012147495187398126914633254395726165962777
43.327073280914999519496122165406782742855460191259
48.005150881167159727942472749427516041686223281277
49.773832477672302181916784678563750745793044449321
52.970321477714460644147602580871262072190579322479
56.446247697063804881189356414275627551328203266708
59.347044002602353079653648674992219565108194603459
60.831778524609809844259171303647976559290381911155
65.112544048081606660875054253197472235841271393420
67.079810529494173714478539942303912864668962193310
69.546401711173979252926733596683667338517489381159
72.067157674481907582522107969826118361062994329191
75.704690699083933168326889784494934814858105158484
77.144840068874805372682664856304637010871228692001""".split()]

A = []
for g in gam:
    rho = mp.mpf('0.5') + 1j*g
    A.append(mp.im(mp.diff(mp.zeta, rho, 2) / mp.diff(mp.zeta, rho)))

xs = [mp.log(g/(2*PI)) for g in gam]
n = len(gam)
sx, sy = sum(xs), sum(A)
sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, A))
a = (n*sxy - sx*sy) / (n*sxx - sx*sx)
b = (sy - a*sx) / n
print(f"  fit A(gamma) = a ln(gamma/2pi) + b: a = {mp.nstr(a,8)}, b = {mp.nstr(b,8)}")
resid = max(abs(y - (a*x + b)) for x, y in zip(xs, A))
print(f"  max residual = {mp.nstr(resid,6)} (over 20 zeros)")

def A_ext(g):
    return a*mp.log(mp.mpf(g)/(2*PI)) + b

# ---- 3) crossing-zone extrapolation table ----
print("\n=== 3) Extrapolation table: radial-vs-angular dominance exchange and B.3 coverage window ===")
print(f"{'gamma':>10} | {'cone=1/2g':>13} | {'lad':>13} | {'N=cone/lad':>10} | {'|A|':>9} | {'u_max':>13} | {'u_max/cone':>13}")
for gs in ['1e2','1e4','1.79e6','1e8','1e12','1e30']:
    g = mp.mpf(gs)
    dc = 1/(2*g); dl = 2*PI/(g*mp.log(g/(2*PI)))
    N  = dc/dl
    Aa = abs(A_ext(g)); um = PI/(2*Aa)
    print(f"{gs:>10} | {mp.nstr(dc,6):>13} | {mp.nstr(dl,6):>13} | {mp.nstr(N,8):>10} | {mp.nstr(Aa,6):>9} | {mp.nstr(um,6):>13} | {mp.nstr(um/dc,6):>13}")

ok = abs(N_cross - 1) < mp.mpf('1e-12')
for gs in ['1e6','1e12','1e30']:
    g = mp.mpf(gs)
    ratio = (PI/(2*abs(A_ext(g)))) / (1/(2*g))
    print(f"  coverage ratio u_max/cone @ gamma={gs} = {mp.nstr(ratio,6)}")
    ok &= ratio > 1000
print("\n" + ("ALL PASS -- crossing-zone extrapolation complete (extrapolation, not a theorem)" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
