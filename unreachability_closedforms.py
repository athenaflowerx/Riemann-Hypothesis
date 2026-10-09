# unreachability_closedforms.py · Unreachability transformation — the four closed forms (five archive operations)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Covers: generating series G_n = Σ (iπ/2)^k/k! vs closed form; remainder R_n = i·γ(n+1, iπ/2)/Γ(n+1)
#         (|R_32| = 3.41e−31); recursion rate λ = π/(4 ln 2) − 1/2; basin rates.
# Usage: python unreachability_closedforms.py  (local only, mpmath required; ~1 second)
# 不可达转化 + 四闭式 验证代码（零联网，numpy/mpmath）
import numpy as np
import mpmath as mp
from mpmath import mpf, gamma, pi, zetazero
mp.mp.dps = 30
ALPHA = np.pi/(2*np.log(2))

def carg(z): return float(mp.atan2(z.imag, z.real))
def E(s): return mpf('0.5')*s*(s-1)*mp.pi**(-s/2)*gamma(s/2)
def hfun(s): return E(s)*mp.diff(mp.zeta, s)          # h = E·ζ′；零点处 = ξ′

def interlock():                                       # ① 互锁（无×无 = +0i）
    A = np.array([np.pi/2, -np.pi/2, np.pi/2, -np.pi/2])
    B = np.array([0.0, np.pi, 0.0, np.pi])
    return (A + B) % (2*np.pi)

TOWER_Z = 0.438283 + 0.360592j                         # ② 凝聚（幂塔 → z*）
def tower(z0, steps=3000, tol=1e-12):
    z = complex(z0)
    for _ in range(steps):
        zn = 1j**z
        if abs(zn-z) < tol or abs(zn) > 1e6: return zn
        z = zn
    return z

def fbar(p):                                           # ④ 迭代收敛（线程归档）
    z = 1 + np.exp(1j*p)
    return (np.arctan2(z.imag, z.real) + ALPHA*np.log(abs(z))) % np.pi
def archive(p0, steps=50):
    p = p0
    for _ in range(steps): p = fbar(p)
    return abs(p - np.pi/2)

def R_closed(n):                                       # 闭式 1：余项（下不完全伽马）
    return 1j*mp.gammainc(n+1, 0, 1j*pi/2, regularized=True)

LAMBDA = np.pi/(4*np.log(2)) - 0.5                     # 闭式 2：递归率
def budget(d0): return abs(d0)/(1-abs(LAMBDA))         # 闭式 3：预算

def dev_derived(k, N, c=0.845):                        # 闭式 4：数据墙关系形
    gk = float(mp.im(zetazero(k))); gN = float(mp.im(zetazero(N)))
    return c*gk*np.log(gN)/(2*np.pi*gN)

def dspec_point(s):                                    # 谱形：缺陷谱
    p = carg(hfun(s)) % float(np.pi)
    return complex(np.exp(1j*(np.pi/2 - p)))

if __name__ == "__main__":
    assert np.allclose(interlock(), np.pi/2)
    for z0 in [0.5+0.5j, 1j, 2, 1+2j]:
        assert abs(tower(z0) - TOWER_Z) < 1e-5  # float64 停滞地板
    for p0 in [1e-9, 0.3, 0.7, 1.0, 1.4, np.pi-1e-9]:
        assert archive(p0) < 1e-9  # float64 下 0.633^50≈4e-11 量级
    n = 8
    assert abs(R_closed(n) - sum((1j*pi/2)**k/mp.factorial(k) for k in range(n+1,300))) < 1e-12
    assert abs(abs(LAMBDA) - 0.633090) < 1e-5
    assert abs(budget(1.0) - 1/(1-0.633090)) < 1e-3
    for k in [1, 2, 3]:
        assert abs(dspec_point(zetazero(k)) - 1) < 1e-8
    print("UNREACHABLE + CLOSED FORMS: ALL PASS")
