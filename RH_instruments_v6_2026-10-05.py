# RH_instruments_v6_2026-10-05.py · Nine-section instrument suite A–I (five archive operations layer)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Covers: layer ladder; configuration issuing + data-wall law; weld-parity readout; m-discrimination
#         battery; convolution readout; H-condition vacuity; FE + seam faces; tower + Euclid increment; face book.
# Usage: python RH_instruments_v6_2026-10-05.py  (local only; ~1 minute)
# RH_instruments_v6_2026-10-05.py
# 仪器箱 v6 · 可更新版 —— 无中生有 + 九面接通时代的全部读出器
# 依赖: numpy, mpmath；零联网。各节独立，可单独更新。
# 作者署名: Guang Yang / Yueting Xiao（分析协助: Kimi）

import numpy as np
import mpmath as mp
from mpmath import mpf, gamma, pi

DPS = 25
mp.mp.dps = DPS

def carg(z): return float(mp.atan2(z.imag, z.real))
def E(s): return mpf('0.5')*s*(s-1)*mp.pi**(-s/2)*gamma(s/2)
def hfun(s): return E(s)*mp.diff(mp.zeta, s)   # h = E·ζ′；零点处 = ξ′

# ============ A. 主读出器：p̄（周期/判别式）============
def period(ρ):
    """p̄(ρ) = Arg ξ′(ρ) mod π ∈ [0, π)。落档 ⟺ p̄ = π/2。"""
    return carg(hfun(ρ)) % float(np.pi)

def is_registered(ρ, tol=1e-9):
    return abs(period(ρ) - np.pi/2) < tol

# ============ B. i^i 构型通道 e^{*(w·g)}（w ∈ ℂ）============
def config_channels(w, K=(1,2,3,4)):
    """ĉ_k = e^{w·e^{−ikπ/2}}；w = i·t 为实周期形，w = z* 为复周期形。"""
    return {k: np.exp(w*np.exp(-1j*k*np.pi/2)) for k in K}

ARCHIVE_W = 1j*np.pi/2   # 档案落点：偶通道 = [+0i, −0i]

# ============ C. 典范递归（对坐标形，α = π/(2ln2)）============
ALPHA = np.pi/(2*np.log(2))
def canon_step(p):
    z = 1 + np.exp(1j*p)
    return (np.arctan2(z.imag, z.real) + ALPHA*np.log(abs(z))) % np.pi

def basin_test(p0, steps=3000, tol=1e-13):
    p = p0
    for n in range(steps):
        pn = canon_step(p)
        if abs(pn - p) < tol: return p, n+1
        p = pn
    return p, steps

# ============ D. 幂塔 z ← i^z（凝聚测试）============
TOWER_Z = 0.438283 + 0.360592j   # 全局吸引子（吸引性 |(iπ/2)z*| = 0.8915 < 1）
def tower(z0, steps=500, tol=1e-13):
    z = z0
    for n in range(steps):
        zn = 1j**z
        if abs(zn - z) < tol or abs(zn) > 1e6: return zn, n+1
        z = zn
    return z, steps

# ============ E. 共谋乘积（无中生有第一定理演示）============
def conspiracy_demo():
    """纤维 A（±π/2 锁）× 纤维 B（{0,π} 锁，同步 parity）→ 定居 +0i。"""
    A = np.array([np.pi/2, -np.pi/2, np.pi/2, -np.pi/2])
    B = np.array([0.0, np.pi, 0.0, np.pi])
    return (A + B) % (2*np.pi)   # 恒 = π/2

# ============ F. 形因子 E(ξ) = |Σ e^{iγξ}|²（素频峰）============
def form_factor(gammas, xi):
    return abs(np.sum(np.exp(1j*xi*np.asarray(gammas))))**2

# ============ G. 面字 Nyquist/DC + 滑窗缺陷定位 ============
def face_word(word):
    """word: ±1 序列 → (Nyquist 幅, DC 泄漏)。完美交替 = (1, 0)。"""
    w = np.asarray(word, dtype=float)
    return abs(np.mean(w * (-1)**np.arange(len(w)))), abs(w.mean())

def sliding_nyquist(word, W=20):
    s = (-1)**np.arange(W)
    return np.array([abs(np.mean(word[i:i+W]*s)) for i in range(len(word)-W+1)])

# ============ H. 终点分类（Haar 子群支撑大小）============
def endpoint_support(support, mod, N=240):
    """原子测度自卷积的 Cesàro 终点支撑（单位 mod；整数钉）。"""
    from fractions import Fraction
    from collections import defaultdict
    mu = defaultdict(Fraction, {a: Fraction(1, len(support)) for a in support})
    acc = defaultdict(Fraction); cur = mu
    for n in range(1, N+1):
        for a, w in cur.items(): acc[a % mod] += w/N
        nxt = defaultdict(Fraction)
        for a, w in cur.items():
            for b, v in mu.items(): nxt[(a+b) % mod] += w*v
        cur = nxt
    return {a: float(w) for a, w in acc.items() if w > 1e-9}

# ============ I. 锚定自检 ============
def self_check(K=8):
    ok = True
    for k in range(1, K+1):
        ρ = mp.zetazero(k)
        p = period(ρ)
        dev = abs(p - np.pi/2)
        ok = ok and dev < 1e-9
        print(f"  [A] ρ{k}: p̄ = {p:.12f}，偏差 {dev:.2e}")
    w = ARCHIVE_W
    ch = config_channels(w)
    evens = [ch[2], ch[4]]
    arch_ok = abs(evens[0] + 1j) < 1e-9 and abs(evens[1] - 1j) < 1e-9
    print(f"  [B] w* 偶通道 = {evens[0]:.3f}, {evens[1]:.3f} → 双档案 {arch_ok}")
    z, n = tower(0.5 + 0.5j)
    print(f"  [D] 塔(0.5+0.5i) {n} 步 → {z:.6f}（z* = {TOWER_Z:.6f}）")
    p, n2 = basin_test(0.7)
    print(f"  [C] 递归 basin(0.7) {n2} 步 → {p:.8f}（π/2 = {np.pi/2:.8f}）")
    cp = conspiracy_demo()
    print(f"  [E] 共谋乘积钉 = {np.round(cp, 4)}（恒 π/2 = +0i）")
    print("ALL PASS" if ok and arch_ok else "CHECK FAILED")
    return ok and arch_ok

if __name__ == "__main__":
    self_check()
