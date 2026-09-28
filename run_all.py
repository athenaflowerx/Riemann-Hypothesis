# run_all.py · One-click regression for the RH discrimination package v3.0 (local, no network)
# Authors: Guang Yang (ORCID: 0000-0003-0599-2881), Yueting Xiao (ORCID: 0009-0002-0268-0655)
# Runs the main instruments sequentially plus legacy regression; any nonzero exit is reported.
import subprocess, sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
JOBS = [
    ("theta_RS_guarded.py",      "theta_RS three-guard self-test"),
    ("blindness_limits.py",      "cross-layer blindness-limit experiment"),
    ("rotation_scan.py",         "rotation-scan discriminant (best expression of T-LOCK)"),
    ("weld_detector.py",         "weld detector (unwrapped phase / argument principle)"),
    ("extrapolation_crossing.py","quantum crossing-zone extrapolation (gamma_x = 2pi e^{4pi})"),
    ("rh_discriminator.py",      "five-layer unified discriminator"),
    ("instrument_discipline.py", "instrument-layer discipline automation"),
    ("verify_0i_program.py",     "RC 0i program numerical verification"),
    ("verify_V1_V6.py",          "checklist V1-V6 verification"),
    ("legacy_v1.9/verify_merged_v5.py",        "[legacy] merged v5"),
    ("legacy_v1.9/rc_T126_T130_experiments.py","[legacy] T126-T130 experiments"),
]
fails = []
for f, desc in JOBS:
    t0 = time.time()
    r = subprocess.run([sys.executable, os.path.join(HERE, f)],
                       capture_output=True, text=True, timeout=900)
    dt = time.time() - t0
    status = "PASS" if r.returncode == 0 else "FAIL"
    print(f"[{status}] {desc}  ({f}, {dt:.0f}s)")
    if r.returncode != 0:
        fails.append(f)
        print(r.stdout[-800:]); print("STDERR:", r.stderr[-800:])
if fails:
    print("FAILED:", fails); sys.exit(1)
print("== all regression passed ==")
