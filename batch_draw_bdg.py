"""Batch draw BDG cross-link spectra with timing."""
import sys, os, time, subprocess

PYTHON_EXE = "python3"  # Python 3.13.14
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BDG_DIR    = r"D:\MSdata\260729-AdK\20260512\BDG"
RAW_DIR    = r"D:\MSdata\260729-AdK\20260512\raw"
OUT_BASE   = os.path.join(BDG_DIR, "cross-link_png")
MAIN_PY    = os.path.join(SCRIPT_DIR, "main.py")

plabel_files = [
    f for f in os.listdir(BDG_DIR)
    if f.endswith(".cross-linked.BDG-H.plabel")
]
plabel_files.sort()

print(f"Found {len(plabel_files)} cross-link .plabel files")
print(f"Python: {PYTHON_EXE}")
print(f"{'='*70}")

total_spectra_drawn = 0
total_time = 0.0
results = []

for i, plabel_name in enumerate(plabel_files, 1):
    base_name = plabel_name.replace(".cross-linked.BDG-H.plabel", "")
    mgf_name = base_name + ".mgf"
    plabel_path = os.path.join(BDG_DIR, plabel_name)
    mgf_path = os.path.join(RAW_DIR, mgf_name)
    out_path = os.path.join(OUT_BASE, base_name)

    if not os.path.exists(mgf_path):
        print(f"[{i}/{len(plabel_files)}] SKIP {base_name}: MGF not found")
        continue

    mgf_size_mb = os.path.getsize(mgf_path) / (1024 * 1024)

    t0 = time.time()
    try:
        result = subprocess.run(
            [PYTHON_EXE, MAIN_PY,
             "--mgf", mgf_path,
             "--ident", plabel_path,
             "--parser", "plink",
             "--out", out_path],
            capture_output=True, text=True, timeout=1800,
            cwd=SCRIPT_DIR,
        )
        elapsed = time.time() - t0
    except subprocess.TimeoutExpired:
        print(f"[{i}/{len(plabel_files)}] TIMEOUT {base_name}")
        continue

    n_drawn = 0
    for line in result.stdout.splitlines():
        if "spectra drawn" in line:
            parts = line.strip().split()
            n_drawn = int(parts[1])

    total_spectra_drawn += n_drawn
    total_time += elapsed

    speed = n_drawn / elapsed if elapsed > 0 else 0
    status = "OK" if result.returncode == 0 else f"ERR({result.returncode})"
    results.append((base_name, mgf_size_mb, n_drawn, elapsed, speed, status))

    print(f"[{i:2d}/{len(plabel_files)}] {status} {base_name}: "
          f"MGF={mgf_size_mb:.0f}MB, drawn={n_drawn}, "
          f"{elapsed:.1f}s ({speed:.1f} spectra/s)")

    if result.returncode != 0:
        stderr_preview = result.stderr.strip()
        if stderr_preview:
            lines = stderr_preview.splitlines()
            print(f"  stderr: {' | '.join(lines[:3])}")

print(f"\n{'='*70}")
print(f"SUMMARY")
print(f"{'='*70}")
print(f"{'Name':<35s} {'MGF(MB)':>8s} {'Drawn':>6s} {'Time(s)':>8s} {'Sp/s':>7s} {'Status':>6s}")
print(f"{'-'*75}")
for name, size, drawn, elapsed, speed, status in results:
    print(f"{name:<35s} {size:>8.0f} {drawn:>6d} {elapsed:>8.1f} {speed:>7.1f} {status:>6s}")

avg_speed = total_spectra_drawn / total_time if total_time > 0 else 0
print(f"{'-'*75}")
print(f"{'TOTAL':<35s} {'':>8s} {total_spectra_drawn:>6d} {total_time:>8.1f} {avg_speed:>7.1f}")
print(f"\nAverage speed: {avg_speed:.1f} spectra/second")
print(f"Total time: {total_time/60:.1f} minutes")
print(f"Total spectra drawn: {total_spectra_drawn}")
