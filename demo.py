"""Tile on cement board, scanned by hand: where each method puts what's behind it, and a drill check."""
import numpy as np
from wallscan import sim, image, detect

case = sim.make(2)                       # dev wall: tile on cement board, 3 studs, PEX pipe
w = case["wall"]
n_true = sum(d * np.sqrt(e) for d, e in w["layers"]) / w["thickness"]
print(f"wall: {w['name']}, {w['thickness'] * 1000:.1f} mm, effective n {n_true:.2f}; "
      f"hand wobble {np.ptp(case['scan']['h_true']) * 1000:.1f} mm peak-to-peak over a 1 m sweep\n")
pipe = next(t for t in w["targets"] if t["kind"] != "stud")
names = {"free_space": "free space, nominal standoff", "nominal_eps": "standoff fixed, drywall eps", "this": "standoff fixed + eps autofocus"}
for m, label in names.items():
    img = image.form(case["scan"], m)
    found = detect.peaks(img)
    near = min(found, key=lambda p: abs(p["x"] - pipe["x"]), default=None)
    where = f"x {near['x'] * 1000:6.1f} mm, depth {near['z'] * 1000:5.1f} mm" if near and abs(near["x"] - pipe["x"]) < 0.03 else "not found"
    x_drill, depth = pipe["x"] + 0.012, 0.045
    call, _ = detect.drill_check(found, x_drill, depth)
    print(f"  {label:32s} n={img['n']:.2f}  {pipe["kind"]} -> {where}   45 mm hole 12 mm beside it: {call}")
print(f"\n  truth: {pipe["kind"]} at x {pipe['x'] * 1000:.1f} mm, front face {pipe['z'] * 1000:.1f} mm deep -> a 45 mm hole there hits it")
