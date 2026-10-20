"""docs/result.png: the demo wall imaged three ways, with true target positions."""
import sys
sys.path.insert(0, ".")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wallscan import sim, image

case = sim.make(2)
w = case["wall"]
fig, ax = plt.subplots(3, 1, figsize=(11, 7.5), sharex=True)
for a, (m, title) in zip(ax, [("free_space", "free space, nominal standoff"), ("nominal_eps", "standoff fixed, drywall eps"),
                              ("this", "standoff fixed + eps autofocus (this)")]):
    img = image.form(case["scan"], m)
    I = np.abs(img["I"]); I /= I.max()
    a.imshow(20 * np.log10(I + 1e-3), extent=[0, 1000, 120, 0], aspect="auto", vmin=-25, vmax=0, cmap="magma")
    for t in w["targets"]:
        a.plot(t["x"] * 1000, t["z"] * 1000, "c+", ms=12, mew=2)
    a.axhline(w["thickness"] * 1000, color="w", lw=0.6, ls="--")
    a.set_title(f"{title}  (n = {img['n']:.2f}; true {sum(d * np.sqrt(e) for d, e in w['layers']) / w['thickness']:.2f})", fontsize=10)
    a.set_ylabel("depth (mm)")
ax[-1].set_xlabel("position along wall (mm)   + = true stud / pipe front face, dashed = back of tile + cement board")
plt.tight_layout(); plt.savefig("docs/result.png", dpi=100)
