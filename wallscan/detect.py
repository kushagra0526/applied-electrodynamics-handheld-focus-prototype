"""Image -> reflectors (x, depth) -> "safe to drill here?" with a stated clearance."""
import numpy as np

CLEAR = 0.010          # a drill needs this much side clearance from anything behind the wall (m)
MARGIN = 0.020         # half-width we assume for a detected reflector we can't size (pipe radius ~11 mm + slack)
DEPTH_MARGIN = 0.003   # stated depth uncertainty added to every reflector
DEPTHS = (0.025, 0.035, 0.045, 0.06)   # hole depths checked at every drill point (anchor, box, deep anchor, through-stud bit)


def peaks(img, k_noise=2.5, sep_x=0.02):
    """Strongest column maxima (over depth), at least sep_x apart, above k_noise x the median column (the image's floor)."""
    I = np.abs(img["I"])
    I = I * (img["z"][:, None] > (img["D"] or 0) + 0.002)      # behind the wall layer only
    col = I.max(axis=0)
    zi = I.argmax(axis=0)
    thr = k_noise * np.median(col)
    out = []
    for j in np.argsort(col)[::-1]:
        if col[j] < thr:
            break
        x = img["x"][j]
        if any(abs(x - p["x"]) < sep_x for p in out):
            continue
        out.append(dict(x=float(x), z=float(img["z"][zi[j]]), amp=float(col[j] / col.max())))
    return out


def drill_check(found, x, depth):
    """NOT_CLEAR if any reflector could be within CLEAR of the drill point and shallower than the hole (+ margin)."""
    near = [p for p in found if abs(p["x"] - x) <= CLEAR + MARGIN and p["z"] < depth + DEPTH_MARGIN]
    return ("NOT_CLEAR", near) if near else ("CLEAR", [])


def score(wall, found, length):
    """Per-target localisation error, misses, false reflectors, and drill-point outcomes."""
    T = wall["targets"]
    used, errs = set(), []
    for t in T:
        c = [(abs(p["x"] - t["x"]), i) for i, p in enumerate(found) if abs(p["x"] - t["x"]) <= t["half"] + 0.02 and i not in used]
        if not c:
            errs.append(dict(kind=t["kind"], found=False))
            continue
        _, i = min(c)
        used.add(i)
        p = found[i]
        errs.append(dict(kind=t["kind"], found=True, dx=max(0.0, abs(p["x"] - t["x"]) - (t["half"] if t["kind"] == "stud" else 0)),
                         dz=p["z"] - t["z"]))
    false = len(found) - len(used)
    danger = false_stop = n = unsafe = 0
    for x in np.arange(0.03, length - 0.03, 0.005):
        for depth in DEPTHS:
            truth = any(abs(x - t["x"]) <= t["half"] + CLEAR and t["z"] < depth for t in T)
            call, _ = drill_check(found, x, depth)
            n += 1; unsafe += truth
            danger += truth and call == "CLEAR"
            false_stop += (not truth) and call == "NOT_CLEAR"
    return dict(targets=errs, false=false, drill_points=n, dangerous_clear=danger, false_not_clear=false_stop, unsafe_points=unsafe)
