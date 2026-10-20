"""Localisation and drill-check outcomes per imaging method.  python evaluate.py [dev|heldout] [out.json]"""
import json, sys, time
import numpy as np
from wallscan import sim, image, detect

METHODS = {"free_space": "free space, nominal standoff", "nominal_eps": "standoff-corrected, drywall eps",
           "this": "standoff-corrected + eps autofocus (this)"}


def run(which):
    if which == "heldout":
        h = json.load(open("heldout/config.json"))
        seeds, cfg = h["seeds"], h["cfg"]
    else:
        seeds, cfg = list(range(1, 16)), None
    R = {m: dict(walls=0, targets=0, found=0, dz=[], dx=[], false=0, drill=0, unsafe=0, danger=0, false_stop=0) for m in METHODS}
    per_wall = []
    for s in seeds:
        case = sim.make(s, cfg)
        w = case["wall"]
        row = dict(seed=s, wall=w["name"], n_true=sum(d * np.sqrt(e) for d, e in w["layers"]) / w["thickness"])
        for m in METHODS:
            img = image.form(case["scan"], m)
            found = detect.peaks(img)
            sc = detect.score(w, found, case["cfg"]["length"])
            r = R[m]
            r["walls"] += 1; r["targets"] += len(sc["targets"]); r["false"] += sc["false"]
            for t in sc["targets"]:
                if t["found"]:
                    r["found"] += 1; r["dz"].append(abs(t["dz"])); r["dx"].append(t["dx"])
            r["drill"] += sc["drill_points"]; r["unsafe"] += sc["unsafe_points"]
            r["danger"] += sc["dangerous_clear"]; r["false_stop"] += sc["false_not_clear"]
            row[m] = dict(n=img["n"], found=sum(t["found"] for t in sc["targets"]), targets=len(sc["targets"]),
                          danger=sc["dangerous_clear"], false=sc["false"])
        per_wall.append(row)
    for r in R.values():
        r["dz_median_mm"] = 1000 * float(np.median(r["dz"])) if r["dz"] else None
        r["dz_p90_mm"] = 1000 * float(np.percentile(r["dz"], 90)) if r["dz"] else None
        r["dx_median_mm"] = 1000 * float(np.median(r["dx"])) if r["dx"] else None
        del r["dz"], r["dx"]
    return dict(set=which, results=R, walls=per_wall)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "dev"
    t0 = time.time()
    r = run(which)
    print(f"{which} ({time.time() - t0:.0f}s)")
    for m, x in r["results"].items():
        print(f"  {METHODS[m]:42s} found {x['found']}/{x['targets']}  depth err median {x['dz_median_mm']:.1f} mm, p90 {x['dz_p90_mm']:.1f}  "
              f"lateral {x['dx_median_mm']:.1f} mm  false {x['false']}   drill checks called CLEAR over something {x['danger']}/{x['unsafe']}  "
              f"false NOT CLEAR {x['false_stop']}/{x['drill'] - x['unsafe']}")
    if len(sys.argv) > 2:
        json.dump(r, open(sys.argv[2], "w"), indent=1, default=float)
