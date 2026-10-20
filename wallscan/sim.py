"""2-D wall cross-section scanned by a handheld stepped-frequency radar.

Coordinates in metres: x along the wall, z into the wall (z = 0 is the wall's front face).
The antenna sits at (x_k, -h_k): x_k from the hand's path, h_k the standoff, which wobbles.
Delay model (straight ray): two-way delay = 2/c * (R + sum_i (n_i - 1) * L_i), L_i = the ray's length inside layer i.
"""
import numpy as np

C = 299_792_458.0
FREQS = np.linspace(4e9, 20e9, 161)               # 16 GHz bandwidth, 100 MHz steps
WALLS = {   # name: [(thickness m, eps_r)], front to back
    "drywall_1/2": [(0.0127, 2.8)],
    "drywall_5/8": [(0.0159, 2.6)],
    "plaster_lath": [(0.019, 4.0)],
    "tile_on_cement_board": [(0.008, 6.5), (0.0127, 7.5)],
    "drywall_double": [(0.0254, 2.7)],
}
DEV = dict(walls=["drywall_1/2", "drywall_5/8", "tile_on_cement_board"], eps_jitter=0.15, h0=0.025, wobble_mm=4.0,
           length=1.0, step=0.004, snr_db=25.0)


def path_extra(xa, ha, x, z, layers):
    """Sum over layers of (n-1) * length inside the layer, straight ray from (xa, -ha) to (x, z)."""
    R = np.hypot(x - xa, z + ha)
    sec = R / (z + ha)
    extra, top = 0.0, 0.0
    for d, eps in layers:
        inside = np.clip(z - top, 0, d)
        extra = extra + (np.sqrt(eps) - 1) * inside * sec
        top += d
    return R, extra


def wall(rng, cfg):
    name = rng.choice(cfg["walls"])
    layers = [(d * rng.uniform(0.95, 1.05), e * (1 + rng.uniform(-cfg["eps_jitter"], cfg["eps_jitter"]))) for d, e in WALLS[name]]
    D = sum(d for d, _ in layers)
    T = []          # targets: kind, x, z of the front surface, half-width, amplitude
    x = rng.uniform(0.05, 0.25)
    while x < cfg["length"] - 0.03:
        T.append(dict(kind="stud", x=x, z=D, half=0.019, amp=0.25))
        x += 0.406
    for _ in range(rng.integers(1, 4)):
        kind = rng.choice(["copper", "pex", "cable"])
        xx = rng.uniform(0.04, cfg["length"] - 0.04)
        if any(abs(xx - t["x"]) < t["half"] + 0.03 for t in T):
            continue
        r = {"copper": 0.0111, "pex": 0.0111, "cable": 0.006}[kind]
        zf = D + rng.uniform(0.004, 0.06)
        T.append(dict(kind=kind, x=xx, z=zf, half=r, amp={"copper": 1.0, "pex": 0.18, "cable": 0.35}[kind]))
    return dict(name=name, layers=layers, thickness=D, targets=T)


def scan(rng, w, cfg):
    xs = np.arange(0, cfg["length"] + 1e-9, cfg["step"])
    t = xs / cfg["length"]
    h = cfg["h0"] + cfg["wobble_mm"] / 1000 * sum(rng.normal(0, 0.6) * np.sin(2 * np.pi * (k * t + rng.random())) for k in (1, 2, 3, 5))
    S = np.zeros((len(xs), len(FREQS)), complex)
    k = 2 * np.pi * FREQS / C
    n1 = np.sqrt(w["layers"][0][1])
    G1 = (1 - n1) / (1 + n1)
    opt = sum(d * np.sqrt(e) for d, e in w["layers"])
    for i, (xa, ha) in enumerate(zip(xs, h)):
        S[i] += G1 * np.exp(-2j * k * ha)                               # front face (specular)
        S[i] += -0.6 * G1 * np.exp(-2j * k * (ha + opt))               # back face
        for T in w["targets"]:
            pts = [T["x"]] if T["kind"] != "stud" else np.linspace(T["x"] - T["half"], T["x"] + T["half"], 7)
            for px in pts:
                R, extra = path_extra(xa, ha, px, T["z"], w["layers"])
                cos = (T["z"] + ha) / R
                a = T["amp"] / len(pts) * cos ** 4 * (0.05 / R)       # antenna pattern + spreading
                S[i] += a * np.exp(-2j * k * (R + extra))
    sig = np.sqrt(np.mean(np.abs(S) ** 2)) * 10 ** (-cfg["snr_db"] / 20)
    S += sig * (rng.normal(size=S.shape) + 1j * rng.normal(size=S.shape)) / np.sqrt(2)
    return dict(x=xs, h_true=h, S=S)


def make(seed, cfg=None):
    c = {**DEV, **(cfg or {})}
    rng = np.random.default_rng(seed)
    w = wall(rng, c)
    return dict(wall=w, scan=scan(rng, w, c), cfg=c)
