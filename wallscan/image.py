"""Range profiles -> standoff per position -> wall-echo removal -> backprojection through the wall layer -> eps autofocus."""
import numpy as np
from .sim import C, FREQS, path_extra

PAD = 16
DF = FREQS[1] - FREQS[0]
K = 2 * np.pi * FREQS / C


def profiles(S):
    P = np.fft.ifft(S * np.hanning(len(FREQS)), n=len(FREQS) * PAD, axis=1)
    tau = np.arange(P.shape[1]) / (len(FREQS) * PAD * DF)
    return P, tau


def _peak(mag, i):
    """Parabolic sub-sample peak around index i."""
    if 0 < i < len(mag) - 1:
        a, b, c = mag[i - 1], mag[i], mag[i + 1]
        den = a - 2 * b + c
        return i + (0.5 * (a - c) / den if den else 0.0)
    return float(i)


def standoff(S, h_nominal):
    """Per-position standoff from the front-face echo: the first strong peak near the median profile's front echo."""
    P, tau = profiles(S)
    mag = np.abs(P)
    med = np.median(mag, axis=0)
    lim = int(2 * (h_nominal + 0.02) / C / tau[1])
    i0 = int(np.argmax(med[:lim]))
    w = int(2 * 0.012 / C / tau[1])                     # +-12 mm of travel around it
    h = np.empty(len(S))
    for k in range(len(S)):
        lo = max(0, i0 - w)
        j = lo + int(np.argmax(mag[k, lo:i0 + w]))
        h[k] = _peak(mag[k], j) * tau[1] * C / 2
    return h


def optical_thickness(S, h):
    """sum(n_i d_i) from the delay between the front- and back-face echoes of the aligned mean profile."""
    A = (S * np.exp(2j * K[None, :] * (h[:, None] - h.mean()))).mean(axis=0)
    P, tau = profiles(A[None, :])
    mag = np.abs(P[0])
    i_f = int(np.argmax(mag[: int(2 * (h.mean() + 0.02) / C / tau[1])]))
    lo, hi = i_f + int(2 * 0.006 / C / tau[1]), i_f + int(2 * 0.12 / C / tau[1])
    m = mag[lo:hi]
    peaks = [i for i in range(1, len(m) - 1) if m[i] > m[i - 1] and m[i] >= m[i + 1]]   # skip the front echo's skirt
    i_b = lo + max(peaks, key=lambda i: m[i])
    return (_peak(mag, i_b) - _peak(mag, i_f)) * tau[1] * C / 2


def clean(S, h):
    """Remove the wall's own echoes: align every position on its front face, subtract the mean, undo the shift."""
    ph = np.exp(2j * K[None, :] * (h[:, None] - h.mean()))
    A = S * ph
    return (A - A.mean(axis=0)) / ph


def backproject(S, xs, h, layers, X, Z):
    P, tau = profiles(S)
    I = np.zeros(X.shape, complex)
    t1 = tau[1]
    for k in range(len(xs)):
        R, extra = path_extra(xs[k], h[k], X, Z, layers)
        t = 2 * (R + extra) / C
        idx = t / t1
        i0 = np.clip(idx.astype(int), 0, P.shape[1] - 2)
        f = idx - i0
        v = P[k, i0] * (1 - f) + P[k, i0 + 1] * f
        I += v * np.exp(2j * np.pi * FREQS[0] * t)
    return I


def sharpness(I):
    a = np.abs(I) ** 2
    return float((a ** 2).sum() / a.sum() ** 2)


def form(scan, method, h_nominal=0.025, n_nominal=np.sqrt(2.8), grid=(0.002, 0.001), zmax=0.12):
    """method: 'free_space' | 'nominal_eps' | 'this'. Returns image, axes, and what was estimated."""
    xs, S = scan["x"], scan["S"]
    X1 = np.arange(xs[0], xs[-1] + 1e-9, grid[0]); Z1 = np.arange(0.0, zmax, grid[1])
    X, Z = np.meshgrid(X1, Z1)
    if method == "free_space":
        h = np.full(len(xs), h_nominal)
        Sc = S - S.mean(axis=0)
        return dict(I=backproject(Sc, xs, h, [], X, Z), x=X1, z=Z1, h=h, n=1.0, D=None)
    h = standoff(S, h_nominal)
    Sc = clean(S, h)
    opt = optical_thickness(S, h)
    if method == "nominal_eps":
        n = n_nominal
    else:                                             # autofocus on a coarse grid behind the wall, then refine
        Xc, Zc = np.meshgrid(X1[::3], Z1[::2])
        def score(n):
            return sharpness(backproject(Sc, xs, h, [(opt / n, n * n)], Xc, Zc) * (Zc > opt / n + 0.002))
        cands = np.arange(1.3, 3.21, 0.15)
        best = max(cands, key=score)
        fine = np.arange(best - 0.12, best + 0.121, 0.04)
        n = float(max(fine, key=score))
    D = opt / n
    return dict(I=backproject(Sc, xs, h, [(D, n * n)], X, Z), x=X1, z=Z1, h=h, n=n, D=D)
