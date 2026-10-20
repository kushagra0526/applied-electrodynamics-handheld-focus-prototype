import numpy as np
from wallscan import sim, image, detect


def test_layer_delay_adds_optical_path():
    R, extra = sim.path_extra(0.0, 0.02, 0.0, 0.05, [(0.0127, 4.0)])
    assert abs(R - 0.07) < 1e-12 and abs(extra - 0.0127) < 1e-12        # n - 1 = 1, normal incidence


def test_standoff_tracks_the_hand_to_sub_mm():
    c = sim.make(3)
    h = image.standoff(c["scan"]["S"], 0.025)
    assert np.abs(h - c["scan"]["h_true"]).max() < 0.001


def test_optical_thickness_from_wall_echoes():
    c = sim.make(2)
    h = image.standoff(c["scan"]["S"], 0.025)
    true = sum(d * np.sqrt(e) for d, e in c["wall"]["layers"])
    assert abs(image.optical_thickness(c["scan"]["S"], h) - true) < 0.002


def test_autofocus_finds_tile_velocity_and_depth():
    c = sim.make(2)                                     # tile on cement board
    img = image.form(c["scan"], "this")
    n_true = sum(d * np.sqrt(e) for d, e in c["wall"]["layers"]) / c["wall"]["thickness"]
    assert c["wall"]["name"] == "tile_on_cement_board" and abs(img["n"] - n_true) < 0.15
    sc = detect.score(c["wall"], detect.peaks(img), c["cfg"]["length"])
    assert all(t["found"] and abs(t["dz"]) < 0.006 for t in sc["targets"])


def test_drill_check_uses_depth():
    found = [dict(x=0.5, z=0.040, amp=1.0)]
    assert detect.drill_check(found, 0.51, 0.030)[0] == "CLEAR"            # 30 mm hole stops before a pipe at 40 mm
    assert detect.drill_check(found, 0.51, 0.045)[0] == "NOT_CLEAR"
    assert detect.drill_check(found, 0.56, 0.060)[0] == "CLEAR"            # 60 mm to the side


def test_free_space_puts_targets_too_deep_behind_tile():
    c = sim.make(2)
    sc = detect.score(c["wall"], detect.peaks(image.form(c["scan"], "nominal_eps")), c["cfg"]["length"])
    assert np.median([t["dz"] for t in sc["targets"] if t["found"]]) > 0.008
