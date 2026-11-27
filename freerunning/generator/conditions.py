"""The conditions of the Run 4 readout per channel, from the tables of geometry/data/atlas_conditions_run4/ (sources in its NOTICE.md).

Per channel and gain: counts per MeV of visible energy, pedestal, top of the ADC, pulse (peak 1 at t = 0, ns), noise, sampling phase (ns).
- LAr: counts = E_visible / fsampl / (adc2mev_hg_mc x ramp_r4_gain / ramp_hg_mc), the MeV per count of the Phase-II front end; low gain:
  the medium-gain rows of the pulse and noise-correlation tables; noise: white Gaussian through an autoregressive filter of order 4 that
  reproduces the correlation of the ring at lags 1 to 4 (Yule-Walker), with the sigma of the front end; 15-bit ADC, pedestal 8192.
- Tile: counts_hg = E_visible x sampling fraction / (MeV per count of the Run 4 simulation x 0.525), low gain 40 times fewer counts;
  pulse CR-RC^4, tau = 10 ns, optionally with a first-order AC coupling (long tail); noise: one Gaussian, the sigma that keeps the variance
  of the double Gaussian of the table; photoelectrons per MeV of visible energy for the photoelectron statistics; 12-bit ADC, pedestal 90.
"""
import csv
import math
import os

import numpy as np

import cxx

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "geometry", "data", "atlas_conditions_run4")

TS, NBC = 25.0, 3564
LAR_PEDESTAL, LAR_TOP = 8192.0, 32767
TILE_PEDESTAL, TILE_TOP = 90.0, 4095
TILE_FENICS_SCALE = (26.25 / 12.5) * (1023.0 / 4095.0)  # MeV per count of the FENICS high gain over that of the present electronics
TILE_GAIN_RATIO = 40.0
TILE_N, TILE_TAU = 4, 10.0  # CR-RC^4, tau in ns
GAINS = ("hg", "lg")


def rows(name):
    with open(os.path.join(DATA, name)) as f:
        return list(csv.DictReader(f))


def tile_pulse(t, n=TILE_N, tau=TILE_TAU):
    """CR-RC^n with its peak, 1, at t = 0 (ns)"""
    x = (np.asarray(t, float) + n * tau) / (n * tau)
    out = np.zeros_like(x)
    m = x > 0
    out[m] = x[m] ** n * np.exp(n * (1 - x[m]))
    return out


def tile_pulse_end(eps=1e-6):
    t = np.arange(0, 2000, 0.5)
    return float(t[np.argmax(tile_pulse(t) < eps)])


class TileTail:
    """the CR-RC^4 pulse after a first-order AC coupling of time constant tau (ns): h(t) - u(t), u(t) = (1/tau) int h(t') e^{-(t-t')/tau}
    dt', zero total area. u is tabulated on 0.5 ns up to T1 = 1 us after the peak; beyond, u(t) = u(T1) e^{-(t - T1)/tau} exactly, since h is
    zero there."""

    def __init__(self, tau, t1=1000.0, dt=0.5):
        self.tau, self.t1, self.dt = float(tau), t1, dt
        t0 = -TILE_N * TILE_TAU
        self.grid = np.arange(t0, t1 + dt / 2, dt)
        h = tile_pulse(self.grid)
        a = math.exp(-dt / self.tau)
        u = np.zeros_like(h)
        for i in range(1, len(h)):  # exact for a piecewise-linear h
            u[i] = u[i - 1] * a + (dt / self.tau) * 0.5 * (h[i - 1] * a + h[i])
        self.u = u
        self.u1 = float(u[-1])

    def short(self, t):
        """h - u on the tabulated range (zero outside it)"""
        t = np.asarray(t, float)
        out = tile_pulse(t) - np.interp(t, self.grid, self.u, left=0.0, right=0.0)
        out[t > self.t1] = 0.0
        return out


class Conditions:
    def __init__(self):
        self.lar_const = {(r["layer"], int(r["ring"])): r for r in rows("lar_run4_constants_by_ring.csv")}
        self.lar_pulse = {}
        for r in rows("lar_run4_pulses_by_ring.csv"):
            w = np.array(r["wave"].split(), float)
            t = float(r["t_start_ns"]) + float(r["dt_ns"]) * np.arange(len(w))
            self.lar_pulse[(r["layer"], int(r["ring"]), r["gain"])] = (t, w)
        self.lar_rho = {(r["layer"], int(r["ring"]), r["gain"]): [float(r[f"rho{k}"]) for k in (1, 2, 3, 4)]
                        for r in rows("lar_run4_noise_autocorr_by_ring.csv")}
        self.lar_phase = {(r["layer"], int(r["ring"])): float(r["phase_ns"]) for r in rows("lar_phase_by_ring.csv")}
        key = lambda r: (int(r["section"]), int(r["sampling"]), int(r["tower"]))  # noqa: E731
        self.tile_sf = {key(r): float(r["sampling_fraction"]) for r in rows("tile_run4_sampling_fraction_by_ring.csv")}
        self.tile_scale = {key(r): r for r in rows("tile_run4_scale_by_ring.csv")}
        self.tile_noise = {}
        for r in rows("tile_run4_noise_by_ring.csv"):
            s1, s2, w = float(r["sigma1_counts"]), float(r["sigma2_counts"]), float(r["w_sigma1"])
            self.tile_noise[key(r) + ({"high": "hg", "low": "lg"}[r["gain"]],)] = math.sqrt(w * s1 * s1 + (1 - w) * s2 * s2)
        self.tile_phase = {key(r): float(r["phase_ns"]) for r in rows("tile_phase_by_ring.csv")}
        self.train = np.zeros(NBC, bool)
        for r in rows("hllhc_bunch_train_ip1.csv"):
            if int(r["colliding_ip1"]):
                self.train[int(r["first_slot"]):int(r["first_slot"]) + int(r["n_slots"])] = True

    # ---------------------------------------------------------------- LAr
    def lar(self, layer, ring, gain):
        c = self.lar_const[(layer, ring)]
        g = {"hg": "high", "lg": "medium"}[gain]
        ramp = float(c["ramp_r4_hg" if gain == "hg" else "ramp_r4_lg"])
        mev_per_count = float(c["adc2mev_hg_mc"]) * ramp / float(c["ramp_hg_mc"])
        t, w = self.lar_pulse[(layer, ring, g)]
        return dict(counts_per_visible_mev=1.0 / (float(c["fsampl_mc"]) * mev_per_count), visible_to_total=1.0 / float(c["fsampl_mc"]),
                    pedestal=LAR_PEDESTAL, top=LAR_TOP, pulse_t=t, pulse_w=w,
                    sigma=float(c["noise_adc_r4_hg" if gain == "hg" else "noise_adc_r4_lg"]), rho=self.lar_rho[(layer, ring, g)],
                    phase=self.lar_phase[(layer, ring)])

    # ---------------------------------------------------------------- Tile
    def tile(self, ring, gain):
        sf = self.tile_sf[ring]
        sc = self.tile_scale[ring]
        mev_hg = float(sc["mev_per_count_legacy_hg"]) * TILE_FENICS_SCALE
        # TileDigitsMaker rounds the sampling fraction to 3 decimals; TileHitVecToCntTool rounds the photoelectrons per MeV of visible energy
        counts = round(sf * 1000) / 1000 / mev_hg / (1.0 if gain == "hg" else TILE_GAIN_RATIO)
        npe_mev = round(float(sc["photoelectrons_per_gev"]) / 1000.0 * sf * 1000) / 1000
        return dict(counts_per_visible_mev=counts, visible_to_total=sf, pedestal=TILE_PEDESTAL, top=TILE_TOP,
                    sigma=self.tile_noise[ring + (gain,)], phase=self.tile_phase[ring], npe_per_visible_mev=npe_mev)


def ar4(rho):
    """coefficients a1..a4 and innovation variance (for unit variance) of the AR(4) process with autocorrelation rho1..rho4 (Yule-Walker);
    raises if the process is not stationary"""
    r = np.r_[1.0, rho]
    R = np.array([[r[abs(i - j)] for j in range(4)] for i in range(4)])
    a = np.linalg.solve(R, r[1:])
    var = 1.0 - float(a @ r[1:])
    roots = np.roots(np.r_[1.0, -a])
    if var <= 0 or np.any(np.abs(roots) >= 1):
        raise ValueError(f"the correlation {rho} has no stationary AR(4) process")
    return a, var


def lar_noise(rng, n, sigma, rho, burn=200):
    """n samples of the correlated LAr noise: white Gaussian through the AR(4) filter, scaled to sigma"""
    a, var = ar4(rho)
    e = rng.normal(0.0, math.sqrt(var), n + burn)
    y = np.zeros(n + burn)
    a1, a2, a3, a4 = (float(x) for x in a)
    # y[i] = a1 y[i-1] + a2 y[i-2] + a3 y[i-3] + a4 y[i-4] + e[i] from i = 4, in C++ (the same operations as the loop)
    cxx.load().lzt_ar4(e, n + burn, a1, a2, a3, a4, y)
    return sigma * y[burn:]
