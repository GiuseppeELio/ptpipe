"""Scale bridging: from optical quantities to heat sources (SI units).

Nano -> meso:  sigma_abs (Mie/FEM) + areal or volume density -> effective absorbance /
               attenuation, or Maxwell-Garnett layer for the TMM.
Meso -> macro: absorbance profile a(z) (TMM / FEM) x beam profile I(r, t) -> q(r, z, t).
"""
import numpy as np


def gauss_I(r, P, w):
    """Intensity (W/m2) of a Gaussian beam of power P and 1/e^2 radius w."""
    return 2 * P / (np.pi * w ** 2) * np.exp(-2 * r ** 2 / w ** 2)


def smooth_rect(t, t_on, t_off, tau=0.05):
    """Smoothed on/off window (logistic edges of width tau, seconds)."""
    a = np.clip((t - t_on) / tau, -60, 60); b = np.clip((t - t_off) / tau, -60, 60)
    return 1 / (1 + np.exp(-a)) - 1 / (1 + np.exp(-b))


def monolayer_absorbance(n_s, sigma_abs):
    """Independent-particle (dilute) absorbance of a monolayer: n_s [1/m2], sigma [m2].
    Valid while n_s*sigma_ext << 1 and particles are not near-field coupled;
    beyond that use FEM/RCWA on the periodic cell (Floquet) instead."""
    return n_s * sigma_abs


def volume_attenuation(N_v, sigma_ext):
    """Beer-Lambert attenuation coefficient (1/m) of a dilute dispersion."""
    return N_v * sigma_ext


def source_from_profile(z_edges_layer, a_per_m, z0, P, w, window=None):
    """Return a callable q(r, z, t) [W/m3] = I(r) * a(z - z0) * window(t).
    z_edges_layer: positions (m, relative to the layer start) where a(z) [1/m] is given."""
    def q(r, z, t):
        zz = z - z0
        inside = (zz >= z_edges_layer[0]) & (zz <= z_edges_layer[-1])
        a = np.where(inside, np.interp(zz, z_edges_layer, a_per_m), 0.0)
        f = 1.0 if window is None else window(t)
        return gauss_I(r, P, w) * a * f
    return q
