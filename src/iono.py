import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from .geometry import earth_radius_at_lat, radius_along_path, slant_segment

HMIN = 50e3
HMAX = 1049e3
F_L1 = 1.57542e9  # GPS L1 [Hz]

def read_iri_txt(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Load IRI text data with columns 'km' and 'Ne/cm-3'."""
    df = pd.read_csv(path, sep=r"\s+", header=0)
    if 'km' not in df.columns or 'Ne/cm-3' not in df.columns:
        raise ValueError("File must contain 'km' and 'Ne/cm-3' columns.")
    h = df['km'].astype(float).to_numpy() * 1e3
    ne = df['Ne/cm-3'].astype(float).to_numpy() * 1e6
    return h, ne

def slant_tec(h: np.ndarray, ne: np.ndarray, elevation_deg: float, lat_deg: float) -> float:
    """Integrate the slant TEC along the line of sight."""
    R0 = earth_radius_at_lat(lat_deg)
    Rh = np.array([radius_along_path(hi, elevation_deg, lat_deg) for hi in h])
    s = slant_segment(Rh, h, R0, elevation_deg)
    ds = np.diff(s)
    ne_avg = 0.5 * (ne[:-1] + ne[1:])
    return float(np.sum(ne_avg * ds))

def iono_delay_m(slant_tec_e_m2: float, freq_hz: float = F_L1) -> float:
    """Compute ionospheric delay [m] from TEC and frequency."""
    return (40.3 * slant_tec_e_m2 / (freq_hz**2)) / 0.163

def elevation_perturbed(elev_deg: float, orbital_err_m: float, lat_deg: float) -> float:
    """Estimate the new elevation after an orbital position error."""
    R0 = earth_radius_at_lat(lat_deg)
    Rh_min = radius_along_path(HMIN, elev_deg, lat_deg)
    Rh_max = radius_along_path(HMAX, elev_deg, lat_deg)

    under = slant_segment(np.array([Rh_min]), np.array([HMIN]), R0, elev_deg)[0]
    through = slant_segment(np.array([Rh_max]), np.array([HMAX]), R0, elev_deg)[0] - under
    total = max(under + through, 1e-9)

    delta_el_rad = np.arcsin(np.clip(orbital_err_m / total, -1.0, 1.0))
    return float(elev_deg - np.degrees(delta_el_rad))

def sweep_year(
    iri_path: str,
    elevation_range_deg: np.ndarray,
    lat_deg: float,
    orbital_errors_m=(2.5e3, 10e3),
) -> Dict[str, List[float]]:
    """Sweep through elevation angles for one IRI dataset."""
    h, ne = read_iri_txt(iri_path)
    delays_nom, delays_err1, delays_err2 = [], [], []
    err_pct_1, err_pct_2, diff_abs_1, diff_abs_2 = [], [], [], []

    for elev in elevation_range_deg:
        tec = slant_tec(h, ne, elev, lat_deg)
        d_nom = iono_delay_m(tec)

        e1 = elevation_perturbed(elev, orbital_errors_m[0], lat_deg)
        e2 = elevation_perturbed(elev, orbital_errors_m[1], lat_deg)

        d1 = iono_delay_m(slant_tec(h, ne, e1, lat_deg))
        d2 = iono_delay_m(slant_tec(h, ne, e2, lat_deg))

        pct1 = abs((d1 - d_nom) / d_nom) * 100 if d_nom != 0 else 0.0
        pct2 = abs((d2 - d_nom) / d_nom) * 100 if d_nom != 0 else 0.0
        diff1 = abs(d1 - d_nom)
        diff2 = abs(d2 - d_nom)

        delays_nom.append(d_nom)
        delays_err1.append(d1)
        delays_err2.append(d2)
        err_pct_1.append(pct1)
        err_pct_2.append(pct2)
        diff_abs_1.append(diff1)
        diff_abs_2.append(diff2)

    return dict(
        delays_nom=delays_nom,
        delays_err1=delays_err1,
        delays_err2=delays_err2,
        err_pct_1=err_pct_1,
        err_pct_2=err_pct_2,
        diff_abs_1=diff_abs_1,
        diff_abs_2=diff_abs_2,
        h=h,
        ne=ne,
    )

def vertical_tec_profile(iri_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the vertical cumulative TEC profile."""
    df = pd.read_csv(iri_path, sep=r"\s+", header=0)
    alt_km = df['km'].values
    ne = df['Ne/cm-3'].values * 1e6
    delta_h = np.diff(alt_km) * 1e3
    tecu_cum = np.cumsum(ne[:-1] * delta_h) / 1e16
    tecu_cum = np.insert(tecu_cum, 0, 0.0)
    return alt_km, tecu_cum
