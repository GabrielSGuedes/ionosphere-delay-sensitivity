import numpy as np
from scipy.optimize import minimize_scalar


A = 6378136.3  # semi-major axis [m]
E = 0.081819190842622  # eccentricity
B = A * np.sqrt(1 - E**2)  # semi-minor axis [m]

def earth_radius_at_lat(latitude_deg: float) -> float:
    """Compute the local Earth radius at a given geodetic latitude."""
    lat = np.radians(latitude_deg)
    num = (A**2 * np.cos(lat))**2 + (B**2 * np.sin(lat))**2
    den = (A * np.cos(lat))**2 + (B * np.sin(lat))**2
    return float(np.sqrt(num / den))

def radius_along_path(altitude_m: float, elevation_deg: float, latitude_receiver: float) -> float:
    """Find the shell radius along the line of sight at a given altitude."""
    def objective(phi_deg: float) -> float:
        lat_point = latitude_receiver + phi_deg
        R_point = earth_radius_at_lat(lat_point)
        num = earth_radius_at_lat(latitude_receiver) * np.sin(np.radians(90 + elevation_deg))
        den = R_point + altitude_m
        rhs = num / max(den, 1e-9)
        lhs = np.cos(np.radians(phi_deg + elevation_deg))
        return abs(lhs - rhs)

    result = minimize_scalar(objective, bounds=(0.0, 89.0), method='bounded')
    if not result.success:
        raise ValueError("Could not find φ for this altitude.")
    phi_opt = result.x
    return earth_radius_at_lat(latitude_receiver + phi_opt)

def slant_segment(Rh: np.ndarray, h: np.ndarray, R0: float, elev_deg: float) -> np.ndarray:
    """Compute the slant segment length for each height sample."""
    elev = np.radians(elev_deg)
    return np.sqrt((Rh + h)**2 + (R0 * np.sin(elev))**2 - Rh**2) - R0 * np.sin(elev)
