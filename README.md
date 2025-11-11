# Ionospheric Delay and TEC Sensitivity Analysis

This repository contains a numerical model developed to quantify how satellite orbital errors affect ionospheric signal delay estimation for Low Earth Orbit (LEO) satellites.  
The simulation computes the **Slant Total Electron Content (TEC)** and corresponding **ionospheric delay** using electron density profiles from the **International Reference Ionosphere (IRI)** model for different solar activity periods.

---
## Motivation

The Sun’s 11-year magnetic cycle plays a key role in shaping the Total Electron Content (TEC) of Earth's ionosphere.  
Increased solar activity intensifies radiation and solar wind, increasing ionization in the upper atmosphere and, consequently, signal delay.

With the rapid expansion of **LEO constellations** for navigation, communications, and Earth observation, understanding how **satellite ephemeris errors** affect ionospheric delay is critical.  
Orbit propagation models like **Two-Line Elements (TLE)** can exhibit initial position errors of several kilometers, which propagate over time due to unmodeled forces — primarily atmospheric drag.

This project quantifies how those orbit errors translate into variations in the estimated ionospheric delay, providing a physically grounded model that can be extended to any frequency or orbital regime.

---

## Overview

Simulations were conducted using IRI-generated electron density profiles near **Chicago (42° N, 88° W)** at **14:00 local solar time**, for two key periods:

- **December 1, 2019:** Solar Minimum  
- **August 1, 2024:** Solar Maximum  

Vertical profiles were computed at 1 km resolution between **50 km** and **1049 km** altitude.  
For each case, the model evaluated slant paths from **0°–89° elevation** and assessed orbital altitude errors of **2.5 km** and **10 km** across a full LEO altitude range.

---

## Features

- Computes **Slant TEC** from IRI electron density data  
- Models **ionospheric delay** in a **frequency-independent** form  
- Uses **ellipsoidal Earth geometry** and iterative path radius estimation  
- Evaluates **absolute** and **percentage delay errors** for orbital offsets  
- Generates 8 analytical plots showing geometric and plasma effects  
- Compares **solar minimum vs. solar maximum** conditions  
- Integrates via the **trapezoidal rule** for accurate path discretization  

---

## Key Equations

The Slant Total Electron Content (TEC) is defined by integrating the electron density \( Nₑ(h) \) along the signal path \( s \):

\[
TEC = \int Nₑ(h) \, ds
\]

The corresponding ionospheric delay for any radio frequency \( f \) is:

\[
\Delta s = \frac{40.3 \times TEC}{f^2}
\]

Since the TEC is calculated first, the delay can be evaluated for **any frequency band** (VHF, UHF, S-band, X-band, GPS L1, etc.) simply by changing \( f \).

---

## Repository Structure

ionosphere-delay-sensitivity
│
├── main.py # Main simulation script
├── geometry.py # Earth geometry and line-of-sight calculations
├── iono.py # TEC and ionospheric delay computations
├── plotting.py # Visualization and graph generation
├── requirements.txt # Python dependencies
└── data/
├── IRI_2019_Chicago.txt
└── IRI_2024_Chicago.txt

## Plots

| # | Title | X-Axis | Y-Axis | Purpose / Insight |
|---|--------|--------|--------|------------------|
| 1 | **Slant TEC vs Elevation (Solar Max)** | Elevation (°) | TEC (TECU) | Shows how TEC increases with decreasing elevation and how orbital errors distort it. |
| 2 | **Absolute TEC Error (Solar Max)** | Elevation (°) | Absolute Error (TECU) | Quantifies deviation from true TEC caused by 2.5 km / 10 km orbital errors. |
| 3 | **TEC Error Comparison (Min vs Max Solar)** | Elevation (°) | TEC Difference (TECU) | Compares ionospheric error magnitude between 2019 and 2024. |
| 4 | **TEC Error Percentage (Solar Max)** | Elevation (°) | Error (%) | Relative percentage error due to orbital uncertainty. |
| 5 | **TEC Error Percentage (Solar Min)** | Elevation (°) | Error (%) | Same analysis as #4, for solar minimum. |
| 6 | **TEC Error Percentage – Min vs Max Solar** | Elevation (°) | Error (%) | Combined view of both years and error levels. |
| 7 | **Delay Error vs Satellite Altitude (Solar Max)** | Satellite Altitude (m) | Error (%) | Evaluates how delay uncertainty changes with orbital height. |
| 8 | **Vertical TEC vs Altitude – Min vs Max Solar** | Altitude (km) | Vertical TEC (TECU) | Compares integrated TEC profiles for 2019 and 2024. |

**Additional Diagnostic Plots:**
- Orbital error in degrees vs. elevation (2.5 km case)  
- Slant path through the ionosphere vs. elevation  

---

## Key Results

- At **500 km altitude**, the **absolute error** in Slant TEC peaks near **20° elevation**,  
  while the **percentage error** peaks near **30°**.  
- During **solar maximum (2024)**, TEC and error magnitudes increase significantly compared to 2019.  
- For a **2.5 km** orbital uncertainty:  
  - Delay errors can reach **≈ 11 ns** for **VHF (137 MHz)** satellites.  
  - At **Ku-band (10.7 GHz)**, the delay reaches **≈ 1.8 ps** — small but relevant in precision timing applications.  
- Both **absolute and percentage errors** decrease with increasing satellite altitude.  
- The results confirm that **orbital geometry** dominates ionospheric delay sensitivity at mid-elevations.

---

## Parameters and Constants

| Symbol | Description | Value / Unit |
|---------|-------------|--------------|
| \( a \) | Earth equatorial radius | 6378.136 km |
| \( e \) | Earth eccentricity | 0.08181919 |
| \( h_{min} \) | Lower ionosphere boundary | 50 km |
| \( h_{max} \) | Upper ionosphere boundary | 1049 km |
| \( f \) | Radio frequency (variable) | e.g., 150 MHz – 8.4 GHz |
| \( c \) | Speed of light | 299,792,458 m/s |
| \( \Delta h \) | Orbital magnitude of errors in the direction that magnifies the elevation angle difference | 2.5 km, 10 km |

---

## Running the Simulation

1. Clone this repository:
   ```bash
   git clone https://github.com/<username>/ionosphere-delay-sensitivity.git
   cd ionosphere-delay-sensitivity
2. Install the dependencies:
pip install -r requirements.txt

3. Place the IRI data files in the directory:
IRI_2019_Chicago.txt
IRI_2024_Chicago.txt

4. Run the simulation:
python main.py

5. View the results:
Plots are generated automatically
Each figure corresponds to one of the eight analytical comparisons.

---

## Applications

GNSS and LEO signal propagation modeling.
Error sensitivity for orbit determination.
Space weather impact assessment.
Satellite communication and navigation system design.
Validation of ionospheric correction algorithms.

---

## Future Work

Extend geometry to non-meridional azimuths (az ≠ 0° or 180°) to include full 3D propagation paths.
Incorporate time-varying IRI profiles to simulate diurnal and seasonal changes.
Compare modeled orbital uncertainties with real satellite ephemerides to validate sensitivity predictions.
Integrate multi-frequency compensation models for advanced GNSS or radar missions.

---

## Scientific Context

The ionosphere is a dispersive plasma whose refractive index depends on electron density and radio frequency. Accurate modeling of its total electron content (TEC) is critical for satellite-based systems, where propagation paths vary dynamically with orbit geometry.
This simulation framework combines:
- Orbital geometry (Earth curvature and elevation)
- Plasma physics (from IRI profiles)
- Numerical integration (slant path TEC estimation)

It provides a quantitative framework to evaluate how orbital and environmental parameters influence signal delay.

---

## Author

Gabriel Souza Guedes
Aerospace Engineering Student – Illinois Institute of Technology
Research focus: orbital dynamics, ionospheric modeling, and satellite systems

---

## LinkedIn
https://www.linkedin.com/in/gabrielsouzaguedes/

---

## License

Creative Commons Attribution–NonCommercial 4.0 International (CC BY-NC 4.0)
Free for academic and research use with attribution. Commercial use prohibited.

---

## Acknowledgments

International Reference Ionosphere (IRI 2020) model for electron density data
Illinois Institute of Technology – Space Weather Lab
Prof. Seebany Datta-Barua for mentorship in ionospheric research
Supported by the Armour R&D Research Grant, Illinois Institute of Technology

If used in academic or technical work, please cite as:

Guedes, G. S. (2025). Ionospheric Delay and TEC Sensitivity Analysis.
Illinois Institute of Technology, Space Weather Laboratory.
GitHub Repository
