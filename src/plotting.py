import numpy as np
import matplotlib.pyplot as plt

def plot_slant_tec_vs_elev(elev_range, data_2024):
    plt.figure(figsize=(10, 6))
    plt.plot(elev_range, data_2024["delays_nom"], label='Correct trajectory', color='black', linewidth=2)
    plt.plot(elev_range, data_2024["delays_err1"], '--', label='2.5 km error', color='blue')
    plt.plot(elev_range, data_2024["delays_err2"], ':', label='10 km error', color='red')
    plt.xlabel('Elevation Angle (°)', fontsize=14)
    plt.ylabel('Ionospheric Delay (m)', fontsize=14)
    plt.title('Slant Delay vs Elevation — Max Solar (2024)', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_error_abs(elev_range, data_2024):
    plt.figure(figsize=(10, 6))
    plt.plot(elev_range, np.zeros(len(elev_range)), label='Reference (no error)', color='black', linewidth=2)
    plt.plot(elev_range, data_2024['diff_abs_1'], '--', label='2.5 km error', color='red')
    plt.plot(elev_range, data_2024['diff_abs_2'], ':', label='10 km error', color='blue')
    plt.xlabel('Elevation Angle (°)', fontsize=14)
    plt.ylabel('Absolute Error (m)', fontsize=14)
    plt.title('Slant Delay Absolute Error — 2024', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_compare_years(elev_range, data_2019, data_2024):
    plt.figure(figsize=(12, 7))
    plt.plot(elev_range, np.zeros(len(elev_range)), label='2019 Reference', color='black', linewidth=2)
    plt.plot(elev_range, np.array(data_2024['delays_nom']) - np.array(data_2019['delays_nom']),
             '-', label='2024 vs 2019 (Nominal)', color='green')
    plt.plot(elev_range, np.array(data_2024['delays_err1']) - np.array(data_2019['delays_nom']),
             '--', label='2024 2.5 km vs 2019', color='orange')
    plt.plot(elev_range, np.array(data_2024['delays_err2']) - np.array(data_2019['delays_nom']),
             ':', label='2024 10 km vs 2019', color='purple')
    plt.xlabel('Elevation Angle (°)', fontsize=14)
    plt.ylabel('Difference (m)', fontsize=14)
    plt.title('Slant Delay Differences — Min vs Max Solar', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
    plt.tight_layout()
    plt.show()

def plot_error_percent(elev_range, data, label):
    plt.figure(figsize=(10, 6))
    plt.plot(elev_range, np.zeros(len(elev_range)), label=f'{label} Reference', color='black', linewidth=2)
    plt.plot(elev_range, data['err_pct_1'], '--', label=f'{label} — 2.5 km Error', color='blue')
    plt.plot(elev_range, data['err_pct_2'], ':', label=f'{label} — 10 km Error', color='red')
    plt.xlabel('Elevation Angle (°)', fontsize=14)
    plt.ylabel('Error Percentage (%)', fontsize=14)
    plt.title(f'Slant Delay Error (%) — {label}', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_vertical_tec(alt19, tec19, alt24, tec24):
    plt.figure(figsize=(10, 6))
    plt.plot(alt19, tec19, 'b-', label='2019')
    plt.plot(alt24, tec24, 'r-', label='2024')
    plt.xlabel('Altitude (km)', fontsize=14)
    plt.ylabel('Vertical TEC (TECU)', fontsize=14)
    plt.title('Vertical TEC vs Altitude — Min vs Max Solar', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()
