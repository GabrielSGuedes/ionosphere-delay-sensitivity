import numpy as np
from .iono import sweep_year, vertical_tec_profile
from .plotting import (
    plot_slant_tec_vs_elev,
    plot_error_abs,
    plot_compare_years,
    plot_error_percent,
    plot_vertical_tec,
)
# Fixed receiver location: Chicago, IL
LATITUDE_RECEIVER = 41.8781
LONGITUDE_RECEIVER = -87.6298
ALTITUDE_RECEIVER_M = 180

IRI_2019 = r"data/IRI_2019_Chicago.txt"
IRI_2024 = r"data/IRI_2024_Chicago.txt"

def main():
    elevation_range = np.linspace(1, 89, 89)

    print("Processing 2019...")
    data_2019 = sweep_year(IRI_2019, elevation_range, LATITUDE_RECEIVER)
    print("Processing 2024...")
    data_2024 = sweep_year(IRI_2024, elevation_range, LATITUDE_RECEIVER)

    plot_slant_tec_vs_elev(elevation_range, data_2024)
    plot_error_abs(elevation_range, data_2024)
    plot_compare_years(elevation_range, data_2019, data_2024)
    plot_error_percent(elevation_range, data_2024, "2024")
    plot_error_percent(elevation_range, data_2019, "2019")

    alt19, tec19 = vertical_tec_profile(IRI_2019)
    alt24, tec24 = vertical_tec_profile(IRI_2024)
    plot_vertical_tec(alt19, tec19, alt24, tec24)

if __name__ == "__main__":
    main()
