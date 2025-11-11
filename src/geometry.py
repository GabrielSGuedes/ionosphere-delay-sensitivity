import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize_scalar
from skyfield.api import EarthSatellite, load, wgs84
from datetime import timedelta

#Constants
a = 6378136.3  # Radius of the Earth in meters
ex = 0.081819190842622  # Eccentricity of the Earth
b = a*np.sqrt(1-ex**2)  # Semi-minor axis of the Earth in meters
hmin = 50e3  # Minimum height of the ionosphere in meters
hmax = 1049e3  # Maximum height of the ionosphere in meters
f1 = 1.57542e9  # Frequency of L1 in Hz
fvhf: 150e6
fuhf: 400e6
fS: 2.2e9
fx: 8.4e9
latitude_receiver = 41.8781  # Latitude of the receiver in degrees
longitude_receiver = -87.6298 # Longitude of the receiver in degrees
altitude_receiver = 180  # Altitude of the receiver in meters
c = 299,792,458

satellite_altitude = 500e3
elevation_range = np.linspace(1,89,89) # Zenith angle
iono_delays = []
iono_delays_err = []
error_percent = []

# Function to calculate the radius of the Earth at a given latitude
def R_Lat(latitude_deg):
    lat_rad = np.radians(latitude_deg)
    numerator = (a**2 * np.cos(lat_rad))**2 + (b**2 * np.sin(lat_rad))**2
    denominator = (a * np.cos(lat_rad))**2 + (b * np.sin(lat_rad))**2
    r = np.sqrt(numerator / denominator)
    
    return r

# Função para encontrar raio ao longo da linha de visada
def get_radius_on_path(altitude_point, elevation_deg, latitude_receiver):
    def objective(phi_deg):
        lat_point = latitude_receiver + phi_deg
        R_point = R_Lat(lat_point)
        num = R_Lat(latitude_receiver) * np.sin(np.radians(90+elevation_deg))
        den = R_point + altitude_point
        rhs = num / den
        lhs = np.cos(np.radians(phi_deg + elevation_deg))
        return abs(lhs - rhs)
    
    result = minimize_scalar(objective, bounds=(0, 89), method='bounded')
    if result.success:
        phi_opt = result.x
        lat_target = latitude_receiver + phi_opt
        R_target = R_Lat(lat_target)
        return R_target
    else:
        raise ValueError("Não foi possível encontrar φ para a altitude fornecida.")

R_receiver = get_radius_on_path(0,0,latitude_receiver)  # Radius of the Earth at the receiver's latitude
# R_hmin, R_hmax, and R_satellite will be computed per elevation inside the loop

def normalize(arr):
    arr = np.array(arr, dtype=float)
    min_val = np.min(arr)
    max_val = np.max(arr)
    if max_val == min_val:
        return np.zeros_like(arr)  # or np.ones_like(arr)
    return (arr - min_val) / (max_val - min_val)

all_iono_delays = {'2019': [], '2024': []}
all_iono_delays_new1 = {'2019': [], '2024': []}
all_iono_delays_new2 = {'2019': [], '2024': []}
all_error_percent1 = {'2019': [], '2024': []}
all_error_percent2 = {'2019': [], '2024': []}
all_diff_iono1 = {'2019': [], '2024':[]}
all_diff_iono2 = {'2019': [], '2024':[]} 
combined_impact_by_year = {'2019': [], '2024':[]}
normalized_TEC_by_year = {'2019': [], '2024':[]}
normalized_error_by_year = {'2019': [], '2024':[]}


# TEC Calculation
IRI_2024 = r'C:\Users\gabri\Downloads\IRI_2024_Chicago.txt'
IRI_2019 = r'C:\Users\gabri\Downloads\IRI_2019_Chicago.txt'

for year, file_path in [('2019', IRI_2019), ('2024', IRI_2024)]:
    try:
        print(f"\nProcessing {year} file...")
        df = pd.read_csv(file_path, sep=r'\s+', header=0)
        print(f"File read successfully. First few rows:\n{df.head()}")
        
        # Check critical columns exist
        if 'km' not in df.columns or 'Ne/cm-3' not in df.columns:
            print(f"Error: Required columns missing in {year} data")
            continue
            
        df['km'] = df['km'].astype(float)
        df['Ne/cm-3'] = df['Ne/cm-3'].astype(float)
        print(f"Data types converted. km range: {df['km'].min()} to {df['km'].max()} km")
        
        height = df['km'].values * 1e3
        ne = df['Ne/cm-3'].values * 1e6
        print(f"Height array shape: {height.shape}, Ne array shape: {ne.shape}")
        
        iono_delays = []
        iono_delays_new1 = []
        iono_delays_new2 = []
        error_percent1 = []
        error_percent2 = []
        slant_TEC_values = []
        difference_iono_delay1 = []
        difference_iono_delay2 = []
        orbital_errors1_deg = []
        orbital_errors2_deg = []
        new_elevation1 = []
        new_elevation2 = []
        slant_ionosphere_g = []



        for elevation in elevation_range:
            # Compute R_hmin, R_hmax, and R_satellite for the current elevation
            R_hmin = get_radius_on_path(hmin, elevation, latitude_receiver)
            R_hmax = get_radius_on_path(hmax, elevation, latitude_receiver)
            R_satellite = get_radius_on_path(satellite_altitude, elevation, latitude_receiver)

            # Slant Path beneath the Ionosphere
            slant_under_ionosphere = np.sqrt((R_hmin+hmin)**2 +(R_receiver* np.sin(np.radians(elevation)))**2-R_hmin**2) - R_receiver*np.sin(np.radians(elevation))

            #Slant Ionosphere Path Calculation
            if satellite_altitude< hmax:
                slant_ionosphere =np.sqrt((R_satellite+satellite_altitude)**2 +(R_receiver* np.sin(np.radians(elevation)))**2-R_satellite**2) - R_receiver*np.sin(np.radians(elevation)) - slant_under_ionosphere
            else:
                slant_ionosphere = np.sqrt((R_hmax+hmax)**2 +(R_receiver* np.sin(np.radians(elevation)))**2-R_hmax**2) - R_receiver*np.sin(np.radians(elevation)) - slant_under_ionosphere
            total_slant_path = slant_ionosphere + slant_under_ionosphere
            # Angle Errors from the TLE
            orbital_error_1_m = 2.5e3  # Random error in meters
            orbital_error_2_m = 10e3
            delta_elevation_1 = orbital_error_1_m / (total_slant_path)
            delta_elevation_2 = orbital_error_2_m / (total_slant_path)
            delta_elevation_1_rad = np.arcsin(orbital_error_1_m / (total_slant_path)) # Convert radians to degrees
            delta_elevation_2_rad = np.arcsin(orbital_error_2_m/total_slant_path)
            new_elevation_1 = elevation - np.degrees(delta_elevation_1_rad)  # New elevation angle after error
            new_elevation_2 = elevation - np.degrees(delta_elevation_2_rad)
            # Calculate the slant path lengths
            R_height = np.array([get_radius_on_path(h, elevation, latitude_receiver) for h in height])
            R_new_height_1 = np.array([get_radius_on_path(h, new_elevation_1, latitude_receiver) for h in height])
            R_new_height_2 = np.array([get_radius_on_path(h, new_elevation_2, latitude_receiver) for h in height])

            s = np.sqrt((R_height+height)**2 +(R_receiver* np.sin(np.radians(elevation)))**2-R_height**2) - R_receiver*np.sin(np.radians(elevation))
            s_new_1 = np.sqrt((R_new_height_1+height)**2 +(R_receiver* np.sin(np.radians(new_elevation_1)))**2-R_new_height_1**2) - R_receiver*np.sin(np.radians(new_elevation_1))
            s_new_2 = np.sqrt((R_new_height_2+height)**2 +(R_receiver* np.sin(np.radians(new_elevation_2)))**2-R_new_height_2**2) - R_receiver*np.sin(np.radians(new_elevation_2))     
            delta_s = np.diff(s)  # Calculate the differential path length
            delta_s_new_1 = np.diff(s_new_1)  # Calculate the differential path length for new elevation
            delta_s_new_2 = np.diff(s_new_2)

            ne_avg = 0.5*(ne[:-1] + ne[1:])  # Average Ne for trapezoidal integration (length matches delta_s)

            slant_TEC = np.sum(ne_avg * delta_s)  # Slant TEC in /m^2
            slant_TEC_new_1 = np.sum(ne_avg * delta_s_new_1)
            slant_TEC_new_2 = np.sum(ne_avg*delta_s_new_2)
            slant_TEC_TECU = slant_TEC / 1e16  # Convert to TECU
            
            slant_TEC_new_TECU = slant_TEC_new_1 / 1e16
            slant_TECU_new_1 = slant_TEC_new_1/1e16
            slant_TECU_new_2 = slant_TEC_new_2/1e16

            # After calculating slant_TEC
            if np.isnan(slant_TEC) or np.isinf(slant_TEC):
                print(f"Warning: Invalid TEC at elevation {elevation}°")
                slant_TEC = 0  # Fallback value

            ionosphere_delay = (40.3 * slant_TEC / f1**2)/0.163  # Ionosphere delay in meters
            ionosphere_delay_new_1 = (40.3 * slant_TEC_new_1 / f1**2)/0.163  # Error Ionosphere delay in meters
            ionosphere_delay_new_2 = (40.3*slant_TEC_new_2 /f1**2)/0.163
            error_perc_iono_delay_1 = abs((ionosphere_delay_new_1 - ionosphere_delay) / ionosphere_delay) * 100
            error_iono_delay_1 = abs(ionosphere_delay_new_1 - ionosphere_delay)
            error_perc_iono_delay_2 = abs((ionosphere_delay_new_2 - ionosphere_delay) / ionosphere_delay) * 100
            error_iono_delay_2 = abs(ionosphere_delay_new_2 - ionosphere_delay)
            
            iono_delays.append(ionosphere_delay)
            iono_delays_new1.append(ionosphere_delay_new_1)
            iono_delays_new2.append(ionosphere_delay_new_2)
            error_percent1.append(error_perc_iono_delay_1)
            error_percent2.append(error_perc_iono_delay_2)  
            difference_iono_delay1.append(error_iono_delay_1)
            difference_iono_delay2.append(error_iono_delay_2)
            orbital_errors1_deg.append(delta_elevation_1)
            orbital_errors2_deg.append(delta_elevation_2)
            slant_ionosphere_g.append(slant_ionosphere)


        # Armazenar resultados para este arquivo
        all_iono_delays[year] = iono_delays
        all_iono_delays_new1[year] = iono_delays_new1
        all_iono_delays_new2[year] = iono_delays_new2
        all_error_percent1[year] = error_percent1
        all_error_percent2[year] = error_percent2
        all_diff_iono1[year] = difference_iono_delay1
        all_diff_iono2[year] = difference_iono_delay2

    except Exception as e:
        print(f"Error reading or processing the file: {e}")
        TEC = 0  # Default value if file reading fails

print(f"First 5 elevation values: {elevation_range[:5]}")
print(f"First 5 delay values (2024): {all_iono_delays['2024'][:5]}")
print(f"First 5 delay values (2024): {all_error_percent1['2024'][:5]}")
print(f"First 5 delay values (2024): {all_diff_iono1['2024'][:5]}")
print("\nPrimeiros 35 erros percentuais (2024 - 2.5 km):")
print(all_error_percent1['2024'][:35])

print("\nPrimeiros 35 erros percentuais (2024 - 10 km):")
print(all_error_percent2['2024'][:35])



# Dicionários para armazenar os dados
altitudes = {'2019': [], '2024': []}
tec_values = {'2019': [], '2024': []}

for year, file_path in [('2019', IRI_2019), ('2024', IRI_2024)]:
    try:
        df = pd.read_csv(file_path, sep=r'\s+', header=0)
        
        # Extrair altitudes (convertendo para km)
        altitudes[year] = df['km'].values
        
        # Calcular TEC acumulado até cada altitude (integral da densidade eletrônica)
        ne = df['Ne/cm-3'].values * 1e6  # Converter de cm⁻³ para m⁻³
        delta_h = np.diff(df['km'].values) * 1e3  # Diferença de altura em metros
        tec_cumulative = np.cumsum(ne[:-1] * delta_h) / 1e16  # TECU (1 TECU = 1e16 elétrons/m²)
        
        # Adicionar um zero no início para que o TEC comece em 0
        tec_values[year] = np.insert(tec_cumulative, 0, 0)
        
    except Exception as e:
        print(f"Erro ao processar {year}: {e}")

satellite_altitude_range = np.linspace(160e3, 2000e3, 46)
elevations = [15,30,45,60,75]
error_by_elevation = {elev: [] for elev in elevations}  # Initialize the dictionary to store errors by elevation
# Read and parse the file once before the loop
try:
    df = pd.read_csv(IRI_2024, sep=r'\s+', header=0)
    df['km'] = df['km'].astype(float)  # Ensure km is float
    df['Ne/cm-3'] = df['Ne/cm-3'].astype(float)  # Ensure Ne/cm-3 is float
    height = df['km'].values * 1e3  # Convert km to meters
    ne = df['Ne/cm-3'].values * 1e6  # Convert Ne from cm^-3 to m^-3
except Exception as e:
    print(f"Error reading or processing the file: {e}")
    TEC = 0  # Default value if file reading fails
    height = np.array([])
    ne = np.array([])

for elevation in elevations:
    for satellite_altitude in satellite_altitude_range:
        # Compute R_hmin, R_hmax, and R_satellite for the current elevation
        R_hmin = get_radius_on_path(hmin, elevation, latitude_receiver)
        R_hmax = get_radius_on_path(hmax, elevation, latitude_receiver)
        R_satellite = get_radius_on_path(satellite_altitude, elevation, latitude_receiver)

        slant_under_ionosphere = np.sqrt((R_hmin+hmin)**2 +(R_receiver* np.sin(np.radians(elevation)))**2-R_hmin**2) - R_receiver*np.sin(np.radians(elevation))

        if satellite_altitude< hmax:
            slant_ionosphere =np.sqrt((R_satellite+satellite_altitude)**2 +(R_receiver* np.sin(np.radians(elevation)))**2-R_satellite**2) - R_receiver*np.sin(np.radians(elevation)) - slant_under_ionosphere
        else:
            slant_ionosphere = np.sqrt((R_hmax+hmax)**2 +(R_receiver* np.sin(np.radians(elevation)))**2-R_hmax**2) - R_receiver*np.sin(np.radians(elevation)) - slant_under_ionosphere

        total_slant_path = slant_ionosphere + slant_under_ionosphere
        orbital_error_m = 2.5e3  # Random error in meters
        delta_elevation_rad = np.arcsin(orbital_error_m / (total_slant_path)) # Convert radians to degrees
        new_elevation = elevation - np.degrees(delta_elevation_rad)  # New elevation angle after error

        R_height = np.array([get_radius_on_path(h, elevation, latitude_receiver) for h in height])
        R_new_height = np.array([get_radius_on_path(h, new_elevation, latitude_receiver) for h in height])
        s = np.sqrt((R_height+height)**2 +(R_receiver* np.sin(np.radians(elevation)))**2-R_height**2) - R_receiver*np.sin(np.radians(elevation))
        s_new = np.sqrt((R_new_height+height)**2 +(R_receiver* np.sin(np.radians(new_elevation)))**2-R_new_height**2) - R_receiver*np.sin(np.radians(new_elevation))
        delta_s = np.diff(s)  # Calculate the differential path length
        delta_s_new = np.diff(s_new)  # Calculate the differential path length for new elevation
        ne_avg = 0.5*(ne[:-1] + ne[1:])  # Average Ne for trapezoidal integration (length matches delta_s)
        slant_TEC = np.sum(ne_avg * delta_s)  # Slant TEC in /m^2
        slant_TEC_new = np.sum(ne_avg * delta_s_new)
        slant_TEC_TECU = slant_TEC / 1e16  # Convert to TECU
        slant_TEC_new_TECU = slant_TEC_new / 1e16
        ionosphere_delay = (40.3 * slant_TEC / f1**2)/0.163  # Ionosphere delay in mevters
        ionosphere_delay_new = (40.3 * slant_TEC_new / f1**2)/0.163  # Error Ionosphere delay in meters
        error_iono_delay = abs((ionosphere_delay_new - ionosphere_delay) / ionosphere_delay) * 100
        error_by_elevation[elevation].append(error_iono_delay)

# ------------------------------------------------------------
# Dicionários para armazenar, para cada ano, a elevação de máximo erro percentual
# e a elevação de máximo erro absoluto
elev_max_percent = {'2019': [], '2024': []}
elev_max_abs     = {'2019': [], '2024': []}

# ------------------------------------------------------------
# 2) Para cada ano, recalcule o erro para cada altitude e encontre a elevação que maximiza
for year, file_path in [('2019', IRI_2019), ('2024', IRI_2024)]:
    # --- leia o IRI e monte os arrays de height, ne como já faz ---
    df = pd.read_csv(file_path, sep=r'\s+', header=0)
    height = df['km'].values * 1e3
    ne     = df['Ne/cm-3'].values * 1e6

    for sat_alt in satellite_altitude_range:
        pct_err_list = []
        abs_err_list = []

        # para cada elevação, reproduza o cálculo de delay e erro (usando orbital_error_1 de 2.5 km)
        for elev in elevation_range:
            # --- calcule R_satellite com a função que você já tem ---
            R_hmin      = get_radius_on_path(hmin, elev, latitude_receiver)
            R_hmax      = get_radius_on_path(hmax, elev, latitude_receiver)
            R_satellite = get_radius_on_path(sat_alt, elev, latitude_receiver)

            # --- calcule o slant path igual ao seu código acima ---
            slant_under = np.sqrt((R_hmin+hmin)**2 + (R_receiver*np.sin(np.radians(elev)))**2 - R_hmin**2) \
                          - R_receiver*np.sin(np.radians(elev))
            if sat_alt < hmax:
                slant_iono = np.sqrt((R_satellite+sat_alt)**2 + (R_receiver*np.sin(np.radians(elev)))**2 - R_satellite**2) \
                             - R_receiver*np.sin(np.radians(elev)) - slant_under
            else:
                slant_iono = np.sqrt((R_hmax+hmax)**2 + (R_receiver*np.sin(np.radians(elev)))**2 - R_hmax**2) \
                             - R_receiver*np.sin(np.radians(elev)) - slant_under
            total_path = slant_under + slant_iono

            # --- erro orbital de 2.5 km ---
            orbital_error_m = 2.5e3
            delta_el_rad   = np.arcsin(orbital_error_m / total_path)
            new_elev       = elev - np.degrees(delta_el_rad)

            # --- calcule o TEC e o delay original e com erro ---
            R_h_arr       = np.array([get_radius_on_path(h, elev, latitude_receiver) for h in height])
            R_h_new_arr   = np.array([get_radius_on_path(h, new_elev, latitude_receiver) for h in height])
            s             = np.sqrt((R_h_arr+height)**2 + (R_receiver*np.sin(np.radians(elev)))**2 - R_h_arr**2) \
                            - R_receiver*np.sin(np.radians(elev))
            s_new         = np.sqrt((R_h_new_arr+height)**2 + (R_receiver*np.sin(np.radians(new_elev)))**2 - R_h_new_arr**2) \
                            - R_receiver*np.sin(np.radians(new_elev))
            delta_s       = np.diff(s)
            delta_s_new   = np.diff(s_new)
            ne_avg        = 0.5*(ne[:-1] + ne[1:])
            tec           = np.sum(ne_avg * delta_s)
            tec_new       = np.sum(ne_avg * delta_s_new)
            delay         = (40.3 * tec  / f1**2) / 0.163
            delay_new     = (40.3 * tec_new / f1**2) / 0.163

            # --- erros percentuais e absolutos ---
            pct_err_list.append(abs((delay_new - delay) / delay) * 100)
            abs_err_list.append(abs(delay_new - delay))

        # --- escolha a elevação de máximo erro ---
        idx_pct = int(np.nanargmax(pct_err_list))
        idx_abs = int(np.nanargmax(abs_err_list))
        elev_max_percent[year].append(elevation_range[idx_pct])
        elev_max_abs[year].append(elevation_range[idx_abs])

# Plotting the results

#1st Graph
plt.figure(figsize=(10, 6))

# Correct trajectory
plt.plot(
    elevation_range,
    all_iono_delays['2024'],
    label='Correct trajectory',
    color='black',
    linewidth=2
)

# Trajectories with orbital errors
plt.plot(
    elevation_range,
    all_iono_delays_new1['2024'],
    '--',
    label='2.5 km error',
    color='blue'
)

plt.plot(
    elevation_range,
    all_iono_delays_new2['2024'],
    ':',
    label='10 km error',
    color='red'
)

plt.xlabel('Elevation Angle (°)', fontsize=14)
plt.ylabel('Total Electron Content (TECU)', fontsize=14)
plt.title('Slant TEC vs Elevation Angle - Maximum Solar', fontsize=18)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()


#2nd Graph
if len(all_iono_delays['2024']) == len(elevation_range):
    plt.figure(figsize=(10, 6))
    # --- Referência (0 m) ---
    plt.plot(
        elevation_range,
        np.zeros(len(elevation_range)),  # Linha de referência em 0
        label='2024 Correct (Reference)',
        color='black',
        linewidth=2
    )

    # --- Erros em 2024 ---
    # 1. Erro de 2.5 km (dados de all_diff_iono1['2024'])
    plt.plot(
        elevation_range,
        all_diff_iono1['2024'],
        '--',
        label='2.5 km error',
        color='red'
    )

    # 2. Erro de 10 km (dados de all_diff_iono2['2024'])
    plt.plot(
        elevation_range,
        all_diff_iono2['2024'],
        ':',
        label='10 km error',
        color='blue'
    )

    # Formatação
    plt.xlabel('Elevation Angle (°)', fontsize=14)
    plt.ylabel('Absolute Error (TECU)', fontsize=14)
    plt.title('Slant TEC Errors - Maximum Solar', fontsize=18)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.show()
else:
    print(f"Cannot plot: Data length mismatch. Elevation: {len(elevation_range)}, Delays: {len(all_iono_delays['2024'])}")

#3rd Graph

plt.figure(figsize=(12, 7))

# --- Referência (0 m) ---
plt.plot(
    elevation_range,
    np.zeros(len(elevation_range)),  # Linha de referência em 0
    label='2019 Correct Path (Reference)',
    color='black',
    linewidth=2
)

# --- Erros em 2019 ---
# 1. 2019 - Erro de 2.5 km (já calculado em all_diff_iono1['2019'])
plt.plot(
    elevation_range,
    all_diff_iono1['2019'],
    '--',
    label='2019 - 2.5 km error',
    color='blue'
)

# 2. 2019 - Erro de 10 km (já calculado em all_diff_iono2['2019'])
plt.plot(
    elevation_range,
    all_diff_iono2['2019'],
    ':',
    label='2019 - 10 km error',
    color='red'
)

# --- Erros em 2024 ---
# 3. 2024 - Trajetória correta vs 2019 correta
plt.plot(
    elevation_range,
    np.array(all_iono_delays['2024']) - np.array(all_iono_delays['2019']),
    '-',
    label='2024 Correct Path vs 2019',
    color='green'
)

# 4. 2024 - Erro de 2.5 km vs 2019 correta
plt.plot(
    elevation_range,
    np.array(all_iono_delays_new1['2024']) - np.array(all_iono_delays['2019']),
    '--',
    label='2024 - 2.5 km error vs 2019',
    color='orange'
)

# 5. 2024 - Erro de 10 km vs 2019 correta
plt.plot(
    elevation_range,
    np.array(all_iono_delays_new2['2024']) - np.array(all_iono_delays['2019']),
    ':',
    label='2024 - 10 km error vs 2019',
    color='purple'
)

# Formatação
plt.xlabel('Elevation Angle (°)', fontsize=14)
plt.ylabel('Absolute Error (TECU)', fontsize=14)
plt.title('Slant TEC Errors - Minimum vs Maximum Solar', fontsize=18)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
plt.tight_layout()
plt.show()

#4th Graph

plt.figure(figsize=(10, 6))

# Trajetória correta (0% de erro)
plt.plot(
    elevation_range,
    np.zeros(len(elevation_range)),
    label='2024 Correct Path',
    color='black',
    linewidth=2
)

# Erros percentuais em 2024 (2.5 km e 10 km)
plt.plot(
    elevation_range,
    all_error_percent1['2024'],  # % erro para 2.5 km
    '--',
    label='2024 - 2.5 km Error',
    color='blue'
)

plt.plot(
    elevation_range,
    all_error_percent2['2024'],  # % erro para 10 km
    ':',
    label='2024 - 10 km Error',
    color='red'
)

plt.xlabel('Elevation Angle (°)', fontsize=14)
plt.ylabel('Error Percentage (%)', fontsize=14)
plt.title('Slant TEC Error (%) - Maximum Solar', fontsize=18)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()

#5th Graph

plt.figure(figsize=(10, 6))

# Trajetória correta (0% de erro)
plt.plot(
    elevation_range,
    np.zeros(len(elevation_range)),
    label='2019 Correct Path',
    color='black',
    linewidth=2
)

# Erros percentuais em 2024 (2.5 km e 10 km)
plt.plot(
    elevation_range,
    all_error_percent1['2019'],  # % erro para 2.5 km
    '--',
    label='2019 - 2.5 km Error',
    color='blue'
)

plt.plot(
    elevation_range,
    all_error_percent2['2019'],  # % erro para 10 km
    ':',
    label='2019 - 10 km Error',
    color='red'
)

plt.xlabel('Elevation Angle (°)', fontsize=14)
plt.ylabel('Error Percentage (%)', fontsize=14)
plt.title('Slant TEC Error (%) - Minimum Solar', fontsize=18)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()


#6th Graph

plt.figure(figsize=(12, 7))

# --- 2024 ---
plt.plot(
    elevation_range,
    all_error_percent1['2024'],
    '--',
    label='2024 - 2.5 km Error',
    color='red'
)

plt.plot(
    elevation_range,
    all_error_percent2['2024'],
    ':',
    label='2024 - 10 km Error',
    color='blue'
)

# --- 2019 ---
plt.plot(
    elevation_range,
    all_error_percent1['2019'],
    '--',
    label='2019 - 2.5 km Error',
    color='green'
)

plt.plot(
    elevation_range,
    all_error_percent2['2019'],
    ':',
    label='2019 - 10 km Error',
    color='purple'
)

plt.xlabel('Elevation Angle (°)', fontsize=14)
plt.ylabel('Error Percentage (%)', fontsize=14)
plt.title('Slant TEC Error (%) - Minimum vs Maximum Solar', fontsize=18)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
plt.tight_layout()
plt.show()

#7th Graph

plt.figure(figsize=(10, 6))
for elev in elevations:
    plt.plot(satellite_altitude_range, error_by_elevation[elev], label=f'Elevation {elev}°')
plt.ylabel('Percentage Error in Ionosphere Delay (%)', fontsize=14)
plt.xlabel('Satellite Altitude (m)', fontsize=14)
plt.title('2.5 km Error (%) in Slant TEC vs Altitude - Maximum Solar Period', fontsize=18)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#8th Graph

plt.figure(figsize=(10, 6))
plt.plot(altitudes['2019'], tec_values['2019'], 'b-', label='2019')
plt.plot(altitudes['2024'], tec_values['2024'], 'r-', label='2024')
plt.xlabel('Altitude (km)', fontsize=14)
plt.ylabel('Total Electron Content (TECU)', fontsize=14)
plt.title('Vertical TEC vs Altitude - Minimum vs Maximum Solar', fontsize=18)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()

plt.plot(elevation_range, orbital_errors1_deg)
plt.xlabel("Elevation angle (°)", fontsize=14)
plt.ylabel("Elevation angle error due to orbital uncertainty (°)", fontsize=14)
plt.title("Orbital Error (2.5 km) vs. Elevation", fontsize=18)
plt.grid(True)
plt.show()

plt.plot (elevation_range, slant_ionosphere_g)
plt.xlabel("Elevation angle (°)", fontsize=14)
plt.ylabel("Ionosphere path", fontsize=14)
plt.title('Path through Ionosphere vs Elevation', fontsize=18)
plt.grid(True)
plt.show()

# ------------------------------------------------------------
# 3) Plote: Elevação de máximo erro percentual vs altitude
plt.figure(figsize=(8,5))
plt.plot(satellite_altitude_range/1e3, elev_max_percent['2019'], '-o', label='2019')
plt.plot(satellite_altitude_range/1e3, elev_max_percent['2024'], '-o', label='2024')
plt.xlabel('Altitude do Satélite (km)', fontsize=12)
plt.ylabel('Elevação de Máximo Erro % (°)', fontsize=12)
plt.title('Elevação de Máximo Erro Percentual vs Altitude', fontsize=14)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

# 4) Plote: Elevação de máximo erro absoluto vs altitude
plt.figure(figsize=(8,5))
plt.plot(satellite_altitude_range/1e3, elev_max_abs['2019'], '-s', label='2019')
plt.plot(satellite_altitude_range/1e3, elev_max_abs['2024'], '-s', label='2024')
plt.xlabel('Altitude do Satélite (km)', fontsize=12)
plt.ylabel('Elevação de Máximo Erro Absoluto (°)', fontsize=12)
plt.title('Elevação de Máximo Erro Absoluto vs Altitude', fontsize=14)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
# ------------------------------------------------------------