import os
import sys
PACKAGE_TOP_LEVEL = sys.path[0]
assert PACKAGE_TOP_LEVEL.endswith("Julia_ELFIN_Tools"), f"Script directory {PACKAGE_TOP_LEVEL} is not in /Julia_ELFIN_Tools"

# Library includes
import requests # for web downloading
import spacepy # for L-shell and B-field calculations
from spacepy import pycdf, coordinates, irbempy
import numpy as np
import re # regex
from multiprocessing import Pool # allows multithreading
from wakepy import keep # provides functionality to prevent computer sleeping

N_THREADS = 8 # Number of threds to run L-shell calculations on

"""
elfin_data_download.py
This script automatically downloads all available level 2 ELFIN data for both satellites. Call
the script with command line argument `-delete` to delete all existing data and start fresh.

Written by Julia Claxton (julia.claxton@colorado.edu)
Released under MIT License (see LICENSE.txt for full license)
"""

def download_all_elfin_data(delete_prior = False):
    print("")
    # Create data directories if they don't already exist
    if os.path.isdir(f"{PACKAGE_TOP_LEVEL}/data/raw") == False:
        os.mkdir(f"{PACKAGE_TOP_LEVEL}/data/raw")

    if os.path.isdir(f"{PACKAGE_TOP_LEVEL}/data/processed_scientific_data") == False:
        os.mkdir(f"{PACKAGE_TOP_LEVEL}/data/processed_scientific_data")

    if os.path.isdir(f"{PACKAGE_TOP_LEVEL}/data/processed_position_data") == False:
        os.mkdir(f"{PACKAGE_TOP_LEVEL}/data/processed_position_data")

    if delete_prior == True:
        answer = input("\033[31;1;4mDelete all existing data and download all data from ELFIN's entire lifetime? (Yy/Nn)\033[0m ")
        if answer.lower() != "y":
            print("Cancelled, exiting...")
            return
        print("")
        answer = input("\033[31;1;4mARE YOU SURE? Type \"I WANT TO DELETE 9 HOURS OF PROCESSOR TIME\" to proceed, anything else to cancel.\033[0m\n")
        if answer != "I WANT TO DELETE 9 HOURS OF PROCESSOR TIME":
            print("Cancelled, exiting...")
            return
        print("")

        os.system(f"rm -f {PACKAGE_TOP_LEVEL}/data/raw/*");                       print(f"rm -f {PACKAGE_TOP_LEVEL}/data/raw/*")
        os.system(f"rm -f {PACKAGE_TOP_LEVEL}/data/processed_scientific_data/*"); print(f"rm -f {PACKAGE_TOP_LEVEL}/data/processed_scientific_data/*")
        os.system(f"rm -f {PACKAGE_TOP_LEVEL}/data/processed_position_data/*");   print(f"rm -f {PACKAGE_TOP_LEVEL}/data/processed_position_data/*")



    # For each satellite and each year
    print("")
    print(f"------------------- DOWNLOADING ALL DATA -------------------")
    print("")
    print(" Date         Data Type   Satellite   Downloaded   Processed")
    print("════════════════════════════════════════════════════════════")

    
    for year in ["2020", "2021", "2022"]:
        for sat in ["a", "b"]:
            url = f"https://data.elfin.ucla.edu/el{sat}/l2/epd/fast/electron/{year}/"
            html = str(requests.get(url).content)

            all_dates = re.findall('href="el[A-Za-z]_l2_epdef_[0-9]+_v01\.cdf"', html)
            for i in range(len(all_dates)):
                all_dates[i] = all_dates[i][19:-9] # Trim to only the yyyymmdd part
            N = len(all_dates)

            # Download data

            for date in all_dates:
                # Science data
                _get_science_data(date, sat)
                _get_position_data(date, sat)
                print("")

    
    print("\nDownloading metadata...")
    download_metadata()


def _get_science_data(date, sat):
# Checks data availability for a date, downloads it if available, and cleans it.
    # Parse input time
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    print(f"{year}/{month}/{day}    ", end = "")
    print("Science     ", end = "")
    print(f"ELFIN-{sat.upper()}     ", end = "")


    url = f"https://data.elfin.ucla.edu/el{sat}/l2/epd/fast/electron/{year}/el{sat}_l2_epdef_{date}_v01.cdf"
    raw_path = f"{PACKAGE_TOP_LEVEL}/data/raw/science_{date}_el{sat}.cdf" # Path to save the raw .cdf to
    processed_path = f"{PACKAGE_TOP_LEVEL}/data/processed_scientific_data/{date}_el{sat}.npz" # Path to save the cleaned .npz to

    # Don't redownload existing files
    if os.path.isfile(processed_path): 
        print("\033[92mYes          \033[0m", end = "")
        print("\033[92mYes\033[0m")
        return

    # Go ahead with downloading data
    r = requests.get(url)
    # If url doesn't work, warn
    if r.status_code != 200:
        print(f"\033[91m█████████ NO █████████\033[0m")
        return

    # Create raw file
    raw_file = open(raw_path, 'wb')
    raw_file.write(r.content)
    raw_file.close()
    print("\033[92mYes          \033[0m", end = "")

    # Create cleaned file
    _save_science_npz(raw_path, processed_path, sat)

def _save_science_npz(source_path, destination_path, sat):
# Packages a .cdf file into a .npz file for analysis in Julia.
    sat = "el" + sat
    file = pycdf.CDF(source_path)

    # Set flag indicating what satellite the data is from
    if sat == "ela":
        ela = True
        t = "T" # This is because (for some godforsaken reason), the "Tspin" field is "Tspin" on ELFIN-A and "tspin" on ELFIN-B, so we need to vary capitalization of that field based on satellite.
    else:
        ela = False
        t = "t"

    # Save processed .cdf as .npz for analysis in Julia
    np.savez(destination_path, 
        ela = ela,
        et_time          = _DateTime_to_string(file[f"{sat}_pef_et_time"][:]),
        hs_time          = _DateTime_to_string(file[f"{sat}_pef_hs_time"][:]),
        fs_time          = _DateTime_to_string(file[f"{sat}_pef_fs_time"][:]),
        Et_nflux         = file[f"{sat}_pef_Et_nflux"][:],
        Et_eflux         = file[f"{sat}_pef_Et_eflux"][:],
        Et_dfovf         = file[f"{sat}_pef_Et_dfovf"][:],
        energy_bins_mean = file[f"{sat}_pef_energies_mean"][:],
        energy_bins_min  = file[f"{sat}_pef_energies_min"][:],
        energy_bins_max  = file[f"{sat}_pef_energies_max"][:],
        pa               = file[f"{sat}_pef_pa"][:],
        spinphase        = file[f"{sat}_pef_spinphase"][:],
        sectnum          = file[f"{sat}_pef_sectnum"][:],
        Tspin            = file[f"{sat}_pef_{t}spin"][:],
        hs_Epat_nflux    = file[f"{sat}_pef_hs_Epat_nflux"][:],
        hs_Epat_eflux    = file[f"{sat}_pef_hs_Epat_eflux"][:],
        hs_Epat_dfovf    = file[f"{sat}_pef_hs_Epat_dfovf"][:],
        hs_LCdeg         = file[f"{sat}_pef_hs_LCdeg"][:],
        hs_antiLCdeg     = file[f"{sat}_pef_hs_antiLCdeg"][:],
        hs_epa_spec      = file[f"{sat}_pef_hs_epa_spec"][:],
        fs_Epat_nflux    = file[f"{sat}_pef_fs_Epat_nflux"][:],
        fs_Epat_eflux    = file[f"{sat}_pef_fs_Epat_eflux"][:],
        fs_Epat_dfovf    = file[f"{sat}_pef_fs_Epat_dfovf"][:],
        fs_LCdeg         = file[f"{sat}_pef_fs_LCdeg"][:],
        fs_antiLCdeg     = file[f"{sat}_pef_fs_antiLCdeg"][:],
        fs_epa_spec      = file[f"{sat}_pef_fs_epa_spec"][:],
        nspinsinsum      = file[f"{sat}_pef_nspinsinsum"][:],
        nsectors         = file[f"{sat}_pef_nsectors"][:],
        sect2add         = file[f"{sat}_pef_sect2add"][:],
        spinph2add       = file[f"{sat}_pef_spinph2add"][:]
    )
    print("\033[92mYes\033[0m")
    file.close()
    os.system(f"rm -f {source_path}") # Delete raw file

def _get_position_data(date, sat):
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    print(f"{year}/{month}/{day}    ", end = "")
    print("Position    ", end = "")
    print(f"ELFIN-{sat.upper()}     ", end = "")

    url = f"https://data.elfin.ucla.edu/el{sat}/l1/state/defn/{year}/el{sat}_l1_state_defn_{date}_v02.cdf"    
    raw_path = f"{PACKAGE_TOP_LEVEL}/data/raw/{url[-29:]}" # Path to save the raw .cdf to
    processed_path = f"{PACKAGE_TOP_LEVEL}/data/processed_position_data/{date}_el{sat}.npz" # Path to save the cleaned .npz to

    # Don't redownload existing files
    if os.path.isfile(processed_path): 
        print("\033[92mYes          \033[0m", end = "")
        print("\033[92mYes\033[0m")
        return
    
    # Go ahead with downloading data
    r = requests.get(url)
    # If url doesn't work, warn
    if r.status_code != 200:
        print(f"\033[91m█████████ NO █████████\033[0m")
        return

    # Create raw file
    raw_file = open(raw_path, 'wb')
    raw_file.write(r.content)
    raw_file.close()
    print("\033[92mYes          \033[0m", end = "")

    # Create cleaned file
    _save_position_npz(raw_path, processed_path, sat)

def _save_position_npz(source_path, destination_path, sat):
    sat = "el" + sat
    file = pycdf.CDF(source_path)
    print("...\b\b\b", end = "")
    
    # Convert position to L & MLT
    Re = 6378 # Earth radius, km
    time = file[f"{sat}_state_time"][::10] # Take every 10th index (L-shell readings every 10 seconds)
    position = file[f"{sat}_pos_gei"][::10]
    altitude = (position[:,0]**2 + position[:,1]**2 + position[:,2]**2)**(1/2) - Re # km
    position = position / Re # Divide by radius of Earth to convert km to Re

    # Start threads and process data
    threads = Pool(N_THREADS)
    L_MLT_results = threads.starmap(_single_thread_calculate_L_and_MLT, zip(time, position))
    L   = [i[0] for i in L_MLT_results]
    MLT = [i[1] for i in L_MLT_results]

    # Save results
    np.savez(destination_path,
        state_time = _DateTime_to_string(time),
        L          = L,
        MLT        = MLT,
        altitude   = altitude,
        position   = position
    )

    print("\033[92mYes\033[0m")
    file.close()
    os.system(f"rm -f {source_path}") # Delete raw file


def _single_thread_calculate_L_and_MLT(time, position):
    time = spacepy.time.Ticktock(time, 'UTC')
    position = coordinates.Coords(position, 'GEI', 'car', units = ['Re', 'Re', 'Re'], use_irbem = False)
    result = irbempy.get_Lm(time, position, [90], extMag = 'T87LONG') # https://spacepy.github.io/irbempy.html
    
    L = result['Lm'][0]
    MLT = result['MLT'][0]

    return [L, MLT]

def download_metadata():
    for sat in ["a", "b"]:
        url = f"https://data.elfin.ucla.edu/el{sat}/data_availability/el{sat}_epd_all.csv"
        r = requests.get(url)

        if r.status_code != 200:
            print(f"{url} did not respond correctly. Please check url or try again later.\033[0m")
            sys.exit(1)
            
        file = open(f"{PACKAGE_TOP_LEVEL}/data/el{sat}_epd_all.csv", 'wb')
        file.write(r.content)
        file.close()

    geomag_indices_urls = ["https://data.elfin.ucla.edu/kp/elfin_kp.csv", "https://data.elfin.ucla.edu/dst/elfin_dst.csv"]
    geomag_names = ["kp", "dst"]

    for i in range(len(geomag_indices_urls)):
        url = geomag_indices_urls[i]
        name = geomag_names[i]
        r = requests.get(url)

        if r.status_code != 200:
            print(f"{url} did not respond correctly. Please check url or try again later.\033[0m")
            sys.exit(1)

        file = open(f"{PACKAGE_TOP_LEVEL}/data/{name}.csv", 'wb')
        file.write(r.content)
        file.close()



def _url_exists(url):
# Check if a url returns HTTP 200 code (OK)
    r = requests.head(url)
    if r.status_code == 200: # URL was found & request was successful
        return True
    else:
        return False # Anything else indicates an error


def _DateTime_to_string(raw_time):
# Converts a DateTime object to a string. Useful since Julia doesn't understand Python DateTime objects.
    clean_time = []
    for i in range(0, len(raw_time)):
        datetime_string = raw_time[i].strftime("%Y.%m.%d.%H.%M.%S.%f")
        clean_time.append(datetime_string[0:-3]) # Save only up to millisecond precision, as Julia DateTime doesn't support finer time resolutions.
    return clean_time


if __name__ == "__main__": # For thread safety
# Executed at runtime. Downloads all ELFIN data.
    delete_prior = False
    command_line_args = sys.argv
    if len(command_line_args) > 1:
        if command_line_args[1] == "-delete":
            delete_prior = True

    with keep.presenting(): # keeps screen awake
        download_all_elfin_data(delete_prior)