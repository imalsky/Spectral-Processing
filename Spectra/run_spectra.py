#!/usr/bin/env python

from shutil import copyfile
from sys import exit
import os
import sys
import pandas as pd
import numpy as np
import run_grid
import altitude_regridding
import convert_fort_files


# Phases in degrees, inclination in radians (sorry)
# An inclination of 0 corresponds to edge on
phases = [22.5]
inclinations = [0.0]
system_obliquity = 0

# I recommend leaving these as is
# The NLAT and NLON can be changed, but these values work well
INITIAL_NTAU = 65
NTAU = 250

# Please don't touch these
NLAT = 48
NLON = 96

# Whethere there are clouds
# 0 is no clouds, 1 is clouds
# This is also important for filling in the correct number of 0s for the input files
CLOUDS = 0

# 0 is off
# 1 is everything
# 2 is Wind only
# 3 is rotation only
dopplers = [0,1]

# If you only need to change the phase you can use this knob
# It skips a lot of steps for the regridding
ONLY_PHASE = True

# If you only have the fort files use this
# Please don't only have the fort files
# It requires that the fort files are named according to something particular
# In the correct directory
USE_FORT_FILES = True

# There are low resolution spectra and high resolution spectra that can be created
# There are somethings that need to be changed in the template inputs file to make this happen
# If you change the underlying data files these might need to be changed
high_res = False

# These are the planet files that you need to run the code
# So These should be in New_Jups/Planets
# They should be pretty big files, and don't include the .txt with the names here
planet_name = 'Hayley-Tests'


# This is specifically for the regridding
# You can change the grid density
# I wouldn't mess with this though
grid_lat_min = -87.16
grid_lat_max = 87.16
grid_lon_min = 0.0
grid_lon_max = 356.25

def add_zeros(planet_name):
    data_file = '../Planets/' + planet_name + '.txt'
    df = pd.read_csv(data_file, delim_whitespace=True, names=('lat', 'lon', 'level', 'alt', 'pres', 'temp', 'u', 'v', 'w',
                                                                       'total aerosol tau_lw', 'mean_aerosol_g0_lw', 'mean_aerosol_pi0_lw'))

    df.drop(columns=['total aerosol tau_lw', 'mean_aerosol_g0_lw', 'mean_aerosol_pi0_lw'])

    df['aero_sw_tau_1'] = 0
    df['aero_sw_tau_2'] = 0
    df['aero_sw_tau_3'] = 0
    df['aero_sw_tau_4'] = 0

    df['sw_asym_1'] = 0
    df['sw_asym_2'] = 0
    df['sw_asym_3'] = 0
    df['sw_asym_4'] = 0

    df['sw_pi0_1'] = 0
    df['sw_pi0_2'] = 0
    df['sw_pi0_3'] = 0
    df['sw_pi0_4'] = 0

    df = df[['lat', 'lon', 'level',
             'alt', 'pres', 'temp',
             'u', 'v', 'w',
             'aero_sw_tau_1', 'sw_asym_1', 'sw_pi0_1',
             'aero_sw_tau_2', 'sw_asym_2', 'sw_pi0_2',
             'aero_sw_tau_3', 'sw_asym_3', 'sw_pi0_3',
             'aero_sw_tau_4', 'sw_asym_4', 'sw_pi0_4']]

    np.savetxt(data_file, df,
            fmt='%12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E\t')
    return None


def get_run_lists(phases, inclinations):
    for phase in phases:
        for inc in inclinations:
            phase = str(phase)
            inc = str(inc)

            input_paths.append('DATA/init_' + planet_name + '_phase_{}_inc_{}.txt'.format(phase, inc))
            inclination_strs.append(inc)
            phase_strs.append(phase)

    return input_paths, inclination_strs, phase_strs


def run_exo(input_paths, inclination_strs, phase_strs, doppler_val):
    """
    This runs Eliza's code
    """
    inputs_file = 'input.h'
    output_paths = []

    # The output paths should be similar to the input paths
    # Minus the .dat file extension and saved to OUT/
    for file_path in input_paths:
        output_paths.append('OUT/Spec_' + str(doppler_val) + '_' + file_path[10:-4])

    # Each Run needs to have a specific input.h file
    # With the correct input and output paths
    for i in range(len(input_paths)):
        output_temp = output_paths[i]
        input_temp  = input_paths[i]
        
        # Copy the template for inputs
        try:
            copyfile('template_inputs.h', inputs_file)
        except IOError as e:
            print("Unable to copy file. %s" % e)
            exit(1)
        except:
            print("Unexpected error:", sys.exc_info())
            exit(1)
        
        # Read in the file
        with open(inputs_file, 'r') as file :
            filedata = file.read()

        # Replace the input and output paths
        filedata = filedata.replace("<<output_file>>", "\"" + output_temp + str(W0_VAL) + str(G0_VAL) +"\"")
        filedata = filedata.replace("<<input_file>>", "\"" + input_temp + "\"")
        filedata = filedata.replace("<<doppler>>", str(doppler_val))
        filedata = filedata.replace("<<inclination>>", inclination_strs[i])
        filedata = filedata.replace("<<phase>>", phase_strs[i])

        filedata = filedata.replace("<<CLOUDS>>", str(CLOUDS))
        filedata = filedata.replace("<<NTAU>>", str(NTAU))
        filedata = filedata.replace("<<NLAT>>", str(NLAT))
        filedata = filedata.replace("<<NLON>>", str(NLON))

        filedata = filedata.replace("<<W0_VAL>>", str(W0_VAL))
        filedata = filedata.replace("<<G0_VAL>>", str(G0_VAL))

        # This is the part that changes the low-res vs high res
        if high_res == True:
            filedata = filedata.replace("<<num_pressure_points>>", "17")

            filedata = filedata.replace("<<CHEM_FILE>>", "\"DATA/eos_solar_doppler.dat\"")
            filedata = filedata.replace("<<CH4_FILE>>",  "\"DATA/opacCH4_hires.dat\"")
            filedata = filedata.replace("<<CO2_FILE>>",  "\"DATA/opacCO2_hires.dat\"")
            filedata = filedata.replace("<<CO_FILE>>",   "\"DATA/opacCO_hires.dat\"")
            filedata = filedata.replace("<<H2O_FILE>>",  "\"DATA/opacH2O_hires.dat\"")
            filedata = filedata.replace("<<NH3_FILE>>",  "\"DATA/opacNH3_hires.dat\"")
            filedata = filedata.replace("<<O2_FILE >>",  "\"DATA/opacO2_hires.dat\"")
            filedata = filedata.replace("<<O3_FILE>>",   "\"DATA/opacO3_hires.dat\"")
        else:
            filedata = filedata.replace("<<num_pressure_points>>", "13")

            filedata = filedata.replace("<<CHEM_FILE>>", "\"DATA/eos_solar_doppler_2016_cond.dat\"")
            filedata = filedata.replace("<<CH4_FILE>>",  "\"DATA/opacCH4.dat\"")
            filedata = filedata.replace("<<CO2_FILE>>",  "\"DATA/opacCO2.dat\"")
            filedata = filedata.replace("<<CO_FILE>>",   "\"DATA/opacCO.dat\"")
            filedata = filedata.replace("<<H2O_FILE>>",  "\"DATA/opacH2O.dat\"")
            filedata = filedata.replace("<<NH3_FILE>>",  "\"DATA/opacNH3.dat\"")
            filedata = filedata.replace("<<O2_FILE >>",  "\"DATA/opacO2.dat\"")
            filedata = filedata.replace("<<O3_FILE>>",   "\"DATA/opacO3.dat\"")

        # Write the file out again
        with open(inputs_file, 'w') as file:
            file.write(filedata)
        
        # Run Eliza's code
        os.system('make clean')
        os.system('make rt_emission_aerosols.exe') 
        os.system('./rt_emission_aerosols.exe')

    return None


input_paths = []              
output_paths = []
inclination_strs = []
phase_strs = []


# Convert the fort files to the correct format
if USE_FORT_FILES == True:
    convert_fort_files.convert_to_correct_format('', planet_name)
    print ("Converted the fort files to the new format")
else:
    pass


# Regrid the file to constant altitude and the correct number of layers
altitude_regridding.regrid_gcm_to_constant_alt(CLOUDS, planet_name, NLAT, NLON, INITIAL_NTAU, NLON, NTAU)
print ("Regridded the planet to constant altitude")


# Add the zeros if the initial file doesn't have them
if CLOUDS == 0:
    add_zeros(planet_name)
    print ("Added a bunch of zeros because the code assumes the initial file has no clouds")
    print ("The spectral processing needs the number of columns to be correct")
else:
    pass


# If you already have the Final planet file creates you can commend out run_grid and double planet file
run_grid.run_all_grid(planet_name, phases, inclinations, system_obliquity, NTAU, NLAT, NLON, grid_lat_min, grid_lat_max, grid_lon_min, grid_lon_max, ONLY_PHASE)


# Get all the files that you want to run
input_paths, inclination_strs, phase_strs = get_run_lists(phases, inclinations)


print ("here!!!!!")
print (input_paths)
# If you want to manually set these values you can leave them here
# Normally they will not affect it, unless you manually set them in two_stream.h
W0_VALS = [0.0]
G0_VALS = [0.0]

for G0_VAL in G0_VALS:
    for W0_VAL in W0_VALS:
        for doppler_val in dopplers:
            run_exo(input_paths, inclination_strs, phase_strs, doppler_val)