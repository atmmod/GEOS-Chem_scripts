#!/usr/bin/env python3

# Complete Automated testing framework script for the GEOS-Chem-hyd manuscript

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from matplotlib.ticker import ScalarFormatter
from scipy import stats
import warnings
import gcpy.plot as gcplot
from cmcrameri import cm
# from gcpy.units import check_units, data_unit_is_mol_per_mol
import os
from joblib import Parallel, delayed
from matplotlib import colors
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap

# Ignore warnings
warnings.filterwarnings("ignore")

# Reading Files

root_path = '/glade/work/sakinjole/rundirs/gc_4x5_merra2_fullchem'
species = 'NO'
specie1 = 'NO'
specie2 = 'O3'
duration = '20min'
mechanism = 'fullchem_hemco_1'
mechanism2 = 'fullchem_quad'
pert_layer = '/allcells'  # else ground layer
ptype = 'H'

h1 = 0.01
hc = 0.01
h2 = 0.0001

print(mechanism)
# A dictionary for species variables and indices to plot (This serves as the numerator)
# Only for the species concentration diagnostic
# plot_species = {'CO':251, 'CO2':86, 'SO4':94, 'SO2':95, 'NO2':131, 'NO':132, 'O3':127, 'NH4': 135, 'NH3': 136}
plot_species = {'CO': 251, 'SO4': 94, 'SO2': 95, 'NO2': 131, 'NO': 132, 'O3': 127,
                'NH4': 135, 'NH3': 136, 'NIT': 134, 'HNO3': 217, 'HCl': 224, 'ISOP': 173}

#hemco_species = {'EmisSO4_Total': 3, 'EmisNO2_Total': 10, 'EmisNO_Total': 11,
#                 'EmisNH3_Total': 12, 'EmisCO_Total': 20}
                 
hemco_species = {'EmisO3_Total': 10, 'EmisNO2_Total': 11, 'EmisNO_Total': 12, 'EmisHNO3_Total': 19}



# Avoid using gcpy for this:


def data_unit_is_mol_per_mol(da):
    """
    Check if the units of an xarray DataArray are mol/mol based on a set
    list of unit strings mol/mol may be.

    Args:
        da: xarray DataArray
            Data array containing a units attribute

    Returns:
        is_molmol: bool
            Whether input units are mol/mol
    """
    conc_units = ["mol mol-1 dry", "mol/mol", "mol mol-1"]
    is_molmol = False
    if da.units.strip() in conc_units:
        is_molmol = True
    return is_molmol


def default_path(mechanism):
    # Folder paths default run
    real_path1 = root_path + '/default/' + mechanism2 + '/' + \
        duration + '/GEOSChem.SpeciesConc.20190701_0000z.nc4'
    hyd_path1 = root_path + '/default_hyd/' + mechanism2 + '/' + \
        duration + '/GEOSChem.SpeciesConc.20190701_0000z.nc4'
    real_path2 = root_path + '/default/' + mechanism2 + '/' + \
        duration + '/GEOSChem.AerosolMass.20190701_0000z.nc4'
    hyd_path2 = root_path + '/default_hyd/' + mechanism2 + '/' + \
        duration + '/GEOSChem.AerosolMass.20190701_0000z.nc4'
    # Read the arrays:
    real1 = xr.open_dataset(real_path1)
    hyd1 = xr.open_dataset(hyd_path1)
    real2 = xr.open_dataset(real_path2)
    hyd2 = xr.open_dataset(hyd_path2)
    return real1, hyd1, real2, hyd2
    
def def_path(mechanism):
    # Folder paths default run
    real_path = root_path + '/first_output/' + mechanism + '/' + \
        duration + pert_layer + '/def' + specie1 + '/HEMCO_diagnostics.201907010000.nc'
    # Read the arrays:
    real = xr.open_dataset(real_path)
    return real


def main_path(ty, layer=''):  # where ty stands for the order of operation
    if ty == 'first':
        if layer == '/allcells':
            # Folder paths first order - updated to central difference
            real_path = '/first_output/' + mechanism + '/' + \
                duration + layer + '/dec' + specie1 + '/'

            # pertreal_path = '/first_output/incNO/'
            pertreal_path = '/first_output/' + mechanism + \
                '/' + duration + layer + '/inc' + specie1 + '/'
            # perthyd_path = '/second_output_hyd/perturbedNO_boxN/'
            perthyd_path = '/first_output_hyd/' + mechanism + \
                '/' + duration + layer + '/hyd' + specie1 + '/'
        else:
            # Folder paths first order - updated to central difference
            real_path = '/first_output/' + mechanism + '/' + \
                duration + '/ground/dec' + specie1 + '/'

            # pertreal_path = '/first_output/incNO/'
            pertreal_path = '/first_output/' + mechanism + \
                '/' + duration + '/ground/inc' + specie1 + '/'
            # perthyd_path = '/second_output_hyd/perturbedNO_boxN/'
            perthyd_path = '/first_output_hyd/' + mechanism + \
                '/' + duration + '/ground/hyd' + specie1 + '/'
    elif ty == 'second':
        if layer == '/allcells':

            # Folder paths second order
            real_path = '/second_output/' + mechanism + '/' + \
                duration + layer + '/dec' + specie1 + '/'

            pertreal_path = '/second_output/' + mechanism + \
                '/' + duration + layer + '/inc' + specie1 + '/'
            perthyd_path = '/second_output_hyd/' + mechanism + \
                '/' + duration + layer + '/hyd' + specie1 + '/'
        else:
            # Folder paths second order
            real_path = '/second_output/' + mechanism + \
                '/' + duration + '/ground/dec' + specie1 + '/'

            pertreal_path = '/second_output/' + mechanism + \
                '/' + duration + '/gound/inc' + specie1 + '/'
            perthyd_path = '/second_output_hyd/' + mechanism + \
                '/' + duration + '/ground/hyd' + specie1 + '/'

    elif ty == 'cross':
        # Folder paths hybrid
        if layer == '/allcells':
            real_path = '/cross_output/' + mechanism + '/' + \
                duration + layer + '/dec' + specie1 + '_' + specie2 + '/'
        # hyd_path = '/hybrid_output_hyd/default'

            pertreal_path = '/cross_output/' + mechanism + '/' + \
                duration + layer + '/inc' + specie1 + '_' + specie2 + '/'
            perthyd_path = '/cross_output_hyd/' + mechanism + '/' + \
                duration + layer + '/hyd' + specie1 + '_' + specie2 + '/'
    return real_path, pertreal_path, perthyd_path


# File paths
def file(ty, pt):  # where ty is the diagnostic type, pt is the required path
    if ty == 'A':
        if pt == 'real':
            return 'GEOSChem.AerosolMass.20190701_0000z.nc4'
        elif pt == 'dx1':
            return 'GEOSChem.AerosolMassdx1.20190701_0000z.nc4'
        elif pt == 'dx2':
            return 'GEOSChem.AerosolMassdx2.20190701_0000z.nc4'
        elif pt == 'dx1x2':
            return 'GEOSChem.AerosolMassdx1x2.20190701_0000z.nc4'
        elif pt == 'isens':
            return 'GEOSChem.Sensitivity.20190701_0000z.nc4'

    elif ty == 'S':
        if pt == 'real':
            return 'GEOSChem.SpeciesConc.20190701_0000z.nc4'
        elif pt == 'dx1':
            return 'GEOSChem.SpeciesConcdx1.20190701_0000z.nc4'
        elif pt == 'dx2':
            return 'GEOSChem.SpeciesConcdx2.20190701_0000z.nc4'
        elif pt == 'dx1x2':
            return 'GEOSChem.SpeciesConcdx1x2.20190701_0000z.nc4'
        elif pt == 'isens':
            return 'GEOSChem.Sensitivity.20190701_0000z.nc4'

    elif ty == 'H':
        if pt == 'real':
            return 'HEMCO_diagnostics.201907010000.nc'
        elif pt == 'dx1':
            return 'HEMCOdx1_diagnostics.201907010000.nc'
        elif pt == 'dx2':
            return 'HEMCOdx2_diagnostics.201907010000.nc'
        elif pt == 'dx1x2':
            return 'HEMCOdx1x2_diagnostics.201907010000.nc'


def read_file(loc):
    '''Read the real files with xarray'''
    return xr.open_dataset(loc)


# gets all file paths for a particular diagnostic
def get_file(array, order, real_path, pertreal_path, perthyd_path):
    if array == 'A':
        if order == 'first':
            ref_real = read_file(root_path + real_path + file('A', 'real'))
            ref_pert = read_file(root_path + pertreal_path + file('A', 'real'))
            hyd_dx1 = read_file(root_path + perthyd_path + file('A', 'dx1'))
            hyd_dx2 = read_file(root_path + perthyd_path + file('A', 'dx2'))
            hyd_dx1x2 = read_file(
                root_path + perthyd_path + file('A', 'dx1x2'))
#            ipercent = read_file(root_path + perthyd_path + file('A', 'isens'))
        elif order == 'second':
            ref_real = read_file(root_path + real_path + file('A', 'dx2'))
            ref_pert = read_file(root_path + pertreal_path + file('A', 'dx2'))
            hyd_dx1 = read_file(root_path + perthyd_path + file('A', 'dx1'))
            hyd_dx2 = read_file(root_path + perthyd_path + file('A', 'dx2'))
            hyd_dx1x2 = read_file(
                root_path + perthyd_path + file('A', 'dx1x2'))
            # ipercent = read_file(root_path + perthyd_path + file('A', 'isens'))

    elif array == 'S':
        if order == 'first':
            ref_real = read_file(root_path + real_path + file('S', 'real'))
            ref_pert = read_file(root_path + pertreal_path + file('S', 'real'))
            hyd_dx1 = read_file(root_path + perthyd_path + file('S', 'dx1'))
            hyd_dx2 = read_file(root_path + perthyd_path + file('S', 'dx2'))
            hyd_dx1x2 = read_file(
                root_path + perthyd_path + file('S', 'dx1x2'))
#            ipercent = read_file(root_path + perthyd_path + file('S', 'isens'))
#            ipercentb = read_file(root_path + real_path + file('S', 'isens'))
#            ipercentf = read_file(root_path + pertreal_path + file('S', 'isens'))
        elif order == 'second':
            ref_real = read_file(root_path + real_path + file('S', 'dx2'))
            ref_pert = read_file(root_path + pertreal_path + file('S', 'dx2'))
            hyd_dx1 = read_file(root_path + perthyd_path + file('S', 'dx1'))
            hyd_dx2 = read_file(root_path + perthyd_path + file('S', 'dx2'))
            hyd_dx1x2 = read_file(
                root_path + perthyd_path + file('S', 'dx1x2'))
#            ipercent = read_file(root_path + perthyd_path + file('S', 'isens'))
#            ipercentb = read_file(root_path + real_path + file('S', 'isens'))
#            ipercentf = read_file(root_path + pertreal_path + file('S', 'isens'))
        elif order == 'cross':
            ref_real = read_file(root_path + real_path + file('S', 'dx2'))
            ref_pert = read_file(root_path + pertreal_path + file('S', 'dx2'))
            hyd_dx1 = read_file(root_path + perthyd_path + file('S', 'dx1'))
            hyd_dx2 = read_file(root_path + perthyd_path + file('S', 'dx2'))
            hyd_dx1x2 = read_file(
                root_path + perthyd_path + file('S', 'dx1x2'))
#            ipercent = read_file(root_path + perthyd_path + file('S', 'isens'))
#            ipercentb = read_file(root_path + real_path + file('S', 'isens'))
#            ipercentf = read_file(root_path + pertreal_path + file('S', 'isens'))

    elif array == 'H':
        if order == 'first':
            ref_real = read_file(root_path + real_path + file('H', 'real'))
            ref_pert = read_file(root_path + pertreal_path + file('H', 'real'))
            hyd_dx1 = read_file(root_path + perthyd_path + file('H', 'dx1'))
            hyd_dx2 = read_file(root_path + perthyd_path + file('H', 'dx2'))
            hyd_dx1x2 = read_file(
                root_path + perthyd_path + file('H', 'dx1x2'))

        elif order == 'second':
            ref_real = read_file(root_path + real_path + file('H', 'dx2'))
            ref_pert = read_file(root_path + pertreal_path + file('H', 'dx2'))
            hyd_dx1 = read_file(root_path + perthyd_path + file('H', 'dx1'))
            hyd_dx2 = read_file(root_path + perthyd_path + file('H', 'dx2'))
            hyd_dx1x2 = read_file(
                root_path + perthyd_path + file('H', 'dx1x2'))

    return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2
    # return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, ipercent, ipercentb, ipercentf


def full_path(order, array):
    if order == 'first':
        real_path, pertreal_path, perthyd_path = main_path('first', pert_layer)
    elif order == 'second':
        real_path, pertreal_path, perthyd_path = main_path(
            'second', pert_layer)
    elif order == 'cross':
        real_path, pertreal_path, perthyd_path = main_path('cross', pert_layer)

    # save the paths to be used later for the one to one plots
    ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2 = get_file(
        array, order, real_path, pertreal_path, perthyd_path)
#    ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2 = get_file(array, order, real_path, pertreal_path, perthyd_path)
#    ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, ipercent, ipercentb, ipercentf = get_file(array, order, real_path, pertreal_path, perthyd_path)
    # return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, real_path, pertreal_path, perthyd_path
    return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, perthyd_path
    # return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, perthyd_path, ipercent, ipercentb, ipercentf

# Run the main path here:


# Running all plots sequentially:
if ptype == 'S':
    # ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, real_pathf, pertreal_pathf, perthyd_pathf = full_path('first', 'S')
    #    ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, perthyd_pathf, ipercent, ipercentb, ipercentf = full_path('first', 'S')
    ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, perthyd_pathf = full_path(
        'first', 'S')
# ref_realfa, ref_pertfa, hyd_dx1fa, hyd_dx2fa, hyd_dx1x2fa, perthyd_pathfa, ipercenta = full_path('first', 'A')
#    ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, perthyd_paths, ipercent2, ipercent2b, ipercent2f = full_path('second', 'S')
    ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, perthyd_paths = full_path(
        'second', 'S')
# ref_realc, ref_pertc, hyd_dx1c, hyd_dx2c, hyd_dx1x2c, real_pathc, pertreal_pathc, perthyd_pathc = full_path('cross', 'S')
elif ptype == 'A':
    # ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, perthyd_pathf, ipercent = full_path('first', 'A')
    ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, perthyd_pathf = full_path(
        'first', 'A')
    # ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, perthyd_paths, ipercent2 = full_path('second', 'A')
    ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, perthyd_paths = full_path(
        'second', 'A')

elif ptype == 'H':
    ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, perthyd_pathf = full_path(
        'first', 'H')
    ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, perthyd_paths = full_path(
        'second', 'H')
        
    #Run again with the new denominator, only needed once:
    old = specie1
    specie1 = specie2
    mechanism = 'fullchem_hemco_0.01'
    ref_realf2, ref_pertf2, hyd_dx1f2, hyd_dx2f2, hyd_dx1x2f2, perthyd_pathf2 = full_path(
        'first', 'H')
    ref_reals2, ref_perts2, hyd_dx1s2, hyd_dx2s2, hyd_dx1x2s2, perthyd_paths2 = full_path(
        'second', 'H')
    
    #Reassign specie1 to old name
    specie1 = old
    
        
    #Default file path
    #ref_reald = def_path(mechanism)
# full_path('second', 'A')
# full_path('hybrid', 'A')

# Get the default files
# default_file1, hyd_file1, default_file2, hyd_file2 = default_path(mechanism)

# add unit conversions to pbb from mol/mol


def unit_convert(ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, useless):
    for i in ref_real:
        if i in useless:
            continue
        else:
            if data_unit_is_mol_per_mol(ref_real[i]):
                ref_real[i].values = ref_real[i].values * 1e9
                ref_real[i].attrs["units"] = "ppb"
    for i in ref_pert:
        if i in useless:
            continue
        else:
            if data_unit_is_mol_per_mol(ref_pert[i]):
                ref_pert[i].values = ref_pert[i].values * 1e9
                ref_pert[i].attrs["units"] = "ppb"
    for i in hyd_dx1:
        if i in useless:
            continue
        else:
            if data_unit_is_mol_per_mol(hyd_dx1[i]):
                hyd_dx1[i].values = hyd_dx1[i].values * 1e9
                hyd_dx1[i].attrs["units"] = "ppb"
    for i in hyd_dx2:
        if i in useless:
            continue
        else:
            if data_unit_is_mol_per_mol(hyd_dx2[i]):
                hyd_dx2[i].values = hyd_dx2[i].values * 1e9
                hyd_dx2[i].attrs["units"] = "ppb"
    for i in hyd_dx1x2:
        if i in useless:
            continue
        else:
            if data_unit_is_mol_per_mol(hyd_dx1x2[i]):
                hyd_dx1x2[i].values = hyd_dx1x2[i].values * 1e9
                hyd_dx1x2[i].attrs["units"] = "ppb"
    return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2

# add unit conversions to pbb from mol/mol for one argument


def units_convert(array, useless):
    for i in array:
        if i in useless:
            continue
        else:
            if data_unit_is_mol_per_mol(array[i]):
                array[i].values = array[i].values * 1e9
                array[i].attrs["units"] = "ppb"
    return array

# Calculate finite differences

# Calculating the defined perturbation size for a variable


# To create an array of 0s to map each variable and iterate through
lev = 72
lat = 46
lon = 72

# Applies for both Species_Conc and aerosolMass
# useless variables in the array:
useless = ['lat_bnds', 'lon_bnds', 'hyam',
           'hybm', 'hyai', 'hybi', 'P0', 'AREA']
num_var = len(ref_realf) - 8
# num_var = len(ref_reals) - 8


# comment out if not needed
# ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2 = unit_convert(ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, useless)

# Running all sequentially:
ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f = unit_convert(
    ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, useless)
ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s = unit_convert(
    ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, useless)
# ref_realc, ref_pertc, hyd_dx1c, hyd_dx2c, hyd_dx1x2c = unit_convert(ref_realc, ref_pertc, hyd_dx1c, hyd_dx2c, hyd_dx1x2c, useless)

# Convert the default files to ppb:
# default_file1 = units_convert(default_file1, useless)
# default_file2 = units_convert(default_file2, useless)
# hyd_file1 = units_convert(hyd_file1, useless)
# hyd_file2 = units_convert(hyd_file2, useless)


def first_order_finite(base, pert, h, pert_specie, ptype, useless, mode):
    # Index to the species of interest
    # where h is the perturbation, pert_specie is the perturbed variable
    if ptype == 'H':
        # number of variables, ignore useless variables
        species_num = len(base) - 4
    else:
        # number of variables, ignore useless variables
        species_num = len(base) - 8
    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape

    ct = 0  # counter
    specie = []  # list of species names

    for i in base:
        if i in useless:
            continue
        # if mode == 'A': #post process for units to ug/m3
        else:

            if mode == 'semi':
                array[ct, :, :, :, :] = (
                    pert[i] - base[i]) / (h-1)  # multiplicative
            elif mode == 'normal':
                array[ct, :, :, :, :] = (pert[i] - base[i]) / h
            # add variable name to list:
            specie.append(i)

        ct += 1

    return array, shape, specie


def central_first_order(dec, inc, h, pert_specie, ptype, useless, mode):
    # Index to the species of interest
    # where h is the perturbation, pert_specie is the perturbed variable

    species_num = 0
    if ptype == 'H':
        # First, count the number of variables that satisfy the constraints
        for i in dec:
            # Skip useless variables
            if i in useless:
                continue

            # For 'H' type, only consider variables with 'Total' in the name, or consider only ship
            #if ptype == 'H' and 'Total' not in i:
#            if ptype == 'H' and 'Ship' not in i:
#                continue

#            # Count only 4D variables
#            if base[i].ndim == 4:
            # Count only 3D variables
            if dec[i].ndim == 4:
                species_num += 1
    else:
#        # number of variables, ignore useless variables
#    if ptype == 'H':
#        species_num = len(base) - 4
#    else:
        species_num = len(dec) - 8
    time, lev, lat, lon = dec[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
#    time, lat, lon = base[pert_specie].shape
#    # array holding all variables
#    array = np.zeros([species_num, time, lat, lon])
    shape = array.shape

    ct = 0  # counter
    specie = []  # list of species names

    for i in dec:
        if i in useless:
            continue

        elif dec[i].ndim == 4:
            if mode == 'semi':
                array[ct, :, :, :, :] = (
                    inc[i] - dec[i]) / (2 * h)  # multiplicative
            elif mode == 'normal':
                array[ct, :, :, :, :] = (inc[i] - dec[i]) / h
            # add variable name to list:
            specie.append(i)
            ct += 1

    return array, shape, specie


# Calculating Hyperdual sensitivities

def hyd_first_order(base, pert_specie, ptype, mode, useless, hyd_pert=1):

    species_num = 0
    if ptype == 'H':
        # First, count the number of variables that satisfy the constraints
        for i in base:
            # Skip useless variables
            if i in useless:
                continue

#            # Count only 4D variables
#            if base[i].ndim == 4:

            if base[i].ndim == 4:
                species_num += 1
    else:

        species_num = len(base) - 8
    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])

    shape = array.shape

    ct = 0  # counter
    specie = []  # list of species names
    for i in base:
        if i in useless:
            continue

        elif base[i].ndim == 4:
            array[ct, :, :, :, :] = base[i][:, :, :, :] / \
                hyd_pert  # semi normalization (multiplicative)
            # add variable name to list:
            specie.append(i)
            ct += 1

    return array, shape, specie


def hybrid_second_order(base, pert, pert_specie, ptype, mode, useless, hc, h=1):
    # where pert represents the diagnostic with both real and dual perturbation, while base is the
    # only dual perturbation

    if ptype == 'H':
        # number of variables, ignore useless variables
        species_num = len(base) - 4
    else:
        # number of variables, ignore useless variables
        species_num = len(base) - 8
    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape

    # perturbation size

    ct = 0  # counter
    specie = []  # list of species names
    for i in base:
        if i in useless:
            continue

        else:
            if mode == 'semi':
                # array[ct, :, :, :, :] = (pert[i][:,:,:,:] - base[i][:,:,:,:]) /((hc-1)*h)
                array[ct, :, :, :, :] = (
                    pert[i][:, :, :, :] / (hc*h*(hc-1))) - (base[i][:, :, :, :] / ((hc-1)*h))
            elif mode == 'normal additive':
                array[ct, :, :, :, :] = (
                    (pert[i][:, :, :, :] / h) - (base[i][:, :, :, :] / h)) / h1
            # add variable name to list:
            specie.append(i)
        ct += 1

    return array, shape, specie

#dec: the netcdf array corresponding to the dx2 component of the decreased case
#inc: the netcdf array corresponding to the dx2 component of the increased case
#pert_specie: a variable name of the array
def central_hybrid_second_order(dec, inc, pert_specie, ptype, mode, useless, hc, h=1):
    # where pert represents the diagnostic with both real and dual perturbation, while base is the
    # only dual perturbation

    species_num = 0
    if ptype == 'H':
        # First, count the number of variables that satisfy the constraints
        for i in dec:
            # Skip useless variables
            if i in useless:
                continue


            if dec[i].ndim == 4:
                species_num += 1
    else:

        species_num = len(dec) - 8
    time, lev, lat, lon = dec[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
#    time, lat, lon = base[pert_specie].shape
#    # array holding all variables
#    array = np.zeros([species_num, time, lat, lon])
    shape = array.shape

    # perturbation size

    ct = 0  # counter
    specie = []  # list of species names
    for i in dec:
        if i in useless:
            continue

        elif dec[i].ndim == 4:
            if mode == 'semi':
                array[ct, :, :, :, :] = (
                    inc[i][:, :, :, :] / (2*hc*h*(1+hc))) - (dec[i][:, :, :, :] / (2*h*hc*(1-hc)))
            # elif mode == 'normal additive':
#                 array[ct, :, :, :, :] = ((pert[i][:,:,:,:] / h) - (base[i][:,:,:,:] /h)) / h1
            # add variable name to list:
            specie.append(i)
            ct += 1

    return array, shape, specie


def central_finite_second_order(dec, base, inc, pert_specie, ptype, mode, useless, hc, h=1):
    # where pert represents the diagnostic with both real and dual perturbation, while base is the
    # only dual perturbation

    species_num = 0
    if ptype == 'H':
        # First, count the number of variables that satisfy the constraints
        for i in dec:
            # Skip useless variables
            if i in useless:
                continue
#            # Count only 4D variables
#            if base[i].ndim == 4:

            if dec[i].ndim == 4:
                species_num += 1
    else:
        # number of variables, ignore useless variables
        species_num = len(dec) - 8
    time, lev, lat, lon = dec[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape

    # perturbation size

    ct = 0  # counter
    specie = []  # list of species names
    for i in dec:
        if i in useless:
            continue

        elif dec[i].ndim == 4:
            if mode == 'semi':
                # array[ct, :, :, :, :] = (pert[i][:,:,:,:] - base[i][:,:,:,:]) /((hc-1)*h)
                array[ct, :, :, :, :] = (inc[i][:,:,:,:] - (2*base[i][:,:,:,:]) + dec[i][:,:,:,:]) / (hc*hc)
            # elif mode == 'normal additive':
#                 array[ct, :, :, :, :] = ((pert[i][:,:,:,:] / h) - (base[i][:,:,:,:] /h)) / h1
            # add variable name to list:
            specie.append(i)
            ct += 1

    return array, shape, specie


def hyd_second_order(base, pert_specie, ptype, mode, useless, hyd_pert=1):
    # Index to the species of interest
    # bspecies = base['PM25dx1x2']

    species_num = 0
    if ptype == 'H':
        # First, count the number of variables that satisfy the constraints
        for i in base:
            # Skip useless variables
            if i in useless:
                continue

#            # Count only 4D variables

            if base[i].ndim == 4:
                species_num += 1
    else:
        # number of variables, ignore useless variables
        species_num = len(base) - 8
    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
#    time, lat, lon = base[pert_specie].shape
#    # array holding all variables
#    array = np.zeros([species_num, time, lat, lon])
    shape = array.shape

    ct = 0  # counter
    specie = []  # list of species names
    for i in base:
        if i in useless:
            continue
        elif base[i].ndim == 4:
            if mode == 'semi':
                array[ct, :, :, :, :] = base[i][:, :, :, :] / (hyd_pert**2.0)
            elif mode == 'normal additive':
                array[ct, :, :, :, :] = base[i][:, :, :, :] / (hyd_pert**2.0)
            # add variable name to list:
            specie.append(i)
            ct += 1

    return array, shape, specie


def hybrid_cross_sensitivity(base, pert, pert_specie, ptype, mode, useless, hc, h=1):
    # Index to the species of interest

    if ptype == 'H':
        # number of variables, ignore useless variables
        species_num = len(base) - 4
    else:
        # number of variables, ignore useless variables
        species_num = len(base) - 8
    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape

    # perturbation size

    ct = 0  # counter
    specie = []  # list of species names
    for i in base:
        if i in useless:
            continue
        else:
            array[ct, :, :, :, :] = (
                (pert[i][:, :, :, :] / h) - (base[i][:, :, :, :] / h)) / ((hc - 1)*h)
            # add variable name to list:
            specie.append(i)

        ct += 1

    return array, shape, specie


def central_hybrid_cross_sensitivity(base, pert, pert_specie, ptype, mode, useless, hc, h=1):
    # where pert represents the diagnostic with both real and dual perturbation, while base is the
    # only dual perturbation

    if ptype == 'H':
        # number of variables, ignore useless variables
        species_num = len(base) - 4
        # Use number of 'total' variables instead
        species_num = 32
    else:
        # number of variables, ignore useless variables
        species_num = len(base) - 8
    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape

    # perturbation size

    ct = 0  # counter
    specie = []  # list of species names
    for i in base:
        if i in useless:
            continue

        else:
            if mode == 'semi':
                # array[ct, :, :, :, :] = (pert[i][:,:,:,:] - base[i][:,:,:,:]) /((hc-1)*h)
                array[ct, :, :, :, :] = (
                    pert[i][:, :, :, :] / (2*hc*h)) - (base[i][:, :, :, :] / (2*h*hc))
            # elif mode == 'normal additive':
#                 array[ct, :, :, :, :] = ((pert[i][:,:,:,:] / h) - (base[i][:,:,:,:] /h)) / h1
            # add variable name to list:
            specie.append(i)
        ct += 1

    return array, shape, specie


def hyd_cross_sensitivity(base, pert_specie, ptype, useless, hyd_pert=1):
    # Index to the species of interest
    # bspecies = base['PM25dx1x2']

    if ptype == 'H':
        # number of variables, ignore useless variables
        species_num = len(base) - 4
        # Use number of 'total' variables instead
        species_num = 32
    else:
        # number of variables, ignore useless variables
        species_num = len(base) - 8
    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape

    ct = 0  # counter
    specie = []  # list of species names
    for i in base:
        if i in useless:
            continue
        else:
            array[ct, :, :, :, :] = base[i][:, :, :, :] / (hyd_pert**2.0)
            # add variable name to list:
            specie.append(i)
        ct += 1

    return array, shape, specie


def analytical_second_order(base, pert_specie, useless):
    # number of variables, ignore useless variables
    species_num = len(base) - 8
    time, lev, lat, lon = base[pert_specie].shape

    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape

    ct = 0  # counter
    specie = []  # list of species names
    for i in base:
        if i in useless:
            continue

        else:
            array[ct, :, :, :, :] = base[i][:, :, :, :] * 6.0

            # Then convert to ppb:
            array[ct, :, :, :, :] = array[ct, :, :, :, :] * 1e9
            # semi-normalize:
            array[ct, :, :, :, :] = array[ct, :, :, :, :] * \
                (base[i][:, :, :, :] ** 2.0)
            specie.append(i)
        ct += 1

    return array, shape, specie


# def one_to_one(fd, hyd, mean, meanb, meanf, utype, fname, pert, t, sensitivity, x, y, z='', base='conc'):  #x
# x and y are strings for bottom and up, t is simulation time
def one_to_one1(ax, fd, hyd, utype, fname, pert, t, sensitivity, x, y, z='', base='conc', color='blue', lab=None):

    # Flatten the arrays:
    fd = fd.flatten()
    hyd = hyd.flatten()
#    if data is not None:
#        data = data.flatten()

    # t = '20min'


    rangemin = min(np.nanmin(fd), np.nanmin(hyd))
    rangemax = max(np.nanmax(fd), np.nanmax(hyd))
    v = [rangemin, rangemax, rangemin, rangemax]
    print(f'The maximum is {rangemax:1.6f}.')
    print(f'The minimum is {rangemin:1.6f}.')

    ax.axis(v)

    # to plot the 1 to 1 line
    ax.plot([rangemin, rangemax], [rangemin, rangemax], 'k-')

    # retrieving the axis to ax
    #ax = plt.gca()
    
    #Scatterplot with specified colors and markers
    ax.scatter(fd, hyd, c=color, marker='o', alpha=0.6)
    

    # setting axis conditions
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

#    # Create the scatterplot with the colors based on the 'data' array
#    #cplot = ax.scatter(fd, hyd, c=colors, cmap=cm.batlow, norm=norm)
#    cplot = ax.scatter(fd, hyd, c=colors.reshape(-1,4))
    # Scatter plot with specified color and markers
#    ax.scatter(fd[~mask_c], hyd[~mask_c], c='blue', label='Ship emission', marker='o', alpha=0.6)
#    ax.scatter(fd[mask_c], hyd[mask_c], c='white', label='No ship emission', marker='o', alpha=0.6)

    # ax.legend(shadow=True, loc="upper left")
    # Mask of nan numbers and infinite numbers
    # Using isfinte to mask for the two scenarios
    mask = ~np.isnan(fd) & ~np.isnan(hyd)
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        fd[mask], hyd[mask])
    # slope, intercept, r_value, p_value, std_err = stats.linregress(fd, hyd)

    # defining axis and notations
    ax.minorticks_on()
    ax.yaxis.grid(True, which='major', color='black',
                  linestyle='--', linewidth='0.4')
    ax.yaxis.grid(True, which='minor', color='gray',
                  linestyle=':', linewidth='0.3',)
    ax.xaxis.grid(True, which='major', color='black',
                  linestyle='--', linewidth='0.4')
    ax.xaxis.grid(True, which='minor', color='gray',
                  linestyle=':', linewidth='0.3')

    if base == 'conc':
        if sensitivity == 'first':
            if pert == '/allcells':
                # ax.set_xlabel(r'finite difference $\frac{{\partial ' + y.replace("_", "\_") + '}_{(x,y,z,t=1)}}{{\partial ' + x.replace("_", "\_") + '}_{(x,y,z,t=1)}}$') #r treats as raw strings, escape code ignored
                # ax.set_ylabel(r'Hyd $\frac{{\partial ' + y.replace("_", "\_") + '}_{(x,y,z,t=1)}}{{\partial ' + x.replace("_", "\_") + '}_{(x,y,z,t=1)}}$')
                xlabel = (r'$\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    # r treats as raw strings, escape code ignored
                    "SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$')
                ylabel = (r'$\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$')
            else:  # ground layer perturbation
                xlabel = (r'$\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    # r treats as raw strings, escape code ignored
                    "SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$')
                ylabel = (r'$\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$')
        elif sensitivity == 'second':
            if pert == '/allcells':
                if utype == 'A':
                    xlabel = (r'$\frac{{{\partial}^2 [' + y.replace("dx2", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                        "SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}  {[' + x.replace("SpeciesConc_", "") + ']}^2$')
                    ylabel = (r'$\frac{{{\partial}^2 [' + y.replace("dx2", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                        "SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$')
                else:
                    # ax.set_xlabel(r'Hybrid finite difference $\frac{{\partial}^2 ' + y.replace("_", "\_").replace("dx2", "") + '}{{\partial ' + x.replace("_", "\_") + '}^2}$')
                    #         ax.set_ylabel(r'Hyd $\frac{{\partial}^2 ' + y.replace("_", "\_").replace("dx2", "") + '}{{\partial ' + x.replace("_", "\_") + '}^2}$')
                    xlabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                        "SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}  {[' + x.replace("SpeciesConc_", "") + ']}^2$')
                    ylabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                        "SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$')
            else:
                xlabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}^2}_{(x,y,z=1,t=0hr)}}  {[' + x.replace("SpeciesConc_", "") + ']}^2$')
                ylabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}^2}_{(x,y,z=1,t=0hr)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$')
        elif sensitivity == 'cross':
            if pert == '/allcells':
                # ax.set_xlabel(r'finite difference Sensitivities $\frac{{\partial}^2 ' + y.replace("_", "\_") + '}{{\partial ' + x.replace("_", "\_") + '}\partial ' + z + '}$')

                xlabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}{\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
                ylabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}{\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
            else:
                xlabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}{\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
                ylabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}{\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
            # ax.set_ylabel(r'Hyd $\frac{{\partial}^2 ' + y.replace("SpeciesConcdx2_", "") + '}{{\partial ' + x.replace("SpeciesConc_", "") + '}\partial ' + z.replace("SpeciesConc_", "") + '}$')
    elif base == 'emis':
        if sensitivity == 'first':
            xlabel = (r'$\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                          r'{{\partial (E_{' + x.replace("SpeciesConc_",
                                                         "") + '})}_{(x,y,z,t=0)}} '
                          r'(E_{' + x.replace("SpeciesConc_", "") + '})$')
            ylabel = (r'$\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                          r'{{\partial (E_{' + x.replace("SpeciesConc_",
                                                         "") + '})}_{(x,y,z,t=0)}} '
                          r'(E_{' + x.replace("SpeciesConc_", "") + '})$')
        elif sensitivity == 'second':
            xlabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                          r'{{{\partial (E_{' + x.replace("SpeciesConc_",
                                                          "") + '})}^2}_{(x,y,z,t=0)}}   '
                          r'{(E_{' + x.replace("SpeciesConc_", "") + '})}^2$')
            ylabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                          r'{{{\partial (E_{' + x.replace("SpeciesConc_",
                                                          "") + '})}^2}_{(x,y,z,t=0)}}   '
                          r'{(E_{' + x.replace("SpeciesConc_", "") + '})}^2$')

    else:
        if sensitivity == 'first':
            xlabel = (r'$\frac{{\partial (E_{' + y.replace("Emis", "").replace("_","\_") + '})}_{(x,y,z,t=' + t + ')}}'
                          r'{{\partial [' +
                          x.replace("SpeciesConc_", "") +
                          ']}_{(x,y,z,t=0)}} '
                          r'[' + x.replace("SpeciesConc_", "") + ']$')
            ylabel = (r'$\frac{{\partial (E_{' + y.replace("Emis", "").replace("_","\_") + '})}_{(x,y,z,t=' + t + ')}}'
                          r'{{\partial [' +
                          x.replace("SpeciesConc_", "") +
                          ']}_{(x,y,z,t=0)}} '
                          r'[' + x.replace("SpeciesConc_", "") + ']$')
        elif sensitivity == 'second':
            xlabel = (r'$\frac{{{\partial}^2 (E_{' + y.replace("Emis", "").replace("_","\_") + '})}_{(x,y,z,t=' + t + ')}}'
                          r'{{{\partial [' + x.replace("SpeciesConc_", "") +
                          ']}^2}_{(x,y,z,t=0)}}   '
                          r'{[' + x.replace("SpeciesConc_", "") + ']}^2$')
            ylabel = (r'$\frac{{{\partial}^2 (E_{' + y.replace("Emis", "").replace("_","\_") + '})}_{(x,y,z,t=' + t + ')}}'
                          r'{{{\partial [' + x.replace("SpeciesConc_", "") +
                          ']}^2}_{(x,y,z,t=0)}}   '
                          r'{[' + x.replace("SpeciesConc_", "") + ']}^2$')

    # ax.set_title(r'Graph of $\frac{\partial b}{\partial x}$')
    # ax.set_title(r'Graph of ${CO}^3$')
    #ax.legend()
    ax.set_title(ylabel)
    # Add legend at the bottom right
    #ax.legend(loc='lower right')
    # ax.set_title(r'Hyd sensitivities against finite sensitivities')
    # previous: xy=(0.05, 0.7), fontsize=14
#    ax.annotate('R$^2$ = {:0.3F}\nSlope = {:0.2F}\nIntercept = {:0.2E}\nBlack Line = 1:1\nIgnored Calls = {:0.3F}%\nForward calls = {:0.3F}%\nBackward calls = {:0.3F}%'ignored_calls, forward_calls, backward_calls)
    ax.annotate('R$^2$ = {:0.3F}\nSlope = {:0.2F}\nIntercept = {:0.2E}\nBlack Line = 1:1'.format(r_value**2., slope, intercept), xy=(0.05, 0.7), fontsize=14, fontweight='bold', xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='bottom')
    #add another annotation for labelling
    ax.annotate(lab, xy=(0.90, 0.05), fontsize=16, fontweight='bold', xycoords='axes fraction')
#    plt.tight_layout()
#    pdfFile.savefig()
#    pdfFile.close()
#    plt.close()


def drawmap(shape, dataset, mark, ty, y, x, fname, utype):
    # Convert data back to xarray

    # units = r'$\frac{PM_{2.5} \mu g / m^3}{NH4 \mu g / m^3}$'
    base = xr.DataArray(dataset, coords=shape.coords,
                        dims=shape.dims, attrs=shape.attrs)
    if mark == "hyd":
        if ty == 'first':
            if utype == 'A':
                units = r'$\frac{(' + y.replace("_", "\_") + ' \mu g / m^3)}{(' + \
                    x.replace("_", "\_") + ' \mu g / m^3)}$'
            elif utype == 'S':
                units = 'ppb'
            # testing
            # units = r'$\frac{( {y} \mu g / m^3)}{( {x} \mu g / m^3)}$'
            # funits = units.format(y=y, x=x)
        # gcplot.single_panel(base, title=r'HYD Sensitivities $\frac{\partial PM2.5}{\partial NH4}$', pdfname="./hydPM25", unit=units)
            gcplot.single_panel(base, title=r'HYD Sensitivities $\frac{\partial ' + y.replace(
                "_", "\_") + '}{\partial ' + x.replace("_", "\_") + '}$', pdfname="./map/" + fname, unit=units, comap=cm.batlow)
            plt.close()
        elif ty == 'second':
            units = 'ppb'
            # units = r'${ppb}^2$'
            gcplot.single_panel(base, title=r'Hyd sensitivities $\frac{{{\partial}^2 ' + y.replace("SpeciesConcdx2_", "") + '}_{(x,y,z,t=20min)}}{{{\partial ' + x.replace(
                "SpeciesConc_", "") + '}^2}_{(x,y,z,t=0min)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$', pdfname="./map/" + fname, unit=units, comap=cm.batlow)
            plt.close()
    else:
        if ty == 'first':
            if utype == 'A':
                units = r'$\frac{(' + y.replace("_", "\_") + ' \mu g / m^3)}{(' + \
                    x.replace("_", "\_") + ' \mu g / m^3)}$'
            elif utype == 'S':
                units = 'ppb'
            gcplot.single_panel(base, title=r'Finite Difference Sensitivities $\frac{\partial ' + y.replace(
                "_", "\_") + '}{\partial ' + x.replace("_", "\_") + '}$', pdfname="./map/" + fname, unit=units, comap=cm.batlow)
            plt.close()
        elif ty == 'second':
            units = 'ppb'
            # units = r'${ppb}^2$'
            gcplot.single_panel(base, title=r'Hybrid Sensitivities $\frac{{\partial}^2 ' + y.replace(
                "_", "\_") + '}{{\partial ' + x.replace("_", "\_") + '}^2}$', pdfname="./map/" + fname, unit=units)
            plt.close()
        # gcplot.single_panel(base, title=r'Finite difference Sensitivities $\frac{\partial PM2.5}{\partial NH4}$', pdfname="./finitePM25", unit=units)
        # gcplot.single_panel(base, title=r'Finite difference Sensitivities $\frac{{\partial}^2 PM_{2.5}}{{\partial siagrowth}\partial NH_4}$', pdfname="./hybridPM25", unit=units)




def main():
    # make the one-to-one plot: for Aerosolmass
    colors = ['red', 'blue', 'green', 'purple']  # Different colors for each panel
    labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'] #label list for 8-panel plot
    
    # First-order sensitivity plots
#    fig1, axs1 = plt.subplots(1, 3, figsize=(18, 6))
#    axs1 = axs1.flatten()
#    fig1, axs1 = plt.subplots(2, 4, figsize=(24, 12))
#    axs1 = axs1.flatten()
    
    #Run a separate four-panel plot for the second-order plot:
    fig1, axs1 = plt.subplots(2, 2, figsize=(12, 12))
    axs1 = axs1.flatten()

    # Second-order sensitivity plots
#    fig2, axs2 = plt.subplots(2, 4, figsize=(24, 12))
#    axs2 = axs2.flatten()

    if ptype == 'S':
        for idx, (spec, i) in enumerate(plot_species.items()):
            print(i)
            print(specief[i])
            # First-order plot
            fd = fd_sensf[i, 0, :, :, :]
            hyd = hyd_sensf[i, 0, :, :, :]
            one_to_one1(axs1[idx], fd, hyd, 'S', 'one-to-one_' + duration + '_first_' + mechanism + '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + specie1, y=specief[i], base='conc', color=colors[idx], lab=labels[idx])
            
            # Second-order plot
#            fd = fd_senss[i, 0, :, :, :]
#            hyd = hyd_senss[i, 0, :, :, :]
#            one_to_one1(axs2[idx], fd, hyd, 'S', 'one-to-one_' + duration + '_second_' + mechanism + '_' + specief[i], pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie1, y=speciess[i], base='conc' color=colors[idx])
            
    elif ptype == 'H':
#        for idx, (spec, i) in enumerate(hemco_species.items()):
#            print(i)
#            print(specief[i])
#            # First-order plot
#            fd = fd_sensf[i, 0, :, :, :]
#            hyd = hyd_sensf[i, 0, :, :, :]
#            one_to_one1(axs1[idx], fd, hyd, 'S', 'one-to-one_' + duration + '_first_' + mechanism + '_' + species, pert_layer, duration, sensitivity='first', x='SpeciesConc_' + specie1, y=specief[i], base='hemco', color=colors[idx])
        #Using a loop instead for the for panel plot for first order hemco
#        idx = 0 #Initialie the outside loop to keep track of it across loops
#        for spec, i in hemco_species.items(): #first loop, first subplot row
#            print(i)
#            print(specief[i])
#            # First-order plot
#            fd = fd_sensf[i, 0, :, :, :]
#            hyd = hyd_sensf[i, 0, :, :, :]
#            one_to_one1(axs1[idx], fd, hyd, 'H', 'one-to-one_' + duration + '_first_' + mechanism + '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + specie1, y=specief[i], base='hemco', color=colors[idx], lab=labels[idx])
#            #Increment idx to keep track of the subplot number across loops
#            idx += 1
#            
#        for id, (spec, i) in enumerate(hemco_species.items()): #Second loop, second subplot row
#            print(i)
#            print(specief[i])
#            # First-order plot
#            fd = fd_sensf2[i, 0, :, :, :]
#            hyd = hyd_sensf2[i, 0, :, :, :]
#            one_to_one1(axs1[idx], fd, hyd, 'H', 'one-to-one_' + duration + '_first_' + mechanism + '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + specie2, y=specief[i], base='hemco', color=colors[id], lab=labels[idx])
#            #Increment idx to keep track of the subplot number across loops
#            idx += 1
            # Second-order plot
#            fd = fd_senss[i, 0, :, :, :]
#            hyd = hyd_senss[i, 0, :, :, :]
#            one_to_one1(axs2[idx], fd, hyd, mean, meanb, meanf, 'S', 'one-to-one_' + duration + '_second_' + mechanism + '_' + species, pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie1, y=speciess[i], color=colors[idx])
        #Using the loop for second order hemco
        idx = 0 #Initialie the outside loop to keep track of it across loops
        for spec, i in hemco_species.items(): #first loop, first subplot row
            print(i)
            print(specief[i])
            # Second-order plot
            fd = fd_senss[i, 0, :, :, :]
            hyd = hyd_senss[i, 0, :, :, :]
            one_to_one1(axs1[idx], fd, hyd, 'H', 'one-to-one_' + duration + '_second_' + mechanism + '_' + specief[i], pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie1, y=speciess[i], base='hemco', color=colors[idx], lab=labels[idx])
            #Increment idx to keep track of the subplot number across loops
            idx += 1
            
            #Exit the loop if idx becomes 3
            if idx == 3:
                break
            
        for id, (spec, i) in enumerate(hemco_species.items()): #Second loop, second subplot row
            print(i)
            print(specief[i])
            # Second-order plot
            #Only do O3 to O3:
            if spec == 'EmisO3_Total':
                fd = fd_senss2[i, 0, :, :, :]
                hyd = hyd_senss2[i, 0, :, :, :]
                one_to_one1(axs1[idx], fd, hyd, 'H', 'one-to-one_' + duration + '_second_' + mechanism + '_' + specief[i], pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie2, y=speciess[i], base='hemco', color=colors[idx], lab=labels[idx])
                #Increment idx to keep track of the subplot number across loops
                idx += 1


    # Second-order plot labels
#    fig2.text(0.5, 0.04, 'Hybrid [kg/m2/s]', ha='center', fontsize=14)
#    fig2.text(0.04, 0.5, 'Hyd [kg/m2/s]', va='center', rotation='vertical', fontsize=14)
    
    fig1.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, hspace=0.3)
    fig1.text(0.5, 0.01, 'Hybrid [kg/m2/s]', ha='center', fontsize=18)  # Lowered the bottom text
    fig1.text(0.01, 0.5, 'Hyd [kg/m2/s]', va='center', rotation='vertical', fontsize=18)  # Shifted the left text
    #fig1.subplots_adjust(left=0.20, bottom=0.20, right=0.95, top=0.95, hspace=0.4, wspace=0.4)
    #fig1.subplots_adjust(left=0.1, bottom=0.12, right=0.95, top=0.95, hspace=0.3) #for hemco first order

	# Adjust margins and spacing between subplots
    # fig1.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15, wspace=0.4, hspace=0.4)
    # fig2.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15, wspace=0.4, hspace=0.4)
    #plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
    
    # Save both figures to PDF
    #pdfFile1 = PdfPages('four_panel_plot_first_order_NOx_shipemis_1hr_third.pdf')
    #pdfFile1 = PdfPages('8_panel_plot_first_order_shipemis.pdf')
    pdfFile1 = PdfPages('4_panel_plot_second_order_shipemis.pdf')
    #plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
    #plt.tight_layout()
    pdfFile1.savefig(fig1)
    pdfFile1.close()
    
    #pdfFile2 = PdfPages('four_panel_plot_second_order_NOx_shipemis_1hr_third.pdf')
#    pdfFile2 = PdfPages('three_panel_plot_second_order_NOx_emis_1hr_final20%.pdf')
#    # plt.tight_layout(rect=[0.05, 0.05, 1, 1])
#    #plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
#    plt.tight_layout()
#    pdfFile2.savefig(fig2)
#    pdfFile2.close()
    
    
    plt.close(fig1)
#    plt.close(fig2)
    
    # plot_var = {'O3': 127} #For NOx
#     for i in plot_var.values():
#         #plot first order
#         drawmap(shaped, hyd_sensf[i, 0, 0, :, :], "hyd", "first", y=specief[i], x=specie1, fname=specief[i] +'_NOxemis_1hr.pdf', utype='S', t='1hr')
#         #plot second order
#         drawmap(shaped, hyd_senss[i, 0, 0, :, :], "hyd", "second", y=speciess[i], x=specie1, fname=specie_hyds[i] + '_NOxemis_1hr.pdf', utype='S', t='1hr')
    
    #Do the mapped plots
    #plot_var = {'ISOP': 173, 'O3': 127} For ISOP
    
    #shaped = ref_realf['SpeciesConc_CO'].isel(time=0, lev=0)

#Get the reversed colormap:
#rev_WhGrYlRd = WhGrYlRd.reversed()
white_cmap = ListedColormap(['white']) #white colormap

colors = [(0, 0, 0.5), (1, 1, 1), (0.5, 0, 0)]  # Deep red, white, deep blue
positions = [0, 0.5, 1]  # Positions of the colors

# Create the colormap
custom_cmap = LinearSegmentedColormap.from_list("deep_red_white_deep_blue", list(zip(positions, colors)))

    


# Running all plots at once:
# Write a loop to handle S or A:
if ptype == "S":
    hyd_sensf, shape, specie_hydf = hyd_first_order(
        hyd_dx2f, 'SpeciesConcdx2_' + species, ptype, 'S', useless, hyd_pert=1)
    fd_sensf, shape, specief = central_first_order(
        ref_realf, ref_pertf, h1, 'SpeciesConc_' + species, ptype, useless, 'semi')
# hyd_sensf, shape, specie_hydf = hyd_first_order(hyd_dx2f, 'SpeciesConcdx2_' + specie1, 'S', useless, hyd_pert=1)
# fd_sensf, shape, specief = central_first_order(ref_realf, ref_pertf, h1, 'SpeciesConc_' + specie1, useless, 'semi')

    hyd_senss, shape, specie_hyds = hyd_second_order(
        hyd_dx1x2s, 'SpeciesConcdx1x2_' + species, ptype, 'semi', useless, hyd_pert=1)
    fd_senss, shape, speciess = central_hybrid_second_order(
        ref_reals, ref_perts, 'SpeciesConcdx2_' + species, ptype, 'semi', useless, hc, h=1)

# Aerosol mass case:
elif ptype == 'A':

    # First order
    hyd_sensf, shape, specie_hydf = hyd_first_order(
        hyd_dx2f, 'AerMassNH4dx2', ptype, 'A', useless, hyd_pert=1)
    fd_sensf, shape, specief = central_first_order(
        ref_realf, ref_pertf, h1, 'AerMassNH4', ptype, useless, 'semi')
# Second order
    hyd_senss, shape, specie_hyds = hyd_second_order(
        hyd_dx1x2s, 'AerMassNH4dx1x2', ptype, 'semi', useless, hyd_pert=1)
    fd_senss, shape, speciess = central_hybrid_second_order(
        ref_reals, ref_perts, 'AerMassNH4dx2', ptype, 'semi', useless, hc, h=1)

elif ptype == 'H':

#    hyd_sensf, shape, specie_hydf = hyd_first_order(
#        hyd_dx2f, 'Emis' + species + '_Total', ptype, 'semi', useless, hyd_pert=1)
#    fd_sensf, shape, specief, species_num = central_first_order(
#        ref_realf, ref_pertf, h1, 'Emis' + species + '_Total', ptype, useless, 'semi')
#
#    hyd_senss, shape, specie_hyds = hyd_second_order(
#        hyd_dx1x2s, 'Emis' + species + '_Total', ptype, 'semi', useless, hyd_pert=1)
#    fd_senss, shape, speciess = central_hybrid_second_order(
#        ref_reals, ref_perts, 'Emis' + species + '_Total', ptype, 'semi', useless, hc, h=1)
        
    hyd_sensf, shape, specie_hydf = hyd_first_order(
        hyd_dx2f, 'Emis' + species + '_Total', ptype, 'semi', useless, hyd_pert=1)
    fd_sensf, shape, specief = central_first_order(
        ref_realf, ref_pertf, h1, 'Emis' + species + '_Total', ptype, useless, 'semi')

    hyd_senss, shape, specie_hyds = hyd_second_order(
        hyd_dx1x2s, 'Emis' + species + '_Total', ptype, 'semi', useless, hyd_pert=1)
    fd_senss, shape, speciess = central_hybrid_second_order(
        ref_reals, ref_perts, 'Emis' + species + '_Total', ptype, 'semi', useless, hc, h=1)
        
#    fd_senssf, shape, speciess = central_finite_second_order(
#        ref_reals, ref_reald, ref_perts, 'Emis' + species + '_Total', ptype, 'semi', useless, hc, h=1)
   
   #Second forward order case:
#    hyd_sensf2, shape, specie_hydf2 = hyd_first_order(
#        hyd_dx2f2, 'Emis' + species + '_Total', ptype, 'semi', useless, hyd_pert=1)
#    fd_sensf2, shape, specief2 = central_first_order(
#        ref_realf2, ref_pertf2, h1, 'Emis' + species + '_Total', ptype, useless, 'semi')
        
    #Second second order case:
    hyd_senss2, shape, specie_hyds2 = hyd_second_order(
        hyd_dx1x2s2, 'Emis' + species + '_Total', ptype, 'semi', useless, hyd_pert=1)
    fd_senss2, shape, speciess2 = central_hybrid_second_order(
        ref_reals2, ref_perts2, 'Emis' + species + '_Total', ptype, 'semi', useless, h2, h=1)






main()
