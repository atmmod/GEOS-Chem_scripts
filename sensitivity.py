#!/usr/bin/env python3

# Complete Automated testing framework script

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from matplotlib.ticker import ScalarFormatter
from scipy import stats
import warnings
# import gcpy.plot as gcplot
from cmcrameri import cm
# from gcpy.units import check_units, data_unit_is_mol_per_mol
import os
from joblib import Parallel, delayed
from matplotlib import colors
from matplotlib.colors import BoundaryNorm

# Ignore warnings
warnings.filterwarnings("ignore")

# Reading Files

root_path = '/glade/work/sakinjole/rundirs/gc_4x5_merra2_fullchem'
species = 'NO'
specie1 = 'O3'
specie2 = 'NO2'
duration = '20min'
mechanism = 'fullchem_hemco_0.01'
mechanism2 = 'fullchem_quad'
pert_layer = '/allcells'  # else ground layer
ptype = 'H'

h1 = 0.0001
hc = 0.0001

print(mechanism)
# A dictionary for species variables and indices to plot (This serves as the numerator)
# Only for the species concentration diagnostic
# plot_species = {'CO':251, 'CO2':86, 'SO4':94, 'SO2':95, 'NO2':131, 'NO':132, 'O3':127, 'NH4': 135, 'NH3': 136}
plot_species = {'CO': 251, 'SO4': 94, 'SO2': 95, 'NO2': 131, 'NO': 132, 'O3': 127,
                'NH4': 135, 'NH3': 136, 'NIT': 134, 'HNO3': 217, 'HCl': 224, 'ISOP': 173}

#hemco_species = {'EmisSO4_Total': 3, 'EmisNO2_Total': 10, 'EmisNO_Total': 11,
#                 'EmisNH3_Total': 12, 'EmisCO_Total': 20}
                 
hemco_species = {'EmisO3_Total': 10, 'EmisNO2_Total': 11, 'EmisNO_Total': 12, 'EmisHNO3_Total': 19}

# emis_species = ?

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

# A list of species indices to plot

# Extending by adding sum of absolute difference code for both SpeciesConc and AerosolMass


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

#            # Count only 4D variables
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

            # Count only 4D variables

            if base[i].ndim == 4:
                species_num += 1
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
def one_to_one(fd, hyd, utype, fname, pert, t, sensitivity, x, y, z='', base='conc',data=None):

    # Flatten the arrays:
    fd = fd.flatten()
    hyd = hyd.flatten()
    if data is not None:
        data = data.flatten()

    # t = '20min'
#    ignored_calls = np.mean(mean)
#    backward_calls = np.mean(meanb)
#    forward_calls = np.mean(meanf)

    # Save files to the appropriate concentration file folder. In this case saving only at the hyd path
    if sensitivity == 'first':
        # pdfFile = PdfPages('/glade/work/sakinjole/rundirs/gc_4x5_merra2_fullchem/1hr_plot/'+ fname +'.pdf')
        pdfFile = PdfPages(root_path + perthyd_pathf + fname + '.pdf')
    elif sensitivity == 'second':
        # pdfFile = PdfPages('/glade/work/sakinjole/rundirs/gc_4x5_merra2_fullchem/plot_second/'+ fname +'.pdf')
        pdfFile = PdfPages(root_path + perthyd_paths + fname + '.pdf')
    elif sensitivity == 'cross':
        # pdfFile = PdfPages('/glade/work/sakinjole/rundirs/gc_4x5_merra2_fullchem/plot_cross/'+ fname +'.pdf')
        pdfFile = PdfPages(root_path + perthyd_pathc + fname + '.pdf')
    rangemin = min(np.nanmin(fd), np.nanmin(hyd))
    rangemax = max(np.nanmax(fd), np.nanmax(hyd))
    v = [rangemin, rangemax, rangemin, rangemax]
    print(f'The maximum is {rangemax:1.2f}.')
    print(f'The minimum is {rangemin:1.2f}.')

    plt.axis(v)

    # to plot the 1 to 1 line
    plt.plot([rangemin, rangemax], [rangemin, rangemax], 'k-')

    # retrieving the axis to ax
    ax = plt.gca()

    # setting axis conditions
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    # Adopting Ed's code here for axis conditions
    # ax.set_xticks(np.linspace(rangemin, rangemax, 6))
#     ax.set_yticks(np.linspace(rangemin, rangemax, 6))
#
#     # Copied from https://stackoverflow.com/questions/42142144/displaying-first-decimal-digit-in-scientific-notation-in-matplotlib
#     class ScalarFormatterForceFormat(ScalarFormatter):
#         def _set_format(self):  # Override function that finds format to use.
#             self.format = "%1.2f"  # Give format here
#     yfmt = ScalarFormatterForceFormat()
#     yfmt.set_powerlimits((0,0))
#     ax.yaxis.set_major_formatter(yfmt)
#     ax.xaxis.set_major_formatter(yfmt)


    # Assuming your data range is [vmin, vmax]
    mask_c = data == 0 #Identify where data is zero
    #Initialize all to white (RGBA [1,1,1,1])
#    colors = np.ones((*data.shape, 4)) # Creating RGBA color array filled with white
#    
#    #Set non-zero values to blue
#    blue_color = [0,0,1,1] #You can adjust shade of blue
#    colors[~mask_c] = blue_color
#    
#
#    # Create the scatterplot with the colors based on the 'data' array
#    #cplot = ax.scatter(fd, hyd, c=colors, cmap=cm.batlow, norm=norm)
#    cplot = ax.scatter(fd, hyd, c=colors.reshape(-1,4))
    # Scatter plot with specified color and markers
#    ax.scatter(fd[~mask_c], hyd[~mask_c], c='blue', label='Ship emission', marker='o', alpha=0.6)
#    ax.scatter(fd[mask_c], hyd[mask_c], c='white', label='No ship emission', marker='o', alpha=0.6)

    #plot with colormap
    cplot = ax.scatter(fd, hyd, c=data, cmap='binary')
    if utype == 'S':
        cbar = plt.colorbar(cplot, ax=ax, label=r'ppb')
    elif utype == 'H':
        cbar = plt.colorbar(cplot, ax=ax, label=r'kg/m2/s')
    else:
        cbar = plt.colorbar(cplot, ax=ax, label=r'$\mu g / m^3$')
#    cplot = ax.scatter(fd, hyd, c=hyd, cmap=cm.batlow)
#    if utype == 'S':
#        cbar = plt.colorbar(cplot, ax=ax, label=r'ppb')
#    elif utype == 'H':
#        cbar = plt.colorbar(cplot, ax=ax, label=r'kg/m2/s')
#    else:
#        cbar = plt.colorbar(cplot, ax=ax, label=r'$\mu g / m^3$')
    

        # cbar = plt.colorbar(cplot, ax=ax, label=r'${ppb}^2$')
    # cbar = plt.colorbar(cplot, ax=ax, label=r'$\frac{{ppb}^2}{{ppb}^2}')
    # cbar = plt.colorbar(cplot, ax=ax)
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
                ax.set_xlabel(r'finite difference $\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    # r treats as raw strings, escape code ignored
                    "SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$')
                ax.set_ylabel(r'Hyd $\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$')
            else:  # ground layer perturbation
                ax.set_xlabel(r'finite difference $\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    # r treats as raw strings, escape code ignored
                    "SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$')
                ax.set_ylabel(r'Hyd $\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$')
        elif sensitivity == 'second':
            if pert == '/allcells':
                if utype == 'A':
                    ax.set_xlabel(r'Hybrid $\frac{{{\partial}^2 [' + y.replace("dx2", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                        "SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}  {[' + x.replace("SpeciesConc_", "") + ']}^2$')
                    ax.set_ylabel(r'Hyd $\frac{{{\partial}^2 [' + y.replace("dx2", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                        "SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$')
                else:
                    # ax.set_xlabel(r'Hybrid finite difference $\frac{{\partial}^2 ' + y.replace("_", "\_").replace("dx2", "") + '}{{\partial ' + x.replace("_", "\_") + '}^2}$')
                    #         ax.set_ylabel(r'Hyd $\frac{{\partial}^2 ' + y.replace("_", "\_").replace("dx2", "") + '}{{\partial ' + x.replace("_", "\_") + '}^2}$')
                    ax.set_xlabel(r'Hybrid $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                        "SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}  {[' + x.replace("SpeciesConc_", "") + ']}^2$')
                    ax.set_ylabel(r'Hyd $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                        "SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$')
            else:
                ax.set_xlabel(r'Hybrid $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}^2}_{(x,y,z=1,t=0hr)}}  {[' + x.replace("SpeciesConc_", "") + ']}^2$')
                ax.set_ylabel(r'Hyd $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}^2}_{(x,y,z=1,t=0hr)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$')
        elif sensitivity == 'cross':
            if pert == '/allcells':
                # ax.set_xlabel(r'finite difference Sensitivities $\frac{{\partial}^2 ' + y.replace("_", "\_") + '}{{\partial ' + x.replace("_", "\_") + '}\partial ' + z + '}$')

                ax.set_xlabel(r'Hybrid $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}{\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
                ax.set_ylabel(r'Hyd $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}{\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
            else:
                ax.set_xlabel(r'Hybrid $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}{\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
                ax.set_ylabel(r'Hyd $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace(
                    "SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}{\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
            # ax.set_ylabel(r'Hyd $\frac{{\partial}^2 ' + y.replace("SpeciesConcdx2_", "") + '}{{\partial ' + x.replace("SpeciesConc_", "") + '}\partial ' + z.replace("SpeciesConc_", "") + '}$')
    elif base == 'emis':
        if sensitivity == 'first':
            ax.set_xlabel(r'$\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                          r'{{\partial (E_{' + x.replace("SpeciesConc_",
                                                         "") + '})}_{(x,y,z,t=0hr)}} '
                          r'(E_{' + x.replace("SpeciesConc_", "") + '})$')
            ax.set_ylabel(r'$\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                          r'{{\partial (E_{' + x.replace("SpeciesConc_",
                                                         "") + '})}_{(x,y,z,t=0hr)}} '
                          r'(E_{' + x.replace("SpeciesConc_", "") + '})$')
        elif sensitivity == 'second':
            ax.set_xlabel(r'Hybrid $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                          r'{{{\partial (E_{' + x.replace("SpeciesConc_",
                                                          "") + '})}^2}_{(x,y,z,t=0hr)}}   '
                          r'{(E_{' + x.replace("SpeciesConc_", "") + '})}^2$')
            ax.set_ylabel(r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                          r'{{{\partial (E_{' + x.replace("SpeciesConc_",
                                                          "") + '})}^2}_{(x,y,z,t=0hr)}}   '
                          r'{(E_{' + x.replace("SpeciesConc_", "") + '})}^2$')

    else:
        if sensitivity == 'first':
            ax.set_xlabel(r'FD $\frac{{\partial (E_{' + y.replace("Emis", "").replace("_","\_") + '})}_{(x,y,z,t=' + t + ')}}'
                          r'{{\partial [' +
                          x.replace("SpeciesConc_", "") +
                          ']}_{(x,y,z,t=0hr)}} '
                          r'[' + x.replace("SpeciesConc_", "") + ']$')
            ax.set_ylabel(r'Hyd $\frac{{\partial (E_{' + y.replace("Emis", "").replace("_","\_") + '})}_{(x,y,z,t=' + t + ')}}'
                          r'{{\partial [' +
                          x.replace("SpeciesConc_", "") +
                          ']}_{(x,y,z,t=0hr)}} '
                          r'[' + x.replace("SpeciesConc_", "") + ']$')
        elif sensitivity == 'second':
            ax.set_xlabel(r'Hybrid $\frac{{{\partial}^2 (E_{' + y.replace("Emis", "").replace("_","\_") + '})}_{(x,y,z,t=' + t + ')}}'
                          r'{{{\partial [' + x.replace("SpeciesConc_", "") +
                          ']}^2}_{(x,y,z,t=0hr)}}   '
                          r'{[' + x.replace("SpeciesConc_", "") + ']}^2$')
            ax.set_ylabel(r'Hyd $\frac{{{\partial}^2 (E_{' + y.replace("Emis", "").replace("_","\_") + '})}_{(x,y,z,t=' + t + ')}}'
                          r'{{{\partial [' + x.replace("SpeciesConc_", "") +
                          ']}^2}_{(x,y,z,t=0hr)}}   '
                          r'{[' + x.replace("SpeciesConc_", "") + ']}^2$')

    # ax.set_title(r'Graph of $\frac{\partial b}{\partial x}$')
    # ax.set_title(r'Graph of ${CO}^3$')
    ax.legend()
    ax.set_title('Hemco simulation, 0.01% perturbation')
    # Add legend at the bottom right
    ax.legend(loc='lower right')
    # ax.set_title(r'Hyd sensitivities against finite sensitivities')
    # previous: xy=(0.05, 0.7), fontsize=14
#    ax.annotate('R$^2$ = {:0.3F}\nSlope = {:0.2F}\nIntercept = {:0.2E}\nBlack Line = 1:1\nIgnored Calls = {:0.3F}%\nForward calls = {:0.3F}%\nBackward calls = {:0.3F}%'ignored_calls, forward_calls, backward_calls)
    ax.annotate('R$^2$ = {:0.3F}\nSlope = {:0.2F}\nIntercept = {:0.2E}\nBlack Line = 1:1'.format(r_value**2., slope, intercept), xy=(0.05, 0.65), fontsize=12, fontweight='bold', xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='bottom')
    plt.tight_layout()
    pdfFile.savefig()
    pdfFile.close()
    plt.close()


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


def sum_abs(base, pert, pert_specie, useless):

    # number of variables, ignore useless variables
    species_num = len(base) - 8

    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    # shape = array.shape

    ct = 0  # counter
    # specie = [] #list of species names
    for i in base:
        if i in useless:
            continue

        else:
            array[ct, :, :, :, :] = base[i][:, :, :, :] - \
                pert[i][:, :, :, :]  # take the difference
            # add variable name to list:
            # specie.append(i)
        ct += 1

    return array


def absrel(base, pert, pert_specie, useless):

    # number of variables, ignore useless variables
    species_num = len(base) - 8

    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    # shape = array.shape

    ct = 0  # counter
    # specie = [] #list of species names
    for i in base:
        if i in useless:
            continue

        else:
            # take the relative difference
            array[ct, :, :, :, :] = (
                base[i][:, :, :, :] - pert[i][:, :, :, :]) / base[i][:, :, :, :]
            # add variable name to list:
            # specie.append(i)
        ct += 1

    return array


def sum_rel(base, pert, pert_specie, useless):

    # number of variables, ignore useless variables
    species_num = len(base) - 8

    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    # shape = array.shape

    ct = 0  # counter
    # specie = [] #list of species names
    for i in base:
        if i in useless:
            continue

        else:
            # Create a mask to avoid division by zero
            mask = base[i][:, :, :, :] != 0
            # base_mask = base[i][:,:,:,:][mask]
            # array[ct, :, :, :, :] =(pert[i][:,:,:,:] - base[i][:,:,:,:][mask]) /base[i][:,:,:,:][mask]
            # array[ct] = (pert[i][:,:,:,:] - denominator) / denominator  # relative difference
            # array[ct][~mask] = 0 # Set relative difference to zero where the denominator is zero
            array[ct, :, :, :, :] = np.where(
                # relative difference
                mask, (pert[i][:, :, :, :] - base[i][:, :, :, :]) / base[i][:, :, :, :], 0)
            # array[ct, :, :, :, :][:, ~mask] = 0
        ct += 1

        # array[ct, :, :, :, :] =(pert[i][:,:,:,:] - base[i][:,:,:,:]) /base[i][:,:,:,:]  #relative difference
#             #add variable name to list:
#             #specie.append(i)
#         ct += 1

    return array


def sum_reln(base, pert, pert_specie, useless):

    # number of variables, ignore useless variables
    species_num = len(base) - 8

    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    # shape = array.shape

    ct = 0  # counter
    # specie = [] #list of species names
    for i in base:
        if i in useless:
            continue

        else:
            # Create a mask to avoid division by zero
            base_mask = base[i][:, :, :, :] != 0
            base_masked = np.where(base_mask, np.nan, base[i][:, :, :, :])
            pert_masked = np.where(base_mask, np.nan, pert[i][:, :, :, :])
            # base_mask = base[i][:,:,:,:][mask]
            array[ct, :, :, :, :] = (pert_masked - base_masked) / base_masked
            # array[ct] = (pert[i][:,:,:,:] - denominator) / denominator  # relative difference
            # array[ct][~mask] = 0 # Set relative difference to zero where the denominator is zero
            # array[ct, :, :, :, :] = np.where(mask, (pert[i][:, :, :, :] - base[i][:,:,:,:]) / base[i][:,:,:,:], 0)  # relative difference
            # array[ct, :, :, :, :][:, ~mask] = 0
        ct += 1

        # array[ct, :, :, :, :] =(pert[i][:,:,:,:] - base[i][:,:,:,:]) /base[i][:,:,:,:]  #relative difference
#             #add variable name to list:
#             #specie.append(i)
#         ct += 1

    return array

# Check for NIT and NH4


def sum_check(base, pert, pert_specie, useless):

    # number of variables, ignore useless variables
    species_num = len(base) - 8

    time, lev, lat, lon = base[pert_specie].shape
    # array holding all variables
    array = np.zeros([species_num, time, lev, lat, lon])
    # shape = array.shape

    ct = 0  # counter
    # specie = [] #list of species names
    for i in base:
        if i in useless:
            continue

        elif i == 'SpeciesConc_NIT' or i == 'SpeciesConc_NH4':
            # elif i == 'SpeciesConc_CO':
            array[ct, :, :, :, :] = base[i][:, :, :, :] - \
                pert[i][:, :, :, :]  # take the difference
            # add variable name to list:
            # specie.append(i)
        else:
            continue
        ct += 1

    return array


def main():

    if ptype == 'S':
        for i in plot_species.values():

            # for i in range(19):
            #
            #
            print(i)
            print(specief[i])
        # print(f"Type of specief[i]: {type(specief[i])}")
        # print(f"Value of specief[i]: {specief[i]}")
        # first order
        # if specief[i] == 'SpeciesConc_HCl':
#              one_to_one(fd_senss[i, -1, :, :, :], hyd_senss[i, -1, :, :, :], 'one-to-one_' + duration + '_second_' + mechanism + '_' + speciess[i], pert_layer, duration, sensitivity='second', x='SpeciesConc_' + species, y=speciess[i])
#          else:
#              continue
            one_to_one(fd_sensf[i, 1, :, :, :], hyd_sensf[i, 1, :, :, :], 'S', 'one-to-one_' + duration + '_first_' + mechanism +
                       '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + specie1, y=specief[i], base='conc')
        # one_to_one(fd_sensf[i, -1, :, :, :], hyd_sensf[i, -1, :, :, :], isens[-1, :, :, :], 'S', 'one-to-one_' + duration + '_first_' + mechanism + '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + species, y=specief[i])
        # Second order
            one_to_one(fd_senss[i, 1, :, :, :], hyd_senss[i, 1, :, :, :], 'S', 'one-to-one_' + duration + '_second_' + mechanism +
                       '_' + speciess[i], pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie1, y=speciess[i], base='conc')
        # one_to_one(fd_sensc[i, -1, :, :, :], hyd_sensc[i, -1, :, :, :], 'one-to-one_' + duration + '_cross_' + mechanism + '_' + speciec[i], pert_layer, duration, sensitivity='cross', x='SpeciesConc_' + species, y=speciec[i], z='SpeciesConc_' + specie2)

        # one_to_one(fd_senss[i, -1, :, :, :], hyd_senss[i, -1, :, :, :], 'one-to-one_' + duration + '_second_' + mechanism + '_' + species[i], pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie1, y=species[i])

        # one_to_one(fd_sensc[i, -1, :, :, :], hyd_sensc[i, -1, :, :, :], 'one-to-one_' + duration + '_cross_' + mechanism + '_' + speciec[i], pert_layer, duration, sensitivity='cross', x='SpeciesConc_' + species, y=speciec[i], z='SpeciesConc_' + specie2)
        # else:
    elif ptype == 'A':
        for i in range(19):

            print(i)
            print(specief[i])

            one_to_one(fd_sensf[i, -1, :, :, :], hyd_sensf[i, -1, :, :, :], 'A', 'one-to-one_' + duration + '_first_' +
                       mechanism + '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + specie1, y=specief[i])
            one_to_one(fd_senss[i, -1, :, :, :], hyd_senss[i, -1, :, :, :], 'A', 'one-to-one_' + duration + '_second_' +
                       mechanism + '_' + speciess[i], pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie1, y=speciess[i])

    elif ptype == 'H':
        # for i in range(171): #Hemco dx1
        #for i in range(36):  # Hemco default
        #for i in hemco_species.values():
        for i in hemco_species.values():

            print(i)
            print(specief[i])

            one_to_one(fd_sensf[i, -1, :, :, :], hyd_sensf[i, -1, :, :, :], 'H', 'one-to-one_' + duration + '_first_' + mechanism +
                       '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + specie1, y=specief[i], base='hemco', data=ref_pertf['EmisHNO3_Total'][-1,:,:,:].values)
        # one_to_one(fd_sensf[i, -1, :, :, :], hyd_sensf[i, -1, :, :, :], isens[-1, :, :, :], 'S', 'one-to-one_' + duration + '_first_' + mechanism + '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + species, y=specief[i])
        # Second order
            one_to_one(fd_senss[i, -1, :, :, :], hyd_senss[i, -1, :, :, :], 'H', 'one-to-one_' + duration + '_second_' + mechanism +
                       '_' + speciess[i], pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie1, y=speciess[i], base='hemco', data=ref_pertf['EmisHNO3_Total'][-1,:,:,:].values)

#            one_to_one(fd_sensf[i, 0, :, :], hyd_sensf[i, 0, :, :], 'H', 'one-to-one_' + duration + '_first_' + mechanism +
#                       '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + specie1, y=specief[i], base='hemco')
#        # one_to_one(fd_sensf[i, -1, :, :, :], hyd_sensf[i, -1, :, :, :], isens[-1, :, :, :], 'S', 'one-to-one_' + duration + '_first_' + mechanism + '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + species, y=specief[i])
#        # Second order
#            one_to_one(fd_senss[i, 0, :, :], hyd_senss[i, 0, :, :], 'H', 'one-to-one_' + duration + '_second_' + mechanism +
#                       '_' + speciess[i], pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie1, y=speciess[i], base='hemco')


#             continue
    # hyd sensitivity map plot
    # for species_conc:
    # for i in range(307):
#         if i == 86 or i == 251: #CO and CO2
#              drawmap(shaped, hyd_senss[i, 1, 0, :, :], "hyd", "second", y=species[i], x='SpeciesConc_CO', fname=specie_hyds[i], utype='S')
        # drawmap(shaped, hyd_sens[i, 1, 0, :, :], "hyd", "second", y=specie[i], x='SpeciesConc_CO', fname=specie_hydf[i], utype='S')
#
#     #finite difference:
    # for i in range(307):
#           drawmap(shaped, fd_sens[i, 0, 0, :, :], "finite", "first", y=specie[i], x='SpeciesConc_CO', fname=specie[i], utype='S')
#
    # finite difference map plot
    # drawmap(shape, fd_sens, mark='')


# fd_sens, shape, specie = first_order_finite(ref_real, ref_pert, h, 'AerMassNH4', useless)
# hyd_sens, shape, specie_hyd = hyd_first_order(hyd_dx2, 'AerMassNH4dx2', 'A', useless, hyd_pert=1)
# hyd_sens, shape, specie_hyd = hyd_first_order(hyd_dx2, 'SpeciesConcdx2_CO', 'S', useless, hyd_pert=1)
# fd_sens, shape, specie = first_order_finite(ref_real, ref_pert, h, 'SpeciesConc_CO', useless, 'semi')



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
        


# hyd_sensc, shape, specie_hydc = hyd_cross_sensitivity(hyd_dx1x2c, 'SpeciesConcdx1x2_' + species, useless, hyd_pert=1)
# fd_sensc, shape, speciec = central_hybrid_cross_sensitivity(ref_realc, ref_pertc, 'SpeciesConcdx2_' + species, 'semi', useless, hc, h=1)

# Find the indicies of negative values
# negative_indices = np.where(hyd_sensf[251, -1, :, :, :] < 0)
# negative_indices = np.where(hyd_sensf[251, -1, :, :, :] < 0) #Second order sensitivity

# Calculate the absolute difference between hyd_sensf and fd_sensf
abs_diff = np.abs(hyd_sensf - fd_sensf)
# abs_diff = np.abs(hyd_senss - fd_senss)
diff = hyd_sensf - fd_sensf

# Find the indices where the absolute difference is greater than 1.0
# diff_indices = np.where(abs_diff[127, 0, :, :, :] > 10) #NH3 O3
# # diff_indices2 = np.where(diff[136, -1, :, :, :] < -0.7)
# # #diff_indices = np.where(abs_diff[224, -1, :, :, :] > 3.5e-12) #find the difference for second order
# #
# print("Indices with absolute difference greater than 10 at {} timestep:".format(duration))
# for dim, indices in enumerate(diff_indices):
#     print(f"Dimension {dim}: {indices}")
#
# print("Indices with absolute difference less than -0.7 at second transport timestep:")
# for dim, indices in enumerate(diff_indices2):
#     print(f"Dimension {dim}: {indices}")

# Loop through the second index of abs_diff
# found = False
# for i in range(abs_diff.shape[1]):
#    diff_indices = np.where(abs_diff[134, i, :, :, :] > 1)  # NH3 O3- 127, NIT - 134
#    if any(len(indices) > 0 for indices in diff_indices):
#        print(f"Indices with absolute difference greater than 1 NIT at {i} index timestep:")
#        for dim, indices in enumerate(diff_indices):
#            print(f"Dimension {dim}: {indices}")
#        found = True
#        break
#
# if not found:
#    print("No indices with absolute difference greater than 1 found in any index timestep.")
#
# second = False
# for i in range(abs_diff.shape[1]):
#    diff_indices2 = np.where(abs_diff[217, i, :, :, :] > 1)  # HNO3
#    if any(len(indices) > 0 for indices in diff_indices2):
#        print(f"Indices with absolute difference greater than 2.5 HNO3 at {i} index timestep:")
#        for dim, indices in enumerate(diff_indices2):
#            print(f"Dimension {dim}: {indices}")
#        second = True
#        break
#
# if not second:
#    print("No indices with absolute difference greater than 1 found in any index timestep.")


# #
# # Print the indices
# print("Indices with negative values:")
# print(negative_indices)
#
# print("Elements which are negative values")
#
# array_new = hyd_sensf[251, 2, :, :, :]
# print(array_new[negative_indices])
#
# print("Indices with negative values:")
# for dim, indices in enumerate(negative_indices):
#     print(f"Dimension {dim}: {indices}")

# Test to be sure you got the write array
# print("Testing the index 19, 71")
# print(hyd_sensf[251, 1, 26, 19, 71])

# print("finite sensitivity of CO at layer 0 and index 19, 71")
# print(fd_sensf[251, -1, 0, 19, 71])
#
# print("Hyd sensitivity of CO at layer 0 and index 19, 71")
# print(hyd_sensf[251, -1, 0, 19, 71])
#
# print("finite sensitivity of CO at layer 26 and index 19, 71")
# print(fd_sensf[251, -1, 26, 19, 71])
#
# print("Hyd sensitivity of CO at layer 26 and index 19, 71")
# print(hyd_sensf[251, -1, 26, 19, 71])
#

main()
