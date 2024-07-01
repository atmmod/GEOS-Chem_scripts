#!/usr/bin/env python3

# Complete Automated testing framework script for IGC

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
#from gcpy.units import check_units, data_unit_is_mol_per_mol
import os
from joblib import Parallel, delayed
from gcpy import WhGrYlRd
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap

# Ignore warnings
warnings.filterwarnings("ignore")

# Reading Files

root_path = '/glade/work/sakinjole/rundirs/gc_4x5_merra2_fullchem'
species = 'NO'
specie1 = 'NOx'
specie2 = 'NO'
duration = '24hr'
mechanism = 'fullchem_quad_1e3_20emispercent'
mechanism2 = 'fullchem_quad'
pert_layer = '/allcells'  # else ground layer
ptype = 'S'

print(mechanism)
#A dictionary for species variables and indices to plot (This serves as the numerator)
#Only for the species concentration diagnostics
#plot_species = {'CO':251, 'CO2':86, 'SO4':94, 'SO2':95, 'NO2':131, 'NO':132, 'O3':127, 'NH4': 135, 'NH3': 136}
#plot_species = {'CO':251, 'SO4':94, 'SO2':95, 'NO2':131, 'NO':132, 'O3':127, 'NH4': 135, 'NH3': 136, 'NIT': 134, 'HNO3': 217, 'HCl': 224, 'ISOP': 173}
plot_species = {'ISOP': 173, 'NO2':131, 'O3':127, 'NH4': 135}

#Avoid using gcpy for this:
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

#A list of species indices to plot

#Extending by adding sum of absolute difference code for both SpeciesConc and AerosolMass
def default_path(mechanism):
	#Folder paths default run
	real_path1 = root_path + '/default/' + mechanism2 + '/' + duration + '/GEOSChem.SpeciesConc.20190701_0000z.nc4'
	hyd_path1 = root_path + '/default_hyd/' + mechanism2 + '/' + duration + '/GEOSChem.SpeciesConc.20190701_0000z.nc4'
	real_path2 = root_path + '/default/' + mechanism2 + '/' + duration + '/GEOSChem.AerosolMass.20190701_0000z.nc4'
	hyd_path2 = root_path + '/default_hyd/' + mechanism2 + '/' + duration + '/GEOSChem.AerosolMass.20190701_0000z.nc4'
	#Read the arrays:
	real1 = xr.open_dataset(real_path1)
	hyd1 = xr.open_dataset(hyd_path1)
	real2 = xr.open_dataset(real_path2)
	hyd2 = xr.open_dataset(hyd_path2)
	return real1, hyd1, real2, hyd2
	
def main_path(ty, layer =''):     #where ty stands for the order of operation
	if ty == 'first':
		if layer == '/allcells':
			#Folder paths first order - updated to central difference
			real_path = '/first_output/' + mechanism + '/' + duration + layer + '/dec' + specie1 + '/'

			#pertreal_path = '/first_output/incNO/'
			pertreal_path = '/first_output/' + mechanism + '/' + duration + layer + '/inc' + specie1 + '/'
			#perthyd_path = '/second_output_hyd/perturbedNO_boxN/'
			perthyd_path = '/first_output_hyd/' + mechanism + '/' + duration + layer + '/hyd' + specie1 + '/'
		else:
			#Folder paths first order - updated to central difference
			real_path = '/first_output/' + mechanism + '/' + duration + '/ground/dec' + specie1 + '/'

			#pertreal_path = '/first_output/incNO/'
			pertreal_path = '/first_output/' + mechanism + '/' + duration + '/ground/inc' + specie1 + '/'
			#perthyd_path = '/second_output_hyd/perturbedNO_boxN/'
			perthyd_path = '/first_output_hyd/' + mechanism + '/' + duration + '/ground/hyd' + specie1 + '/'
	elif ty == 'second':
		if layer == '/allcells':
	
			#Folder paths second order
			real_path = '/second_output/' + mechanism + '/' + duration + layer + '/dec' + specie1 + '/'

			pertreal_path = '/second_output/' + mechanism + '/' + duration + layer + '/inc' + specie1 + '/'
			perthyd_path = '/second_output_hyd/' + mechanism + '/' + duration + layer + '/hyd' + specie1 + '/'
		else:
			#Folder paths second order
			real_path = '/second_output/' + mechanism + '/' + duration + '/ground/dec' + specie1 + '/'

			pertreal_path = '/second_output/' + mechanism + '/' + duration + '/gound/inc' + specie1 + '/'
			perthyd_path = '/second_output_hyd/' + mechanism + '/' + duration + '/ground/hyd' + specie1 + '/'

	elif ty == 'cross':
		#Folder paths hybrid
		if layer == '/allcells':
			real_path = '/cross_output/' + mechanism + '/' + duration + layer + '/dec' + specie1 + '_' + specie2 + '/'
		#hyd_path = '/hybrid_output_hyd/default'

			pertreal_path = '/cross_output/' + mechanism + '/' + duration + layer + '/inc' + specie1 + '_' + specie2 + '/'
			perthyd_path = '/cross_output_hyd/' + mechanism + '/' + duration + layer + '/hyd' + specie1 + '_' + specie2 + '/'
	return real_path, pertreal_path, perthyd_path


#File paths
def file(ty, pt):   #where ty is the diagnostic type, pt is the required path
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
 

def read_file(loc):
    '''Read the real files with xarray'''
    return xr.open_dataset(loc)

def get_file(array, order, real_path, pertreal_path, perthyd_path):  #gets all file paths for a particular diagnostic
	if array == 'A':
		if order == 'first':
			ref_real = read_file(root_path + real_path + file('A', 'real'))
			ref_pert = read_file(root_path + pertreal_path + file('A', 'real'))
			hyd_dx1 = read_file(root_path + perthyd_path + file('A', 'dx1'))
			hyd_dx2 = read_file(root_path + perthyd_path + file('A', 'dx2'))
			hyd_dx1x2 = read_file(root_path + perthyd_path + file('A', 'dx1x2'))
			ipercent = read_file(root_path + perthyd_path + file('A', 'isens'))
		elif order == 'second':
			ref_real = read_file(root_path + real_path + file('A', 'dx2'))
			ref_pert = read_file(root_path + pertreal_path + file('A', 'dx2'))
			hyd_dx1 = read_file(root_path + perthyd_path + file('A', 'dx1'))
			hyd_dx2 = read_file(root_path + perthyd_path + file('A', 'dx2'))
			hyd_dx1x2 = read_file(root_path + perthyd_path + file('A', 'dx1x2'))
			ipercent = read_file(root_path + perthyd_path + file('A', 'isens'))
			
	elif array == 'S':
		if order == 'first':
			ref_real = read_file(root_path + real_path + file('S', 'real'))
			ref_pert = read_file(root_path + pertreal_path + file('S', 'real'))
			hyd_dx1 = read_file(root_path + perthyd_path + file('S', 'dx1'))
			hyd_dx2 = read_file(root_path + perthyd_path + file('S', 'dx2'))
			hyd_dx1x2 = read_file(root_path + perthyd_path + file('S', 'dx1x2'))
			ipercent = read_file(root_path + perthyd_path + file('S', 'isens'))
			ipercentb = read_file(root_path + real_path + file('S', 'isens'))
			ipercentf = read_file(root_path + pertreal_path + file('S', 'isens')) 
		elif order == 'second':
			ref_real = read_file(root_path + real_path + file('S', 'dx2'))
			ref_pert = read_file(root_path + pertreal_path + file('S', 'dx2'))
			hyd_dx1 = read_file(root_path + perthyd_path + file('S', 'dx1'))
			hyd_dx2 = read_file(root_path + perthyd_path + file('S', 'dx2'))
			hyd_dx1x2 = read_file(root_path + perthyd_path + file('S', 'dx1x2'))
			ipercent = read_file(root_path + perthyd_path + file('S', 'isens'))
			ipercentb = read_file(root_path + real_path + file('S', 'isens'))
			ipercentf = read_file(root_path + pertreal_path + file('S', 'isens'))
		elif order == 'cross':
			ref_real = read_file(root_path + real_path + file('S', 'dx2'))
			ref_pert = read_file(root_path + pertreal_path + file('S', 'dx2'))
			hyd_dx1 = read_file(root_path + perthyd_path + file('S', 'dx1'))
			hyd_dx2 = read_file(root_path + perthyd_path + file('S', 'dx2'))
			hyd_dx1x2 = read_file(root_path + perthyd_path + file('S', 'dx1x2'))
			ipercent = read_file(root_path + perthyd_path + file('S', 'isens'))
			ipercentb = read_file(root_path + real_path + file('S', 'isens'))
			ipercentf = read_file(root_path + pertreal_path + file('S', 'isens')) 
		
	return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, ipercent, ipercentb, ipercentf
		
def full_path(order, array):
	if order == 'first':
		real_path, pertreal_path, perthyd_path = main_path('first', pert_layer)
	elif order == 'second':
		real_path, pertreal_path, perthyd_path = main_path('second', pert_layer)	
	elif order == 'cross':
		real_path, pertreal_path, perthyd_path = main_path('cross', pert_layer)
		
	#save the paths to be used later for the one to one plots
	ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, ipercent, ipercentb, ipercentf = get_file(array, order, real_path, pertreal_path, perthyd_path)
	#return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, real_path, pertreal_path, perthyd_path
	return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, perthyd_path, ipercent, ipercentb, ipercentf
	
#Run the main path here:

#Running all plots sequentially:
if ptype == 'S':
#ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, real_pathf, pertreal_pathf, perthyd_pathf = full_path('first', 'S')
	ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, perthyd_pathf, ipercent, ipercentb, ipercentf = full_path('first', 'S')
#ref_realfa, ref_pertfa, hyd_dx1fa, hyd_dx2fa, hyd_dx1x2fa, perthyd_pathfa, ipercenta = full_path('first', 'A')
	ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, perthyd_paths, ipercent2, ipercent2b, ipercent2f = full_path('second', 'S')
#ref_realc, ref_pertc, hyd_dx1c, hyd_dx2c, hyd_dx1x2c, real_pathc, pertreal_pathc, perthyd_pathc = full_path('cross', 'S')
else:
	ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, perthyd_pathf, ipercent = full_path('first', 'A')
	ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, perthyd_paths, ipercent2 = full_path('second', 'A')
#full_path('second', 'A')
#full_path('hybrid', 'A')

#Get the default files
#default_file1, hyd_file1, default_file2, hyd_file2 = default_path(mechanism)

#add unit conversions to pbb from mol/mol
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

#add unit conversions to pbb from mol/mol for one argument
def units_convert(array, useless):
	for i in array:
		if i in useless:
			continue
		else:
			if data_unit_is_mol_per_mol(array[i]):
				array[i].values = array[i].values * 1e9
				array[i].attrs["units"] = "ppb"
	return array


#To create an array of 0s to map each variable and iterate through
lev = 72
lat = 46
lon = 72

#Applies for both Species_Conc and aerosolMass
#useless variables in the array:
useless = ['lat_bnds', 'lon_bnds', 'hyam', 'hybm', 'hyai', 'hybi', 'P0', 'AREA']
num_var = len(ref_realf) - 8
#num_var = len(ref_reals) - 8


#comment out if not needed
#ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2 = unit_convert(ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, useless)

#Running all sequentially:
ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f = unit_convert(ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, useless)
ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s = unit_convert(ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, useless)
#ref_realc, ref_pertc, hyd_dx1c, hyd_dx2c, hyd_dx1x2c = unit_convert(ref_realc, ref_pertc, hyd_dx1c, hyd_dx2c, hyd_dx1x2c, useless)
    		

def first_order_finite(base, pert, h, pert_specie, useless, mode):
    #Index to the species of interest
    #where h is the perturbation, pert_specie is the perturbed variable
    
    species_num = len(base) - 8  #number of variables, ignore useless variables
    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape
    
    ct = 0 #counter
    specie = [] #list of species names
   
    for i in base:
    	if i in useless:
    		continue
    	#if mode == 'A': #post process for units to ug/m3
    	else:
    		
    		if mode == 'semi':
    			array[ct, :, :, :, :] = (pert[i] - base[i]) / (h-1)  #multiplicative
    		elif mode == 'normal':
    		    array[ct, :, :, :, :] = (pert[i] - base[i]) / h
    		#add variable name to list:
    		specie.append(i)
    		
    	ct += 1
    	
    return array, shape, specie

def central_first_order(base, pert, h, pert_specie, useless, mode):
    #Index to the species of interest
    #where h is the perturbation, pert_specie is the perturbed variable
    
    species_num = len(base) - 8  #number of variables, ignore useless variables
    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape
    
    ct = 0 #counter
    specie = [] #list of species names
   
    for i in base:
    	if i in useless:
    		continue
    	#if mode == 'A': #post process for units to ug/m3
    	else:
    		
    		if mode == 'semi':
    			array[ct, :, :, :, :] = (pert[i] - base[i]) / (2 * h)  #multiplicative
    		elif mode == 'normal':
    		    array[ct, :, :, :, :] = (pert[i] - base[i]) / h
    		#add variable name to list:
    		specie.append(i)
    		
    	ct += 1
    	
    return array, shape, specie

    	
# Calculating Hyperdual sensitivities    
    
def hyd_first_order(base, pert_specie, mode, useless, hyd_pert=1):
    
    species_num = len(base) - 8  #number of variables, ignore useless variables
    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape
    
    ct = 0 #counter
    specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	else:
    		array[ct, :, :, :, :] = base[i][:,:,:,:] / hyd_pert #semi normalization (multiplicative)
    		#add variable name to list:
    		specie.append(i)
    	ct += 1
    
    return array, shape, specie  

    
def hybrid_second_order(base, pert, pert_specie, mode, useless, hc, h=1):
	#where pert represents the diagnostic with both real and dual perturbation, while base is the
	#only dual perturbation
  
    
    species_num = len(base) - 8  #number of variables, ignore useless variables
    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape
    
    #perturbation size
    
    
    ct = 0 #counter
    specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	else:  	
    	    if mode == 'semi':
    	        #array[ct, :, :, :, :] = (pert[i][:,:,:,:] - base[i][:,:,:,:]) /((hc-1)*h)
    	        array[ct, :, :, :, :] = (pert[i][:,:,:,:] /(hc*h*(hc-1))) - (base[i][:,:,:,:] /((hc-1)*h))
    	    elif mode == 'normal additive':
    		    array[ct, :, :, :, :] = ((pert[i][:,:,:,:] / h) - (base[i][:,:,:,:] /h)) / h1
    		#add variable name to list:
    	    specie.append(i)
    	ct += 1
    
    return array, shape, specie 
    
def central_hybrid_second_order(base, pert, pert_specie, mode, useless, hc, h=1):
	#where pert represents the diagnostic with both real and dual perturbation, while base is the
	#only dual perturbation
  
    
    species_num = len(base) - 8  #number of variables, ignore useless variables
    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape
    
    #perturbation size
    
    
    ct = 0 #counter
    specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	else:  	
    	    if mode == 'semi':
    	        #array[ct, :, :, :, :] = (pert[i][:,:,:,:] - base[i][:,:,:,:]) /((hc-1)*h)
    	        array[ct, :, :, :, :] = (pert[i][:,:,:,:] /(2*hc*h*(1+hc))) - (base[i][:,:,:,:] /(2*h*hc*(1-hc)))
    	    # elif mode == 'normal additive':
#     		    array[ct, :, :, :, :] = ((pert[i][:,:,:,:] / h) - (base[i][:,:,:,:] /h)) / h1
    		#add variable name to list:
    	    specie.append(i)
    	ct += 1
    
    return array, shape, specie 
    
    
def hyd_second_order(base, pert_specie, mode, useless, hyd_pert=1):
    #Index to the species of interest
    #bspecies = base['PM25dx1x2']
    
    
    species_num = len(base) - 8  #number of variables, ignore useless variables
    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape
    
    ct = 0 #counter
    specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	else:
    	    if mode == 'semi':
    	    	array[ct, :, :, :, :] = base[i][:,:,:,:] / (hyd_pert**2.0)
    	    elif mode == 'normal additive':
    		    array[ct, :, :, :, :] = base[i][:,:,:,:] / (hyd_pert**2.0)
    		#add variable name to list:
    	    specie.append(i)
    	ct += 1
    
    return array, shape, specie
    
def hybrid_cross_sensitivity(base, pert, pert_specie, mode, useless, hc, h=1):
     #Index to the species of interest
    
    
    species_num = len(base) - 8  #number of variables, ignore useless variables
    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape
    
    #perturbation size
    
    ct = 0 #counter
    specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	else:
    		array[ct, :, :, :, :] = ((pert[i][:,:,:,:] / h) - (base[i][:,:,:,:] /h)) / ((hc - 1)*h)
    		#add variable name to list:
    		specie.append(i)
            
    	ct += 1
    
    return array, shape, specie
    
def central_hybrid_cross_sensitivity(base, pert, pert_specie, mode, useless, hc, h=1):
	#where pert represents the diagnostic with both real and dual perturbation, while base is the
	#only dual perturbation
  
    
    species_num = len(base) - 8  #number of variables, ignore useless variables
    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape
    
    #perturbation size
    
    
    ct = 0 #counter
    specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	else:  	
    	    if mode == 'semi':
    	        #array[ct, :, :, :, :] = (pert[i][:,:,:,:] - base[i][:,:,:,:]) /((hc-1)*h)
    	        array[ct, :, :, :, :] = (pert[i][:,:,:,:] /(2*hc*h)) - (base[i][:,:,:,:] /(2*h*hc))
    	    # elif mode == 'normal additive':
#     		    array[ct, :, :, :, :] = ((pert[i][:,:,:,:] / h) - (base[i][:,:,:,:] /h)) / h1
    		#add variable name to list:
    	    specie.append(i)
    	ct += 1
    
    return array, shape, specie 
    
    
def hyd_cross_sensitivity(base, pert_specie, useless, hyd_pert=1):
    #Index to the species of interest
    #bspecies = base['PM25dx1x2']
    
    
    species_num = len(base) - 8  #number of variables, ignore useless variables
    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    shape = array.shape
    
    ct = 0 #counter
    specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue
    	# elif mode == 'A': #post process for units to ug/m3
#     		base[i] = base[i] / 1e+9
#     		array[ct, :, :, :, :] = base[i][:,:,:,:] / (hyd_pert**2.0)
    	else:
    		array[ct, :, :, :, :] = base[i][:,:,:,:] / (hyd_pert**2.0)
    		#add variable name to list:
    		specie.append(i)
    	ct += 1
    
    return array, shape, specie
    
def analytical_second_order(base, pert_specie, useless):
	species_num = len(base) - 8  #number of variables, ignore useless variables
	time, lev, lat, lon = base[pert_specie].shape
	
	array = np.zeros([species_num, time, lev, lat, lon])
	shape = array.shape
	
	ct = 0 #counter
	specie = [] #list of species names
	for i in  base:
		if i in useless:
			continue
			
		else:
			array[ct, :, :, :, :] = base[i][:,:,:,:] * 6.0
			
			#Then convert to ppb:
			array[ct, :, :, :, :] = array[ct, :, :, :, :] * 1e9
			#semi-normalize:
			array[ct, :, :, :, :] = array[ct, :, :, :, :] * (base[i][:,:,:,:] ** 2.0)
			specie.append(i)
		ct += 1
		
	return array, shape, specie
	
    
#use_cmap_RdBu=True
		

def drawmap(shape, dataset, mark, ty, y, x, fname, utype, t):
	#Convert data back to xarray
	#t = '6hr'
	#units = r'$\frac{PM_{2.5} \mu g / m^3}{NH4 \mu g / m^3}$'
	base = xr.DataArray(dataset, coords=shape.coords, dims=shape.dims, attrs=shape.attrs)
	if mark == "hyd":
		if ty == 'first':
			if utype == 'A':
				units = r'$\frac{(' + y.replace("_", "\_") + ' \mu g / m^3)}{(' + x.replace("_", "\_") + ' \mu g / m^3)}$'
			elif utype == 'S':
				units ='ppb'

			gcplot.single_panel(base, title=r'Hyd sensitivities $\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                  r'{{\partial (E_{' + x.replace("SpeciesConc_", "") + '})}_{(x,y,z,t=0hr)}} '
                  r'(E_{' + x.replace("SpeciesConc_", "") + '})$', pdfname="./igcmap/" + fname, comap=custom_cmap, unit=units, vmin=-0.06, vmax=0.06)
			plt.close()
		elif ty == 'second':
			#units = r'${ppb}^2$'
			units = 'ppb'
			gcplot.single_panel(base, title=r'Hyd sensitivities $\frac{{{\partial}^2 (E_{' + y.replace("SpeciesConcdx2_", "") + '})}_{(x,y,z,t=' + t + ')}}'
                  r'{{{\partial (E_{' + x.replace("SpeciesConc_", "") + '})}^2}_{(x,y,z,t=0hr)}}   '
                  r'{(E_{' + x.replace("SpeciesConc_", "") + '})}^2$', pdfname="./igcmap/" + fname, comap=custom_cmap, unit=units, vmin=-0.002, vmax=0.002)
			plt.close()
	else:
		if ty == 'first':
			if utype == 'A':
				units = r'$\frac{(' + y.replace("_", "\_") + ' \mu g / m^3)}{(' + x.replace("_", "\_") + ' \mu g / m^3)}$'
			elif utype == 'S':
				units ='ppb'
			gcplot.single_panel(base, title=r'Finite Difference Sensitivities $\frac{\partial ' + y.replace("_", "\_") + '}{\partial ' + x.replace("_", "\_") + '}$', pdfname="./map/" + fname, unit=units, comap=cm.batlow)
			plt.close()
		elif ty == 'second':
			units = r'ppb'
			gcplot.single_panel(base, title=r'Hybrid Sensitivities $\frac{{\partial}^2 ' + y.replace("_", "\_") + '}{{\partial ' + x.replace("_", "\_") + '}^2}$', pdfname="./map/" + fname, unit=units)
			plt.close()
		#gcplot.single_panel(base, title=r'Finite difference Sensitivities $\frac{\partial PM2.5}{\partial NH4}$', pdfname="./finitePM25", unit=units)
		#gcplot.single_panel(base, title=r'Finite difference Sensitivities $\frac{{\partial}^2 PM_{2.5}}{{\partial siagrowth}\partial NH_4}$', pdfname="./hybridPM25", unit=units)


def sum_abs(base, pert, pert_specie, useless):
    
    species_num = len(base) - 8  #number of variables, ignore useless variables

    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    #shape = array.shape
    
    ct = 0 #counter
    #specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	else:
    		array[ct, :, :, :, :] =base[i][:,:,:,:] - pert[i][:,:,:,:] #take the difference
    		#add variable name to list:
    		#specie.append(i)
    	ct += 1
    
    return array 
    
def absrel(base, pert, pert_specie, useless):
    
    species_num = len(base) - 8  #number of variables, ignore useless variables

    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    #shape = array.shape
    
    ct = 0 #counter
    #specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	else:
    		array[ct, :, :, :, :] =(base[i][:,:,:,:] - pert[i][:,:,:,:]) / base[i][:,:,:,:] #take the relative difference
    		#add variable name to list:
    		#specie.append(i)
    	ct += 1
    
    return array 

def sum_rel(base, pert, pert_specie, useless):
    
    species_num = len(base) - 8  #number of variables, ignore useless variables

    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    #shape = array.shape
    
    ct = 0 #counter
    #specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	else:
    		mask = base[i][:,:,:,:] != 0  # Create a mask to avoid division by zero
    		#base_mask = base[i][:,:,:,:][mask]
    		#array[ct, :, :, :, :] =(pert[i][:,:,:,:] - base[i][:,:,:,:][mask]) /base[i][:,:,:,:][mask]
    		#array[ct] = (pert[i][:,:,:,:] - denominator) / denominator  # relative difference
    		#array[ct][~mask] = 0 # Set relative difference to zero where the denominator is zero
    		array[ct, :, :, :, :] = np.where(mask, (pert[i][:, :, :, :] - base[i][:,:,:,:]) / base[i][:,:,:,:], 0)  # relative difference
    		#array[ct, :, :, :, :][:, ~mask] = 0  
    	ct += 1
        	
       		
        	
        
    		# array[ct, :, :, :, :] =(pert[i][:,:,:,:] - base[i][:,:,:,:]) /base[i][:,:,:,:]  #relative difference
#     		#add variable name to list:
#     		#specie.append(i)
#     	ct += 1
    
    return array 
    
def sum_reln(base, pert, pert_specie, useless):
    
    species_num = len(base) - 8  #number of variables, ignore useless variables

    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    #shape = array.shape
    
    ct = 0 #counter
    #specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	else:
    		base_mask = base[i][:,:,:,:] != 0  # Create a mask to avoid division by zero
    		base_masked = np.where(base_mask, np.nan, base[i][:,:,:,:])
    		pert_masked = np.where(base_mask, np.nan, pert[i][:, :, :, :])
    		#base_mask = base[i][:,:,:,:][mask]
    		array[ct, :, :, :, :] =(pert_masked - base_masked) /base_masked
    		#array[ct] = (pert[i][:,:,:,:] - denominator) / denominator  # relative difference
    		#array[ct][~mask] = 0 # Set relative difference to zero where the denominator is zero
    		#array[ct, :, :, :, :] = np.where(mask, (pert[i][:, :, :, :] - base[i][:,:,:,:]) / base[i][:,:,:,:], 0)  # relative difference
    		#array[ct, :, :, :, :][:, ~mask] = 0  
    	ct += 1
        	
       		
        	
        
    		# array[ct, :, :, :, :] =(pert[i][:,:,:,:] - base[i][:,:,:,:]) /base[i][:,:,:,:]  #relative difference
#     		#add variable name to list:
#     		#specie.append(i)
#     	ct += 1
    
    return array 

#Check for NIT and NH4
def sum_check(base, pert, pert_specie, useless):
    
    species_num = len(base) - 8  #number of variables, ignore useless variables

    time, lev, lat, lon = base[pert_specie].shape
    #array holding all variables 
    array = np.zeros([species_num, time, lev, lat, lon])
    #shape = array.shape
    
    ct = 0 #counter
    #specie = [] #list of species names
    for i in  base:
    	if i in useless:
    		continue

    	elif i == 'SpeciesConc_NIT' or i == 'SpeciesConc_NH4':
    	#elif i == 'SpeciesConc_CO':
    		array[ct, :, :, :, :] =base[i][:,:,:,:] - pert[i][:,:,:,:] #take the difference
    		#add variable name to list:
    		#specie.append(i)
    	else:
    		continue
    	ct += 1
    
    return array
    


def one_to_one1(ax, fd, hyd, mean, meanb, meanf, utype, fname, pert, t, sensitivity, x, y, z='', color='blue', finemask=None, coarsemask=None, finemask2b=None, finemask2f=None, coarsemask2b=None, coarsemask2f=None):
    # Flatten the arrays:
    fd = fd.flatten()
    hyd = hyd.flatten()
    
    t = '1hr'
    ignored_calls = np.mean(mean)
    backward_calls = np.mean(meanb)
    forward_calls = np.mean(meanf)

    # Combine masks for masked and non-masked data points
    if sensitivity == 'first':
        mask = (finemask.flatten() > 0) | (coarsemask.flatten() > 0)
    elif sensitivity == 'second':
        mask = ((finemask.flatten() > 0) | (finemask2b.flatten() > 0) | (finemask2f.flatten() > 0) |
                (coarsemask.flatten() > 0) | (coarsemask2b.flatten() > 0) | (coarsemask2f.flatten() > 0))

    rangemin = min(np.nanmin(fd), np.nanmin(hyd))
    rangemax = max(np.nanmax(fd), np.nanmax(hyd))
    v = [rangemin, rangemax, rangemin, rangemax]
    print(f'The maximum is {rangemax:1.6f}.')
    print(f'The minimum is {rangemin:1.6f}.')

    ax.axis(v)

    # Plot the 1 to 1 line
    ax.plot([rangemin, rangemax], [rangemin, rangemax], 'k-')

    # Scatter plot with specified color and markers
    ax.scatter(fd[~mask], hyd[~mask], c=color, label='No Mask', marker='o', alpha=0.6)
    ax.scatter(fd[mask], hyd[mask], c='gold', label='Mask', marker='x', alpha=0.6)

    # Set axis conditions
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    # Mask of nan numbers and infinite numbers
    valid_mask = ~np.isnan(fd) & ~np.isnan(hyd)
    slope, intercept, r_value, p_value, std_err = stats.linregress(fd[valid_mask], hyd[valid_mask])


    # Defining axis and notations
    ax.minorticks_on()
    ax.yaxis.grid(True, which='major', color='black', linestyle='--', linewidth='0.4')
    ax.yaxis.grid(True, which='minor', color='gray', linestyle=':', linewidth='0.3')
    ax.xaxis.grid(True, which='major', color='black', linestyle='--', linewidth='0.4')
    ax.xaxis.grid(True, which='minor', color='gray', linestyle=':', linewidth='0.3')
    
    if sensitivity == 'first':
        xlabel = (r'$\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                  r'{{\partial [' + x.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}} '
                  r'[' + x.replace("SpeciesConc_", "") + ']$')
        ylabel = (r'$\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                  r'{{\partial (E_{' + x.replace("SpeciesConc_", "") + '})}_{(x,y,z,t=0hr)}} '
                  r'(E_{' + x.replace("SpeciesConc_", "") + '})$')
    elif sensitivity == 'second':
        xlabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                  r'{{{\partial [' + x.replace("SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}  '
                  r'{[' + x.replace("SpeciesConc_", "") + ']}^2$')
        ylabel = (r'$\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}'
                  r'{{{\partial (E_{' + x.replace("SpeciesConc_", "") + '})}^2}_{(x,y,z,t=0hr)}}   '
                  r'{(E_{' + x.replace("SpeciesConc_", "") + '})}^2$')

    # ax.set_xlabel(r'Hybrid $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace("SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}  {[' + x.replace("SpeciesConc_", "") + ']}^2$')
    # ax.set_ylabel(r'Hyd $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace("SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$')
		
    # ax.set_xlabel(xlabel)
    # ax.set_ylabel(ylabel)

    ax.legend()

    ax.set_title(ylabel)

    # Add legend at the bottom right
    ax.legend(loc='lower right')
    
   #  ax.annotate('R$^2$ = {:0.3F}\nSlope = {:0.2F}\nIntercept = {:0.2E}\nBlack Line = 1:1\nIgnored Calls = {:0.3F}%\nForward calls = {:0.3F}%\nBackward calls = {:0.3F}%'.format(r_value**2., slope, intercept, ignored_calls, forward_calls, backward_calls),
#                 xy=(0.05, 0.65), fontsize=12, fontweight='bold', xycoords='axes fraction', 
#                 horizontalalignment='left', verticalalignment='bottom')
    ax.annotate('R$^2$ = {:0.3F}\nSlope = {:0.2F}\nIntercept = {:0.2E}'.format(r_value**2., slope, intercept, ignored_calls, forward_calls, backward_calls),
                xy=(0.05, 0.7), fontsize=14, fontweight='bold', xycoords='axes fraction', 
                horizontalalignment='left', verticalalignment='bottom')

#plot_species = {'ISOP': 173, 'NO2': 131, 'O3': 127, 'NO': 132}
#plot_species = {'ISOP': 173, 'O3': 127, 'NO2': 131} # For ISOP
#plot_species = {'MGLY': 156, 'O3': 127, 'NO2': 131} # For ISOP
plot_species = {'NO2': 131, 'O3': 127, 'HNO3': 217} # For NOx.
ptype = 'S'

def main():
    colors = ['red', 'blue', 'green']  # Different colors for each panel
    #colors = [cm.batlow, cm.lajolla, cm.cork, cm.vik]
    
    # First-order sensitivity plots
    fig1, axs1 = plt.subplots(1, 3, figsize=(18, 6))
    axs1 = axs1.flatten()

    # Second-order sensitivity plots
    fig2, axs2 = plt.subplots(1, 3, figsize=(18, 6))
    axs2 = axs2.flatten()

    if ptype == 'S':
        for idx, (spec, i) in enumerate(plot_species.items()):
            print(i)
            print(specief[i])
            # First-order plot
            fd = fd_sensf[i, 0, :, :, :]
            hyd = hyd_sensf[i, 0, :, :, :]
            mean = isens[0, :, :, :]
            meanb = isensb[0, :, :, :]
            meanf = isensf[0, :, :, :]
            fine = finemask[0, :, :, :]
            coarse = coarsemask[0, :, :, :]
            one_to_one1(axs1[idx], fd, hyd, mean, meanb, meanf, 'S', 'one-to-one_' + duration + '_first_' + mechanism + '_' + species, pert_layer, duration, sensitivity='first', x='SpeciesConc_' + specie1, y=specief[i], color=colors[idx], finemask=fine, coarsemask=coarse)
            
            # Second-order plot
            fd = fd_senss[i, 0, :, :, :]
            hyd = hyd_senss[i, 0, :, :, :]
            mean = isens2[0, :, :, :]
            meanb = isens2b[0, :, :, :]
            meanf = isens2f[0, :, :, :]
            fine2 = finemask2[0, :, :, :]
            coarse2 = coarsemask2[0, :, :, :]
            fine2b = finemask2b[0, :, :, :]
            fine2f = finemask2f[0, :, :, :]
            coarse2b = coarsemask2b[0, :, :, :]
            coarse2f = coarsemask2f[0, :, :, :]
            one_to_one1(axs2[idx], fd, hyd, mean, meanb, meanf, 'S', 'one-to-one_' + duration + '_second_' + mechanism + '_' + species, pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie1, y=speciess[i], color=colors[idx], finemask=fine2, coarsemask=coarse2, finemask2b=fine2b, finemask2f=fine2f, coarsemask2b=coarse2b, coarsemask2f=coarse2f)

    fig1.text(0.5, 0.04, 'Finite difference [ppb]', ha='center', fontsize=14)
    fig1.text(0.04, 0.5, 'Hyd [ppb]', va='center', rotation='vertical', fontsize=14)
    
    # Second-order plot labels
    fig2.text(0.5, 0.04, 'Hybrid [ppb]', ha='center', fontsize=14)
    fig2.text(0.04, 0.5, 'Hyd [ppb]', va='center', rotation='vertical', fontsize=14)
    
    # fig1.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.3, hspace=0.3)
    # fig2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.3, hspace=0.3)
    fig1.subplots_adjust(left=0.1, bottom=0.2)
    fig2.subplots_adjust(left=0.1, bottom=0.2)
    # fig1.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15, wspace=0.4, hspace=0.4)
    # fig2.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15, wspace=0.4, hspace=0.4)
    
    # Save both figures to PDF
    #pdfFile1 = PdfPages('four_panel_plot_first_order_NOx_shipemis_1hr_third.pdf')
    pdfFile1 = PdfPages('three_panel_plot_first_order_NOx_emis_1hr_final20%.pdf')
    #plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
    plt.tight_layout()
    pdfFile1.savefig(fig1)
    pdfFile1.close()
    
    #pdfFile2 = PdfPages('four_panel_plot_second_order_NOx_shipemis_1hr_third.pdf')
    pdfFile2 = PdfPages('three_panel_plot_second_order_NOx_emis_1hr_final20%.pdf')
    # plt.tight_layout(rect=[0.05, 0.05, 1, 1])
    #plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
    plt.tight_layout()
    pdfFile2.savefig(fig2)
    pdfFile2.close()
    
    
    plt.close(fig1)
    plt.close(fig2)
    
    # plot_var = {'O3': 127} #For NOx
#     for i in plot_var.values():
#     	#plot first order
#     	drawmap(shaped, hyd_sensf[i, 0, 0, :, :], "hyd", "first", y=specief[i], x=specie1, fname=specief[i] +'_NOxemis_1hr.pdf', utype='S', t='1hr')
#     	#plot second order
#     	drawmap(shaped, hyd_senss[i, 0, 0, :, :], "hyd", "second", y=speciess[i], x=specie1, fname=specie_hyds[i] + '_NOxemis_1hr.pdf', utype='S', t='1hr')
    
    #Do the mapped plots
	#plot_var = {'ISOP': 173, 'O3': 127} For ISOP
	
	#shaped = ref_realf['SpeciesConc_CO'].isel(time=0, lev=0)

#Get the reversed colormap:
rev_WhGrYlRd = WhGrYlRd.reversed()
white_cmap = ListedColormap(['white']) #white colormap

colors = [(0, 0, 0.5), (1, 1, 1), (0.5, 0, 0)]  # Deep red, white, deep blue
positions = [0, 0.5, 1]  # Positions of the colors

# Create the colormap
custom_cmap = LinearSegmentedColormap.from_list("deep_red_white_deep_blue", list(zip(positions, colors)))
	
	

#fd_sens, shape, specie = first_order_finite(ref_real, ref_pert, h, 'AerMassNH4', useless)
#hyd_sens, shape, specie_hyd = hyd_first_order(hyd_dx2, 'AerMassNH4dx2', 'A', useless, hyd_pert=1)

# hyd_sens, shape, specie_hyd = hyd_first_order(hyd_dx2, 'SpeciesConcdx2_CO', 'S', useless, hyd_pert=1)
# fd_sens, shape, specie = first_order_finite(ref_real, ref_pert, h, 'SpeciesConc_CO', useless, 'semi')
h1 = 0.2
hc = 0.2


#Running all plots at once:
#Write a loop to handle S or A:
if ptype == "S":
	hyd_sensf, shape, specie_hydf = hyd_first_order(hyd_dx2f, 'SpeciesConcdx2_' + species, 'S', useless, hyd_pert=1)
	fd_sensf, shape, specief = central_first_order(ref_realf, ref_pertf, h1, 'SpeciesConc_' + species, useless, 'semi')
#hyd_sensf, shape, specie_hydf = hyd_first_order(hyd_dx2f, 'SpeciesConcdx2_' + specie1, 'S', useless, hyd_pert=1)
#fd_sensf, shape, specief = central_first_order(ref_realf, ref_pertf, h1, 'SpeciesConc_' + specie1, useless, 'semi')

	hyd_senss, shape, specie_hyds = hyd_second_order(hyd_dx1x2s, 'SpeciesConcdx1x2_' + species, 'semi', useless, hyd_pert=1)
	fd_senss, shape, speciess = central_hybrid_second_order(ref_reals, ref_perts, 'SpeciesConcdx2_' + species, 'semi', useless, hc, h=1)

#Aerosol mass case:
else:

#First order
	hyd_sensf, shape, specie_hydf = hyd_first_order(hyd_dx2f, 'AerMassNH4dx2', 'A', useless, hyd_pert=1)
	fd_sensf, shape, specief = central_first_order(ref_realf, ref_pertf, h1, 'AerMassNH4', useless, 'semi')
#Second order
	hyd_senss, shape, specie_hyds = hyd_second_order(hyd_dx1x2s, 'AerMassNH4dx1x2', 'semi', useless, hyd_pert=1)
	fd_senss, shape, speciess = central_hybrid_second_order(ref_reals, ref_perts, 'AerMassNH4dx2', 'semi', useless, hc, h=1)




#hyd_sensc, shape, specie_hydc = hyd_cross_sensitivity(hyd_dx1x2c, 'SpeciesConcdx1x2_' + species, useless, hyd_pert=1)
#fd_sensc, shape, speciec = central_hybrid_cross_sensitivity(ref_realc, ref_pertc, 'SpeciesConcdx2_' + species, 'semi', useless, hc, h=1)

#Find the indicies of negative values
#negative_indices = np.where(hyd_sensf[251, -1, :, :, :] < 0)
#negative_indices = np.where(hyd_sensf[251, -1, :, :, :] < 0) #Second order sensitivity

# Calculate the absolute difference between hyd_sensf and fd_sensf
abs_diff = np.abs(hyd_sensf - fd_sensf)
#abs_diff = np.abs(hyd_senss - fd_senss)
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
found = False
for i in range(abs_diff.shape[1]):
    diff_indices = np.where(abs_diff[127, i, :, :, :] > 10)  # NH3 O3- 127
    if any(len(indices) > 0 for indices in diff_indices):
        print(f"Indices with absolute difference greater than 10 O3 at {i} index timestep:")
        for dim, indices in enumerate(diff_indices):
            print(f"Dimension {dim}: {indices}")
        found = True
        break

if not found:
    print("No indices with absolute difference greater than 10 found in any second index timestep.")
    
second = False
for i in range(abs_diff.shape[1]):
    diff_indices2 = np.where(abs_diff[217, i, :, :, :] > 2.5)  # HNO3 
    if any(len(indices) > 0 for indices in diff_indices2):
        print(f"Indices with absolute difference greater than 2.5 HNO3 at {i} index timestep:")
        for dim, indices in enumerate(diff_indices2):
            print(f"Dimension {dim}: {indices}")
        second = True
        break

if not second:
    print("No indices with absolute difference greater than 10 found in any second index timestep.")

    
#Obtain the sensitivity mean
#The percentage will be the same across all species, just one variable = ipercent
isens = ipercent['Ipercent'][:,:,:,:].values
isensb = ipercentb['Ipercent'][:,:,:,:].values
isensf = ipercentf['Ipercent'][:,:,:,:].values
isens2 = ipercent2['Ipercent'][:,:,:,:].values
isens2b = ipercent2b['Ipercent'][:,:,:,:].values
isens2f = ipercent2f['Ipercent'][:,:,:,:].values

#Obtain the mask values
finemask = ipercent['Finemask'][:,:,:,:].values
finemaskb = ipercentb['Finemask'][:,:,:,:].values
finemaskf = ipercentf['Finemask'][:,:,:,:].values
finemask2 = ipercent2['Finemask'][:,:,:,:].values
finemask2b = ipercent2b['Finemask'][:,:,:,:].values
finemask2f = ipercent2f['Finemask'][:,:,:,:].values

coarsemask = ipercent['Coarsemask'][:,:,:,:].values
coarsemaskb = ipercentb['Coarsemask'][:,:,:,:].values
coarsemaskf = ipercentf['Coarsemask'][:,:,:,:].values
coarsemask2 = ipercent2['Coarsemask'][:,:,:,:].values
coarsemask2b = ipercent2b['Coarsemask'][:,:,:,:].values
coarsemask2f = ipercent2f['Coarsemask'][:,:,:,:].values

#Define mask variables:
#fine_mask_hyd1 = ipercent['Finemask'] 
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

#Test to be sure you got the write array
#print("Testing the index 19, 71")
#print(hyd_sensf[251, 1, 26, 19, 71])

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
# hyd_sensc, shape, specie_hydc = hyd_cross_sensitivity(hyd_dx1x2c, 'SpeciesConcdx1x2_' +specie1, useless, hyd_pert=1)
# fd_sensc, shape, speciec = hybrid_cross_sensitivity(ref_realc, ref_pertc, 'SpeciesConcdx2_' + specie1, 'semi', useless, hc, h=1)

#anal_sens, shape, speciea = analytical_second_order(ref_realf, 'SpeciesConc_CO', useless)

#fd_sens = hybrid_cross_sensitivity(ref_real, ref_pert, 1)
#hyd_sens = hyd_cross_sensitivity(hyd_pert, 1)

#to store dimensions!
shaped = ref_realf['SpeciesConc_CO'].isel(time=0, lev=0)

#For the sum plots:
#Below lines not needed any more:
#First Calculate the difference species conc
# array1 = sum_abs(default_file1, hyd_file1, 'SpeciesConc_' + specie1, useless)
# 
# #Calculate the difference aerosolmass
# array2 = sum_abs(default_file2, hyd_file2, 'PM10', useless)
# 
# #Sum check for NH4 and NIT
# array3 = sum_check(default_file1, hyd_file1, 'SpeciesConc_' + specie1, useless)
# 
# #Calculate the relative difference for species conc
# array4 = absrel(default_file1, hyd_file1, 'SpeciesConc_' + specie1, useless)
# 
# #Calculate the relative difference aerosolmass
# array5 = absrel(default_file2, hyd_file2, 'PM10', useless)
# 
# # Calculate the sum of absolute differences across all data variables
# sum_diff1 = np.sum(np.abs(array1))
# sum_diff2 = np.sum(np.abs(array2))
# sum_diff3 = np.sum(np.abs(array3))
# 
# #Obtain the max absolute relative difference:
# #first obtain finite values:
# finite_values1 = np.isfinite(array4)
# finite_values2 = np.isfinite(array5)
# # max_diff1 = np.max(np.abs(array4.where(finite_values1).values))
# # max_diff2 = np.max(np.abs(array5.where(finite_values2).values))
# max_diff1 = np.max(np.abs(array4[finite_values1]))
# max_diff2 = np.max(np.abs(array5[finite_values2]))
# 
# #Obtain the min absolute relative difference:
# # min_diff1 = np.min(np.abs(array4.where(finite_values1).values))
# # min_diff2 = np.min(np.abs(array5.where(finite_values2).values))
# min_diff1 = np.min(np.abs(array4[finite_values1]))
# min_diff2 = np.min(np.abs(array5[finite_values2]))
# 
# 
# # # Specify the output filename
# # output_filename = 'sum_diff.txt'
# # 
# # # Save the sum_diff to the text file
# # save_sum_diff_to_file(sum_diff, output_filename)
# 
# #Print the absolute difference to the command window
# print('sum of the abs of def - hyd species conc')
# print(sum_diff1)
# print('sum of the abs of def - hyd aerosol mass')
# print(sum_diff2)
# print('sum of the abs of def - hyd for NH4 and NIT in ppb')
# print(sum_diff3)
# 
# #Sum of base 
# 
# #Sum diff of 1 divided by sum of base to give the relative difference
# 
# #To obtain the absolute relative difference:
# print('maximum of the abs relative def - hyd species conc')
# print(max_diff1)
# print('maximum of the abs relative def - hyd aerosol mass')
# print(max_diff1)
# print('minimum of the abs relative def - hyd species conc')
# print(min_diff1)
# print('minimum of the abs relative def - hyd aerosol mass')
# print(min_diff1)

main()