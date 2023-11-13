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
import gcpy.plot as gcplot
from cmcrameri import cm
from gcpy.units import check_units, data_unit_is_mol_per_mol
import os

# Ignore warnings
warnings.filterwarnings("ignore")

# Reading Files

root_path = '/glade/work/sakinjole/rundirs/gc_4x5_merra2_fullchem'
specie1 = 'CO'
specie2 = 'H2'
duration = '1hr'
mechanism = 'fullchem'
pert_layer = '/allcells'  # else ground layer

#A dictionary for species variables and indices to plot (This serves as the numerator)
#Only for the species concentration diagnostic
#plot_species = {'CO':251, 'CO2':86, 'SO4':94, 'SO2':95, 'NO2':131, 'NO':132, 'O3':127, 'NH4': 135, 'NH3': 136}
plot_species = {'CO2':251, 'SO4':94, 'SO2':95, 'NO2':131, 'NO':132, 'O3':127, 'NH4': 135, 'NH3': 136}


#A list of species indices to plot

#Extending by adding sum of absolute difference code for both SpeciesConc and AerosolMass
def default_path(mechanism):
	#Folder paths default run
	real_path1 = root_path + '/default/' + mechanism + '/' + duration + '/GEOSChem.SpeciesConc.20190701_0000z.nc4'
	hyd_path1 = root_path + '/default_hyd/' + mechanism + '/' + duration + '/GEOSChem.SpeciesConc.20190701_0000z.nc4'
	real_path2 = root_path + '/default/' + mechanism + '/' + duration + '/GEOSChem.AerosolMass.20190701_0000z.nc4'
	hyd_path2 = root_path + '/default_hyd/' + mechanism + '/' + duration + '/GEOSChem.AerosolMass.20190701_0000z.nc4'
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
			
	elif ty == 'S':
		if pt == 'real':
			return 'GEOSChem.SpeciesConc.20190701_0000z.nc4'
		elif pt == 'dx1':
			return 'GEOSChem.SpeciesConcdx1.20190701_0000z.nc4'
		elif pt == 'dx2':
			return 'GEOSChem.SpeciesConcdx2.20190701_0000z.nc4'
		elif pt == 'dx1x2':
			return 'GEOSChem.SpeciesConcdx1x2.20190701_0000z.nc4'
 

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
		elif order == 'second':
			ref_real = read_file(root_path + real_path + file('A', 'dx2'))
			ref_pert = read_file(root_path + pertreal_path + file('A', 'dx2'))
			hyd_dx1 = read_file(root_path + perthyd_path + file('A', 'dx1'))
			hyd_dx2 = read_file(root_path + perthyd_path + file('A', 'dx2'))
			hyd_dx1x2 = read_file(root_path + perthyd_path + file('A', 'dx1x2'))
			
	elif array == 'S':
		if order == 'first':
			ref_real = read_file(root_path + real_path + file('S', 'real'))
			ref_pert = read_file(root_path + pertreal_path + file('S', 'real'))
			hyd_dx1 = read_file(root_path + perthyd_path + file('S', 'dx1'))
			hyd_dx2 = read_file(root_path + perthyd_path + file('S', 'dx2'))
			hyd_dx1x2 = read_file(root_path + perthyd_path + file('S', 'dx1x2')) 
		elif order == 'second':
			ref_real = read_file(root_path + real_path + file('S', 'dx2'))
			ref_pert = read_file(root_path + pertreal_path + file('S', 'dx2'))
			hyd_dx1 = read_file(root_path + perthyd_path + file('S', 'dx1'))
			hyd_dx2 = read_file(root_path + perthyd_path + file('S', 'dx2'))
			hyd_dx1x2 = read_file(root_path + perthyd_path + file('S', 'dx1x2'))
		elif order == 'cross':
			ref_real = read_file(root_path + real_path + file('S', 'dx2'))
			ref_pert = read_file(root_path + pertreal_path + file('S', 'dx2'))
			hyd_dx1 = read_file(root_path + perthyd_path + file('S', 'dx1'))
			hyd_dx2 = read_file(root_path + perthyd_path + file('S', 'dx2'))
			hyd_dx1x2 = read_file(root_path + perthyd_path + file('S', 'dx1x2')) 
		
	return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2
		
def full_path(order, array):
	if order == 'first':
		real_path, pertreal_path, perthyd_path = main_path('first', pert_layer)
	elif order == 'second':
		real_path, pertreal_path, perthyd_path = main_path('second', pert_layer)	
	elif order == 'cross':
		real_path, pertreal_path, perthyd_path = main_path('cross', pert_layer)
		
	#save the paths to be used later for the one to one plots
	ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2 = get_file(array, order, real_path, pertreal_path, perthyd_path)
	return ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, real_path, pertreal_path, perthyd_path
	
#Run the main path here:

#Running all plots sequentially:
ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, real_pathf, pertreal_pathf, perthyd_pathf = full_path('first', 'S')
ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, real_paths, pertreal_paths, perthyd_paths = full_path('second', 'S')
ref_realc, ref_pertc, hyd_dx1c, hyd_dx2c, hyd_dx1x2c, real_pathc, pertreal_pathc, perthyd_pathc = full_path('cross', 'S')
#full_path('second', 'A')
#full_path('hybrid', 'A')

#Get the default files
default_file1, hyd_file1, default_file2, hyd_file2 = default_path(mechanism)

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


#Calculate finite differences

#Calculating the defined perturbation size for a variable

#To create an array of 0s to map each variable and iterate through
lev = 72
lat = 46
lon = 72

#Applies for both Species_Conc and aerosolMass
#useless variables in the array:
useless = ['lat_bnds', 'lon_bnds', 'hyam', 'hybm', 'hyai', 'hybi', 'P0', 'AREA']
num_var = len(ref_realf) - 8


#comment out if not needed
#ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2 = unit_convert(ref_real, ref_pert, hyd_dx1, hyd_dx2, hyd_dx1x2, useless)

#Running all sequentially:
ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f = unit_convert(ref_realf, ref_pertf, hyd_dx1f, hyd_dx2f, hyd_dx1x2f, useless)
ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s = unit_convert(ref_reals, ref_perts, hyd_dx1s, hyd_dx2s, hyd_dx1x2s, useless)
ref_realc, ref_pertc, hyd_dx1c, hyd_dx2c, hyd_dx1x2c = unit_convert(ref_realc, ref_pertc, hyd_dx1c, hyd_dx2c, hyd_dx1x2c, useless)

    		

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
    	# elif mode == 'A': #post process for units to ug/m3
#     		base[i] = base[i] / 1e+9
#     		pert[i] = pert[i] / 1e+9
#     		array[ct, :, :, :, :] = ((pert[i][:,:,:,:] / h) - (base[i][:,:,:,:] /h)) / h1
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
	
    


def one_to_one(fd, hyd, fname, pert, t, sensitivity, x, y, z=''):  #x and y are strings for bottom and up, t is simulation time

	#Flatten the arrays:
	fd = fd.flatten()
	hyd = hyd.flatten()
	
	#Save files to the appropriate concentration file folder. In this case saving only at the hyd path
	if sensitivity == 'first':
		#pdfFile = PdfPages('/glade/work/sakinjole/rundirs/gc_4x5_merra2_fullchem/1hr_plot/'+ fname +'.pdf')
		pdfFile = PdfPages(root_path + perthyd_pathf + fname +'.pdf')
	elif sensitivity == 'second':
		#pdfFile = PdfPages('/glade/work/sakinjole/rundirs/gc_4x5_merra2_fullchem/plot_second/'+ fname +'.pdf')
		pdfFile = PdfPages(root_path + perthyd_paths + fname +'.pdf')
	elif sensitivity == 'cross':
		#pdfFile = PdfPages('/glade/work/sakinjole/rundirs/gc_4x5_merra2_fullchem/plot_cross/'+ fname +'.pdf')
		pdfFile = PdfPages(root_path + perthyd_pathc + fname +'.pdf')
	rangemin = min(np.nanmin(fd), np.nanmin(hyd))
	rangemax = max(np.nanmax(fd), np.nanmax(hyd))
	v = [rangemin, rangemax, rangemin, rangemax] 
	print(f'The maximum is {rangemax:1.2f}.')
	print(f'The minimum is {rangemin:1.2f}.')

	plt.axis(v)

	#to plot the 1 to 1 line
	plt.plot([rangemin, rangemax], [rangemin, rangemax], 'k-')

	#retrieving the axis to ax
	ax = plt.gca()

	#setting axis conditions
	# ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
# 	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	
	#Adopting Ed's code here for axis conditions
	ax.set_xticks(np.linspace(rangemin, rangemax, 6))
	ax.set_yticks(np.linspace(rangemin, rangemax, 6))
	
	# Copied from https://stackoverflow.com/questions/42142144/displaying-first-decimal-digit-in-scientific-notation-in-matplotlib
	class ScalarFormatterForceFormat(ScalarFormatter):
		def _set_format(self):  # Override function that finds format to use.
			self.format = "%1.2f"  # Give format here
	yfmt = ScalarFormatterForceFormat()
	yfmt.set_powerlimits((0,0))
	ax.yaxis.set_major_formatter(yfmt)
	ax.xaxis.set_major_formatter(yfmt)
	
	#plot with colormap
	cplot = ax.scatter(fd, hyd, c=hyd, cmap=cm.batlow)
	if sensitivity == 'first':
		cbar = plt.colorbar(cplot, ax=ax, label='ppb')
	else:
		cbar = plt.colorbar(cplot, ax=ax, label=r'${ppb}^2$')
		#cbar = plt.colorbar(cplot, ax=ax, label=r'ppb')
	#cbar = plt.colorbar(cplot, ax=ax, label=r'$\frac{{ppb}^2}{{ppb}^2}')
	#cbar = plt.colorbar(cplot, ax=ax)
	#ax.legend(shadow=True, loc="upper left")
	# Mask of nan numbers and infinite numbers
	# Using isfinte to mask for the two scenarios
	mask = ~np.isnan(fd) & ~np.isnan(hyd)
	slope, intercept, r_value, p_value, std_err = stats.linregress(fd[mask], hyd[mask])
	#slope, intercept, r_value, p_value, std_err = stats.linregress(fd, hyd)

	#defining axis and notations
	ax.minorticks_on()
	ax.yaxis.grid(True, which='major', color='black', linestyle='--', linewidth='0.4')
	ax.yaxis.grid(True, which='minor', color='gray', linestyle=':', linewidth='0.3',)
	ax.xaxis.grid(True, which='major', color='black', linestyle='--', linewidth='0.4')
	ax.xaxis.grid(True, which='minor', color='gray', linestyle=':', linewidth='0.3')
	
	if sensitivity == 'first':
		if pert == '/allcells':
		#ax.set_xlabel(r'finite difference $\frac{{\partial ' + y.replace("_", "\_") + '}_{(x,y,z,t=1)}}{{\partial ' + x.replace("_", "\_") + '}_{(x,y,z,t=1)}}$') #r treats as raw strings, escape code ignored
		#ax.set_ylabel(r'Hyd $\frac{{\partial ' + y.replace("_", "\_") + '}_{(x,y,z,t=1)}}{{\partial ' + x.replace("_", "\_") + '}_{(x,y,z,t=1)}}$')	
			ax.set_xlabel(r'finite difference $\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$') #r treats as raw strings, escape code ignored
			ax.set_ylabel(r'Hyd $\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$')
		else: #ground layer perturbation
			ax.set_xlabel(r'finite difference $\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$') #r treats as raw strings, escape code ignored
			ax.set_ylabel(r'Hyd $\frac{{\partial [' + y.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + ']$')
	elif sensitivity == 'second':
		if pert == '/allcells':
		# ax.set_xlabel(r'Hybrid finite difference $\frac{{\partial}^2 ' + y.replace("_", "\_").replace("dx2", "") + '}{{\partial ' + x.replace("_", "\_") + '}^2}$')
# 		ax.set_ylabel(r'Hyd $\frac{{\partial}^2 ' + y.replace("_", "\_").replace("dx2", "") + '}{{\partial ' + x.replace("_", "\_") + '}^2}$')
			ax.set_xlabel(r'Hybrid $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace("SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}  {[' + x.replace("SpeciesConc_", "") + ']}^2$')
			ax.set_ylabel(r'Hyd $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace("SpeciesConc_", "") + ']}^2}_{(x,y,z,t=0hr)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$')
		else:
			ax.set_xlabel(r'Hybrid $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace("SpeciesConc_", "") + ']}^2}_{(x,y,z=1,t=0hr)}}  {[' + x.replace("SpeciesConc_", "") + ']}^2$')
			ax.set_ylabel(r'Hyd $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{{\partial [' + x.replace("SpeciesConc_", "") + ']}^2}_{(x,y,z=1,t=0hr)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$')
	elif sensitivity == 'cross':
		if pert == '/allcells':
		# ax.set_xlabel(r'finite difference Sensitivities $\frac{{\partial}^2 ' + y.replace("_", "\_") + '}{{\partial ' + x.replace("_", "\_") + '}\partial ' + z + '}$')
# 		ax.set_ylabel(r'Hyd Sensitivities $\frac{{\partial}^2 ' + y.replace("_", "\_") + '}{{\partial ' + x.replace("_", "\_") + '}\partial ' + z.replace("_", "\_") + '}$')
			ax.set_xlabel(r'Hybrid  $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}}\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
			ax.set_ylabel(r'Hyd  $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace("SpeciesConc_", "") + ']}_{(x,y,z,t=0hr)}}\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
		else:
			ax.set_xlabel(r'Hybrid  $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}}\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
			ax.set_ylabel(r'Hyd  $\frac{{{\partial}^2 [' + y.replace("SpeciesConcdx2_", "") + ']}_{(x,y,z,t=' + t + ')}}{{\partial [' + x.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}}\partial [' + z.replace("SpeciesConc_", "") + ']}_{(x,y,z=1,t=0hr)}} [' + x.replace("SpeciesConc_", "") + '][' + z.replace("SpeciesConc_", "") + ']$')
		#ax.set_ylabel(r'Hyd $\frac{{\partial}^2 ' + y.replace("SpeciesConcdx2_", "") + '}{{\partial ' + x.replace("SpeciesConc_", "") + '}\partial ' + z.replace("SpeciesConc_", "") + '}$')
	#ax.set_title(r'Graph of $\frac{\partial b}{\partial x}$')
	#ax.set_title(r'Graph of ${CO}^3$')
	ax.set_title('Fullchem simulation')
	#ax.set_title(r'Hyd sensitivities against finite sensitivities')
	ax.annotate('R$^2$ = {:0.3F}\nSlope = {:0.2F}\nIntercept = {:0.2E}\nBlack Line = 1:1'.format(r_value**2., slope, intercept), xy=(0.05, 0.7), fontsize=14, fontweight='bold', xycoords='axes fraction', \
horizontalalignment='left', verticalalignment='bottom')
	plt.tight_layout()
	pdfFile.savefig()
	pdfFile.close()
	plt.close()
		

def drawmap(shape, dataset, mark, ty, y, x, fname, utype):
	#Convert data back to xarray
	
	#units = r'$\frac{PM_{2.5} \mu g / m^3}{NH4 \mu g / m^3}$'
	base = xr.DataArray(dataset, coords=shape.coords, dims=shape.dims, attrs=shape.attrs)
	if mark == "hyd":
		if ty == 'first':
			if utype == 'A':
				units = r'$\frac{(' + y.replace("_", "\_") + ' \mu g / m^3)}{(' + x.replace("_", "\_") + ' \mu g / m^3)}$'
			elif utype == 'S':
				units ='ppb'
			#testing
			#units = r'$\frac{( {y} \mu g / m^3)}{( {x} \mu g / m^3)}$'
			#funits = units.format(y=y, x=x)
		#gcplot.single_panel(base, title=r'HYD Sensitivities $\frac{\partial PM2.5}{\partial NH4}$', pdfname="./hydPM25", unit=units)
			gcplot.single_panel(base, title=r'HYD Sensitivities $\frac{\partial ' + y.replace("_", "\_") + '}{\partial ' + x.replace("_", "\_") + '}$', pdfname="./map/" + fname, unit=units, comap=cm.batlow)
			plt.close()
		elif ty == 'second':
			units = r'${ppb}^2$'
			gcplot.single_panel(base, title=r'Hyd sensitivities $\frac{{{\partial}^2 ' + y.replace("SpeciesConcdx2_", "") + '}_{(x,y,z,t=20min)}}{{{\partial ' + x.replace("SpeciesConc_", "") + '}^2}_{(x,y,z,t=0min)}}   {[' + x.replace("SpeciesConc_", "") + ']}^2$', pdfname="./map/" + fname, unit=units, comap=cm.batlow)
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
			units = r'${ppb}^2$'
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

def main():
	#make the one-to-one plot: for Aerosolmass
	# for i in range(8):
# 		one_to_one(fd_sens[i, 0, 0, :, :], hyd_sens[i, 0, 0, :, :], 'one-to-one_1hr_aerosol_' + specie[i], sensitivity='first', x='AerMassNH4', y=specie[i])
	

	
	#add a 
	# for i in range(307):
# 		print(i)
# 		print(specie[i])
# 		one_to_one(fd_sens[i, 1, 0:39, :, :], hyd_sens[i, 1, 0:39, :, :], 'one-to-one_20min_fullchem_' + specie[i], sensitivity='cross', x='SpeciesConc_NO2', y=specie[i], z='SpeciesConc_OH')
		
	#Running all plots at once for Species Conc
	for i in plot_species.values():
 			
		
		print(i)
		print(specief[i])
		#First order:
			#one_to_one(fd_sensf[i, 1, :58, :, :], hyd_sensf[i, 1, :58, :, :], 'one-to-one_20min_fullchem_' + specief[i], sensitivity='first', x='SpeciesConc_NO', y=specief[i])
		one_to_one(fd_sensf[i, -1, :, :, :], hyd_sensf[i, -1, :, :, :], 'one-to-one_' + duration + '_first_' + fullchem + '_' + specief[i], pert_layer, duration, sensitivity='first', x='SpeciesConc_' + specie1, y=specief[i])
		#Second order
			#one_to_one(fd_senss[i, 1, :58, :, :], hyd_senss[i, 1, :58, :, :], 'one-to-one_20min_fullchem_' + species[i], sensitivity='second', x='SpeciesConc_NO', y=species[i])
		one_to_one(fd_senss[i, -1, :, :, :], hyd_senss[i, -1, :, :, :], 'one-to-one_' + duration + '_second_' + fullchem + '_' + species[i], pert_layer, duration, sensitivity='second', x='SpeciesConc_' + specie1, y=species[i])
			#one_to_one(fd_senss[i, 1, :, :, :], anal_sens[i, 1, :, :, :], 'one-to-one_20min_fullchem_' + species[i], sensitivity='second', x='SpeciesConc_SO2', y=species[i])
		#Cross order
		one_to_one(fd_sensc[i, -1, :, :, :], hyd_sensc[i, -1, :, :, :], 'one-to-one_' + duration + '_cross_' + fullchem '_' + speciec[i], pert_layer, duration, sensitivity='cross', x='SpeciesConc_' + specie1, y=speciec[i], z='SpeciesConc_' + specie2)
		# else:
# 			continue
	#hyd sensitivity map plot
	#for species_conc:
	# for i in range(307):
# 		if i == 86 or i == 251: #CO and CO2
#  			drawmap(shaped, hyd_senss[i, 1, 0, :, :], "hyd", "second", y=species[i], x='SpeciesConc_CO', fname=specie_hyds[i], utype='S')
 		#drawmap(shaped, hyd_sens[i, 1, 0, :, :], "hyd", "second", y=specie[i], x='SpeciesConc_CO', fname=specie_hydf[i], utype='S')
# 	
# 	#finite difference:
	# for i in range(307):
#   		drawmap(shaped, fd_sens[i, 0, 0, :, :], "finite", "first", y=specie[i], x='SpeciesConc_CO', fname=specie[i], utype='S')
# 	
	#finite difference map plot
	#drawmap(shape, fd_sens, mark='')
	
	

#fd_sens, shape, specie = first_order_finite(ref_real, ref_pert, h, 'AerMassNH4', useless)
#hyd_sens, shape, specie_hyd = hyd_first_order(hyd_dx2, 'AerMassNH4dx2', 'A', useless, hyd_pert=1)

# hyd_sens, shape, specie_hyd = hyd_first_order(hyd_dx2, 'SpeciesConcdx2_CO', 'S', useless, hyd_pert=1)
# fd_sens, shape, specie = first_order_finite(ref_real, ref_pert, h, 'SpeciesConc_CO', useless, 'semi')
h1 = 0.1
hc = 0.1


#Running all plots at once:
hyd_sensf, shape, specie_hydf = hyd_first_order(hyd_dx2f, 'SpeciesConcdx2_' + specie1, 'S', useless, hyd_pert=1)
fd_sensf, shape, specief = central_first_order(ref_realf, ref_pertf, h1, 'SpeciesConc_' + specie1, useless, 'semi')

hyd_senss, shape, specie_hyds = hyd_second_order(hyd_dx1x2s, 'SpeciesConcdx1x2_' + specie1, 'semi', useless, hyd_pert=1)
fd_senss, shape, species = central_hybrid_second_order(ref_reals, ref_perts, 'SpeciesConcdx2_' + specie1, 'semi', useless, hc, h=1)

hyd_sensc, shape, specie_hydc = hyd_cross_sensitivity(hyd_dx1x2c, 'SpeciesConcdx1x2_' +specie1, useless, hyd_pert=1)
fd_sensc, shape, speciec = hybrid_cross_sensitivity(ref_realc, ref_pertc, 'SpeciesConcdx2_' + specie1, 'semi', useless, hc, h=1)

#anal_sens, shape, speciea = analytical_second_order(ref_realf, 'SpeciesConc_CO', useless)

#fd_sens = hybrid_cross_sensitivity(ref_real, ref_pert, 1)
#hyd_sens = hyd_cross_sensitivity(hyd_pert, 1)

#to store dimensions!
shaped = ref_realf['SpeciesConc_CO'].isel(time=0, lev=0)

#For the sum plots:
#First Calculate the difference species conc
array1 = sum_abs(default_file1, hyd_file1, 'SpeciesConc_' + specie1, useless)

#Calculate the difference aerosolmass
array2 = sum_abs(default_file2, hyd_file2, 'PM10', useless)

# Calculate the sum of absolute differences across all data variables
sum_diff1 = np.sum(np.abs(array1))
sum_diff2 = np.sum(np.abs(array2))

# # Specify the output filename
# output_filename = 'sum_diff.txt'
# 
# # Save the sum_diff to the text file
# save_sum_diff_to_file(sum_diff, output_filename)

#Print the absolute difference to the command window
print('sum difference for species conc')
print(sum_diff1)
print('sum difference for aerosol mass')
print(sum_diff2)

main()