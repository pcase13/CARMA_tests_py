"""Reads output from the CARMA test carma_swelltest.F90.

Version History:
Parker Case (10/11/2018) First crack.
"""
import h5py
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# Open netcdf file
nc_file = netCDF4.Dataset('carma_swelltest.nc', 'r')
#print(nc_file.variables)

# Open relevant variables
time = nc_file.variables['time'][:]/365 # (days) years since 0000-01-01 00:00:00
height = nc_file.variables['lev'][:]/100 # (kPa) Layer midpoint heights
bins = nc_file.variables['NONEdrymass'][:] # (g) Bin center masses
fitz_bins = nc_file.variables['NONEdrymass'][:] # (g) Bin center masses
gerb_bins = nc_file.variables['NONEdrymass'][:] # (g) Bin center masses
nd_profile = nc_file.variables['NONEND'][:]
fitz_nd_profile = nc_file.variables['FITZND'][:]
gerb_nd_profile = nc_file.variables['GERBND'][:]
#dz = np.zeros(len(height))+1000 # (m) Height of layers
#dz_cm = dz * 100 # (cm) height of layers in cm

# Create a list of the bin names
bin_names = ['NONE01', 'NONE02', 'NONE03', 'NONE04', \
            'NONE05', 'NONE06', 'NONE07', 'NONE08', \
            'NONE09', 'NONE10', 'NONE11', 'NONE12',
            'NONE13', 'NONE14', 'NONE15', 'NONE16']
fitz_bin_names = ['FITZ01', 'FITZ02', 'FITZ03', 'FITZ04', \
            'FITZ05', 'FITZ06', 'FITZ07', 'FITZ08', \
            'FITZ09', 'FITZ10', 'FITZ11', 'FITZ12',
            'FITZ13', 'FITZ14', 'FITZ15', 'FITZ16']
gerb_bin_names = ['GERB01', 'GERB02', 'GERB03', 'GERB04', \
            'GERB05', 'GERB06', 'GERB07', 'GERB08', \
            'GERB09', 'GERB10', 'GERB11', 'GERB12',
            'GERB13', 'GERB14', 'GERB15', 'GERB16']

# Create numpy array for number density
bins_nd = np.zeros((len(bin_names), len(time), len(height))) # (#/cm3) number density
fitz_bins_nd = np.zeros((len(bin_names), len(time), len(height))) # (#/cm3) number density
gerb_bins_nd = np.zeros((len(bin_names), len(time), len(height))) # (#/cm3) number density
bins_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr
fitz_bins_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr
gerb_bins_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr
# Loop over variables to fill number densities
for i, bin_name in zip(range(len(bin_names)), bin_names):
    bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]
for i, bin_name in zip(range(len(bin_names)), fitz_bin_names):
    fitz_bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]
for i, bin_name in zip(range(len(bin_names)), gerb_bin_names):
    gerb_bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]

# Summing & Averaging
mmr_profile = np.sum(bins_mmr, axis=0)
fitz_mmr_profile = np.sum(fitz_bins_mmr, axis=0)
gerb_mmr_profile = np.sum(gerb_bins_mmr, axis=0)

fig, axs = plt.subplots(1, 3, figsize=(10,5))
axs[0].semilogx(mmr_profile[0,:], height, label='None')
axs[0].semilogx(fitz_mmr_profile[0,:], height, label='Fitz')
axs[0].semilogx(gerb_mmr_profile[0,:], height, label='Gerb')
axs[1].semilogx(mmr_profile[-100,:], height, label='None')
axs[1].semilogx(fitz_mmr_profile[-100,:], height, label='Fitz')
axs[1].semilogx(gerb_mmr_profile[-100,:], height, label='Gerber')
axs[2].semilogx(mmr_profile[-1,:], height, label='None')
axs[2].semilogx(fitz_mmr_profile[-1,:], height, label='Fitz')
axs[2].semilogx(gerb_mmr_profile[-1,:], height, label='Gerber')
axs[0].invert_yaxis()
axs[1].invert_yaxis()
axs[2].invert_yaxis()
axs[0].set_ylabel('hPa')
axs[1].set_ylabel('hPa')
axs[2].set_ylabel('hPa')
axs[0].set_xlabel('Aerosol MMR (kg/kg)')
axs[1].set_xlabel('Aerosol MMR (kg/kg)')
axs[2].set_xlabel('Aerosol MMR (kg/kg)')
axs[0].set_title('t=0')
axs[1].set_title('t=something')
axs[2].set_title('t=final')
axs[0].legend()
axs[1].legend()
axs[2].legend()
axs[0].grid()
axs[1].grid()
axs[2].grid()
axs[0].set_xlim(1e-12, 5e-9)
axs[1].set_xlim(1e-12, 5e-9)
axs[2].set_xlim(1e-12, 5e-9)
plt.tight_layout()
plt.savefig('read_swelltest.png')
plt.show()
