"""Reads output from the CARMA test carma_solubletest.F90.

Version History:
Parker Case (10/11/2018) First crack.
"""
import h5py
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# Open netcdf file
nc_file = netCDF4.Dataset('carma_solubletest.nc', 'r')
#print(nc_file.variables)

# Open relevant variables
time = nc_file.variables['time'][:] # (days) years since 0000-01-01 00:00:00
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
bind_names = ['NONED01', 'NONED02', 'NONED03', 'NONED04', \
            'NONED05', 'NONED06', 'NONED07', 'NONED08', \
            'NONED09', 'NONED10', 'NONED11', 'NONED12',
            'NONED13', 'NONED14', 'NONED15', 'NONED16']
fitz_bin_names = ['FITZ01', 'FITZ02', 'FITZ03', 'FITZ04', \
            'FITZ05', 'FITZ06', 'FITZ07', 'FITZ08', \
            'FITZ09', 'FITZ10', 'FITZ11', 'FITZ12',
            'FITZ13', 'FITZ14', 'FITZ15', 'FITZ16']
fitz_bind_names = ['FITZD01', 'FITZD02', 'FITZD03', 'FITZD04', \
            'FITZD05', 'FITZD06', 'FITZD07', 'FITZD08', \
            'FITZD09', 'FITZD10', 'FITZD11', 'FITZD12',
            'FITZD13', 'FITZD14', 'FITZD15', 'FITZD16']
gerb_bin_names = ['GERB01', 'GERB02', 'GERB03', 'GERB04', \
            'GERB05', 'GERB06', 'GERB07', 'GERB08', \
            'GERB09', 'GERB10', 'GERB11', 'GERB12',
            'GERB13', 'GERB14', 'GERB15', 'GERB16']
gerb_bind_names = ['GERBD01', 'GERBD02', 'GERBD03', 'GERBD04', \
            'GERBD05', 'GERBD06', 'GERBD07', 'GERBD08', \
            'GERBD09', 'GERBD10', 'GERBD11', 'GERBD12',
            'GERBD13', 'GERBD14', 'GERBD15', 'GERBD16']

# Create numpy array for number density
#bins_nd = np.zeros((len(bin_names), len(time), len(height))) # (#/cm3) number density
#fitz_bins_nd = np.zeros((len(bin_names), len(time), len(height))) # (#/cm3) number density
#gerb_bins_nd = np.zeros((len(bin_names), len(time), len(height))) # (#/cm3) number density

bins_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr
binsd_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr

fitz_bins_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr
fitz_binsd_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr

gerb_bins_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr
gerb_binsd_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr

# Loop over variables to fill number densities
for i, bin_name in zip(range(len(bin_names)), bin_names):
    bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]
for i, bin_name in zip(range(len(bind_names)), bind_names):
    binsd_mmr[i,:,:] = nc_file.variables[bin_name][:,:]

for i, bin_name in zip(range(len(bin_names)), fitz_bin_names):
    fitz_bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]
for i, bin_name in zip(range(len(bin_names)), fitz_bind_names):
    fitz_binsd_mmr[i,:,:] = nc_file.variables[bin_name][:,:]

for i, bin_name in zip(range(len(bin_names)), gerb_bin_names):
    gerb_bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]
for i, bin_name in zip(range(len(bin_names)), gerb_bind_names):
    gerb_binsd_mmr[i,:,:] = nc_file.variables[bin_name][:,:]

# Summing & Averaging
mmr_profile = np.sum(bins_mmr, axis=0)
mmrd_profile = np.sum(binsd_mmr, axis=0)

fitz_mmr_profile = np.sum(fitz_bins_mmr, axis=0)
fitz_mmrd_profile = np.sum(fitz_binsd_mmr, axis=0)

gerb_mmr_profile = np.sum(gerb_bins_mmr, axis=0)
gerb_mmrd_profile = np.sum(gerb_binsd_mmr, axis=0)

colors = ['#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58']
fig, axs = plt.subplots(2, 3, figsize=(10,10))
axs[0,0].semilogx(mmr_profile[0,:], height, color=colors[0], linestyle='--', label='None Salt')
axs[1,0].semilogx(mmrd_profile[0,:], height, color=colors[0], label='None Dust')

axs[0,0].semilogx(fitz_mmr_profile[0,:], height, color=colors[2], linestyle='--', label='Fitz Salt')
axs[1,0].semilogx(fitz_mmrd_profile[0,:], height, color=colors[2], label='Fitz Dust')

axs[0,0].semilogx(gerb_mmr_profile[0,:], height, color=colors[-1], linestyle='--', label='Gerb Salt')
axs[1,0].semilogx(gerb_mmrd_profile[0,:], height, color=colors[-1], label='Gerb Dust')

axs[0,1].semilogx(mmr_profile[-200,:], height, color=colors[0], linestyle='--', label='None Salt')
axs[1,1].semilogx(mmrd_profile[-200,:], height, color=colors[0], label='None Dust')

axs[0,1].semilogx(fitz_mmr_profile[-200,:], height, color=colors[2], linestyle='--', label='Fitz Salt')
axs[1,1].semilogx(fitz_mmrd_profile[-200,:], height, color=colors[2], label='Fitz Dust')

axs[0,1].semilogx(gerb_mmr_profile[-200,:], height, color=colors[-1], linestyle='--', label='Gerb Salt')
axs[1,1].semilogx(gerb_mmrd_profile[-200,:], height, color=colors[-1], label='Gerb Dust')

axs[0,2].semilogx(mmr_profile[-1,:], height, color=colors[0], linestyle='--', label='None Salt')
axs[1,2].semilogx(mmrd_profile[-1,:], height, color=colors[0], label='None Dust')

axs[0,2].semilogx(fitz_mmr_profile[-1,:], height, color=colors[2], linestyle='--', label='Fitz Salt')
axs[1,2].semilogx(fitz_mmrd_profile[-1,:], height, color=colors[2], label='Fitz Dust')

axs[0,2].semilogx(gerb_mmr_profile[-1,:], height, color=colors[-1], linestyle='--', label='Gerb Salt')
axs[1,2].semilogx(gerb_mmrd_profile[-1,:], height, color=colors[-1], label='Gerb Dust')

axs[0,0].invert_yaxis()
axs[0,1].invert_yaxis()
axs[0,2].invert_yaxis()
axs[1,0].invert_yaxis()
axs[1,1].invert_yaxis()
axs[1,2].invert_yaxis()
axs[0,0].set_ylabel('hPa')
axs[0,1].set_ylabel('hPa')
axs[0,2].set_ylabel('hPa')
axs[1,0].set_ylabel('hPa')
axs[1,1].set_ylabel('hPa')
axs[1,2].set_ylabel('hPa')
axs[0,0].set_xlabel('Aerosol MMR (kg/kg)')
axs[0,1].set_xlabel('Aerosol MMR (kg/kg)')
axs[0,2].set_xlabel('Aerosol MMR (kg/kg)')
axs[1,0].set_xlabel('Aerosol MMR (kg/kg)')
axs[1,1].set_xlabel('Aerosol MMR (kg/kg)')
axs[1,2].set_xlabel('Aerosol MMR (kg/kg)')
axs[0,0].set_title('t=' + str(round(time[0])) + ' days')
axs[0,1].set_title('t=' + str(round(time[-200])) + ' days')
axs[0,2].set_title('t=' + str(round(time[-1])) + ' days')
axs[1,0].set_title('t=' + str(round(time[0])) + ' days')
axs[1,1].set_title('t=' + str(round(time[-200])) + ' days')
axs[1,2].set_title('t=' + str(round(time[-1])) + ' days')
axs[0,0].legend()
axs[0,1].legend()
axs[0,2].legend()
axs[1,0].legend()
axs[1,1].legend()
axs[1,2].legend()
axs[0,0].grid()
axs[0,1].grid()
axs[0,2].grid()
axs[1,0].grid()
axs[1,1].grid()
axs[1,2].grid()
axs[0,0].set_xlim(1e-12, 5e-9)
axs[0,1].set_xlim(1e-12, 5e-9)
axs[0,2].set_xlim(1e-12, 5e-9)
axs[1,0].set_xlim(1e-12, 5e-9)
axs[1,1].set_xlim(1e-12, 5e-9)
axs[1,2].set_xlim(1e-12, 5e-9)
plt.tight_layout()
plt.savefig('read_solubletest.png')
plt.show()
