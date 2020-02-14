"""Reads output from the CARMA test carma_swelltest.F90.

Version History:
Parker Case (10/11/2018) First crack.
"""
import h5py
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# Open netcdf file
nc_file = netCDF4.Dataset('carma_mixedgrouptest.nc', 'r')
#print(nc_file.variables)

# Open relevant variables
time = nc_file.variables['time'][:]/365 # (days) years since 0000-01-01 00:00:00
height = nc_file.variables['lev'][:]/100 # (kPa) Layer midpoint heights
pc_bins = nc_file.variables['SULFdryr'][:] # (m) bin center radii
pure_bins = nc_file.variables['PUREdryr'][:] # (m) bin center radii
#dz = np.zeros(len(height))+1000 # (m) Height of layers
#dz_cm = dz * 100 # (cm) height of layers in cm

# Create a list of the bin names
pc_bin_names = ['SULF01', 'SULF02', 'SULF03', 'SULF04', \
            'SULF05', 'SULF06', 'SULF07', 'SULF08', \
            'SULF09', 'SULF10', 'SULF11', 'SULF12', \
            'SULF13', 'SULF14', 'SULF15', 'SULF16', \
            'SULF17', 'SULF18', 'SULF19', 'SULF20', \
            'SULF21', 'SULF22', 'SULF23', 'SULF24', \
            'SULF25', 'SULF26', 'SULF27', 'SULF28']
core_bin_names = ['CORE01', 'CORE02', 'CORE03', 'CORE04', \
            'CORE05', 'CORE06', 'CORE07', 'CORE08', \
            'CORE09', 'CORE10', 'CORE11', 'CORE12', \
            'CORE13', 'CORE14', 'CORE15', 'CORE16', \
            'CORE17', 'CORE18', 'CORE19', 'CORE20', \
            'CORE21', 'CORE22', 'CORE23', 'CORE24', \
            'CORE25', 'CORE26', 'CORE27', 'CORE28']
pure_bin_names = ['PURE01', 'PURE02', 'PURE03', 'PURE04', \
            'PURE05', 'PURE06', 'PURE07', 'PURE08', \
            'PURE09', 'PURE10', 'PURE11', 'PURE12', \
            'PURE13', 'PURE14', 'PURE15', 'PURE16', \
            'PURE17', 'PURE18', 'PURE19', 'PURE20', \
            'PURE21', 'PURE22', 'PURE23', 'PURE24', \
            'PURE25', 'PURE26', 'PURE27', 'PURE28']

# Create numpy array for number density
pc_bins_nd = np.zeros((len(pc_bin_names), len(time), len(height))) # (#/cm3) number density
core_bins_nd = np.zeros((len(pc_bin_names), len(time), len(height))) # (#/cm3) number density
pure_bins_nd = np.zeros((len(pc_bin_names), len(time), len(height))) # (#/cm3) number density
pc_bins_mmr = np.zeros((len(pc_bin_names), len(time), len(height))) # (kg/kg) mmr
core_bins_mmr = np.zeros((len(pc_bin_names), len(time), len(height))) # (kg/kg) mmr
pure_bins_mmr = np.zeros((len(pc_bin_names), len(time), len(height))) # (kg/kg) mmr
# Loop over variables to fill number densities
for i, bin_name in zip(range(len(pc_bin_names)), pc_bin_names):
    pc_bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]
for i, bin_name in zip(range(len(core_bin_names)), core_bin_names):
    core_bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]
for i, bin_name in zip(range(len(pure_bin_names)), pure_bin_names):
    pure_bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]

# Summing & Averaging
sulf_mmr = pc_bins_mmr - core_bins_mmr
#mmr_profile = np.sum(bins_mmr, axis=0)
#core_mmr_profile = np.sum(core_bins_mmr, axis=0)
#pure_mmr_profile = np.sum(pure_bins_mmr, axis=0)

# Test for equality
print('pc element')
print(pc_bins_mmr[:,0,0])
print('core element')
print(core_bins_mmr[:,0,0])
print('implied sulfate')
print(sulf_mmr[:,0,0])

fig, ax = plt.subplots(1, 1, figsize=(10,10))
plt.loglog(pc_bins, sulf_mmr[:,0,0], '--', color='red', label='t=0, pc_sulfs')
plt.loglog(pc_bins, core_bins_mmr[:,0,0], '--', color='blue', label='t=0, core')
plt.loglog(pure_bins, pure_bins_mmr[:,0,0], '--', color='purple', label='t=0, pure')
plt.loglog(pc_bins, sulf_mmr[:,-1,0], color='red', label='t=500, pc_sulfs')
plt.loglog(pc_bins, core_bins_mmr[:,-1,0], color='blue', label='t=500, core')
plt.loglog(pure_bins, pure_bins_mmr[:,-1,0], color='purple', label='t=500, pure')
ax.set_ylabel('mmr (kg/kg)')
ax.set_xlabel('radius (cm)')
ax.legend()
ax.grid()
#ax.set_xlim(1e-12, 5e-9)
plt.tight_layout()
plt.savefig('read_mixedgrouptest.png')
plt.show()
