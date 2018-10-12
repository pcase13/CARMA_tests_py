"""Reads output from the CARMA test carma_falltest.F90.

Version History:
Parker Case (10/11/2018) First crack.
"""
import h5py
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# Open netcdf file
nc_file = netCDF4.Dataset('carma_falltest.nc', 'r')

# Open relevant variables
time = nc_file.variables['time'][:]*24*60*60 # (seconds) time since 0000-01-01 00:00:00
height = nc_file.variables['Height'][:] # (m) Layer midpoint heights
bins = nc_file.variables['DUSTdrymass'][:] # (g) Bin center masses
dz = np.zeros(len(height))+100 # (m) Height of layers
dz_cm = dz * 100 # (cm) height of layers in cm

# Create a list of the bin names
bin_names = ['DUST01', 'DUST02', 'DUST03', 'DUST04', \
            'DUST05', 'DUST06', 'DUST07', 'DUST08']

# Create numpy array for number density
bins_nd = np.zeros((len(bin_names), len(time), len(height))) # (#/cm3) number density
bins_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr
# Loop over variables to fill number densities
for i, bin_name in zip(range(len(bin_names)), bin_names):
    bins_nd[i,:,:] = nc_file.variables[bin_name+'nd'][:,:]
    bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]

# Create numpy array for mass
bins_mass = np.zeros((len(bin_names), len(time), len(height))) # (g/cm3) mixing ratio
bins_mass = np.swapaxes(np.swapaxes(bins_nd, 0, 2)*bins, 0, 2) # number density * bin mass

# Summing & averaging...
total_nd = np.sum(bins_nd,axis=0) # (#/cm3) Number density summed across bins
total_mass = np.sum(bins_mass, axis=0) # (g/cm3) mass summed across bins
column_mass = np.sum(total_mass*dz_cm, axis=1) # (g/cm2) mass summed across bins and column
km = height/1000 # (km) Height of layers, in km

# Create figure, axes
fig, axs = plt.subplots(2,1,figsize=(7,8))

# For 4 timesteps evenly spaced through simulation, plot mass column
colors = ['#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58']
for it, i in zip(range(0, len(time), int(len(time)/4)), range(4)):
    axs[0].plot(total_mass[it,:], km, color=colors[i], label=str(time[it]) + ' s')
axs[0].set_title('Falling History')
axs[0].set_xlabel('Mass Concentration ($g\ cm^{-3}$)')
axs[0].set_ylabel('Altitude ($km$)')
axs[0].legend(bbox_to_anchor=(1.05,1.0))
axs[0].grid()

# Plot timeseries of total mass
axs[1].plot(column_mass, color=colors[-1])
axs[1].set_title('Column Mass Evolution')
axs[1].set_xlabel('Timesteps')
axs[1].set_ylabel('Column Mass ($g\ cm^{-2}$)')
axs[1].grid()

plt.tight_layout()
plt.savefig('read_falltest.png')
plt.show()
