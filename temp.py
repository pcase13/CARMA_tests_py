"""Reads output from the CARMA test carma_coagtest.F90.

Version History:
Parker Case (10/11/2018) First crack.
"""
import h5py
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# Open netcdf file
nc_file = netCDF4.Dataset('carma_coagtest.nc', 'r')

# Open relevant variables
time = nc_file.variables['time'][:]*24*60*60 # (seconds) time since 0000-01-01 00:00:00
height = nc_file.variables['Z'][0,:] # (m) Layer midpoint heights
bins = nc_file.variables['DUSTdrymass'][:] # (g) Bin center masses
dz = np.zeros(len(height))+100 # (m) Height of layers
dz_cm = dz * 100 # (cm) height of layers in cm

# Create a list of the bin names
bin_names = ['DUST01', 'DUST02', 'DUST03', 'DUST04', \
            'DUST05', 'DUST06', 'DUST07', 'DUST08', \
            'DUST09', 'DUST10', 'DUST11', 'DUST12', \
            'DUST13', 'DUST14', 'DUST15', 'DUST16', \
            'DUST17', 'DUST18', 'DUST19', 'DUST20']
for i in range(21,77):
    bin_names.append('DUST' + str(i))

# Create numpy array for number density
bins_nd = np.zeros((len(bin_names), len(time), len(height))) # (#/cm3) number density
bins_mmr = np.zeros((len(bin_names), len(time), len(height))) # (kg/kg) mmr
# Loop over variables to fill number densities
for i, bin_name in zip(range(len(bin_names)), bin_names):
    bins_nd[i,:,:] = nc_file.variables[bin_name+'nd'][:,:]
    bins_mmr[i,:,:] = nc_file.variables[bin_name][:,:]

# Sum all vertical levels
bins_nd = np.nansum(bins_nd, axis=2)
bins_mmr = np.nansum(bins_mmr, axis=2)
bins_mass = np.swapaxes(np.swapaxes(bins_nd, 0, 1) * bins, 0, 1)

for i in np.linspace(0,len(time)-1,4):
    i = int(i)
    plt.semilogx(bins, bins_mass[:,i], label=str(int(time[i])) + ' seconds')
plt.legend()
plt.xlabel('Bin Mass (g)')
plt.ylabel('Mass in Bin (g)')
plt.grid()
plt.tight_layout()
plt.savefig('read_coagtest.png')
plt.show()

