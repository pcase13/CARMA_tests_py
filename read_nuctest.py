"""Reads output from the CARMA test carma_nuctest.F90.

Note: There is a version of carma_nuctest.F90 that will cause this script to
fail. When writing output for the drymass and dryr variables for each group,
the loop fails to update 'sname', causing the wrong variables to be written.
Make sure that CARMAGROUP_Get retrieves 'sname' in this loop.

Version History:
Parker Case (10/11/2018) First crack.
"""
import h5py
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# Open netcdf file
nc_file = netCDF4.Dataset('carma_nuctest.nc', 'r')

# Open relevant variables
time = nc_file.variables['time'][:]*24*60*60 # (seconds) time since 0000-01-01 00:00:00
in_radii = nc_file.variables['CRINdryr'][:] # (m) Bin center radius
ice_radii = nc_file.variables['CRICEdryr'][:] # (m) Bin center radius
total_in_mmr = nc_file.variables['CRINMR'][:] # (g) Total mass mixing ratio across bins
total_ice_mmr = nc_file.variables['CRICEMR'][:] # (g) Total mass mixing ratio across bins
total_core_mmr = nc_file.variables['CRCOREMR'][:] # (g) Total mass mixing ratio across bins
h2o_mmr = nc_file.variables['H2O'][:] # (g) Total mass mixing ratio across bins
satice = nc_file.variables['H2Osatice'][:] # Saturation ratio over ice
satliq = nc_file.variables['H2Osatliq'][:] # Saturation ratio over ice

# Create lists of the bin names
in_element_name = 'CRIN'
ice_element_name = 'CRICE'
core_element_name = 'CRCORE'
bin_numbers = ['01','02','03','04','05','06','07','08','09','10','11','12', \
              '13','14','15','16']
in_bin_names = [in_element_name + bin_number for bin_number in bin_numbers]
ice_bin_names = [ice_element_name + bin_number for bin_number in bin_numbers]
core_bin_names = [core_element_name + bin_number for bin_number in bin_numbers]

# Create numpy array for number density
bins_mmr = np.zeros((3, len(bin_numbers), len(time))) # (kg/kg) mmr
# Loop over variables to fill number densities
for j, bin_names in zip(range(3), [in_bin_names, ice_bin_names, core_bin_names]):
    for i, bin_name in zip(range(len(bin_names)), bin_names):
        bins_mmr[j,i,:] = nc_file.variables[bin_name][:]

# Summing & averaging...
# None here!

# Create figure, axes
fig, axs = plt.subplots(4,1,figsize=(7,8))

# For 4 timesteps evenly spaced through simulation, plot mass column
colors = ['#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58']

# Plot the sulfate ice nuclei mmr distribution
# Past timesteps in grey, last timestep in blue
for i in range(len(time)-1):
    axs[0].loglog(in_radii*10e7, bins_mmr[0,:,i], color='#999999')
axs[0].loglog(in_radii*10e7, bins_mmr[0,:,-1], color=colors[-1], label='Sulfate')
axs[0].set_ylim(1e-20, 1e-10)
axs[0].set_ylabel('mmr ($kg/kg$)')
axs[0].set_xlabel('Radius ($nm$)')
axs[0].legend()
axs[0].grid()

# Plot the ice and core mmr distribution
# Past timesteps in grey, last timestep in blue
for i in range(len(time)-1):
    axs[1].loglog(ice_radii*10e7, bins_mmr[1,:,i], color='#999999')
    axs[1].loglog(ice_radii*10e7, bins_mmr[2,:,i], color='#CCCCCC')
axs[1].loglog(ice_radii*10e7, bins_mmr[1,:,-1], color=colors[-1], label='Ice')
axs[1].loglog(ice_radii*10e7, bins_mmr[2,:,-1], color=colors[-4], label='Sulfate Core')
axs[1].set_ylim(1e-18, 1e-3)
axs[1].set_ylabel('mmr ($kg/kg$)')
axs[1].set_xlabel('Radius ($nm$)')
axs[1].legend()
axs[1].grid()

# Plot the integrated mmr for ice and vapor
axs[2].plot(time, total_ice_mmr, color=colors[1], label='Ice MMR')
axs[2].plot(time, h2o_mmr, color=colors[3], label='H2O Gas MMR')
axs[2].plot(time, h2o_mmr + total_ice_mmr, color=colors[5], label='Ice + H2O Gas MMR')
axs[2].legend(loc=4)
axs[2].set_ylabel('mmr ($kg/kg$)')
axs[2].set_xlabel('Time (seconds)')
axs[2].grid()

# Plot the timeseries of ice & liquid saturation ratio
axs[3].axhline(1, color='black')
axs[3].plot(time, satice, color=colors[1], label='Saturation Ratio (ice)')
axs[3].plot(time, satliq, color=colors[5], label='Saturation Ratio (water)')
axs[3].set_ylabel('s')
axs[3].set_xlabel('Time (seconds)')
axs[3].set_ylim(0.,2.)
axs[3].legend()
axs[3].grid()

plt.tight_layout()
plt.savefig('read_nuctest.png')
plt.show()
