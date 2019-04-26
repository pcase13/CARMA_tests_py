"""Reads output from the CARMA test carma_coagtest.F90.

Version History:
Parker Case (10/11/2018) First crack.
"""
import h5py
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import carmabins

nbins = 24
rmrat = 3.7515201
rmin = 2.6686863e-10
rhop = 1923. # kg m-3

# Set up variable names and carma bins
rmass, rmassup, r, rup, dr, rlow, masspart = carmabins.carmabins(nbins, rmrat, rmin, rhop)

# Open netcdf file
nc_file = netCDF4.Dataset('carma_ststest.nc', 'r')
nc_file_ctl = netCDF4.Dataset('carma_sulfatetest.nc', 'r')

# Open relevant variables
time = nc_file.variables['time'][:]*24*60*60 # (seconds) time since 0000-01-01 00:00:00
height = nc_file.variables['Z'][:] # (m) Layer midpoint heights
bins = nc_file.variables['SULFdrymass'][:] # (g) Bin center masses
bins_dryr = nc_file.variables['SULFdryr'][:] * .01 # (m) Bin center radius
satliq = nc_file.variables['H2Osatliq'][:]
dz = np.zeros(len(height))+100 # (m) Height of layers
dz_cm = dz * 100 # (cm) height of layers in cm

# Create a list of the bin names
bin_names = ['SULF01', 'SULF02', 'SULF03', 'SULF04', \
            'SULF05', 'SULF06', 'SULF07', 'SULF08', \
            'SULF09', 'SULF10', 'SULF11', 'SULF12', \
            'SULF13', 'SULF14', 'SULF15', 'SULF16', \
            'SULF17', 'SULF18', 'SULF19', 'SULF20', \
            'SULF21', 'SULF22', 'SULF23', 'SULF24']

# Create numpy array for number density
bins_nd = np.zeros((len(bin_names), len(time))) # (#/cm3) number density
bins_mmr = np.zeros((len(bin_names), len(time))) # (kg/kg) mmr
bins_rwet = np.zeros((len(bin_names), len(time))) # (m) mmr
bins_rhop = np.zeros((len(bin_names), len(time))) # (m) mmr

bins_nd_ctl = np.zeros((len(bin_names), len(time))) # (#/cm3) number density
bins_mmr_ctl = np.zeros((len(bin_names), len(time))) # (kg/kg) mmr
bins_rwet_ctl = np.zeros((len(bin_names), len(time))) # (m) mmr
bins_rhop_ctl = np.zeros((len(bin_names), len(time))) # (m) mmr
# Loop over variables to fill number densities
for i, bin_name in zip(range(len(bin_names)), bin_names):
    bins_nd[i,:] = nc_file.variables[bin_name+'nd'][:]
    bins_mmr[i,:] = nc_file.variables[bin_name][:]
    bins_rwet[i,:] = nc_file.variables[bin_name+'wetr'][:] * .01
    bins_rhop[i,:] = nc_file.variables[bin_name+'rhopwet'][:]

    bins_nd_ctl[i,:] = nc_file_ctl.variables[bin_name+'nd'][:]
    bins_mmr_ctl[i,:] = nc_file_ctl.variables[bin_name][:]
    bins_rwet_ctl[i,:] = nc_file_ctl.variables[bin_name+'wetr'][:] * .01
    bins_rhop_ctl[i,:] = nc_file_ctl.variables[bin_name+'rhopwet'][:]


cdf = np.cumsum(bins_nd[::-1,-1])[::-1]
cdf_ctl = np.cumsum(bins_nd_ctl[::-1,-1])[::-1]

wtn = np.asarray([2.1478113995874924E-009,
    3.6004254127917136E-003,
    1.3142047479136516,
    1.8947429993998743E-002,
    3.5522966725688200E-004,
    9.6259758412181665E-003,
    0.11944385011600316,
    0.50898226643647615,
    1.2806405827483014,
    2.3507882767010284,
    3.5279210909898882,
    4.6361053110177188,
    5.5706721325479140,
    6.2971634389683926,
    6.8299929966664070,
    7.2040789493965010,
    7.4600202627493921,
    7.6309225792906865,
    7.7440168739463875,
    7.8173383421983313,
    7.8661731172420737,
    7.8967364451046054,
    7.9170992399198825,
    7.9296958677713203])
wts = np.asarray([70.183209524083779,
    70.168154563169793,
    68.669305053447147,
    68.387198299494798,
    64.323052503873541,
    58.626148833302544,
    53.556412252877536,
    49.911724008473598,
    47.244386999925602,
    45.083599016915379,
    43.263885564675938,
    41.769932516428177,
    40.601648540121580,
    39.732571488233610,
    39.112092993605700,
    38.683575135511852,
    38.393607470803637,
    38.201201581167680,
    38.074482097922726,
    37.992406490258773,
    37.938082604126571,
    37.903895680389653,
    37.881265627768705,
    37.867159210784024])
hno3_kelv = np.asarray([1.4979227274800483E-013,
    1.7844497941020609E-011,
    3.8687506771473646E-010,
    2.8017316162753846E-009,
    1.0018580509270943E-008,
    2.2748156286037150E-008,
    3.8560989616063768E-008,
    5.4157210495917363E-008,
    6.7388916646450723E-008,
    7.7568292873847533E-008,
    8.4918789070923323E-008,
    9.0013735513076278E-008,
    9.3453244952775701E-008,
    9.5736016536227505E-008,
    9.7234551820705783E-008,
    9.8211350706805548E-008,
    9.8845175320336937E-008,
    9.9255249227171260E-008,
    9.9520060889012138E-008,
    9.9690859972023558E-008,
    9.9800936587944972E-008,
    9.9871843160321495E-008,
    9.9917503331908198E-008,
    9.9946900009062230E-008])
h2o_kelv = np.asarray([1.4018226795180163E-028,
    2.6869709209618734E-022,
    2.9686094470651519E-018,
    1.1876923419129425E-015,
    5.6153532227040240E-014,
    6.7166216732143255E-013,
    3.3172379307428360E-012,
    9.2720373026834634E-012,
    1.7966565406873986E-011,
    2.7501537741190429E-011,
    3.6169987302377820E-011,
    4.3144749521806824E-011,
    4.8329520326466769E-011,
    5.1991298198694805E-011,
    5.4493363235865989E-011,
    5.6166948462215124E-011,
    5.7271107655369835E-011,
    5.7993167593563434E-011,
    5.8462673113980987E-011,
    5.8766842394673959E-011,
    5.8963433873143608E-011,
    5.9090302451404189E-011,
    5.9172095932246239E-011,
    5.9224795848392004E-011])

fig, axs = plt.subplots(1,2,figsize=(12,6))

axs[0].loglog(bins_dryr, cdf, marker='o', color='black', label='Dry radius STS')
axs[0].loglog(bins_rwet[:,-1], cdf, marker='o', color='grey', label='Wet radius STS')

axs[0].loglog(bins_dryr, cdf_ctl, marker='o', color='blue', label='Dry radius H2O')
axs[0].loglog(bins_rwet_ctl[:,-1], cdf_ctl, marker='o', color='cyan', label='Wet radius H2O')

axs[0].set_xlabel('Bin radius ($m$)')
axs[0].set_ylabel('Cumulative number density ($cm^{-3}$)')
axs[0].set_title('Comparison with binary')
axs[0].legend()

axs[1].semilogx(bins_dryr, bins_rhop[:,-1], marker='o', color='black', label='STS')
axs[1].semilogx(bins_dryr, bins_rhop_ctl[:,-1], marker='o', color='blue', label='H2O')

axs[1].set_xlabel('Dry bin radius ($m$)')
axs[1].set_ylabel('Particle wet density ($g\ cm^{3}$)')
axs[1].set_title('Particle density')
axs[1].legend()

plt.tight_layout()
plt.show()

fig, axs = plt.subplots(1,2,figsize=(12,6))

axs[0].semilogx(bins_dryr, wts, marker='o', color='black', label='wt% H2SO4')
axs[0].semilogx(bins_dryr, wtn, marker='o', color='grey', label='wt% HNO3')

axs[0].set_ylim(0,100)
axs[0].set_xlabel('Dry bin radius ($m$)')
axs[0].set_ylabel('Weight percent ($\%$)')
axs[0].set_title('Particle composition')
axs[0].legend()

axs[1].loglog(bins_dryr, h2o_kelv, marker='o', color='black', label='Water')
axs[1].loglog(bins_dryr, hno3_kelv, marker='o', color='grey', label='Nitric acid')

axs[1].set_xlabel('Dry bin radius ($m$)')
axs[1].set_ylabel('Vapor available ($kg\ kg^{-1}$)')
axs[1].set_title('Kelvin effect')
axs[1].legend()

plt.tight_layout()
plt.show()

fig, axs = plt.subplots(1,2,figsize=(12,6))

vol = 4./3. * np.pi * bins_rwet[:,-1] ** 3 * bins_rhop[:,-1]
vol_ctl = 4./3. * np.pi * bins_rwet_ctl[:,-1] ** 3 * bins_rhop_ctl[:,-1]
axs[0].loglog(bins_rwet[:,-1], vol, marker='o', color='black', label='STS')
axs[0].loglog(bins_rwet_ctl[:,-1], vol_ctl, marker='o', color='blue', label='H2SO4/H2O')
axs[0].set_xlabel('Bin radius ($m$)')
axs[0].set_ylabel('Bin mass ($g\ cm^{-3}$)')
axs[0].legend()

axs[1].semilogx(bins_dryr, bins_rhop[:,-1], marker='o', color='black', label='STS')
axs[1].semilogx(bins_dryr, bins_rhop_ctl[:,-1], marker='o', color='blue', label='H2O')
axs[1].set_xlabel('Dry bin radius ($m$)')
axs[1].set_ylabel('Particle wet density ($g\ cm^{3}$)')
axs[1].set_title('Particle density')
axs[1].legend()

plt.tight_layout()
plt.show()

'''
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
'''

