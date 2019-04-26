"""Module to create CARMA like bins for size distribution analysis
Translated from Dr. Peter Colarco's IDL program

Author: Parker A. Case
Changelog: Created 6/21/2017
"""
import numpy as np

def carmabins(nbin, rmrat, rmin, rhop):
    """Create CARMA like size distribution bins

    Arguments:
    nbin -- number of bins
    rmrat -- mass ratio
    rmin -- minimum radius
    rhop -- density of particles
    """
    cpi = 4./3. * np.pi
    # Convert minimum radius to mass
    rmassmin = cpi*rhop*rmin**3.
    # Volume to radius factor
    vrfact = ((3./2. / np.pi / (rmrat+1))**(1./3.))*(rmrat**(1./3.) - 1.)

    rmass = np.zeros(nbin)
    rmassup = np.zeros(nbin)
    r = np.zeros(nbin)
    rup = np.zeros(nbin)
    dr = np.zeros(nbin)
    rlow = np.zeros(nbin)

    for ibin in range(0, nbin):
        rmass[ibin]   = rmassmin*rmrat**ibin # Bin median mass
        rmassup[ibin] = 2.*rmrat/(rmrat+1.)*rmass[ibin] # Bin maximum radius
        r[ibin]       = (rmass[ibin]/rhop/cpi)**(1./3.) # Bin median radius
        rup[ibin]     = (rmassup[ibin]/rhop/cpi)**(1./3.) # Bin maximum radius
        dr[ibin]      = vrfact*(rmass[ibin]/rhop)**(1./3.) # Bin width in radius
        rlow[ibin]    = rup[ibin] - dr[ibin] # Bin minimum radius

    masspart = 4./3. * np.pi * r**3. * rhop # Mass of median radius

    return rmass, rmassup, r, rup, dr, rlow, masspart
