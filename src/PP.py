import numpy as np
from astropy import constants as c


def pebble_predictor(snowline_i=np.zeros(1), **pars):

    """
    Returns estimates of pebble Stokes number and flux from disk.

    **This has been modifed by S. McCloat to work directly with the
    PPOLs Model code.** Main modification incorporates a snow line and
    changes to dust surface density, dust/gas fraction. Also introduces
    tracking of eta and Hgas over time for gas evolution, but not yet
    implemented.

    User inputs disk properties and PP returns estimates of pebble
    Stokes number St(t,r) and pebble flux Mdot_p(t,r) [g/s] at every
    time (tgrid) and location (rgrid).

    Parameters
    ----------
    rgrid : ndarray
        Radial grid positions for disk. [cm]
    tgrid : ndarray
        Time grid for "lifetime" of disk model. [s]
    Mstar : float
        Mass of the central star. [g]
    SigmaGas, SigmaDust : ndarray
        Initial gas, dust surface density at every point of rgrid. Must
        be same size as rgrid. [g/cm^2]
    T : ndarray
        Gas temperature at every point of rgrid. [K]
    alpha : float
        Turbulence strength parameter.
    vfrag : float
        Collisional fragmentation threshold velocity. [cm/s]
    rhop : float
        Bulk density of dust grains. [g/cm^3]

    Returns
    -------
    st : ndarray
        Flux-averaged Stokes value of pebbles at each [t,r].
    flux : ndarray
        Estimated pebble mass flux (i.e. Mdot_p[t,r]).
    etaI : ndarray
        Initial and constant gas pressure gradient. Negative values.
    HgasI : ndarray
        Initial and constant gas scale height. [cm]
    fDG : ndarray
        Dust/Gas fraction, averaged across rgrid for each t.
    hgas_temptest, etas_temptest : ndarray
        Experimental arrays that record eta, hgas, over [t,r] instead of
        just [r]. Intended for future implementation of gas evolution.

    Notes
    -----
    The base and fundamental action of pebble_predictor (PP) comes from
    Drazkowska et al. (2021), and interested readers are directed there
    for full technical details. PP is modified here to incorporate a
    snow line and related changes, and return more parameters.
    """

    rgrid = pars['rgrid']
    tgrid = pars['tgrid']
    Mstar = pars['Mstar']
    SigmaGas = pars['SigmaGas']
    SigmaDust = pars['SigmaDust']
    T = pars['T']
    alpha = pars['alpha']
    vfrag = pars['vfrag']
    rhop = pars['rhop']

    k_b = c.k_B.cgs.value
    m_p = c.m_p.cgs.value
    Grav = c.G.cgs.value
    AH2 = 2.e-15  # This is unused?
    au = c.au.cgs.value

    # Prepare the matrices to fill:
    st = np.zeros((np.size(tgrid), np.size(rgrid)))
    flux = np.zeros((np.size(tgrid), np.size(rgrid)))

    # Solid mass outside of cell, for flux prediction:
    massout = np.zeros((np.size(tgrid), np.size(rgrid)))

    # SM: record radial drift velocities, unused
    rad_vs = np.zeros((np.size(tgrid), np.size(rgrid)))
    fDG = np.zeros((np.size(tgrid)))
    hgas_temptest = np.zeros_like(st)
    etas_temptest = np.zeros_like(st)

    # a0 = starting size of dust [cm]
    a0 = 1.e-4

    # Factor needed for growth faster than drift (Okuzumi et al. )
    fff = 30.

    # cs = sound speed
    # OmegaK = Keplerian frequency
    # rhog = midplane gas density
    # pg = gas pressure
    # Hgas = gas scale height, NOT aspect ratio
    # HgasI = redundant for recording

    Ti = T[0]
    cs = np.sqrt(k_b*Ti / (2.3*m_p))
    OmegaK = np.sqrt(Grav * Mstar / rgrid**3.)
    rhog = SigmaGas*OmegaK / (np.sqrt(2.*np.pi)*cs)
    pg = rhog * cs**2.
    Hgas = cs / OmegaK
    HgasI = Hgas
    hgas_temptest[0] = HgasI

    # Grid cell interfaces estimate (needed for pressure gradient
    # calculation):
    rInt = np.zeros(np.size(rgrid) + 1)
    rInt[0] = 1.5*rgrid[0] - 0.5*rgrid[1]
    rInt[1:-1] = 0.5*(rgrid[1:] + rgrid[:-1])
    rInt[-1] = 1.5*rgrid[-1] - 0.5*rgrid[-2]
    dr = rInt[1:] - rInt[:-1]

    # Gas pressure at the interfaces, pgint:
    pgInt = np.interp(rInt, rgrid, pg)
    eta = (pgInt[1:]-pgInt[:-1]) / dr[:] / (2.*rhog*rgrid*OmegaK**2.)
    etaI = eta
    etas_temptest[0] = etaI

    # SM: lines above set up cell boundaries in rgrid. rInt are the
    # "cell walls", dr is the physical size, between each rgrid cell.
    # For example, dr[0] is the first cell and stretches for 0.00039 AU.

    # Various Stoke number calculations, Eq#s from Drazkowska et al.
    # stfrag = turbulence-induced fragmentation regime, Eq 9
    # stdf = drift-induced fragmentaiton regime, Eq 10
    # st0 = micron-sized monomers, Eq 8?
    stfrag = 0.37*vfrag**2. / (3.*alpha*cs**2.)
    stdf = 0.37*vfrag / (2.*abs(eta)*OmegaK*rgrid)
    st0 = 0.5*np.pi*a0*rhop / SigmaGas

    # Initial condition:
    # SM: This next if-statment properly reduces the SigmaDust with
    # evolving snowline. Otherwise, the updating flux/fDG at bottom will
    # double-cut the fDG into the first time step. Weird but thats it.
    if snowline_i.any() != 0:
        cutsigdust = np.concatenate((SigmaDust[:snowline_i[0]] * 0.5,
                                     SigmaDust[snowline_i[0]:]))
        Z0 = cutsigdust / SigmaGas
    else:
        Z0 = SigmaDust / SigmaGas

    fDG[0] = np.mean(Z0)
    # Growth timescale, tgrowth:
    tgrowth = 1. / ((alpha/1.e-4)**0.3333 * Z0 * OmegaK
                    * (rgrid/au)**(-0.3333))
    for ir in range(0, np.size(rgrid)):
        massout[0, ir] = np.sum(2. * np.pi * rgrid[ir:] * dr[ir:]
                                * SigmaDust[ir:])
    st[0, :] = st0

    # vv = pebble drift velocity
    vv = 2.*abs(eta)*OmegaK*rgrid*st[0, :] / (1.+st[0, :]**2.)
    flux[0, :] = 2 * np.pi * rgrid[:] * vv[:] * SigmaDust[:]

    # Time "integration":
    for it in range(1, np.size(tgrid)):
        # Section between ###/### is the second time-iteration to evolve
        # temperature profile, then recalculate necessary parameters for
        # subsequent pebble formulae.

        #####
        cs = np.sqrt(k_b*T[it] / (2.3*m_p))
        OmegaK = np.sqrt(Grav*Mstar / rgrid**3.)
        rhog = SigmaGas*OmegaK / (np.sqrt(2.*np.pi)*cs)
        pg = rhog * cs**2.

        Hgas = cs / OmegaK
        hgas_temptest[it] = Hgas

        # pgINt and eta uses same rInt grid from above.
        pgInt = np.interp(rInt, rgrid, pg)
        eta = (pgInt[1:]-pgInt[:-1]) / dr[:] / (2.*rhog*rgrid*OmegaK**2.)
        etas_temptest[it] = eta

        # Same Stokes calculations as above.
        stfrag = 0.37*vfrag**2. / (3.*alpha*cs**2.)
        stdf = 0.37*vfrag / (2.*abs(eta)*OmegaK*rgrid)
        st0 = 0.5*np.pi*a0*rhop / SigmaGas
        #####

        # Decrease the mass budget:
        massout[it, :] = (massout[it-1, :] - flux[it-1, :]
                          * (tgrid[it] - tgrid[it-1]))

        # Sigma_d/Sigma_g estimate:
        Z = Z0 * (massout[it, :]/massout[0, :])
        stini = st0 * np.exp(tgrid[it]/tgrowth[:])

        # Compares St values and chooses in following order:
        # stdrift = in the drift dominated regime
        # initial growth vs fragmentation;
        # vs radial drift
        # vs St doesn't fall below its initial value
        stdrift = Z / abs(eta) / fff
        st[it, :] = np.minimum(np.minimum(stfrag[it], stdf[it]), stini)
        st[it, :] = np.minimum(st[it, :], stdrift)
        st[it, :] = np.maximum(st[it, :], st0[:])
        vv = 2.*abs(eta)*OmegaK*rgrid*st[it, :] / (1.+st[it, :]**2.)

        # Flux restriction if tgrowth is too long:
        vv = np.minimum(vv[:], rgrid[:] / tgrowth[:] / fff)
        rad_vs[it, :] = vv

        # Cuts SigmaDust in half as SL evolves. Alternate to single cut
        # at first SigmaDust.
        if snowline_i.any() != 0:
            x = np.concatenate((SigmaDust[:snowline_i[it]] * 0.5,
                                SigmaDust[snowline_i[it]:]))
            # Pebble flux corrected by the remaining mass estimate:
            flux[it, :] = (2 * np.pi * rgrid[:] * vv[:] * x
                           * (massout[it, :]/massout[0, :]))
            fDG[it] = np.mean((x * (massout[it, :]/massout[0, :]))
                              / SigmaGas[:])
        else:
            flux[it, :] = (2*np.pi*rgrid[:]*vv[:]*SigmaDust[:]
                           * (massout[it, :]/massout[0, :]))
            fDG[it] = np.mean((SigmaDust[:] * (massout[it, :]/massout[0, :]))
                              / SigmaGas[:])
    Q = cs * OmegaK / (np.pi * Grav * SigmaGas)
    mdotgas = 3 * np.pi * alpha * SigmaGas * Hgas**2 * OmegaK

    return st, flux, etaI, HgasI, fDG, hgas_temptest, etas_temptest, Q, rgrid
