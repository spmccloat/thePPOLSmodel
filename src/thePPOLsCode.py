import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap, BoundaryNorm
from PP import pebble_predictor
import OL18
from astropy import constants as c


# constants in cgs units for PP:
year = 365.25*24*3600     # in s
au = c.au.cgs.value       # in cm
MS = c.M_sun.cgs.value    # in g
ME = c.M_earth.cgs.value  # in g
RE = c.R_earth.cgs.value  # in cm

# iron = 7.9 g/cc
# earth = 4.0 g/cc (uncompressed)
# silic = 3.0 g/cc
# water = 1.0 g/cc
innercomp = 5.5  # density inside SL g/cm3, 50/50 iron/earth
outercomp = 2.5  # desnity outside SL, g/cm3, 50/50 earth/water

def lumscale(msol=1.0):
    """Produces luminosity based on stellar mass, from
    https://en.wikipedia.org/wiki/Mass%E2%80%93luminosity_relation"""
    #if msol < 0.43:
        #return 0.23 * msol**2.3
    #elif msol <= 2.0:4
    return msol**3

class Disk:
    """
    Creates protoplanetary disk object and pebble parameters. ALTERED.

    Calculates and records disk properties, as well as pebble flux and
    St values based on input parameters. Disk object attributes are
    recorded at each space on AU grid, and time. Here it is possible to
    see how pebbles coagulate and drift without any protoplanets.
    Calculations related to placement and growth of seeds are handled
    via Seeds object.

    Parameters
    ----------
    name : string, default 'Unnamed'
        Name used in some plotting and inspection methods.

    alpha : float
        Unitless viscosity parameter as per Shakura & Sunyaev.
        Default is :math:`10^{-3}`.

    vfrag, vfin, vfout : float
        Fragmentation velocity threshold (cm/s) for pebbles.
        If `snowmode` is ``None``, then `vfrag` is used. Otherwise,
        `vfin` and `vfout` refer to values inside and outside the snow
        line. Defaults are 1000, 100, 1000.

    rcrit : float
        Critical radius (au) of disk gas surface density distribution,
        beyond which the surface density drops exponentially.
        Default is 30.

    snowmode : {'none', 'evol', 'temp', float}
        Method used to determine location of the snow line, and where
        dust, solid composition, and pebble properties differ.

            'none' : str
                There is no snow line.

            'evol' : str
                Sets and evolves position according to disk properties
                like the alpha, :math:`Sigma_0, gas` Sigma Gas, and
                decreasing dust/gas fraction. Scaled from Savvidou et.
                al. 2020 Eq 8.

            'temp' : str
                Sets position based on disk temperature >= 160 K. If
                'movetemp' is True, the disk temperature will decrease
                by (a pseudo-physical) scaling with the dust/gas
                fraction.

            float
                Sets the position directly (au). If 'movetemp' is True,
                the disk temperature will decrease by (a
                pseudo-physical) scaling with the dust/gas fraction.

    movetemp : bool, default False
        Determines if the snow line evolves when `snowmode` is 'temp' or
        float. Experimental setting to "artifically" evolve the snow
        line and compare to models using 'evol'.

    MSol : float
        Stellar mass, in solar mass units. Default is 1.0.

    dmf : float
        Disk mass fraction. If `dmf` < 1.0, then sets the total (gas +
        solid) initial mass of disk as fraction of the stellar mass. If
        `dmf` > 1.0, then it represents the total starting solid mass,
        in Earth masses, and calculates the full disk mass using 'z0'.
        Setting the solid (dust) mass directly may be useful for other
        model comparison. Default is 0.2.

    z0 : float
        Metallicity of the disk, specifically - the initial dust/gas
        fraction. Default is 0.0134.

    tmode : 'simp'
        Obsolete: intended to change temperature profiles between a
        power law scaled with stellar mass, and a viscous/irradation
        regime. Calculations for viscous/irradation are present but
        unused.

    Attributes
    ----------
    st : ndarray
        Flux-averaged Stokes number of pebbles at given [t,r].

    flux : ndarray
        Mass flux of pebbles at given [t,r].

    rgrid, tgrid : ndarray
        Array containing positions and times across disk lifetime.
        Typically define shape for other disk property arrays (i.e. over
        the disk radius and with time). Units cgs.

    SigmaGas, SigmaDust : ndarray
        Describing the gas and dust surface density. Units of g/cm2.

    msolids : float
        Initial total solid dust mass of disk, in grams.

    mdisk : float
        Initial total disk mass (gas + solids), in grams.

    Temp : ndarray
        Temperature profile, units K.

    Mstar : float
        Stellar mass in grams.

    eta : ndarray
        Gas pressure gradient - r only.

    Hgas : ndarray
        Gas scale height - r only.

    fDG : float
        Disk average dust/gas fraction. Set with instance creation,
        decreases if `snowmode` is evol.

    snowline_i, snowline_au : ndarray
        Position of the snow line at each time (i.e. len(tgrid)), as the
        closest index in rgrid, or in au.

    comp : ndarray
        Density of solids at each [t,r], in g/cm3. Used to initialize
        seed compositional density.

    name : str
        Name that can be separate from that of ``Disk``, used in some
        plotting.

    alpha : float
        Alpha-viscosity parameter.

    snowmode : str or float
        Records which setting was used to set the snow line.

    MSol : float
        Stellar mass in solar mass units.

    dmf : float
        Disk mass fraction.

    fullHgas, fulletas : ndarray
        Future: same as ``eta`` and ``Hgas`` above, except for all
        [t,r]. Intended for development of gas evolution.

    Methods
    -------
    inspect(plots=False)
        Prints basic informaiton about Disk, option to create plots.
    """

    def __init__(
            self, name='Unnamed', alpha=1.e-3, vfrag=1000, vfin=100,
            vfout=1000, rcrit=30, snowmode='none', MSol=1.0, dmf=0.2,
            z0=0.0134, movetemp=False, tmode='simp'):

        Nr = 301  # number of grid points
        Rin = 0.01*au  # inner disk edge
        Rout = 1000.*au  # outer disk edge
        rgrid = np.logspace(np.log10(Rin), np.log10(Rout), Nr)

        Nt = 1001
        endtime = 1.e7*year
        tgrid = np.logspace(np.log10(year), np.log10(endtime), Nt)

        Mstar = MSol*MS
        Z0 = z0

        # Overall mass of disk, either as fraction or explicit solids,
        # in Earth masses.
        if dmf >= 1.0:
            msolids = dmf
            mdisk = msolids/Z0 * ME  # grams
        else:
            mdisk = dmf * MSol * MS
            msolids = mdisk * Z0

        # Notation to match PP
        rout = rcrit*au

        SigmaGas = (mdisk / (2.*np.pi*rout**2.) * (rgrid/rout)**(-1.)
                    * np.exp(-1.*(rgrid/rout)))
        SigmaDust = Z0*SigmaGas

        alpha = alpha  # turbulence strength parameter
        # vfrag is turned into array to fit into adjusted PP code i.e.
        # takes vfrag at each timestep.
        vfrag = vfrag * np.ones((len(tgrid), len(rgrid)))  # units cm/s
        comp = np.ones_like(vfrag)
        rhop = 1.25  # internal density of dust grains

        # FUTURE: The following are vestigial and unimplemented
        # alternatives to the power law temperature profile that rely on
        # viscous/irradiation processes in disk.
        # lum = lumscale(msol=MSol)
        if tmode == 'simp':
            T = np.ones_like(vfrag) * (280 * MSol**0.5 * (rgrid/au)**(-1/2))
        
        elif tmode == 'ida2019':
            # from Ida 2019:
            tvis = (130 * (alpha/1e-2)**(-1/5) * MSol**(3/10) * (1)**(2/5)
                    * (rgrid/au)**(-9/10))
            tirr = (130 * (lumscale(MSol))**(2/7) * MSol**(-1/7)
                    * (rgrid/au)**(-3/7))
            T = np.maximum(tvis, tirr)
            # T = sci.interpolate.UnivariateSpline(rgrid/au, T, k=2)(rgrid/au)
            # this smooths the corner where tvis < tirr. overkill?
            T = T * np.ones_like(vfrag) #turns the temp prof into array[t,r].
        
        elif tmode == 'ida2016':
            mdot = 10**-8 * MSol**2
            # from Ida 2016:
            tvis = (200 * (MSol**(3/10) * alpha/1e-3)**(-1/5)
                    * (mdot/(1e-8))**(2/5) * (rgrid/au)**(-9/10))
            tirr = (150 * (lumscale(MSol))**(2/7) * MSol**(-1/7)
                    * (rgrid/au)**(-3/7))
            T = np.maximum(tvis, tirr)
            # T = sci.interpolate.UnivariateSpline(rgrid/au, T, k=2)(rgrid/au)
            # this smooths the corner where tvis < tirr. overkill?
            T = T * np.ones_like(vfrag) #turns the temp prof into array[t,r].

    
        # Given a disk, evolved over certain amount of time, produce the
        # types of pebbles [Stokes value] at a given space, at a given
        # time, in the disk.

        # SM: Originally, output st, flux. I modified pp.py to track and
        # output eta (pressure gradient), Hgas (scale height), and fDG.

        ##########################
        # SNOW MODE CALCULATIONS #
        ##########################
        if snowmode == 'none':
            st, flux, etas, Hgas, fDG, Hgas_temptest, etas_temptest, Q, ignore = \
                pebble_predictor(
                    rgrid=rgrid, tgrid=tgrid, Mstar=Mstar, SigmaGas=SigmaGas,
                    T=T, SigmaDust=SigmaDust, alpha=alpha, vfrag=vfrag,
                    rhop=rhop)
            self.snowline_i = None
            self.snowline_au = None
            self.comp = comp * innercomp

        else:
            fDG = pebble_predictor(
                    rgrid=rgrid, tgrid=tgrid, Mstar=Mstar, SigmaGas=SigmaGas,
                    T=T, SigmaDust=SigmaDust, alpha=alpha, vfrag=vfrag,
                    rhop=rhop)[4]

            # Evolve snow line position, using scaled Sav+2020 Eq 8.
            if snowmode == "evol":
                # Run PP once with some initial settings to get the fDG.
                # This assumes NO snow line at first, so SigmaDust is
                # slightly over-estimated.

                AU1 = rgrid.searchsorted(1*au)
                snowline_au = (9.2 * (alpha/1.e-2)**(0.61)
                               * (SigmaGas[AU1]**0.88/1000)**(0.80)
                               * (fDG/0.01)**(0.37))
                snowline_i = rgrid.searchsorted(snowline_au*au)

            elif snowmode == "temp":
                movein = (fDG / np.max(fDG))
                movein = movein**(1/6)  # CAREFUL: this is arbitrary
                snowline_i = np.zeros(len(tgrid), dtype=int)
                snowline_au = np.zeros(len(tgrid))

                for t in np.arange(len(tgrid)):
                    if movetemp:
                        T[t] = T[t] * movein[t]
                    snoindex = np.max(np.where(T[t] > 160))
                    snowline_i[t] = snoindex
                    snowline_au[t] = rgrid[snoindex] / au

            elif type(snowmode) is float:
                # Current fucntion for "float" is a fixed SL that does
                # not evolve, even with "movetemp".
                movein = (fDG / np.max(fDG))
                movein = movein**(1/6)  # CAREFUL: this is arbitrary
                snowline_au = np.ones(len(tgrid)) * snowmode
                snowline_i = (np.ones(len(tgrid), dtype=int)
                              * np.max(np.where(rgrid < snowmode * au)))

                for t in np.arange(len(tgrid)):
                    if movetemp:
                        T[t] = T[t]*movein[t]
                        snowline_au[t] = snowmode * movein[t]
                        snowline_i[t] = np.max(np.where(rgrid < (snowline_au[t]
                                                                 * au)))
                    # SL remains in fixed position.
            ###########################
            for t in np.arange(len(tgrid)):
                vfrag[t, :snowline_i[t]] = vfin
                vfrag[t, snowline_i[t]:] = vfout
                comp[t, :snowline_i[t]] = innercomp
                comp[t, snowline_i[t]:] = outercomp

            st, flux, etas, Hgas, fDG, Hgas_temptest, etas_temptest, Q, ignore = \
                pebble_predictor(
                    snowline_i=snowline_i, rgrid=rgrid, tgrid=tgrid,
                    Mstar=Mstar, SigmaGas=SigmaGas, T=T, SigmaDust=SigmaDust,
                    alpha=alpha, vfrag=vfrag, rhop=rhop)

            # Need np.ones in case movetemp is False.
            self.snowline_i = snowline_i * np.ones_like(tgrid)
            self.snowline_au = snowline_au * np.ones_like(tgrid)
            self.comp = comp
            
            SigmaDust = np.concatenate((SigmaDust[:snowline_i[0]] * 0.5,
                SigmaDust[snowline_i[0]:]))

        # st     = Stokes numbers, array of size disk v time
        # flux   = total mass flux, array at points disk v time
        # etas   = pressure gradient, constant over time in model
        # Hgases = scale height of the gas disk, in AU

        self.st = st  # unitless, tau
        self.flux = flux  # g/cm2 - cgs
        self.rgrid = rgrid  # cgs
        self.tgrid = tgrid  # cgs
        self.MSol = MSol
        self.Mstar = Mstar  # star mass, cgs
        self.SigmaGas = SigmaGas  # cgs
        self.SigmaDust = SigmaDust  # cgs
        self.alpha = alpha  # unitless
        self.Temp = T  # Kelvin
        self.dmf = dmf  # disk mass fraction, 0.1, 0.2, etc.
        self.eta = etas  # negative
        self.Hgas = Hgas  # in cm, scale height, not aspect ratio (hgas)
        self.fDG = fDG
        self.name = str(name)
        self.msolids = msolids  # cgs
        self.mdisk = mdisk  # cgs
        self.snowmode = snowmode
        self.fullHgas = Hgas_temptest  # changes with time
        self.fulletas = etas_temptest  # changes with time
        self.movetemp = movetemp

    def inspect(self, plots=False):
        """Convenience function to print basic disk properties.

        Quickly prints basic physical properties of the disk. Option
        to create plots of st, flux, and surface density profiles.

        Parameters
        ----------
        plots : bool, default=False
            If True, plots 1) `st`, `flux` every 10 timesteps, 2) Gas
            and dust surface density profiles.
        """

        d = self
        ts = np.arange(len(d.tgrid))[::100]
        tc = pl.cm.coolwarm_r(np.linspace(0, 1, len(d.tgrid)))

        # To display solid mass properly:
        if d.dmf >= 1:
            solidmass = d.dmf
        elif d.dmf < 1:
            solidmass = d.msolids / ME

        print("Mstar = {} MSol\n"
              "Disk mass = {:2.3f} MSol / {:.0f} MEarth\n"
              "Dust mass = {:4.0f} MEarth\n"
              "SigGas/Dust at 1 AU = {:.0f} / {:.0f} [g/cm^2]".format(
                d.Mstar / MS,
                d.mdisk / MS, d.mdisk * d.Mstar / MS / ME,
                solidmass,
                d.SigmaGas[d.rgrid.searchsorted(1 * au)],
                d.SigmaDust[d.rgrid.searchsorted(1 * au)]))

        if d.snowmode == 'none':
            print("Snowline = None")
        elif d.snowmode == 'evol' or d.movetemp:
            print("Snowline Starts/Ends = {:.2f}/{:.2f} AU".format(
                d.snowline_au[0], d.snowline_au[-1]))
        elif type(d.snowmode) is float or str('temp'):
            print("Snowline = {:.2f} au".format(d.snowline_au[0]))

        if plots:
            fig, axs = plt.subplots(2, 1, figsize=(8, 8))
            fig.subplots_adjust(hspace=0.35)
            for t in ts:
                axs[0].plot(d.rgrid/au, d.st[t, :], c=tc[t],
                            label='{:1.1e} yr'.format(d.tgrid[t]/year))
                axs[1].plot(d.rgrid/au, d.flux[t, :]/ME*year, c=tc[t],
                            label='{:.2e} yr'.format(d.tgrid[t]/year))

            axs[0].set(yscale='log', xscale='log',
                       ylabel='St', xlabel='r [AU]',
                       ylim=[1.e-9, 10], xlim=[0.1, 300],
                       title='Stokes number')
            axs[0].legend(fontsize='small', ncols=3)
            axs[1].set(yscale='log', xscale='log',
                       ylabel='Pebble mass flux [{}/yr]'.format(
                        r'$M_{\oplus}$'), xlabel='r [AU]',
                       ylim=[1.e-14, 1.e-1], xlim=[0.1, 300],
                       title='Pebble mass flux')

            fig, ax = plt.subplots()
            ax.plot(d.rgrid/au, d.SigmaGas, label='SigmaGas')
            ax.plot(d.rgrid/au, d.SigmaDust, label='SigmaDust')
            ax.legend()
            ax.set(yscale='log', xscale='log',
                   ylim=[1.e-12, 1.e5], xlim=[0.1, 400],
                   ylabel="g/$cm^2$", xlabel='r [AU]')
            fig.suptitle(
                '{}{} + {} dmf\n Snowline: {}'.format(
                    d.Mstar/MS, r'$M_{\odot}$', d.dmf, d.snowmode)
                        )


class Seeds():
    def __init__(self, disk, seeds_au, mass, tintro=1, scale=False):
        """
        Containing location, mass, formation time of protoplanet seeds.

        This class creates Seeds objects that describe the initial
        location, mass, and formation time of protoplanet seed masses.
        This class requries a Disk object to "place" the seeds in as
        parameter. The Seeds object also contains the Disk object and
        its attributes, i.e. myseeds.mydisk.attribute. The method
        Seeds.grow() performs the seed growth, with outputs contained in
        instance attributes. This allows multiple Seed configurations to
        be used with the same disk parameters without creating a new
        Disk object each time.

        Parameters
        ----------
        disk : obj
            The Disk object in which these protoplanet seed masses are
            placed. All Disk attributes accessbile.
        seeds_au : 1d array, float
            Array containing locations of seed masses, au.
        mass : float or array
            Initial seed masses, in Earth mass. Can be single number,
            which is used for all seeds, or must be ``len(seeds_au)``,
            where each initial mass is set individually.
        tintro : float or array, default=1
            Future: can set the formation time of seeds as index in
            tgrid, represented as time when seed will "appear" in disk.
            Can be single number to set for all seeds, or array of
            ``len(seeds_au)`` to set individually.
        scale : bool, default=False
            Experimental: If True, will scale seed locations with
            stellar mass, such that the new locations will have same
            period as with 1 solar mass star.

        Attributes
        ----------
        disk : obj
            Disk object in which seeds are placed/grown.
        seeds_au : ndarray
            Seed semi-major axis, in au.
        seed_rind : ndarray
            Seed locations as indexed in `rgrid`.
        seeds : ndarray
            Seed semi-major axis, in cm.
        mass : ndarray
            Seed masses, in Earth masses.
        qp : ndarray
            Protoplanet seed mass : stellar mass ratio.
        dens : ndarray
            Initial protoplaneet seed compositional density, in g/cm3.
        it : ndarray
            Protoplanet formation/introduciton time index in `tgrid`.

        Methods
        -------
        grow(mode='isocutfiltfrac', epmode='setf', pref=20)
            Grows protoplanet seed masses by pebble accretion.
        """

        tgrid = disk.tgrid
        rgrid = disk.rgrid
        self.disk = disk
        self.seeds_au = seeds_au
        self.seeds = seeds_au * au
        # If one value, use for all seeds. Earth masses.
        self.mass = mass * np.ones(len(seeds_au))
        self.qp = self.mass * ME / disk.Mstar

        if scale:
            # This option scales mass and loc. of seeds to stellar mass.
            # Mass = mass x star mass (solar)
            # SMA  = scaled by constant period.
            self.mass = mass * disk.Mstar/MS * np.ones(len(seeds_au))
            self.qp = self.mass * ME / disk.Mstar
            self.seeds_au = seeds_au * (disk.Mstar/MS)**(1/3.)
            self.seeds = self.seeds_au * au

        self.seed_rind = rgrid.searchsorted(self.seeds)
        self.dens = disk.comp[0, self.seed_rind]

        # Sets introduction/formation time of seed.
        tintro = tintro * year
        # Index in tgrid:
        self.it = tgrid.searchsorted(tintro*np.ones(len(seeds_au)))

    def grow(self, mode='isocutfiltfrac', epmode='setf', pref=20):
        """
        Grows protoplanet seed masses by pebble accretion.

        Grows the protoplanet seeds described by Seeds in the disk
        described by Disk. See theory and description in paper.
        Protoplanets are grown at each time step, starting from the
        outside in. This facilitates calculation of filtering.

        Parameters
        ----------
        mode : {'isocutfiltfrac', 'filtfrac', 'filtraw', 'isocut'}
            Sets certain growth modes regarding filtering of inward
            pebbles, and whether to truncate protoplanet growth at the
            pebble isolation mass. Strings can be combined in any order
            but must be exact. Ex. to filter the physical amount of
            pebbles and cut off growth at isolation mass, use
            'filtrawisocut'.
        epmode : string, default='setf'
            Accretion mode used in epsilon to specifiy, e.g. 2dset,
            3dbalf.
        pref : float, default=20
            Mathematical prefactor (Earth masses) used to calculate
            pebble isolation mass, as values differ among literature.

        Attributes
        ----------
        massgained : ndarray [s,t]
            New mass accreted per seed per time step, Earth masses.
        cumulmass : ndarray [s,t]
            Total cumulative mass of seed per time step, Earth masses.
        finalmass : ndarray [s]
            Final masses of each seed, Earth masses.
        preflux, fluxpast : ndarray [s,t]
            Represents the time-integrated mass in pebbles drifting past
            a given seed mass location. ``preflux`` corresponds to mass
            before any pebbles are removed from flux (filtering);
            fluxpast describes the actual mass in pebbles that make it
            past a seed, accounting for filtering.
        eff : ndarray [s,t]
            Fraction of drifting pebbles that accrete onto seed.
            Calculated via epsilon.
        isos : ndarray [s,t]
            Pebble isolation mass threshold per seed, Earth masses.
        madeiso, isotime : ndarray [s]
            Indices of seeds_au that reached their pebble isolation mass
            and at which index of tgrid did that occur. These are
            recorded "backwards", but since they refer to the seed_au
            and time index.
        qps : ndarray [s,t]
            Protoplanet seed:stellar mass ratio, per seed per time.
            Quantity used in pebble_predictor.
        rads : ndarray [s,t]
            Future: Protoplanet radii, cm, used for some varitions of
            epmode. Currently stops recording a seeds radius after it
            reaches iso.
        epmode : string
            Records the epmode used to run grow.
        seedcomp : ndarray [s,t]
            Per seed bulk density, g/cm3.
        wmf : ndarray [s,t]
            Per seed water mass fraction, i.e. fraction of total mass
            that is water. Unitless.
        tpmass : ndarray [t]
            Total mass contained in protoplanet seeds per time, Earth
            masses."""

        # Init arrays, each row is a seed mass of len(tgrid).
        massgained = np.zeros((len(self.seeds), len(self.disk.tgrid)))
        efficien = np.zeros_like(massgained)
        cumulmass = np.expand_dims(self.mass, axis=1) * np.ones(len(self.disk.tgrid))
        qps = np.expand_dims(self.qp, axis=1) * np.ones(len(self.disk.tgrid))
        
        rads = np.zeros_like(massgained)
        isomass = np.ones_like(massgained)
        isotime = []
        madeiso = []

        # Should produce array of initial seedmass radii, in cm.
        # Used in some epmodes.
        radi = RfromD(m=self.mass*ME, d=2.8)

        # Reshapes flux array from PP, slicing only at seed locations.
        # Used in filtering.
        # Creates the integrated flux past each seed, cgs.
        fluxpast = np.transpose(self.disk.flux[:, self.seed_rind])
        for t in np.arange(len(self.disk.tgrid)):
            if t == 0:
                pass
            else:
                fluxpast[:, t] = fluxpast[:, t] * (self.disk.tgrid[t]
                                                   - self.disk.tgrid[t-1])

        preflux = np.transpose(self.disk.flux[:, self.seed_rind])
        for t in np.arange(len(self.disk.tgrid)):
            if t == 0:
                pass
            else:
                preflux[:, t] = preflux[:, t] * (self.disk.tgrid[t]
                                                 - self.disk.tgrid[t - 1])

        # Iterate through seeds from outside in.
        tnuoc = np.arange(len(self.seeds))
        for i in tnuoc[::-1]:
            ir = self.seed_rind[i]  # index in rgrid of location
            loc = self.seeds[i]  # location of seed mass in cm
            it = self.it[i]  # intro time of seed. i.e. when it appears.

            # massgained: track new mass accreted each timestep.
            # cumulmass: the "growth track".
            cumulmass[i, it] = self.mass[i]
            qps[i, it] = self.qp[i]
            rads[i, it] = radi[i]

            # Initial accretion efficiency.
            # ep   = eccentricity parameter
            # tau  = pebble Stokes number (from pp, at [time, location])
            # qp   = mass ratio of seed mass:star
            # hgas = aspect ratio (not scale height (Hgas)!)
            # eta  = pressure gradient...but is it?
            # Rp   = radius of planet, in radius:distance ratio (a/Rp)

            # An initial accretion efficiency.
            effic = OL18.epsilon_general(tau=self.disk.st[it, ir],
                                         qp=self.qp[i],
                                         eta=(abs(self.disk.eta[ir])),
                                         mode=epmode, alphaz=self.disk.alpha,
                                         hgas=self.disk.Hgas[ir]/loc,
                                         Rp=radi[i]/loc)

            # Set up a few arrays with the useful info per timestep.
            efficien[i, it] = effic

            # newqp is variable needed below.
            newqp = 0

            # TIME STEP THROUGH EACH SEED MASS:
            # Skips first index to avoid nonexistent time index.
            # it = index of introduction time.

            for t in np.arange(it, len(self.disk.tgrid)):
                # Updates isolation mass, with changing scale height
                # each time step, if temp changes.
                isomass[i, t] = Miso(ms=self.disk.Mstar/MS,
                                     lociso=self.seeds_au[i],
                                     hgasiso=self.disk.fullHgas[t, ir]/au,
                                     pref=pref)

            for t in np.arange(it, len(self.disk.tgrid)):
                if t == it:
                    pass
                else:
                    newmass_cgs = effic * fluxpast[i, t]
                    newmass = newmass_cgs / ME
                    cumulmass[i, t] = cumulmass[i, t-1] + newmass
                    # Filtering for all seed masses except innermost seed.
                    if (i > 0
                            and t < len(self.disk.tgrid) - 2
                            and mode.count('filtraw')):
                        filtmass = fluxpast[:i, t+1] - newmass_cgs
                        for f in np.arange(len(filtmass)):
                            if filtmass[f] <= 0:
                                fluxpast[:i, t+1] = 0
                            else:
                                fluxpast[:i, t+1] = filtmass

                    # Alternatively filtering by a fractional decrease
                    # instead of raw substraction.
                    if (i > 0
                            and t < len(self.disk.tgrid) - 2
                            and mode.count('filtfrac')):
                        fluxpast[:i, t+1] = fluxpast[:i, t+1] * (1 - effic)

                    # Isolation mass cut off.
                    if (mode.count('isocut')
                            and cumulmass[i, t] >= isomass[i, t]):
                        cumulmass[i, t:] = isomass[i, t]
                        massgained[i, t:] = 0
                        # Note pebble reduction is non-zero.
                        fluxpast[:i+1, t+1:] *= 0.001
                        isotime.append(t)
                        madeiso.append(i)
                        break

                    # New mass added per year, Earth masses:
                    massgained[i, t] = newmass
                    newqp = cumulmass[i, t] * ME / self.disk.Mstar
                    rad = RfromD(m=newqp*self.disk.Mstar, d=2.8)
                    effic = OL18.epsilon_general(
                        ep=0, tau=self.disk.st[t, ir], qp=newqp, mode=epmode,
                        eta=(abs(self.disk.fulletas[t, ir])),
                        alphaz=self.disk.alpha,
                        hgas=self.disk.fullHgas[t, ir]/loc,
                        Rp=rad/loc)
                    qps[i, t] = newqp
                    efficien[i, t] = effic
                    rads[i, t] = rad

        finalmass = np.ones_like(self.seeds)
        for i in np.arange(len(self.seeds)):
            finalmass[i] = np.nanmax(cumulmass[i])

        massgained[:, 0] = self.mass
        seedcomp = np.empty_like(cumulmass)
        wmf = np.zeros_like(massgained)

        # Seed composition, in bulk density for each time, t. g/cm3.
        for s in np.arange(len(self.mass)):
            seedcomp[s] = (np.nancumsum(
                                massgained[s] * ME
                                * self.disk.comp[:, self.seed_rind[s]])
                           / np.nancumsum(massgained[s]*ME))
            # Water mass fraction (wmf). Checks if seed is past snow
            # line, then assumes half new mass gained is water. wmf is
            # the cumulative wmf at that time step.

            for t in np.arange(len(self.disk.tgrid)):
                if self.disk.snowmode == 'none':
                    wmf[s] = np.cumsum(wmf[s]) / cumulmass[s]
                    break
                elif self.seed_rind[s] > self.disk.snowline_i[t]:
                    wmf[s, t] = massgained[s, t] * 0.5
            wmf[s] = np.cumsum(wmf[s]) / cumulmass[s]
            # This last line should be outside the t-for loop.
            # Otherwise, madness.

        self.massgained = massgained  # Earth masses, integrated
        self.cumulmass = cumulmass  # Earth masses
        self.finalmass = finalmass  # Earth masses, final masses
        self.fluxpast = fluxpast  # g, integrated
        self.preflux = preflux  # int flux (g) from PP, before filtering
        self.eff = efficien  # from OL2018
        self.isos = isomass  # note these go "backward"
        self.madeiso = madeiso  # ind of seed array reached iso mass
        self.isotime = isotime  # ind in tgrid when reached iso mass
        self.qps = qps
        self.rads = rads
        self.ep = epmode
        self.seedcomp = seedcomp
        self.wmf = wmf  # water mass fraction
        self.tpmass = np.sum(cumulmass, axis=0)

    def gtplot(self, xlim=[0.07, 25], ylim=[0.75e-3, 40], ysnow=0.105, **kwargs):
        """Convenience plotter for single a single run.

        Quick plotting feature to show protoplanet masses, water mass
        fraction, pebble isolation mass at each seed location, and the
        position/evolution of the snow line (if applicable).

        Parameters
        ----------
        xlim, ylim : [min, max]
            Sets axis limits for plot, au.

        ysnow : float
            Fine placement control of the label and arrow for the snow
            line. Arrow will be -0.015 label, in axis units.

        **kwargs
            title : str
                User can specify a title of any string text. Default
                title displays the stellar mass / dmf / number of seeds.
        """

        # Set plot title to default, or use kwarg name if given.
        title = '{}{} / {} dmf / {} seeds'.format(
            self.disk.MSol, r'$M_{\odot}$', self.disk.dmf, len(self.seeds_au))
        if 'title' in kwargs:
            title = kwargs['title']
        figsize = [6,5]
        if 'figsize' in kwargs:
             figsize = kwargs['figsize']

        fig, ax = plt.subplots(constrained_layout=True, figsize=figsize)

        # Colors
        cmap = (ListedColormap(
            ['peru', 'firebrick', 'lightgreen', 'lightblue', 'dodgerblue']
            ).with_extremes(under='peru', over='darkblue'))
        bounds = [0, 1e-6, 1e-3, 1e-2, 0.1, 0.25]
        boundlabs = ['{}'.format('$0$'),
                     '{}'.format('$10^{-6}$'),
                     '{}'.format('$10^{-3}$'),
                     '{}'.format('$10^{-2}$'),
                     '{}'.format('$10^{-1}$'),
                     '{}'.format('$0.25$')]
        norm = BoundaryNorm(bounds, cmap.N, clip=False)
        slc = 'royalblue'  # snow line color

        # Colorbar:
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                            ax=ax, label='water mass fraction', extend='max',
                            shrink=0.8, aspect=30)
        cbar.ax.set_yticklabels(boundlabs)

        # Set axes parameters,labels.
        ax.set(yscale='log', xscale='log',
               ylim=ylim, xlim=xlim,
               ylabel=r'Mass [$M_{\oplus}$]', xlabel='r [AU]',
               title=title)

        # Plot final masses color-coded by water mass fraction.
        ax.scatter(self.seeds_au, self.finalmass, zorder=1, c=self.wmf[:, -1],
                   cmap=cmap, norm=norm, marker='o', s=90, edgecolor='k',
                   linewidth=1)

        # Plot the pebble isolation mass as dashed line, mark seeds that
        # reach isolation mass.
        ax.plot(self.seeds_au, self.isos[:, -1], ls=':', c='darkgrey', zorder=0,
                alpha=0.95, label="pebble isolation mass")
        ax.scatter(self.seeds_au[self.madeiso], self.finalmass[self.madeiso],
                   marker='x', c='red', s=25, label='Made Iso Mass')

        # Snow line conditionals...
        if self.disk.snowmode == 'none':
            ax.text(0.1, 0.9, "No snow line", transform=ax.transAxes)

        elif self.disk.snowmode == 'evol' or self.disk.movetemp:
            # Plot initial and final snow line, as vertical lines.
            ax.axvline(self.disk.snowline_au[-1], c=slc, zorder=1, alpha=0.6)
            ax.axvline(self.disk.snowline_au[0], ls='--', c=slc, zorder=0,
                       alpha=0.6)
            ax.annotate("snow line", color=slc,
                        xytext=(self.disk.snowline_au[-1]+0.18, ysnow),
                        xy=(self.disk.snowline_au[-1], ysnow))
            ax.annotate("", xytext=(self.disk.snowline_au[0], ysnow-0.015),
                        xy=(self.disk.snowline_au[-1], ysnow-0.015),
                        arrowprops=dict(arrowstyle="->", color=slc))

        elif type(self.disk.snowmode) is float or str('temp'):
            ax.axvline(self.disk.snowline_au[0], c=slc, zorder=1, alpha=0.6)
            ax.annotate("snow line", color=slc,
                        xytext=(self.disk.snowline_au[0]+0.18, 0.105),
                        xy=(self.disk.snowline_au[0], 0.105))
        ax.legend()
        #return fig

def RfromD(m, d):
    """Calculate radius from mass and density i.e. d=m/v."""
    return (0.75 * m/d * 1/np.pi)**(1/3.)


def Miso(ms=1.0, lociso=10, hgasiso=0.03, pref=20):
    """Returns pebble isolation mass, Earth masses.

    Prefactor set to 20 Em, but see literature for discussion between 20
    or 40 for prefactor."""
    return pref * ms * ((hgasiso/lociso)/0.05)**3
