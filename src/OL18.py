import numpy as np

#Feel free to use, please cite
# - Liu & Ormel (2018) [Paper I]
# - Ormel & Liu (2018) [Paper II]


def epsilon_general (mode='', **pars):
    """
    Returns the pebble accretion efficiency.

    Accepts a variety of "mode" strings that combine to select between
    pebble accretion modes, e.g., 2d vs 3d, ballistic or settling. This
    code has not been modified from the original of Ormel & Liu.

    Parameters
    ----------
    mode : {[2d, 3d,], [f,], [set, bal,]}
    Used to designate accretion modes, can use in combination.
            **[2d, 3d,]** 2d (planar) limit; 3d limit; or mix
            (default)\n
            **[f,]** inclusion of fset modulation factor\n
            **[set, bal,]** settling only; ballistic only;
            mix (default); (mix always triggers calculation of fset
            factor)
        
        Examples
            ``2dset``  2d (planar), settling regime; EqA4, Paper I\n
            ``3dsetf`` 3d regime, w/ modulation factor; Eq21 Paper II\n
            ``set`` interpolation b/w 2d and 3d limits, Eq38 Paper II\n
            ``2dbalf`` 2d, ballistic regime; EqA15 Paper I\n
            ``3dbalf`` 3d, ballistic regime; EqA16 Paper I\n
            ``''`` all effects

    **pars
        A list of keyword arguments.
            qp
                planet-to-star mass ratio
            tau
                dimensionless stoping time (tstop*Omega)
            alphaz
                turbulence diffusivity parameterization in z-direction
            taucorr
                dimensionless turbulence correlation time (tcorr*Omega)
            hgas
                gas aspect ratio
            eta
                gas dimensionless pressure gradient
            sigvec
                turbulence rms velocity components
            nvec
                relative strengths (x,y,z) turbulent velocity components
            Rp
                planet radius over orbital radius (for ballistic regime)

    Notes
    -----
    SM: The calculations of this module remain unchanged from originals, but
    this documentation has been modified to comply with numpydoc style.
    
    OL: In case of turbulence, please specify the rms-velocities explicitly
    via sigvec OR based on alphaz using nvec method. If none is given
    then nvec = (1,1,1) is assumed.
    """

    # These should be present
    tau = pars['tau']
    qp = pars['qp']

    # Cases when we need to calculate fset modulation factor
    if mode.count('f'):
        doCalcfset = True
    else:
        doCalcfset = False

    # When we mix ballistic/settling regimes, always calculate fset
    if mode.count('set')==0 and mode.count('bal')==0:
        doCalcfset = True

    # Obtain the pebble scaleheight (needed for 3D)
    if mode.count('2d')==0:

        # obtain the pebble scaleheight, if not an argument
        # use Youdin & Lithwick expression
        if 'hP' in pars:
            hP = pars['hP']
        else:
            hP = hp_YL07(**pars)

        # find the effective scaleheight in case of inclined planets
        if 'ip' in pars:
            ip = pars['ip']
            heff = heff_app(hP, ip)
        else:
            ip = 0
            heff = hP

    # the 3D limit, settling regime
    if mode.count('2d')==0 and mode.count('bal')==0:
        eps3D = eps_3D (heff=heff, **pars)

        #velocity in z direction
        delVz = ai *ip
    else:
        delVz = 0.
        eps3D = np.inf


    # calculate delVy (usually, but not always, needed below)
    if mode.count('3d')==0 or mode.count('set')==0 or doCalcfset:
        #circular velocity (Paper I)
        vcir = v_circ(**pars)

        if 'ep' in pars:
            ep = pars['ep']
        else:
            ep = 0.0

        delVy = np.maximum(vcir, ae*ep) 
    else:
        delVy = 0.


    # the 2D limiting expression
    if mode.count('3d')==0 and mode.count('bal')==0:
        eps2D = eps_2D (delVy=delVy, **pars)
    else:
        eps2D = np.inf

    # the non-turbulent component of the velocity
    delVvec = [0., delVy, delVz]


    # the turbulent component of the velocity (can be zero)
    if doCalcfset or mode.count('set')==0:
        # the turbulent components of the velocity
        # get it from input parameters (for anisotropic turbulence)
        if 'sigvec' in pars:
            sigvec = pars['sigvec']

        # or relative to alpha-z
        elif 'alphaz' in pars:

            #relative weights specified or not (isotropy)
            if 'nvec' in pars:
                nvec = np.array(pars['nvec'])
            else:
                nvec = np.ones(3)

            sigvec = [nvec[k] *sig_turb(**pars) for k in range(3)]

        # no turbulence
        else:
            sigvec = np.zeros((3))

        #pebble rms velocity different from gas
        if 'taucorr' in pars:
            taucorr = pars['taucorr']
        else:
            taucorr = 1.0

        # particle rms velocities
        sigPvec = [sigvec[k] *np.sqrt(taucorr/(taucorr+tau)) for k in range(3)]
    else:
        sigPvec = np.zeros((3))

    #calculate the fset reduction factor 
    if doCalcfset:
        vast = v_ast(**pars)

        fset = f_set (delVvec, sigPvec, vast)
        fbal = 1-fset
    else:
        fset = 1.
        fbal = 1.

    if mode.count('bal')==0:
        #interpolation formula
        eps2Dset = eps2D *fset
        eps3Dset = eps3D *fset**2

        if mode.count('2d'):
            epsset = eps2Dset
        elif mode.count('3d'):
            epsset = eps3Dset
        else:
            epsset = eps_23 (eps2Dset, eps3Dset)

    else:
        epsset = 0.


    # The ballistic regime (Paper I)
    if mode.count('set')==0:
        # obtain absolute velocity, including turbulence motions
        delV2 = [delVvec[k]**2 +sigPvec[k]**2 for k in range(3)]
        delV = np.sqrt(delV2[0] +delV2[1] +delV2[2])

        
        # if settling is also calculated, reduce eps2Dbal
        # (this mimics aerodynamic deflection)
        if mode.count('3d')==0:
            eps2Dbal = eps_2D_bal (delV=delV, **pars)
            if doCalcfset: eps2Dbal *= (1-fset)
        if mode.count('2d')==0:
            eps3Dbal = eps_3D_bal (delV=delV, hP=heff, **pars)
            if doCalcfset: eps3Dbal *= (1-fset)**2


        if mode.count('2d'):
            epsbal = eps2Dbal
        elif mode.count('3d'):
            epsbal = eps3Dbal
        else:
            # we assume the same mixing expresssion as in the settling case
            epsbal = eps_23 (eps2Dbal, eps3Dbal)

    else:
        epsbal = 0.

    epsgen = epsset +epsbal
    return epsgen


def eps_set(tau, qp, eta, hgas, alphaz):
    """
    Provides settling efficiency WITHOUT accounting for the f_set
    modulation factor or ballistic interactions. ASSUMING zero
    eccentricities, inclinations and standard.
    """
    
    delv = v_circ (tau, eta, qp)  # relative velocity
    # approximation for pebble scale height
    hP = np.sqrt(alphaz / (alphaz+tau)) * hgas
    eps2D = eps_2D (tau, qp, eta, delv)  # 2D expression
    eps3D = eps_3D (tau, qp, eta, hP)  # 3D exprssion
    epsset =  eps_23 (eps2D, eps3D)  # total epsilon
    return epsset


# all fit constants
A2 =    0.322
A3 =    0.393
ash =   0.515
acir =  5.66
aturb = 0.332
ae =    0.764
ai =    0.677
aset =  0.5

#to avoid 0/0
tiny = np.finfo(np.float64).tiny

def v_circ (tau, eta, qp, **dumargs):
    """
    Circular velocity, Eqs. XXXX of Paper I.
    """
    vhw = eta        
    vsh = ash*(qp*tau)**(1./3)

    qc = eta**3/tau
    vcir = vsh +vhw*(1 +acir*(qp/qc))**-1
    return vcir


def hp_YL07 (tau, hgas, alphaz, taucorr = 1.0, **dumargs):
    """
    Returns the particle scaleheight according to Eq.-21 and -28 of
    Youdin & Lithwick (2007).
    """

    h0 = np.sqrt(alphaz/(alphaz+tau)) *hgas  # the Dubrulle expression
    xi = 1 + tau*taucorr**2 /(tau+taucorr)  # correction term

    return h0/np.sqrt(xi)


def heff_app (hP, ip):
    """
    Approximation to the effective scaleheight in case of planet
    inclination (ip<>0). Eq. 26 of OL18.
    """

    # avoid 0/0
    arg = hP**2 + 0.5*np.pi*ip**2*(1 -np.exp(-0.5*ip / (hP+tiny)))
    return np.sqrt(arg)


def sig_turb (alphaz, hgas, taucorr=1.0, **dumargs):
    sigturb = alphaz**0.5*hgas / np.sqrt(taucorr)
    return sigturb


def v_ast (tau, qp, **dumargs):
    """
    Characteristic velocity for pebble accretion.
    """
    return (qp/tau)**(1./3)


def f_set_i (delVi, sigi, vast):
    """
    Calculates fset for a single direction.
    """

    fset = (np.exp(-aset*delVi**2 / (vast**2+aturb*sigi**2)) * vast
           / np.sqrt(vast**2 + aturb*sigi**2))
    return fset

def f_set (delVvec, sigPvec, vast):
    """
    Exponential decay function to account for transition between
    settling and ballistic regimes.
    """

    fset = 1.
    for k in range(3):
        fset *= f_set_i(delVvec[k], sigPvec[k], vast)
    return fset


def eps_3D(tau, qp, eta, heff, **dumargs):
    """
    Calculate accretion efficiency in 3D regime.
    """
    
    eps3D = A3 * qp / (eta*(heff + tiny))
    return eps3D


def eps_2D(tau, qp, eta, delVy, **dumargs):
    """
    Calculate accretion efficiency in 2D regime.
    """

    eps2D = A2 * np.sqrt(qp/tau/eta**2 * delVy)
    return eps2D


def eps_2D_bal(tau, qp, eta, delV, Rp, **dumargs):
    """
    Equation (A.15) of Paper I
    """
    
    eps2Dbal = Rp / (2*np.pi*tau*eta) * np.sqrt(2*qp/Rp + delV**2)
    return eps2Dbal


def eps_3D_bal (tau, qp, eta, delV, Rp, hP, **dumargs):
    """
    Equation (A.15) of Paper I
    """
    
    eps3Dbal = (1. / (4 * np.sqrt(2*np.pi) * tau * eta * hP)
               * (2*qp*Rp/delV + delV*Rp**2))
    return eps3Dbal


def eps_23 (eps2D, eps3D):
    """
    The mixing expressions to transition between 3D and 2D.
    Essentially, eps23 = min(eps2D, eps3D).
    """
    
    eps23 = (eps2D**-2 + eps3D**-2)**-0.5
    return eps23


