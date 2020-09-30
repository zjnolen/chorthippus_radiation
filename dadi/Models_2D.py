#ZN: First set of models are from the dadi_pipeline repo of Portik et al. 2017. Please cite them if using.

import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

#No Migration

def no_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    """
    nu1, nu2, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

#Continuous Migration

def sym_mig(params, ns, pts):
    """
    Split into two populations, with symmetric migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

#Ancestral Migration

def anc_sym_mig(params, ns, pts):
    """
    Split with symmetric migration followed by isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    nu1, nu2, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

#Secondary Contact

def sec_contact_sym_mig(params, ns, pts):
    """
    Split with no gene flow, followed by period of symmetrical gene flow.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m: Migration between pop 2 and pop 1.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    nu1, nu2, m, T1, T2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

# ZN: Following models are from Rougeux et al 2017 and modified by ZJN, please cite if using.

#Continuous Migration and Growth

def CMG(params, ns, pts):
    nu1, nu2, b1, b2, m, T= params
    """
    Model with migration during the divergence.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b: Population growth coefficient
    m: Migration rate between populations (2*Na*m)
    T: The scaled time between the split (in units of 2*Na generations).
    """
    # Define the grid we'll use
    xx = Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # We start the population size change after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/T)
    bnu2_func = lambda t: nu2 * b2**(t/T)
    phi = Integration.two_pops(phi, xx, T, bnu1_func, bnu2_func, m12=m, m21=m)
    ###
    ## Finally, calculate the spectrum.
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

#Secondary Contact and Growth

def SCG(params, ns, pts):
    nu1, nu2, b1, b2, m, Ts, Tsc= params

    """
    Model with split, complete isolation, followed by secondary contact with exponential growth
    nu1: Size of population 1 at split.
    nu2: Size of population 2 at split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """

    # Define the grid we'll use
    xx = Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population reduction after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)
    phi = Integration.two_pops(phi, xx, Tsc, bnu1_func, bnu2_func, m12=m, m21=m)
    ###
    ## Calculate the spectrum
    # Oriented
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
