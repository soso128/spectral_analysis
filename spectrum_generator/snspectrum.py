from numpy import *
from scipy.integrate import quad, simps
import snrate as sn

"""
    This code is used as an interface with the C++ code snrate.cpp to generate and plot
    DSNB spectra and IBD spectra in SK

"""

# Typical choices of parameters for IMF and Cosmic Star Formation Rate (CSFR)
imfs = {"salpeter": {'xi1': 2.35, 'xi2': 2.35}, "kroupa": {'xi1': 1.3, 'xi2': 2.3}, 'BG': {'xi1': 1.5, 'xi2': 2.15}}
csfr_up = {'rho0': 0.0213, 'alpha': 3.6, 'beta': -0.1, 'gamma': -2.5, 'z1': 1, 'z2': 4}
csfr_fiducial = {'rho0': 0.0178, 'alpha': 3.4, 'beta': -0.3, 'gamma': -3.5, 'z1': 1, 'z2': 4}
csfr_lower = {'rho0': 0.0142, 'alpha': 3.2, 'beta': -0.5, 'gamma': -4.5, 'z1': 1, 'z2': 4}

# CCSN rate as a function of redshift
def snrate(z, imfpars, csfrpars):
    """
    CCSN rate as a function of redshift
    imfpars and csfrpars are dictionaries
    """
    return sn.csfr(z, **csfrpars) * sn.integrated_imf(8, 50, **imfpars)/sn.integrated_mimf(0.1, 100, **imfpars)

# Integrand for DSNB rate as a function of redshift
def SRNflux_integrand(z, enu, tnu, imfpars, csfrpars):
    """
    Integrand for DSNB rate as a function of redshift
    imfpars and csfrpars are dictionaries for the corresponding parameters
    """
    return snrate(z, imfpars, csfrpars) * sn.spectrum_blackbody(enu * (1 + z), {"tnu": tnu}) * sn.cosmology(z)

# DSNB rate as a function of neutrino energy
def SRNflux(enu, tnu, imfpars, csfrpars):
    """
    DSNB rate as a function of the neutrino energy enu
    tnu is the effective neutrino temperature assuming a blackbody emission spectrum
    imfpars and csfrpars are dictionaries for the corresponding parameters
    """
    rate = quad(SRNflux_integrand, 0, 5, args = (enu, tnu, imfpars, csfrpars))
    if rate[1] > 0.1 * rate[0]: 
        rate = quad(SRNflux_integrand, 0, 5, args = (enu, tnu, imfpars, csfrpars), epsabs = rate[0]/100)
    return c/H0 * rate[0], c/H0 * rate[1]

# Integrand for SK DSNB rate as a function of cos(theta)
def ibd_integrand(c, ee, specpar, imfpars, csfrpars):
    """
    Integrand for IBD rate in SK (22.5 kton year) as a function of cos(theta)
    theta is the scattering angle
    ee is the positron energy
    tnu is the effective neutrino temperature assuming a blackbody emission spectrum
    imfpars and csfrpars are dictionaries for the corresponding parameters
    """
    en = sn.enu(ee, c)
    args = {**csfrpars, **imfpars}
    # Use C++ function for DSNB rate calculation (muuuch faster)
    res = sn.snrate(en, specpar, *args.values()) * sn.dsigma(en, 0, c) * 0.047390723e42
    return res

# SK DSNB spectrum as a function of positron energy
def ibd_spectrum(ee, specpar, imfpars, csfrpars):
    """
    IBD rate in SK (22.5 kton year)
    ee is the positron energy
    tnu is the effective neutrino temperature assuming a blackbody emission spectrum
    imfpars and csfrpars are dictionaries for the corresponding parameters
    """
    # Scipy quad much faster than C++ Simpson
    return quad(ibd_integrand, -1, 1, args = (ee, specpar, imfpars, csfrpars))[0]

# Scan
# Here, we scan over the BH fraction fBH and the rest is fixed
def param_scan(pbar, fbhmax, fbhstep):
    # Spectrum parameters
    #       type gives the type of parameterization (see snrate.cpp)
    #       fBH is the BH-forming supernova fraction
    #       emean: mean energies
    #       alpha: pinching parameters
    #       lum: luminosities
    #       pbar: survival probability for antineutrinos
    specpar = {"type": 3, "fBH": 0, "emeanE_NS": 15, "emeanE_BH": 24, "alphaE": 3.5, "emeanX_NS": 18, "emeanX_BH": 25, "alphaX": 2.5, "lumE_BH": 2.6, "lumX_NS": 1, "lumX_BH": 0.4 * 2.6, "lumi": 5e52 * 1e3/0.0016, "pbar": pbar}
    f = arange(0, fbhmax + fbhstep, fbhstep)
    x = arange(0, 30.1, 0.1)
    spectra = []
    for ff in f:
        print(ff)
        specpar['fBH'] = ff
        spectra.append(list(map(lambda ee: ibd_spectrum(ee, specpar, imfs['salpeter'], csfr_up), x)))
    return spectra

def gen_spec():
    specpar = {"type": 0, "tnu": 5, "lumi": 3e53 * 1e3/0.0016}
    x = arange(1, 80.1, 0.1)
    spec = array(list(map(lambda ee: ibd_spectrum(ee, specpar, imfs['salpeter'], csfr_up), x)))
    savetxt("/home/sonia/relic/signal/dsnb_tnu{}_lumi{}_nosmear.txt".format(int(specpar['tnu']), int(specpar['lumi'] * 0.0016/1e56)), column_stack((x, spec)))
    return x,spec

# Scan
def param_scan_full(fbhmax, fbhstep):
    # Spectrum parameters
    #       type gives the type of parameterization (see snrate.cpp)
    #       fBH is the BH-forming supernova fraction
    #       emean: mean energies
    #       alpha: pinching parameters
    #       lum: luminosities
    #       pbar: survival probability for antineutrinos
    specpar = {"type": 3, "fBH": 0, "emeanE_NS": 15, "emeanE_BH": 0, "alphaE": 3.5, "emeanX_NS": 18, "emeanX_BH": 0, "alphaX": 2.5, "lumE_BH": 0, "lumX_NS": 1, "lumX_BH": 0, "lumi": 5e52 * 1e3/0.0016, "pbar": 0}
    f = arange(0, fbhmax + fbhstep, fbhstep)
    ebhe = [15,20,25]
    ebhx = [20,25,33]
    lbhe = array([5,10,15,20])/5.
    x = arange(0, 80.1, 0.1)
    with open("spectra_scan.txt", 'w') as fout:
        for ff in f:
            for ebe in ebhe:
                for ebx in ebhx:
                    for lb in lbhe:
                        for pp in [0,0.68]:
                            specpar['fBH'] = ff
                            specpar['pbar'] = pp
                            specpar['emeanE_BH'] = ebe
                            specpar['emeanX_BH'] = ebx
                            specpar['lumE_BH'] = lb
                            specpar['lumX_BH'] = lb * 0.5
                            for i,c in enumerate([csfr_lower, csfr_fiducial, csfr_up]):
                                print(f"{i} {ff} {pp} {ebe} {ebx} {lb}")
                                spec = list(map(lambda ee: ibd_spectrum(ee, specpar, imfs['salpeter'], c), x))
                                specstr = f"{i} {ff} {pp} {ebe} {ebx} {lb} " + " ".join([f"{ss}" for ss in spec])
                                print(specstr, file = fout)
if __name__ == "__main__":
    param_scan_full(0.4,0.1)
