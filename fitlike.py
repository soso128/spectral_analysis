''' Spectral fits for SRN analysis of SK I-IV '''
from __future__ import division
from sys import path, argv
from ast import literal_eval

from numpy import *
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import fmin
import matplotlib.pyplot as plt
import matplotlib as mpl

path.append("spectrum_generator/")
import snspectrum as sns
import likes
from pdf_sk4 import bg_sk4, relic_sk

mpl.rcParams['font.size'] = 15


# livetimes
livetimes = array([1497.44, 793.71, 562.04, 2790.1])

# energy-independent efficiency sys
efftot = {"lma": [0.7975,0.56703,0.77969],
          "malaney": [0.7903, 0.53273, 0.766262],
          "faild": [0.7997,0.5845,0.7831] ,
          "woosley": [0.7876,0.5524,0.7764],
          "ksw" : [0.7908, 0.54527, 0.77371]}
fluxfac = {"lma": 0.535115, "malaney": 0.5845, "ksw": 0.488413,
           "faild": 0.431, "woosley": 0.47447}
sys_eff = array([0.0254, 0.0404, 0.0253, 0.0214])
sys_eff_sk4_ntag = 0.10
regionid = {"low": 0, "medium": 1, "high": 2}
pdfid = {"nue": 0, "numu": 1, "nc": 2, "mupi": 3, "rel": 4}
modelid = {"lma": 0, "faild": -3, "malaney": -1, "ksw": -2, "woosley": -4}

# Signal Cherenkov angle fractions
s_ch_frac = [9.433e-04, 9.925e-01, 6.525e-03]

# 3rd reduction efficiencies
effs_sk4 = loadtxt("efficiencies/efficiencies_sk4.txt")
effsk4 = interp1d(effs_sk4[:,0], effs_sk4[:,1], bounds_error=False,
                  fill_value = (effs_sk4[0,1], effs_sk4[-1,1]))
effs_sk3 = loadtxt("efficiencies/efficiencies_sk3.txt")
effsk3 = interp1d(effs_sk3[:,0], effs_sk3[:,1], bounds_error=False,
                  fill_value = (effs_sk3[0,1], effs_sk3[-1,1]))
effs_sk2 = loadtxt("efficiencies/efficiencies_sk2.txt")
effsk2 = interp1d(effs_sk2[:,0], effs_sk2[:,1], bounds_error=False,
                  fill_value = (effs_sk2[0,1], effs_sk2[-1,1]))
effs_sk1 = loadtxt("efficiencies/efficiencies_sk1.txt")
effsk1 = interp1d(effs_sk1[:,0], effs_sk1[:,1], bounds_error=False,
                  fill_value = (effs_sk1[0,1], effs_sk1[-1,1]))
effs_3rdred = [effsk1, effsk2, effsk3, effsk4]

spaeff_sk1 = array([[16, 0.818], [18, 0.908], [24, 1.0]])
spaeff_sk2 = array([[17.5, 0.762], [20, 0.882], [26, 1.0]])
spaeff_sk3 = array([[16, 0.818], [18, 0.908], [24, 1.0]])
spaeff_sk4 = array([[16, 0.76], [18, 0.863], [20, 0.898], [24, 1.0]])
spaeff = [spaeff_sk1, spaeff_sk2, spaeff_sk3, spaeff_sk4]

soleff_sk1 = array([[16, 0.738], [17, 0.821], [18, 0.878],
                    [19, 0.965], [20, 1]])
soleff_sk2 = array([[17.5, 0.738], [18.02, 0.821], [19.08, 0.878],
                    [20.14, 0.965], [21.2, 1]])
soleff_sk3 = array([[16, 0.738], [17, 0.821], [18, 0.878],
                    [19, 0.965], [20, 1]])
soleff_sk4 = array([[16, 0.7226], [17, 0.8138], [18, 0.8718],
                    [19, 0.965], [20, 1]])
soleff = [soleff_sk1, soleff_sk2, soleff_sk3, soleff_sk4]

# SK4 ntag efficiencies
ntag_ebins = [16, 90]
bdt_cuts = [0.620]
emin, emax = ntag_ebins[0], ntag_ebins[-1]
bdt_roc = genfromtxt('roc_curve_N10gt5_cut6_nlow1_newsys.roc')
cuts_roc, roc_effs, roc_bg = bdt_roc[:,0], bdt_roc[:,1], bdt_roc[:,2]
ntag_eff = interp1d(cuts_roc, roc_effs)
ntag_bg = interp1d(cuts_roc, roc_bg)
ntag_effs = ntag_eff(bdt_cuts) #[eff(c) for c in bdt_cuts]
ntag_bgs = ntag_bg(bdt_cuts) # [bg(c) for c in bdt_cuts]
ntag_eff_ps = 0.455 # N10 > 5, 18-523 microsecs window
ntag_bg_ps = 10.999 # N10 > 5, 18-523 microsecs window

# Scalings between Cherenkov angle regions (from MC)
mupi_rescale_low = [1.367, 1.75, 1.34, 1.34] # mupi from low to medium
mupi_rescale_high = [0.12777, 0.1, 0.13, 0.13] # mupi from low to high
nc_rescale = [1.16313, 1.42, 1.14, 1.14] # NC from high to medium

# Neutron multiplicity systematcs
alpha = ones(len(pdfid) - 1) * 0.4


def load_signal_pdf(sknum, model, elow, ehigh, elow_1n):
    ''' Load relic signal model pdfs, either from a discrete (named) model
    or a parametric model. '''

    def seff_sk(en, ntag=False):
        ''' Signal efficiencies (solar + spallation + 3rd red + ntag)
        (No solar cut with ntag) '''
        if (not ntag and en < elow) or (ntag and en < elow_1n):
            return 0
        eff3 = effs_3rdred[sknum - 1](en)
        if ntag:
            if sknum < 4:
                raise ValueError(f"SK-{sknum} does not have ntag")
            spa = spaeff_sk4
            effspa = lambda x: 1.0 if x < spa[0,0] else spa[(spa[:, 0] <= x), -1][-1]
            return eff3 * effspa(en)
        else:
            spa = spaeff[sknum - 1]
            sol = soleff[sknum - 1]
            effspa = lambda x: 1.0 if x < spa[0,0] else spa[(spa[:, 0] <= x), -1][-1]
            effsol = lambda x: 1.0 if x < sol[0,0] else sol[(sol[:, 0] <= x), -1][-1]
            return eff3 * effspa(en) * effsol(en)
        return 0

    def relic(en, spec):
        eff_func_nontag = lambda z: seff_sk(z, ntag=False)
        if sknum < 4:
            return relic_sk(sknum, en, spec, s_ch_frac, eff_func_nontag,
                            elow=elow, ehigh=ehigh)
        else:
            eff_func_ntag = lambda z: seff_sk(z, ntag=True)
            return relic_sk(sknum, en, spec, s_ch_frac, eff_func_nontag,
                            eff_func_ntag, elow=elow, ehigh=ehigh,
                            elow_n=elow_1n, ntag_ebins=ntag_ebins,
                            ntag_effs=ntag_effs, ntag_bgs=ntag_bgs,
                            ntag_eff_ps=ntag_eff_ps, ntag_bg_ps=ntag_bg_ps)

    def spec_flux(flux, low, high, smear=True):
        """ Return spectrum given neutrino flux,
        and positron energy bounds """
        en = arange(low, high + 0.1, 0.1)
        spec0 = array([sns.ibd_spectrum_flux(ee, flux) for ee in en])
        spec, en = spec0[~isnan(spec0)], en[~isnan(spec0)] # remove undefined
        if smear:
            _, spec = sns.smear_ibd_spectrum(en, column_stack((en, spec)), sknum)
        return en, spec

    def spec_params(mod, imf, csfr, low, high, smear=True):
        """ Return spectrum given parametrization,
        and positron energy bounds """
        en = arange(low, high + 0.1, 0.1)
        spec0 = array([sns.ibd_spectrum(ee, mod, imf, csfr) for ee in en])
        spec, en = spec0[~isnan(spec0)], en[~isnan(spec0)] # remove undefined
        if smear:
            _, spec = sns.smear_ibd_spectrum(en, column_stack((en, spec)), sknum)
        return en, spec

    def tot_rate(en_full, spec_full, low, high, ntag=None, nfracs=None):
        en = en_full[(en_full >= low) & (en_full <= high)]
        spec = spec_full[(en_full >= low) & (en_full <= high)]
        if ntag is not None:
            totrate = 0
            for e, sp in zip(en, spec):
                ebin = digitize(e, ntag_ebins) - 1
                nfrac = nfracs[ntag][ebin]
                totrate += nfrac * sp
        else:
            totrate = spec.sum()
        return totrate * (en[1] - en[0])


    print(f"Load model {model}")
    # elow_spec = elow
    # if sknum == 4:
    #     elow_spec = min(elow, elow_1n)
    if ":" not in model:
        # These are discrete models
        flux = loadtxt(f"models/flux_cross_{model}.dat")  # Enu, flux
        fflux = interp1d(flux[:, 0], flux[:, 1], bounds_error=True)
        totflux, _ = quad(fflux, 17.3, 100)
        en, spec_full = spec_flux(flux, 0, 100)
    else:
        print(f"Loading {model}")
        imf = sns.imfs['salpeter']
        csfr = sns.csfr_fiducial
        mod = literal_eval(model)
        totflux, _ = quad(sns.snflux, 17.3, 100, args=(mod, imf, csfr))
        en, spec_full = spec_params(mod, imf, csfr, 0, 100)
    rel = relic(en, spec_full)
    totrate = tot_rate(en, spec_full, elow, ehigh)
    if sknum == 4 and elow != elow_1n:
        if elow_1n < elow:
            totrate_low = tot_rate(en, spec_full, elow_1n, elow,
                                   ntag=True, nfracs=rel.nregions_frac)
        else:
            totrate_low = tot_rate(en, spec_full, elow, elow_1n,
                                   ntag=False, nfracs=rel.nregions_frac)
            totrate = tot_rate(en, spec_full, elow_1n, ehigh)
        totrate = totrate + totrate_low
    flux_fac = totflux / totrate
    return rel, flux_fac


def pdf(energy, sknum, model, elow, pdfid, region,
        ntag=None, backgrounds=None, signal=None):
    ''' PDF function for signal or background in given region. '''
    if sknum < 4:
        if pdfid == 4:
            return likes.pdf(energy, sknum, modelid[model], elow, 4, region)
        #if pdfid == 4: return signal.pdf(energy,region)
        if pdfid in range(4):
            return likes.pdf(energy, sknum, 0, elow, pdfid, region)
    elif sknum == 4:
        assert isinstance(ntag, bool)
        if pdfid == 4:
            return signal.pdf(energy, region, ntag) # to do (for specific srn models)
        elif pdfid in range(4):
            return backgrounds[pdfid].pdf(energy, region, ntag)
        else: raise ValueError("Invalid pdfid")
    else: raise ValueError("Invalid sknum")


def systematics_atm(energies, sknum, model, elow, ehigh, elow_1n=None,
                    energies_n=None, backgrounds=None):
    '''
    Compute distortion functions due to systematics (for atmospheric spectra)
    Must provide energies=[energies_mid, energies_hi] arrays if sk1/2/3.
    If sk4, energies=[lo, mid, hi], and energies_n=[lo_n, mid_n, hi_n]
    for (0 | >1) neutron region and 1 neutron region respectively.
    '''

    def pdfnorm(pdf_id, region_id, ntag=None):
        ''' Calculate PDF norm for given region.
        sk4: if None is given as region, sum over regions
        '''
        def pdf_en(region_id, ntag=None):
            ''' PDF function depending only on energy '''
            return lambda en: pdf(en, sknum, model, elow, pdf_id,
                                  region_id, ntag, backgrounds)

        if sknum == 4:
            if region_id is not None and ntag is not None:
                if ntag:
                    return quad(pdf_en(region_id, ntag), elow_1n, ehigh)[0]
                else:
                    return quad(pdf_en(region_id, ntag), elow, ehigh)[0]
            elif region_id is not None:
                norm = quad(pdf_en(region_id, False), elow, ehigh)[0]
                norm += quad(pdf_en(region_id, True), elow_1n, ehigh)[0]
                return norm
            elif ntag is not None:
                if ntag:
                    norm = [quad(pdf_en(rid, ntag), elow_1n, ehigh)[0]
                            for rid in range(len(regionid))]
                else:
                    norm = [quad(pdf_en(rid, ntag), elow, ehigh)[0]
                            for rid in range(len(regionid))]
                return sum(norm)
            else:
                norm_1n = [quad(pdf_en(rid, True), elow_1n, ehigh)[0]
                           for rid in range(len(regionid))]
                norm_other = [quad(pdf_en(rid, False), elow, ehigh)[0]
                              for rid in range(len(regionid))]
                return sum(norm_1n) + sum(norm_other)
        else:
            return quad(pdf_en(region_id), elow, ehigh)[0]

    def pdfmoment(pdf_id, region_id):
        ''' First moment of PDF for given Cherenkov region '''
        def integrand(ntag=None):
            return lambda en: en * pdf(en, sknum, model, elow, pdf_id,
                                        region_id, ntag, backgrounds)
        if sknum < 4:
            return quad(integrand(), elow, ehigh)[0]
        else:
            moment = quad(integrand(ntag=True), elow_1n, ehigh)[0]
            moment += quad(integrand(ntag=False), elow, ehigh)[0]
            return moment

    # CC distortion sigmas
    sigmas = arange(-1, 3.5, 0.5)
    if sknum < 4:
        assert energies_n is None
        assert len(energies) == 2
        energies_med, energies_high = energies

        # NC distortion sigmas
        sigmas2 = arange(-1, 3.5, 0.5)

        # Normalization and correction factors for nue CC
        norm0 = pdfnorm(pdfid["nue"], regionid["medium"])
        norm1 = pdfmoment(pdfid["nue"], regionid["medium"])
        normnue = 1. / (1 + 0.5 * sigmas * (norm1 / norm0 - 16) / 74)
        nuefact = 1 + 0.5 * sigmas[newaxis,:] * (energies_med[:,newaxis] - 16)/74
        nuefact *= normnue # (Nenergies x Nsigmas)

        # Correction factors for NC
        normncmed = pdfnorm(pdfid["nc"], regionid["medium"])
        normnchigh = pdfnorm(pdfid["nc"], regionid["high"])
        ncfact_med = 1 + sigmas2  # (Nsigmas2)
        ncfact_high = 1 - sigmas2 * normncmed/normnchigh #  (Nsigmas2)

        # make systematics tensors (Nenergies x Npdfs x Nsigmas x Nsigmas2)
        sysmatrix_med = ones((len(energies_med), 5, len(sigmas), len(sigmas2)))
        sysmatrix_high = ones((len(energies_high), 5, len(sigmas), len(sigmas2)))

        sysmatrix_med[:,pdfid["nue"],:,:] = nuefact[:, :, newaxis]
        sysmatrix_med[:,pdfid["nc"],:,:] = ncfact_med[newaxis, newaxis, :]
        sysmatrix_high[:,pdfid["nc"],:,:] = ncfact_high[newaxis, newaxis, :]
        return sysmatrix_med, sysmatrix_high

    else:
        assert energies_n is not None
        assert len(energies) == len(energies_n) == 3
        energies_low, energies_med, energies_high = energies
        energies_low_n, energies_med_n, energies_high_n = energies_n

        # NC and neutron multiplicity distrtion sigmas
        sigmas2 = arange(-2, 2.5, 0.5) / 2. # TODO
        sigmas3 = arange(-2, 3.5, 0.5)

        # Normalization and correction factors for nue CC
        norm0 = pdfnorm(pdfid["nue"], regionid["medium"])
        norm1 = pdfmoment(pdfid["nue"], regionid["medium"])
        normnue = 1. / (1 + 0.5 * sigmas * (norm1 / norm0 - 16) / 74)
        nuefact = 1 + 0.5*sigmas[newaxis,:]*(energies_med[:,newaxis]-16)/74
        nuefact *= normnue # (Nenergies x Nsigmas)
        nuefact_1n = 1+0.5*sigmas[newaxis,:]*(energies_med_n[:,newaxis]-16)/74
        nuefact_1n *= normnue

        # Correction factors for NC
        normncmed = pdfnorm(pdfid["nc"], regionid["medium"])
        normnchigh = pdfnorm(pdfid["nc"], regionid["high"])
        ncfact_med = 1 + sigmas2  # (Nsigmas2)
        ncfact_high = 1 - sigmas2 * normncmed / normnchigh  # (Nsigmas2)
        #ncfact_high = where(ncfact_high < 0, 0, ncfact_high)

        # Neutron multiplicity correction factors
        neutnorm_1n = array([pdfnorm(pid, region_id=None, ntag=True)
                             for pid in range(len(pdfid) - 1)])[:, newaxis]
        neutnorm_other = array([pdfnorm(pid, region_id=None, ntag=False)
                                for pid in range(len(pdfid) - 1)])[:, newaxis]
        alpha_sigma3 = alpha[:, newaxis] * sigmas3[newaxis, :]
        nfact_1n = 1 + alpha_sigma3 # (Npdfs-1 x Nsigmas3)
        nfact_other = 1 - alpha_sigma3 * neutnorm_1n / neutnorm_other

        # make systematics tensors (Nen x Npdfs x Nsig x Nsig2 x Nsig3)
        sys_shape_low = (len(energies_low), len(pdfid),
                                 len(sigmas), len(sigmas2), len(sigmas3))
        sys_shape_low_1n = (len(energies_low_n),) + sys_shape_low[1:]
        sys_shape_high = (len(energies_high),) + sys_shape_low[1:]
        sys_shape_high_1n = (len(energies_high_n),) + sys_shape_low[1:]
        sys_shape_med = (len(energies_med),) + sys_shape_low[1:]
        sys_shape_med_1n = (len(energies_med_n),) + sys_shape_low[1:]

        sysmatrix_low = ones(sys_shape_low)
        sysmatrix_low_1n = ones(sys_shape_low_1n)
        sysmatrix_med = ones(sys_shape_med)
        sysmatrix_med_1n = ones(sys_shape_med_1n)
        sysmatrix_high = ones(sys_shape_high)
        sysmatrix_high_1n = ones(sys_shape_high_1n)

        nuefact = nuefact[:, :, newaxis, newaxis]
        nuefact_1n = nuefact_1n[:, :, newaxis, newaxis]
        ncfact_med = ncfact_med[newaxis, newaxis, :, newaxis]
        ncfact_high = ncfact_high[newaxis, newaxis, :, newaxis]
        nfact_1n = nfact_1n[newaxis, :, newaxis, newaxis, :]
        nfact_other = nfact_other[newaxis, :, newaxis, newaxis, :]

        sysmatrix_low[:, range(len(pdfid)-1), :, :, :] = nfact_other

        sysmatrix_low_1n[:, range(len(pdfid)-1), :, :, :] = nfact_1n

        sysmatrix_med[:, pdfid["nue"], :, :, :] = nuefact
        sysmatrix_med[:, pdfid["nc"], :, :, :] = ncfact_med
        sysmatrix_med[:, range(len(pdfid)-1), :, :, :] *= nfact_other

        sysmatrix_med_1n[:, pdfid["nue"], :, :, :] = nuefact_1n
        sysmatrix_med_1n[:, pdfid["nc"], :, :, :] = ncfact_med
        sysmatrix_med_1n[:, range(len(pdfid)-1), :, :, :] *= nfact_1n

        sysmatrix_high[:, pdfid["nc"], :, :, :] = ncfact_high
        sysmatrix_high[:, range(len(pdfid)-1), :, :, :] *= nfact_other

        sysmatrix_high_1n[:, pdfid["nc"], :, :, :] = ncfact_high
        sysmatrix_high_1n[:, range(len(pdfid)-1), :, :, :] *= nfact_1n

        sysm = sysmatrix_low, sysmatrix_med, sysmatrix_high
        sysm_1n = sysmatrix_low_1n, sysmatrix_med_1n, sysmatrix_high_1n
        return sysm + sysm_1n


def asym_gaussian():
    ''' Asymmetric gaussian for atm systematics weighting. '''
    return array([0.1643, 0.2517, 0.2636, 0.1888, 0.09240,
                  0.03092, 0.007076, 0.001107, 0.0001184])


def getmaxlike(nrelic, nback_ini, pdfs_low, pdfs_med, pdfs_high, sknum, sys=0):
    ''' Maximum likelihood iteration '''

    def get_like_nosys(ncce, nnc, nmupi):
        ''' Likelihood without systematics, to initialize of bg numbers. '''
        ncce = ncce[0]
        ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
        nccmu = ntot - ncce - nnc - nmupi - nrelic
        if ncce < 0 or nccmu < 0:
            return -1e10
        nevents = array([ncce, nccmu] + [nnc, nmupi] + [nrelic])
        totlike = log((nevents * pdfs_med).sum(axis = 1)).sum() - nevents.sum()
        return totlike

    def get_like(nbackgrounds):
        ''' Likelihood with systematics '''
        if nbackgrounds.min() < 0:
            return -1e10
        wgauss = asym_gaussian()
        wgauss2 = asym_gaussian() if sknum < 4 else 0.20997 * exp(-arange(-2,2.5,0.5)**2/2.)
        nevents = array(list(nbackgrounds) + [nrelic])
        totlike = (log(einsum("j,ijkl", nevents, pdfs_high)).sum(axis = 0)
                + log(dot(nevents,pdfs_low.T)).sum(axis = 0)
                + log(einsum("j,ijkl", nevents, pdfs_med)).sum(axis = 0)
                - nrelic - nbackgrounds.sum()) # maybe double counting?
        totmax = totlike.max()
        likenew = log((exp(totlike - totmax) * wgauss[:, newaxis] * wgauss2[newaxis,:]).sum()) + totmax
        return likenew

    funclike = None
    if sys:
        funclike = lambda nback: -get_like(nback)
        maxlike = fmin(funclike, nback_ini, full_output = True, disp = 0)
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    else:
        funclike = lambda nback: -get_like_nosys(nback, nback_ini[1],
                                                 nback_ini[2])
        maxlike = fmin(funclike, nback_ini[0], full_output = True, disp = 0)
        ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
        nccmu = ntot - maxlike[0] - nback_ini[1] - nback_ini[2] - nrelic
        return concatenate([array([maxlike[0][0], nccmu]), nback_ini[1:],
                            array([nrelic, -maxlike[1]])])


def getmaxlike_sk4(nrelic, nback_ini, pdfs, pdfs_1n, sys=0):
    ''' Maximum likelihood iteration '''

    def get_like_nosys_sk4(ncce, nnc, nmupi):
        '''Likelihood without systematics.
        Only used for initialization of background numbers.
        '''
        assert len(pdfs) == len(pdfs_1n) == 3
        pdfs_low, pdfs_med, pdfs_high = pdfs
        pdfs_low_1n, pdfs_med_1n, pdfs_high_1n = pdfs_1n
        ncce = ncce[0]

        ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
        ntot += len(pdfs_low_1n) + len(pdfs_high_1n) + len(pdfs_med_1n)
        nccmu = ntot - ncce - nnc - nmupi - nrelic
        if ncce < 0 or nccmu < 0:
            return -1e10
        nevents = array([ncce, nccmu] + [nnc, nmupi] + [nrelic])
        totlike = log((nevents * pdfs_med).sum(axis = 1)).sum() - nevents.sum()
        return totlike

    def get_like_sk4(nbackgrounds):
        ''' Likelihood with systematics '''
        assert len(pdfs) == len(pdfs_1n) == 3

        pdfs_dist_low, pdfs_dist_med, pdfs_dist_high = [clip(p, 1e-10, None) for p in pdfs]
        pdfs_dist_low_1n, pdfs_dist_med_1n, pdfs_dist_high_1n = [clip(p, 1e-10, None) for p in pdfs_1n]

        if nbackgrounds.min() < 0:
            return -1e10
        wgauss = asym_gaussian()
        wgauss2 = 0.40823 * exp(-arange(-2, 2.5, 0.5)**2 / 2.)
        wgauss3 = 0.40379 * exp(-arange(-2, 3.5, 0.5)**2 / 2.)
        nevents = array(list(nbackgrounds) + [nrelic])
        # testing = pdfs_dist_high
        # if sum(testing <= 0.0) > 0:
        #     raise ValueError("nonpositivity of distorted pdf")

        totlike = (log(einsum("j,ijklm", nevents, pdfs_dist_high)).sum(axis=0)
                + log(einsum("j,ijklm", nevents, pdfs_dist_med)).sum(axis=0)
                + log(einsum("j,ijklm", nevents, pdfs_dist_low)).sum(axis=0)
                + log(einsum("j,ijklm", nevents, pdfs_dist_high_1n)).sum(axis=0)
                + log(einsum("j,ijklm", nevents, pdfs_dist_med_1n)).sum(axis=0)
                + log(einsum("j,ijklm", nevents, pdfs_dist_low_1n)).sum(axis=0)
                - nrelic - nbackgrounds.sum())
        totmax = totlike.max()
        likenew = log((exp(totlike - totmax)
                    * wgauss[:, newaxis, newaxis]
                    * wgauss2[newaxis, :, newaxis]
                    * wgauss3[newaxis, newaxis, :]).sum()) + totmax
        return likenew

    if sys:
        def funclike(nback):
            return -get_like_sk4(nback)
        maxlike = fmin(funclike, nback_ini, full_output=True, disp=0)
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    else:
        def funclike(nback):
            return -get_like_nosys_sk4(nback, nback_ini[1], nback_ini[2])
        maxlike = fmin(funclike, nback_ini[0], full_output=True, disp=0)
        ntot = sum([len(p) for p in pdfs])
        ntot += sum([len(p) for p in pdfs_1n])
        nccmu = ntot - maxlike[0] - nback_ini[1] - nback_ini[2] - nrelic
        return concatenate([array([maxlike[0][0], nccmu[0]]), nback_ini[1:],
                            array([nrelic, -maxlike[1]])])


def analyse(likes, final=False):
    """ Extract limits """
    lmax = likes[:, -1].max()
    bestpos = likes[:, -1].argmax()
    fact = 0.5 if final else 1
    rel = likes[:, -2] * fact
    best = rel[bestpos]
    norm = exp(likes[:, -1] - lmax).sum()
    errminus = best - rel[searchsorted(likes[:bestpos, -1], lmax - 0.5)]
    errplus = rel[len(likes) - 1 - searchsorted(likes[bestpos:, -1][::-1], lmax - 0.5)] - best
    l90 = rel[searchsorted(exp(likes[:, -1] - lmax).cumsum(), 0.9 * norm)]
    return lmax, best, errplus, errminus, l90


def plotfit(nnue, nnumu, nnc, nmupi, nrelic, model, sknum, elow, ehigh, elow_1n,
            samples, samples_n=None, signal=None, background=None):
    """ Plot spectral fit """
    def plotregion(region, data, ax, elow=elow, ntag=None):
        #plt.figure()
        step = 2
        en = arange(elow, ehigh, 0.1)
        nuecc = nnue * array([pdf(ee, sknum, model, elow, pdfid["nue"], region,
                                ntag=ntag, backgrounds=background) for ee in en])
        numucc = nnumu * array([pdf(ee, sknum, model, elow, pdfid["numu"], region,
                                ntag=ntag, backgrounds=background) for ee in en])
        nc = nnc * array([pdf(ee, sknum, model, elow, pdfid["nc"], region,
                            ntag=ntag, backgrounds=background) for ee in en])
        mupi = nmupi * array([pdf(ee, sknum, model, elow, pdfid["mupi"], region,
                                ntag=ntag, backgrounds=background) for ee in en])
        relic = nrelic * array([pdf(ee, sknum, model, elow, 4, region,
                                    ntag=ntag,signal=signal) for ee in en])
        ax.plot(en, step*nuecc, label = r"$\nu_e$ CC")
        ax.plot(en, step*nc, label = "NC")
        ax.plot(en, step*numucc, label = r"$Decay e^-$")
        ax.plot(en, step*mupi, label = r"$\mu/\pi$")
        ax.plot(en, step*relic, label = "relic")
        ax.plot(en, step*(mupi + nc + numucc + nuecc), label = "all background")
        h = histogram(data, bins = arange(elow, ehigh,step))
        x = h[1][1:] - 0.5 * (h[1][1:] - h[1][:-1])
        ax.errorbar(x, h[0], xerr = step/2, yerr = sqrt(h[0]), fmt = '.', color = 'black')
        if ax.is_first_col():
            ax.legend(loc='upper left', prop={'size': 12})
            if ntag is None:
                ax.set(ylabel="Number of events after cuts")
            elif ntag:
                ax.set(ylabel="Number of events after cuts\n1 neutron tag")
            else:
                ax.set(ylabel="Number of events after cuts\n0 / >1 neutron tags")
        if ax.is_first_row():
            titles = [r"Low sideband (20 < $\Theta_C$ < 38)",
                      r"Signal region (38 < $\Theta_C$ < 50)",
                      r"High sideband (78 < $\Theta_C$ < 90)"]
            ax.set(title=titles[region])
        ax.set(xlabel = "E$_p$ (MeV)")

    plt.style.use("seaborn")
    if sknum < 4:
        _, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
        samplow, sampmed, samphigh = samples
        plotregion(1, sampmed, ax2)
        plotregion(2, samphigh, ax3)
        plotregion(0, samplow, ax1)
    else:
        _, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharey=True, figsize=(13, 10))
        samplow, sampmed, samphigh = samples
        samplow_1n, sampmed_1n, samphigh_1n = samples_n
        plotregion(1, sampmed, ax2, elow=elow, ntag=False)
        plotregion(2, samphigh, ax3, elow=elow, ntag=False)
        plotregion(0, samplow, ax1, elow=elow, ntag=False)
        plotregion(1, sampmed_1n, ax5, elow=elow_1n, ntag=True)
        plotregion(2, samphigh_1n, ax6, elow=elow_1n, ntag=True)
        plotregion(0, samplow_1n, ax4, elow=elow_1n, ntag=True)
    plt.subplots_adjust(wspace=0)


def maxlike(sknum, model, elow, ehigh=90, elow_1n=16, rmin=-5, rmax=100,
            rstep=0.1, quiet=True, outdir='.'):
    '''
    Main maximum likelihood function
    sknum = 1,2,3 (SK phase)
    model = SRN model (right now, only "ando")
    elow = energy threshold (here, 16MeV)
    rmin, rmax, rstep = range and step of numbers of relic events for likelihood maximization
    '''

    def load_sample(ntag=False):
        ''' Data samples for SK I-IV '''
        if ntag:
            low = loadtxt("sk{}/ntag/samplelow.txt".format(int(sknum)))[:, 1]
            med = loadtxt("sk{}/ntag/samplemed.txt".format(int(sknum)))[:, 1]
            high = loadtxt("sk{}/ntag/samplehigh.txt".format(int(sknum)))[:, 1]
            low = low[(low > elow_1n) & (low < ehigh)]
            med = med[(med > elow_1n) & (med < ehigh)]
            high = high[(high > elow_1n) & (high < ehigh)]
        else:
            low = loadtxt("sk{}/samplelow.txt".format(int(sknum)))[:, 1]
            med = loadtxt("sk{}/samplemed.txt".format(int(sknum)))[:, 1]
            high = loadtxt("sk{}/samplehigh.txt".format(int(sknum)))[:, 1]
            low = low[(low > elow) & (low < ehigh)]
            med = med[(med > elow) & (med < ehigh)]
            high = high[(high > elow) & (high < ehigh)]
        return low,med,high

    def get_spasolbins(ntag=False):
        ''' Get efficiency steps for background pdfs '''
        bins = []
        effs = []
        if sknum < 4 or not ntag:
            spa = spaeff[sknum- 1]
            sol = soleff[sknum - 1]
            spabins = spa[:, 0]
            solbins = sol[:, 0]
            effspa = lambda x: 1.0 if x < spa[0,0] else spa[(spa[:, 0] <= x), -1][-1]
            effsol = lambda x: 1.0 if x < sol[0,0] else sol[(sol[:, 0] <= x), -1][-1]
            bins = array(list(set(sorted(append(spabins, solbins)))))
            effs = vectorize(effspa)(bins) * vectorize(effsol)(bins)
            bins = list(bins) + [90.]
        elif sknum == 4 and ntag:
            bins = list(spaeff_sk4[:, 0]) + [90.]
            effs = list(spaeff_sk4[:, 1])
        else:
            raise ValueError(f"No search for SK-{sknum} and ntag={ntag}")
        return bins, effs

    def get_pdfmatrix(energies, region, ntag, signal=None, backgrounds=None):
        ''' Get pdfs for different energies, regions, types.
        Output is an array of pdf values for each energy
        and each type of signal/bg. '''
        p = []
        for e in energies:
            p += [[pdf(e, sknum, model, elow, i, regionid[region],
                signal=signal, backgrounds=backgrounds, ntag=ntag)
                for i in range(len(pdfid))]]
        return array(p)  # (Nen x Npdf)

    def initialize(low, med, high, rel):
        ''' Initialization '''
        nevlow = len(low)
        nevhigh = len(high)
        nevmed = len(med)
        nback = nevlow + nevhigh + nevmed - rel

        # Estimate mupi and nc backgrounds using sidebands
        # Fractions are taken from MC
        mupi = nevlow * mupi_rescale_low[sknum - 1]
        nc = (nevhigh - mupi_rescale_high[sknum - 1] * mupi) * nc_rescale[sknum - 1]

        # Maximize likelihoods over backgrounds in signal region
        likemax = getmaxlike(rel, array([nback/5.,nc,mupi]), low, med, high, sknum)
        return likemax

    def initialize_sk4(low, med, high, low_1n, med_1n, high_1n, rel):
        ''' Initialization '''
        nevlow = len(low) + len(low_1n)
        nevhigh = len(high) + len(high_1n)
        nevmed = len(med) + len(med_1n)
        nback = nevlow + nevhigh + nevmed - rel

        # Estimate mupi and nc backgrounds using sidebands
        # Fractions are taken from MC
        mupi = nevlow * mupi_rescale_low[sknum - 1]
        nc = (nevhigh - mupi_rescale_high[sknum - 1] * mupi) * nc_rescale[sknum - 1]

        # Maximize likelihoods over backgrounds in signal region
        likemax = getmaxlike_sk4(rel, array([nback/5.,nc,mupi]), [low, med, high],
                                [low_1n, med_1n, high_1n])
        return likemax

    def applysys(likes, eff, rmin, rmax, rstep):
        ''' Apply gaussian systematic efficiency error correction'''
        print(f"Signal efficiency is {eff}")
        syseff = sys_eff[sknum - 1]
        lower = max(eff * (1 - 6*syseff), 1e-10)
        upper = min(eff * (1 + 6*syseff), 0.999)
        step = (upper - lower)/1000.
        epsrange = arange(lower, upper+step, step)
        if len(epsrange) > 1001:
            epsrange = epsrange[:-1]
        pgaus = exp(-0.5 * (epsrange - eff)**2 / (syseff * eff)**2)
        pgaus /= pgaus.sum()
        # Convolution (Simpson integration)
        lmax = likes[:, -1].max()
        flikes = interp1d(likes[:, -2], exp(likes[:, -1] - lmax), # like(relic)
                          bounds_error=False, fill_value = 0)
        rates = arange(rmin, rmax, rstep)
        lconv = flikes(rates[:, newaxis] * epsrange[newaxis, :] * livetimes[sknum - 1]/365.25* 0.5) # TODO
        simpsoncoeff = array([step/3.] + list((1 + (arange(1,1000)%2))*2./3 * step) + [step/3.])
        ltot = (lconv * (pgaus * epsrange * simpsoncoeff)).sum(axis = 1) # TODO
        likenew = log(ltot) + lmax # TODO
        return column_stack((likes[:, :-1], likenew))

    samplow, sampmed, samphigh = load_sample()
    samples_n = None
    if sknum == 4:
        samplow_n, sampmed_n, samphigh_n = load_sample(ntag=True) # SK-IV ntag samples
        samples_n = [samplow_n, sampmed_n, samphigh_n]

    # Get signal and background spectra
    signal = None
    effsignal = 1.0
    signal, flux_fac = load_signal_pdf(sknum, model, elow, ehigh, elow_1n)
    if ':' not in model:
        flux_fac = fluxfac[model]

    bgs_sk4 = None
    bg_sk4_dir = "./pdf_bg_sk4"
    if sknum < 4:
        effsignal = efftot[model][sknum - 1]
    else:
        effsignal = signal.overall_efficiency_16_90()
        # Load background pdfs
        # WARNING!!! 16 MeV threshold is hardcoded there!
        cut_bins_ntag, cut_effs_ntag = get_spasolbins(ntag=True)
        cut_bins, cut_effs = get_spasolbins(ntag=False)
        bgs_sk4 = [bg_sk4(i, cut_bins, cut_effs,
                   cut_bins_ntag, cut_effs_ntag, bg_sk4_dir,
                   elow, ehigh=ehigh, elow_n=elow_1n) for i in range(4)]
    print(f"Efficiency is {effsignal}")

    # Get pdfs. PDF matrices: (Nen x Npdf)
    pdfs_high = get_pdfmatrix(samphigh, "high", ntag=False,
                              signal=signal, backgrounds=bgs_sk4)
    pdfs_med = get_pdfmatrix(sampmed, "medium", ntag=False,
                              signal=signal, backgrounds=bgs_sk4)
    pdfs_low = get_pdfmatrix(samplow, "low", ntag=False,
                              signal=signal, backgrounds=bgs_sk4)

    if sknum < 4:
        # Get systematic error matrices
        sysmatrices = systematics_atm([sampmed, samphigh], sknum, model,
                                       elow, ehigh, backgrounds=bgs_sk4)
        sysmatrix_med, sysmatrix_high = sysmatrices
        # Distort pdfs: (Nen x Npdfs x Nsigma x Nsigma2 x Nsigma3)
        pdfs_syshigh = pdfs_high[...,newaxis,newaxis] * sysmatrix_high
        pdfs_sysmed = pdfs_med[...,newaxis,newaxis] * sysmatrix_med

        # Set backgrounds (preliminary estimates)
        init = initialize(pdfs_low, pdfs_med, pdfs_high, rmin)
        nue, numu, nc, mupi, _, _ = init
        numu=numu[0]

        # Main maximization loop
        likedata = []
        rmin = 0
        for i,rel in enumerate(arange(rmin, rmax, rstep)):
            likeres = getmaxlike(rel, array([nue, numu, nc, mupi]),
                                 pdfs_low, pdfs_sysmed, pdfs_syshigh,
                                 sknum, sys=1)
            likedata.append(likeres)
            # Update initial values
            [nue, numu, nc, mupi] = list(likedata[-1][:4])
            if i % 100 == 0:
                print("Step {}/1000, like = {}".format(i, likedata[-1][-1]), flush=True)

    else:
        # samplow_n, sampmed_n, samphigh_n = low_ntag, med_ntag, high_ntag
        pdfs_high_n = get_pdfmatrix(samphigh_n, "high", ntag=True,
                                    signal=signal, backgrounds=bgs_sk4)
        pdfs_med_n = get_pdfmatrix(sampmed_n, "medium", ntag=True,
                                    signal=signal, backgrounds=bgs_sk4)
        pdfs_low_n = get_pdfmatrix(samplow_n, "low", ntag=True,
                                    signal=signal, backgrounds=bgs_sk4)

        # Get systematic error matrices
        sysmatrices = systematics_atm([samplow, sampmed, samphigh], sknum,
                                model, elow, ehigh, elow_1n=elow_1n,
                                backgrounds=bgs_sk4,
                                energies_n=[samplow_n, sampmed_n, samphigh_n])
        sysmatrix_low, sysmatrix_med, sysmatrix_high = sysmatrices[:3]
        sysmatrix_low_1n, sysmatrix_med_1n, sysmatrix_high_1n = sysmatrices[3:]

        if sum(pdfs_high <= 0.0) > 0:
            raise ValueError("zeros in pdf")

        # Distort pdfs
        pdfs_syshigh = pdfs_high[...,newaxis,newaxis,newaxis] * sysmatrix_high
        pdfs_sysmed = pdfs_med[...,newaxis,newaxis,newaxis] * sysmatrix_med
        pdfs_syslow = pdfs_low[...,newaxis,newaxis,newaxis] * sysmatrix_low
        pdfs_syshigh_1n = pdfs_high_n[...,newaxis,newaxis,newaxis] * sysmatrix_high_1n
        pdfs_sysmed_1n = pdfs_med_n[...,newaxis,newaxis,newaxis] * sysmatrix_med_1n
        pdfs_syslow_1n = pdfs_low_n[...,newaxis,newaxis,newaxis] * sysmatrix_low_1n

        # Set backgrounds (preliminary estimates)
        init = initialize_sk4(pdfs_low, pdfs_med, pdfs_high,
                              pdfs_low_n, pdfs_med_n, pdfs_high_n, rmin)
        nue, numu, nc, mupi, _, _ = init

        # Main maximization loop
        likedata = []
        rmin = 0
        for i,rel in enumerate(arange(rmin, rmax, rstep)):
            likeres = getmaxlike_sk4(rel, array([nue, numu, nc, mupi]),
                            [pdfs_syslow, pdfs_sysmed, pdfs_syshigh],
                            [pdfs_syslow_1n, pdfs_sysmed_1n, pdfs_syshigh_1n],
                            sys=1)
            likedata.append(likeres)
            # Update initial values
            [nue, numu, nc, mupi] = list(likedata[-1][:4])
            if i % 100 == 0:
                print("Step {}/1000, like = {}".format(i, likedata[-1][-1]), flush=True)

    results = column_stack((arange(rmin, rmax, rstep), likedata))
    results = results[results[:, 0] >= 0]   # results[i] = rel, nback, rel, like

    # Systematic efficiency error correction + limits
    _, best2, errplus2, errminus2, limit2 = analyse(results)
    results_sys = applysys(results, effsignal, rmin, rmax, rstep)
    _, best, errplus, errminus, limit = analyse(results_sys, final=True)

    flux_90cl = limit * flux_fac
    lpos = results_sys[:, -1].argmax()

    # Save and display results
    savetxt(f"{outdir}/fit_sk{sknum}.txt", column_stack((results, results_sys[:, -1])))
    print("SK-{}".format(sknum), "Best fit:")
    print("{} +{} -{} relic evts/yr".format(best, errplus, errminus))
    print("{} +{} -{} relic evts".format(best2, errplus2, errminus2))
    print("90% c.l. relic event rate: {} ev/yr {}".format(limit, limit2))
    print("90% c.l. {} /cm^2/s > 17.3 MeV".format(flux_90cl))
    print("nu-e events {}".format(results_sys[lpos, 1]))
    print("nu-mu events {}".format(results_sys[lpos, 2]))
    print("NC elastic events {}".format(results_sys[lpos, 3]))
    print("mu/pi events {}".format(results_sys[lpos, 4]))
    print("Max likelihood {}".format(results_sys[lpos, -1]))
    if sknum < 4:
        ev_nums = len(samplow), len(sampmed), len(samphigh)
        print("Numbers of events: {} {} {}".format(*ev_nums))
    else:
        ev_nums = len(samplow), len(sampmed), len(samphigh)
        ev_nums_1n = len(samplow_n), len(sampmed_n), len(samphigh_n)
        print("Numbers of events (!=1 neutrons): {} {} {}".format(*ev_nums))
        print("Numbers of events (1 neutron): {} {} {}".format(*ev_nums_1n))
    if not quiet:
        plotfit(results_sys[lpos,1], results_sys[lpos,2], results_sys[lpos,3],
                results_sys[lpos,4], results_sys[lpos,5], model, sknum,
                elow, ehigh, elow_1n, samples=[samplow, sampmed, samphigh],
                samples_n=samples_n, signal=signal, background=bgs_sk4)
        plt.savefig(f"{outdir}/fit_sk{sknum}.pdf")
        plt.clf()
    return limit * flux_fac, flux_fac, results_sys


def combine(results):
    """ Combine SKI-IV likelihoods """
    liketot = results[0][:, -1] - results[0][:, -1].max()
    for r in results[1:]:
        liketot += r[:, -1] - r[:, -1].max()
    return analyse(column_stack((results[0][:,:-1], liketot)), final=True)


def combine_fluxes(results, fluxfacts):
    """ Combined flux limits """
    flux_sampling = arange(0, 50.1, 0.1)
    rels, likes = results[0][:, 0], results[0][:, -1]
    fluxes = rels * fluxfacts[0]
    flike = interp1d(fluxes, likes, bounds_error=False, fill_value=1e-10)
    newlikes = flike(flux_sampling)
    liketot = newlikes - newlikes.max()
    for i, r in enumerate(results[1:]):
        rels, likes = r[:, 0], r[:, -1]
        fluxes = rels * fluxfacts[i + 1]
        flike = interp1d(fluxes, likes, bounds_error=False, fill_value=1e-10)
        newlikes = flike(flux_sampling)
        liketot += newlikes - newlikes.max()
    return analyse(column_stack((flux_sampling, liketot)))


def fullike(model, elow, ehigh, elow_sk4=None, ehigh_sk4=None, elow_sk4_1n=None,
            rmin=-5, rmax=100, rstep=0.1, quiet=False, outdir='.'):
    """ Fit SK I-IV data """
    if elow_sk4 is None:
        elow_sk4 = elow
    if elow_sk4_1n is None:
        elow_sk4_1n = elow
    if ehigh_sk4 is None:
        ehigh_sk4 = ehigh
    like1 = maxlike(1, model, elow, ehigh, elow_sk4_1n, rmin, rmax,
                    rstep, quiet=quiet, outdir=outdir)
    like2 = maxlike(2, model, elow, ehigh, elow_sk4_1n, rmin, rmax,
                    rstep, quiet=quiet, outdir=outdir)
    like3 = maxlike(3, model, elow, ehigh, elow_sk4_1n, rmin, rmax,
                    rstep, quiet=quiet, outdir=outdir)
    like4 = maxlike(4, model, elow_sk4, ehigh_sk4, elow_sk4_1n, rmin, rmax,
                    rstep, quiet=quiet, outdir=outdir)
    fluxlims = [like1[0], like2[0], like3[0], like4[0]]
    fluxfacs = [like1[1], like2[1], like3[1], like4[1]]
    results = [like1[2], like2[2], like3[2], like4[2]]
    ratelims = array(fluxlims) / array(fluxfacs)
    res = combine(results)
    res_fluxes = combine_fluxes(results, fluxfacs)
    _, ratebest_comb, ratepl_comb, ratemin_comb, ratelim_comb = res
    _, fluxbest_comb, fluxpl_comb, fluxmin_comb, fluxlim_comb = res_fluxes

    if not quiet:
        plt.style.use("seaborn")
        plt.figure()
        plt.xlabel("SRN events/year")
        plt.ylabel("Likelihood")
        x = like1[-1][:, 0]/2
        plt.plot(x, like1[-1][:, -1] - like1[-1][:,-1].max(),
                 label="SK-I", alpha=0.5)
        plt.plot(x, like2[-1][:, -1] - like2[-1][:,-1].max(),
                 '--', label="SK-II", alpha=0.5)
        plt.plot(x, like3[-1][:, -1] - like3[-1][:,-1].max(),
                 '-.', label="SK-III", alpha=0.5)
        plt.plot(x, like4[-1][:, -1] - like4[-1][:,-1].max(),
                 '--', label="SK-IV",linewidth=2)
        likesum = like1[-1][:, -1] + like2[-1][:, -1] + like3[-1][:, -1] + like4[-1][:, -1]
        likesum -= likesum.max()
        plt.plot(x, likesum, label = "Combined", color = 'black')
        plt.plot([0,20], -0.5 * ones(2), 'r')
        plt.xlim(0,20)
        plt.ylim(-2,0.2)
        # plt.grid()
        plt.legend()
        plt.savefig(outdir + "/full_like.pdf")
        plt.clf()

    print(f"Best fit rate: {ratebest_comb} + {ratepl_comb} - {ratemin_comb}")
    print("90%% c.l. rate limits are: %f %f %f %f > 16 MeV" % tuple(ratelims))
    print(f"90%% c.l. combined rate limit: {ratelim_comb} > 16 MeV")

    print(f"Best fit flux: {fluxbest_comb} + {fluxpl_comb} - {fluxmin_comb}")
    print("90%% c.l. flux limits are: %f %f %f %f /cm^2/s > 17.3 MeV" % tuple(fluxlims))
    print(f"90%% c.l. combined flux limit: {fluxlim_comb} /cm^2/s > 17.3 MeV")


if __name__ == "__main__":
    modelname = argv[1]
    directory = argv[2]

    fullike(modelname, elow=16, ehigh=90, elow_sk4=20,
            ehigh_sk4=80, elow_sk4_1n=16, outdir=directory)
