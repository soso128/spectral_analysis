''' Spectral fits for SRN analysis of SK I-IV '''
from __future__ import division
from sys import path
from ast import literal_eval
from time import time
import argparse
import pickle
import os

from numpy import *
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import fmin
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec

path.append("spectrum_generator/")
import snspectrum as sns
import get_flux_from_mc as gtmc
import likes
from pdf_sk4 import bg_sk4, relic_sk, spall_sk
from esys_scale_res import esys
from esys_pdf_distorsions import get_distorsion_functions

plt.style.use(["seaborn","myplot.mplstyle"])
plt.rcParams["xtick.major.size"]=  5
plt.rcParams["ytick.major.size"]=  5
plt.rcParams["xtick.minor.size"]=  5
plt.rcParams["ytick.minor.size"]=  5
plt.rcParams["legend.fontsize"]=  16
plt.rcParams["legend.markerscale"]=  1.5
plt.rcParams["legend.frameon"]=  False
plt.rcParams["grid.alpha"]= 1.0
plt.rcParams['mathtext.default']='regular'

# livetimes
livetimes = array([1497.44, 793.71, 562.04, 2970.1])
aft_eff = 0.94
livetimes[3] *= aft_eff

# energy-independent efficiency sys
efftot = {"lma": [0.7975,0.56703,0.77969],
          "malaney": [0.7903, 0.53273, 0.766262],
          "faild": [0.7997,0.5845,0.7831] ,
          "woosley": [0.7876,0.5524,0.7764],
          "ksw" : [0.7908, 0.54527, 0.77371]}
fluxfac = {"lma": 0.535115, "malaney": 0.5845, "ksw": 0.488413,
           "faild": 0.431, "woosley": 0.47447}
# Systematics for signal: reduction efficiency, livetime, Strumia-Vissani cross-section, FV
sys_eff = array([0.0254, 0.0404, 0.0253, 0.0220])
sys_livetime = 0.001
sys_xsec = 0.01
sys_fv = array([0.013, 0.011, 0.010, 0.015])
sys_eff = sqrt(sys_eff**2 + sys_livetime**2 + sys_xsec**2 + sys_fv**2)
# Ntag cut systematics from AmBe study for SK-IV
sys_eff_sk4_ntag = 0.125

# Definitions
regionids = {"low": 0, "medium": 1, "high": 2}
pdfids = {"nue": 0, "numu": 1, "nc": 2, "mupi": 3, "spall": 4, "rel": 5}
modelids = {"lma": 0, "faild": -3, "malaney": -1, "ksw": -2, "woosley": -4}

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

# Spallation cut efficiencies
spaeff_sk1 = array([[16, 0.818], [18, 0.908], [24, 1.0]])
spaeff_sk2 = array([[17.5, 0.762], [20, 0.882], [26, 1.0]])
spaeff_sk3 = array([[16, 0.818], [18, 0.908], [24, 1.0]])
spaeff_sk4_nontag = array([[12, 0.8], [16, 0.65], [18, 0.63], [20, 0.918], [24, 0.98]])
spaeff_sk4 = array([[16, 0.826], [18, 0.887], [20, 0.918], [24, 0.98]])
spaeff = [spaeff_sk1, spaeff_sk2, spaeff_sk3, spaeff_sk4_nontag]

soleff_sk1 = array([[16, 0.738], [17, 0.821], [18, 0.878],
                    [19, 0.965], [20, 1]])
soleff_sk2 = array([[17.5, 0.738], [18.02, 0.821], [19.08, 0.878],
                    [20.14, 0.965], [21.2, 1]])
soleff_sk3 = array([[16, 0.738], [17, 0.821], [18, 0.878],
                    [19, 0.965], [20, 1]])
soleff_sk4 = array([[16, 0.731], [17, 0.822], [18, 0.883],
                    [19, 0.966], [20, 1]])
soleff = [soleff_sk1, soleff_sk2, soleff_sk3, soleff_sk4]

# SK4 ntag efficiencies
ntag_ebins = [12,14,16,18,20,22,24,26,28,90]
bdt_cuts = [0.958,0.874,0.964,0.855,0.129,0.129,0.241,0.669,0.620]
emin, emax = ntag_ebins[0], ntag_ebins[-1]
bdt_roc = genfromtxt('ROCs/roc_curve_N10gt5_cut6_nlow1.roc')
cuts_roc, roc_effs, roc_bg = bdt_roc[:,0], bdt_roc[:,1], bdt_roc[:,2]
ntag_eff = interp1d(cuts_roc, roc_effs)
ntag_bg = interp1d(cuts_roc, roc_bg)
ntag_effs = ntag_eff(bdt_cuts) #[eff(c) for c in bdt_cuts]
ntag_bgs = ntag_bg(bdt_cuts) # [bg(c) for c in bdt_cuts]
ntag_eff_ps = 0.447 # N10 > 5, 18-523 microsecs window
ntag_bg_ps = 7.05 # N10 > 5, 18-523 microsecs window

# Distorsion coefficients for spallation systematics
spacoeffs_sk = array([[0.0105935649, -0.495682897, 7.93209842, -43.5523139],
                     [0.0138280665, -0.749631175, 13.6659053, -83.7150281],
                     [0.0438680847, -2.13974596, 35.1046340, -193.623584],
                     [0.0103067195, -0.475915126, 7.52829430, -40.9820538]])
#spacoeffs_sk = array([[0.0123007834, -0.579943766, 9.33104800, -51.3586273],
                     #[0.0122423073, -0.665453703, 12.2007772, -75.3689846],
                     #[0.0401672378, -1.94512194, 31.7340679, -174.369394],
                     #[0.0121024891, -0.565396545, 9.02737234, -49.4163994]])
#spacoeffs_sk = array([[0.0268777, -1.29806656, 21.162808, -116.5129],
                     #[0.0135476, -0.7159959, 12.80619, -77.4774],
                     #[0.0291259, -1.374243, 21.968258, -119.050657],
                     #[0.0197376, -0.932184, 14.950667, -81.556668]])

# Dictionaries for energy scale and resolution systematics
# For SK >= 4 (others are hardcoded, from K. Bays 2012 analysis)
esys_scale_dict = {}
esys_res_dict = {}

# Scalings between Cherenkov angle regions (from MC)
mupi_rescale_low = [1.367, 1.75, 1.34, 1.34] # mupi from low to medium
mupi_rescale_high = [0.12777, 0.1, 0.13, 0.13] # mupi from low to high
nc_rescale = [1.16313, 1.42, 1.14, 1.14] # NC from high to medium

# systematics multipliers
nc_mult = 1
cc_mult = 1
n_mult = 1

# Neutron multiplicity systematcs
alpha = ones(len(pdfids) - 2) * 0.4 * n_mult

def load_signal_pdf(sknum, model, elow, ehigh, elow_1n, ineff_scale=1.0):
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
            tot_seff = eff3 * effspa(en)
            return tot_seff + (1-tot_seff)*(1-ineff_scale)
        else:
            spa = spaeff[sknum - 1]
            sol = soleff[sknum - 1]
            effspa = lambda x: 1.0 if x < spa[0,0] else spa[(spa[:, 0] <= x), -1][-1]
            effsol = lambda x: 1.0 if x < sol[0,0] else sol[(sol[:, 0] <= x), -1][-1]
            tot_seff = eff3 * effspa(en) * effsol(en)
            return tot_seff + (1-tot_seff)*(1-ineff_scale)
        return 0

    def relic(en, spec):
        eff_func_nontag = lambda z: seff_sk(z, ntag=False)
        if sknum < 4:
            return relic_sk(sknum, en, spec, s_ch_frac, eff_func_nontag,
                            elow=elow, ehigh=ehigh)
        else:
            # print("ntag_effs:", ntag_effs * ntag_eff_ps)
            # print("ntag_bgs:", ntag_bgs * ntag_bg_ps)
            # print("ntag_eff_ps:", ntag_eff_ps)
            # print("ntag_bg_ps:", ntag_bg_ps)
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
            specfunc = gtmc.smeared_spec_from_mc(en, spec, sknum, high)
            spec = specfunc(en)
            #_, spec = sns.smear_ibd_spectrum(en, column_stack((en, spec)), sknum)
        return en, spec

    def spec_params(mod, imf, csfr, low, high, smear=True):
        """ Return spectrum given parametrization,
        and positron energy bounds """
        en = arange(low, high + 0.1, 0.1)
        spec0 = array([sns.ibd_spectrum(ee, mod, imf, csfr) for ee in en])
        spec, en = spec0[~isnan(spec0)], en[~isnan(spec0)] # remove undefined
        if smear:
            specfunc = gtmc.smeared_spec_from_mc(en, spec, sknum, high)
            spec = specfunc(en)
            #_, spec = sns.smear_ibd_spectrum(en, column_stack((en, spec)), sknum)
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
    if ":" not in model:
        # These are discrete models
        flux = loadtxt(f"models/flux_cross_{model}.dat")  # Enu, flux
        fflux = interp1d(flux[:, 0], flux[:, 1], bounds_error=False, fill_value=0)
        totflux, _ = quad(fflux, 17.3, 100)
        en, spec_full = spec_flux(flux, 0, 100)
    else:
        print(f"Loading {model}")
        imf = sns.imfs['salpeter']
        csfr = sns.csfr_fiducial
        mod = literal_eval(model)
        print(mod, imf, csfr)
        totflux, _ = quad(sns.snflux, 17.3, 100, args=(mod, imf, csfr))
        en, spec_full = spec_params(mod, imf, csfr, 0, 100)
    rel = relic(en, spec_full)
    totrate = tot_rate(en, spec_full, 16, 90)
    flux_fac = totflux / totrate
    return rel, flux_fac, totrate, totflux


def skgd_params(gd_cap_frac, srn_lo=1.0, srn_hi=535.0):
    """ Update fit parameters for fitting against SK-Gd projection
    Concentration can be 0.1 or 0.01 (%) """
    # Variables to update
    global ntag_eff_ps, ntag_effs, ntag_bg_ps, ntag_bgs

    # SK-Gd ntag parameters
    # if concentration == 0.01:
    #     gd_cap_frac = 0.5
    # elif concentration == 0.1:
    #     gd_cap_frac = 0.9
    # else:
    #     raise ValueError("Gd concentration must be 0.01 or 0.1")
    # gd_livetime = 10 * 365.25
    # atm_eff = 1.0 #0.4
    h2o_eff_ps_all = 0.59 #0.562 # Pure water presel. efficiency before time window
    gd_eff_ps_all = 0.995 # Pure Gd presel. efficiency before time window
    ntag_bg_ps_sk6 = 54.1
    n_tau = 204.8  # neutron capture characteristic time, microsecs
    n_tau_gd = 35.  # Gd ncap time
    # srn_lo, srn_hi = 1., 535. # SK-Gd SRN analysis time window
    
    h2o_scale = exp(-srn_lo / n_tau) - exp(-srn_hi / n_tau)
    gd_scale = exp(-srn_lo / n_tau_gd) - exp(-srn_hi  /  n_tau_gd)
    h_eff_ps = h2o_scale * h2o_eff_ps_all
    gd_eff_ps = gd_scale * gd_eff_ps_all
    # gdmix_eff_h2o = (1 - gd_cap_frac) * h2o_scale * h2o_eff_ps_all
    # gdmix_eff_gd = gd_cap_frac * gd_scale * gd_eff_ps_all
    # gdmix_eff_ps = gdmix_eff_h2o + gdmix_eff_gd
    # Neutron tagging BDT SK-Gd efficiencies
    # bdt_roc_gd = genfromtxt('ROCs/bdt22_n6_skgd_1.roc', delimiter=', ')
    # bdt_roc_gd001 = genfromtxt('ROCs/bdt22_n6_skgd_001.roc', delimiter=', ')
    
    h_rocfile="/disk02/usr6/giampaol/ntag-mva/models/bdt22_skg4_0.013_10M/roc_test_H.csv"
    gd_rocfile="/disk02/usr6/giampaol/ntag-mva/models/bdt22_skg4_0.013_10M/roc_test_Gd.csv"
    bdt_roc_h = genfromtxt(h_rocfile, delimiter=', ')
    bdt_roc_gd = genfromtxt(gd_rocfile, delimiter=', ')
    
    # Get efficiency on H captures
    cuts_roc_h, roc_effs_h, roc_bg_h = bdt_roc_h[:,0], bdt_roc_h[:,1], bdt_roc_h[:,2]
    roc_bg_h *= ntag_bg_ps_sk6
    roc_effs_h *= h_eff_ps
    cut_from_bg_h = interp1d(roc_bg_h, cuts_roc_h)
    effh = interp1d(cuts_roc_h, roc_effs_h)
    bdt_cuts_h = cut_from_bg_h(ntag_bgs * ntag_bg_ps)
    ntag_effs_h = effh(bdt_cuts_h)

    # Get efficiency on Gd captures
    cuts_roc_gd, roc_effs_gd, roc_bg_gd = bdt_roc_gd[:,0], bdt_roc_gd[:,1], bdt_roc_gd[:,2]
    roc_bg_gd *= ntag_bg_ps_sk6
    roc_effs_gd *= gd_eff_ps
    cut_from_bg_gd = interp1d(roc_bg_gd, cuts_roc_gd)
    effgd = interp1d(cuts_roc_gd, roc_effs_gd)
    bdt_cuts_gd = cut_from_bg_gd(ntag_bgs * ntag_bg_ps)
    ntag_effs_gd = effgd(bdt_cuts_gd)

    # # 0.01% Gd BDT
    # cuts_roc_gd001, roc_effs_gd001, roc_bg_gd001 = bdt_roc_gd001[:,0], bdt_roc_gd001[:,1], bdt_roc_gd001[:,2]
    # cut_from_bg_gd001 = interp1d(roc_bg_gd001, cuts_roc_gd001)
    # eff001 = interp1d(cuts_roc_gd001, roc_effs_gd001)
    # bdt_cuts_gd001 = cut_from_bg_gd001(ntag_bgs)
    # ntag_effs_gd001 = eff001(bdt_cuts_gd001)
    # if gd_cap_frac == 0.5: # 0.01% Gd
    #     ntag_effs_gd = ntag_effs_gd001
    # elif gd_cap_frac == 0.9:
    #     ntag_effs_gd = ntag_effs_gd * 0.9
    # else:
    #     raise ValueError("Ncap fraction must be either 0.5 or 0.9")

    ntag_eff_ps_gdmix = (1-gd_cap_frac) * h_eff_ps + gd_cap_frac * gd_eff_ps
    ntag_effs_gdmix = gd_cap_frac * ntag_effs_gd  + (1-gd_cap_frac) * ntag_effs_h
    ntag_effs_gdmix /= ntag_eff_ps_gdmix
    print(f"Using Gd neutron capture fraction of {gd_cap_frac}%")
    # print(f"SK4 BDT cuts (eff, bg): {bdt_cuts} ({ntag_effs * ntag_eff_ps}, {ntag_bgs * ntag_bg_ps})")
    # print(f"SK6 BDT cuts (eff on H): {bdt_cuts_h} ({ntag_effs_h})")
    # print(f"SK6 BDT cuts (eff on Gd): {bdt_cuts_gd} ({ntag_effs_gd})")
    # print(f"SK6 BDT eff on Gd mix: {ntag_effs_gdmix * ntag_eff_ps_gdmix}")

    # Update global variables
    ntag_eff_ps = ntag_eff_ps_gdmix
    ntag_effs = ntag_effs_gdmix
    ntag_bgs = ntag_bgs * ntag_bg_ps / ntag_bg_ps_sk6
    ntag_bg_ps = ntag_bg_ps_sk6
    # livetimes[3] = gd_livetime
    # print(f"Checking new ntag efficiencies: {ntag_effs*ntag_eff_ps}")
    # print(f"Checking new ntag bg (PS): {ntag_bgs*ntag_bg_ps} ({ntag_bg_ps})")


def ntag_rescale(ntag_bins, ntag_effs, ntag_bgs, ntag_bg_ps):
    " Pdf rescalings for SK-Gd ntag efficiency "

    def ntag_weight(E, ncap, ntag=True):
        ebin = digitize(E, ntag_bins)
        ef = ntag_effs[ebin-1]
        bgrate = ntag_bgs[ebin-1]

        # Probability of tagging exactly 1 true neutron
        prob_1n = ncap * ef * (1 - ef)**(ncap - 1)
        # Probability of mistagging exactly 1 accidental
        prob_1b = ntag_bg_ps * bgrate * (1 - bgrate)**(ntag_bg_ps - 1)
        prob_0n = (1 - ef)**(ncap)
        prob_0b = (1 - bgrate)**(ntag_bg_ps)

        weight_1n = prob_1n * prob_0b + prob_1b * prob_0n
        if ntag:  # 1-neutron region
            return weight_1n
        else:  # 0 or multiple neutrons region
            return 1 - weight_1n

    with open("atm_mc_info/energies_ntag.p", "rb") as fl:
        E_mupi, E_NC, E_CC_nue, E_decaye_mc, E_decaye = pickle.load(fl)
    with open("atm_mc_info/wgt_ntag.p", "rb") as fl:
        wgt_mupi, wgt_NC, wgt_CC_nue, wgt_decaye_mc = pickle.load(fl)
    with open("atm_mc_info/ncaps_ntag.p", "rb") as fl:
        n_mupi, n_NC, n_CC_nue, n_decaye_mc = pickle.load(fl)
    enes_arr = [E_CC_nue, E_decaye_mc, E_NC, E_mupi, E_decaye]
    wgts_arr = [wgt_CC_nue, wgt_decaye_mc, wgt_NC, wgt_mupi]
    ncaps_arr = [n_CC_nue, n_decaye_mc, n_NC, n_mupi]

    num_ntag_bins = len(ntag_bins) - 1
    ntag_wgt_h2o = zeros((4, 2, num_ntag_bins)) # Water
    ntag_wgt_rescaled = zeros((4, 2, num_ntag_bins)) # Gd (using nominal BDT effs)
    for bg_pdf in range(4):
        for ntag_id, ntw in enumerate(wgts_arr[bg_pdf]):
            for angle_id, weights in enumerate(ntw):
                wt = array(weights)
                en = array(enes_arr[bg_pdf][ntag_id][angle_id])
                nc = array(ncaps_arr[bg_pdf][ntag_id][angle_id])

                binid_ntag = digitize(en, ntag_bins) - 1
                for bid in range(num_ntag_bins):
                    bin_wt = wt[binid_ntag == bid]
                    bin_en = en[binid_ntag == bid]
                    bin_nc = nc[binid_ntag == bid]

                    ntag_wgt_h2o[bg_pdf, ntag_id, bid] += sum(bin_wt)

                    wgt_nt0 = bin_wt * ntag_weight(bin_en,bin_nc,False)
                    wgt_nt1 = bin_wt * ntag_weight(bin_en,bin_nc,True)
                    ntag_wgt_rescaled[bg_pdf, 0, bid] += sum(wgt_nt0)
                    ntag_wgt_rescaled[bg_pdf, 1, bid] += sum(wgt_nt1)

    ntag_rescaling = ntag_wgt_rescaled / ntag_wgt_h2o

    return ntag_rescaling


def pdf(energy, sknum, model, elow, pdfid, region,
        ntag=None, backgrounds=None, signal=None):
    ''' PDF function for signal or background in given region. '''
    numbkgnospa = 4
    if sknum < 4:
        if pdfid == pdfids['rel']:
            # try:
            # return likes.pdf(energy, sknum, modelids[model], elow, 4, region)
            # except KeyError:  # Models not used by Kirk
            return signal.pdf(energy, region) # Always use our own SRN pdf
        if pdfid in range(numbkgnospa):
            return likes.pdf(energy, sknum, 0, elow, pdfid, region)
    elif sknum == 4:
        assert isinstance(ntag, bool)
        if pdfid == pdfids['rel']:
            return signal.pdf(energy, region, ntag) # to do (for specific srn models)
        elif pdfid in range(numbkgnospa):
            return backgrounds[pdfid].pdf(energy, region, ntag)
        elif pdfid != pdfids['spall']: raise ValueError("Invalid pdfid")
    # Add spallation backgrounds
    if pdfid == pdfids['spall']:
        return backgrounds[pdfids['spall']].pdf(energy,region,ntag)
    else: raise ValueError("Invalid sknum")

def systematics_escale_res(energies, sknum, model, elow, ehigh, elow_1n=None,
                    energies_n=None, backgrounds=None, signal=None, use_spall = False):

    def pdfnorm(pdf_id, region_id, ntag=None):
        ''' Calculate PDF norm for given region.
        sk4: if None is given as region, sum over regions
        '''
        def pdf_en(region_id, ntag=None):
            ''' PDF function depending only on energy '''
            ''' times systematics size '''
            return lambda en: pdf(en, sknum, model, elow, pdf_id,
                                  region_id, ntag, backgrounds, signal) * esys(
                                      en, sknum, region_id, pdf_id, esys_scale_dict, esys_res_dict, ntag)

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
                            for rid in range(len(regionids))]
                else:
                    norm = [quad(pdf_en(rid, ntag), elow, ehigh)[0]
                            for rid in range(len(regionids))]
                return sum(norm)
            else:
                norm_1n = [quad(pdf_en(rid, True), elow_1n, ehigh)[0]
                           for rid in range(len(regionids))]
                norm_other = [quad(pdf_en(rid, False), elow, ehigh)[0]
                              for rid in range(len(regionids))]
                return sum(norm_1n) + sum(norm_other)
        else:
            if region_id is not None:
                return quad(pdf_en(region_id), elow, ehigh)[0]
            else:
                return sum(quad(pdf_en(rid, False), elow, ehigh)[0]
                              for rid in range(len(regionids)))

    # Get energies
    energies_low, energies_med, energies_high = energies
    energies_low_1n, energies_med_1n, energies_high_1n = None, None, None

    # make systematics tensors (Nen x Npdfs x Nsig)
    sigmas = arange(-4,4.5,0.5)
    numbkg = 5 if use_spall else 4
    sys_shape_low = (len(energies_low), numbkg + 1, len(sigmas))
    sys_shape_high = (len(energies_high),) + sys_shape_low[1:]
    sys_shape_med = (len(energies_med),) + sys_shape_low[1:]

    sysmatrix_low = ones(sys_shape_low)
    sysmatrix_med = ones(sys_shape_med)
    sysmatrix_high = ones(sys_shape_high)
    sysmatrix_low_1n = None
    sysmatrix_med_1n = None
    sysmatrix_high_1n = None

    # Distorsion factors
    nuenorm = pdfnorm(pdfids['nue'],None)
    ncnorm = pdfnorm(pdfids['nc'],None)
    mupinorm = pdfnorm(pdfids['mupi'],None)
    relicnorm = pdfnorm(pdfids['rel'],None)
    relicfact_low = (1 + esys(energies_low,sknum,0,pdfids['rel'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + relicnorm * sigmas[newaxis, :])
    relicfact_med = (1 + esys(energies_med,sknum,1,pdfids['rel'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + relicnorm * sigmas[newaxis, :])
    relicfact_high = (1 + esys(energies_high,sknum,2,pdfids['rel'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + relicnorm * sigmas[newaxis, :])
    nuefact_low = (1 + esys(energies_low,sknum,0,pdfids['nue'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + nuenorm * sigmas[newaxis, :])
    nuefact_med = (1 + esys(energies_med,sknum,1,pdfids['nue'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + nuenorm * sigmas[newaxis, :])
    nuefact_high = (1 + esys(energies_high,sknum,2,pdfids['nue'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + nuenorm * sigmas[newaxis, :])
    ncfact_low = (1 + esys(energies_low,sknum,0,pdfids['nc'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + ncnorm * sigmas[newaxis, :])
    ncfact_med = (1 + esys(energies_med,sknum,1,pdfids['nc'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + ncnorm * sigmas[newaxis, :])
    ncfact_high = (1 + esys(energies_high,sknum,2,pdfids['nc'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + ncnorm * sigmas[newaxis, :])
    mupifact_low = (1 + esys(energies_low,sknum,0,pdfids['mupi'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + mupinorm * sigmas[newaxis, :])
    mupifact_med = (1 + esys(energies_med,sknum,1,pdfids['mupi'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + mupinorm * sigmas[newaxis, :])
    mupifact_high = (1 + esys(energies_high,sknum,2,pdfids['mupi'], esys_scale_dict, esys_res_dict)[:, newaxis] * sigmas[newaxis, :])/(1 + mupinorm * sigmas[newaxis, :])

    #relicfact = relicfact[:, :, newaxis]
    #nuefact = nuefact[:, :, newaxis]
    #ncfact_med = ncfact_med[:, :, newaxis]
    #ncfact_high = ncfact_high[:, :, newaxis]
    #mupifact_low = mupifact_low[:, :, newaxis]
    #mupifact_med = mupifact_med[:, :, newaxis]
    #mupifact_high = mupifact_high[:, :, newaxis]

    relic_column = pdfids['rel'] if use_spall else pdfids['rel'] - 1
    sysmatrix_med[:, pdfids["nue"], :] = nuefact_med
    sysmatrix_med[:, pdfids["nc"], :] = ncfact_med
    sysmatrix_med[:, pdfids["mupi"], :] = mupifact_med
    sysmatrix_med[:, relic_column, :] = relicfact_med

    sysmatrix_high[:, pdfids["nue"], :] = nuefact_high
    sysmatrix_high[:, pdfids["nc"], :] = ncfact_high
    sysmatrix_high[:, pdfids["mupi"], :] = mupifact_high
    sysmatrix_high[:, relic_column, :] = relicfact_high

    sysmatrix_low[:, pdfids["nue"], :] = nuefact_low
    sysmatrix_low[:, pdfids["nc"], :] = ncfact_low
    sysmatrix_low[:, pdfids["mupi"], :] = mupifact_low
    sysmatrix_low[:, relic_column, :] = relicfact_low

    if sknum == 4:
        assert energies_n is not None
        assert len(energies) == len(energies_n) == 3
        energies_low_n, energies_med_n, energies_high_n = energies_n  

        sys_shape_low_1n = (len(energies_low_n), numbkg + 1, len(sigmas))
        sys_shape_high_1n = (len(energies_high_n),) + sys_shape_low[1:]
        sys_shape_med_1n = (len(energies_med_n),) + sys_shape_low[1:]

        sysmatrix_low_1n = ones(sys_shape_low_1n)
        sysmatrix_med_1n = ones(sys_shape_med_1n)
        sysmatrix_high_1n = ones(sys_shape_high_1n)

        relicfact_low_1n = (1 + esys(energies_low_n,sknum,0,pdfids['rel'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + relicnorm * sigmas[newaxis, :])
        relicfact_med_1n = (1 + esys(energies_med_n,sknum,1,pdfids['rel'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + relicnorm * sigmas[newaxis, :])
        relicfact_high_1n = (1 + esys(energies_high_n,sknum,2,pdfids['rel'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + relicnorm * sigmas[newaxis, :])
        nuefact_low_1n = (1 + esys(energies_low_n,sknum,0,pdfids['nue'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + nuenorm * sigmas[newaxis, :])
        nuefact_med_1n = (1 + esys(energies_med_n,sknum,1,pdfids['nue'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + nuenorm * sigmas[newaxis, :])
        nuefact_high_1n = (1 + esys(energies_high_n,sknum,2,pdfids['nue'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + nuenorm * sigmas[newaxis, :])
        ncfact_low_1n = (1 + esys(energies_low_n,sknum,0,pdfids['nc'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + ncnorm * sigmas[newaxis, :])
        ncfact_med_1n = (1 + esys(energies_med_n,sknum,1,pdfids['nc'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + ncnorm * sigmas[newaxis, :])
        ncfact_high_1n = (1 + esys(energies_high_n,sknum,2,pdfids['nc'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + ncnorm * sigmas[newaxis, :])
        mupifact_low_1n = (1 + esys(energies_low_n,sknum,0,pdfids['mupi'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + mupinorm * sigmas[newaxis, :])
        mupifact_med_1n = (1 + esys(energies_med_n,sknum,1,pdfids['mupi'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + mupinorm * sigmas[newaxis, :])
        mupifact_high_1n = (1 + esys(energies_high_n,sknum,2,pdfids['mupi'], esys_scale_dict, esys_res_dict, ntag = True)[:, newaxis] * sigmas[newaxis, :])/(1 + mupinorm * sigmas[newaxis, :])

        #relicfact_1n = relicfact_1n[:, :, newaxis]
        #nuefact_1n = nuefact_1n[:, :, newaxis]
        #ncfact_1n_med = ncfact_1n_med[:, :, newaxis]
        #ncfact_1n_high = ncfact_1n_high[:, :, newaxis]
        #mupifact_1n_low = mupifact_1n_low[:, :, newaxis]
        #mupifact_1n_med = mupifact_1n_med[:, :, newaxis]
        #mupifact_1n_high = mupifact_1n_high[:, :, newaxis]

        relic_column = pdfids['rel'] if use_spall else pdfids['rel'] - 1
        sysmatrix_med_1n[:, pdfids["nue"], :] = nuefact_med_1n
        sysmatrix_med_1n[:, pdfids["nc"], :] = ncfact_med_1n
        sysmatrix_med_1n[:, pdfids["mupi"], :] = mupifact_med_1n
        sysmatrix_med_1n[:, relic_column, :] = relicfact_med_1n

        sysmatrix_high_1n[:, pdfids["nue"], :] = nuefact_high_1n
        sysmatrix_high_1n[:, pdfids["nc"], :] = ncfact_high_1n
        sysmatrix_high_1n[:, pdfids["mupi"], :] = mupifact_high_1n
        sysmatrix_high_1n[:, relic_column, :] = relicfact_high_1n

        sysmatrix_low_1n[:, pdfids["nue"], :] = nuefact_low_1n
        sysmatrix_low_1n[:, pdfids["nc"], :] = ncfact_low_1n
        sysmatrix_low_1n[:, pdfids["mupi"], :] = mupifact_low_1n
        sysmatrix_low_1n[:, relic_column, :] = relicfact_low_1n

    sysm = sysmatrix_low, sysmatrix_med, sysmatrix_high
    sysm_1n = sysmatrix_low_1n, sysmatrix_med_1n, sysmatrix_high_1n
    return sysm + sysm_1n

def systematics_mupi(energies, sknum, model, elow, ehigh, elow_1n=None,
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
                            for rid in range(len(regionids))]
                else:
                    norm = [quad(pdf_en(rid, ntag), elow, ehigh)[0]
                            for rid in range(len(regionids))]
                return sum(norm)
            else:
                norm_1n = [quad(pdf_en(rid, True), elow_1n, ehigh)[0]
                           for rid in range(len(regionids))]
                norm_other = [quad(pdf_en(rid, False), elow, ehigh)[0]
                              for rid in range(len(regionids))]
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
    assert energies_n is not None
    assert len(energies) == len(energies_n) == 3
    energies_low, energies_med, energies_high = energies
    energies_low_n, energies_med_n, energies_high_n = energies_n

    # Normalization and correction factors for mupi
    norm0 = pdfnorm(pdfids["mupi"], regionids["low"])
    norm1 = pdfmoment(pdfids["mupi"], regionids["low"])
    normmupi = 1. / (1 + 0.5 * sigmas * cc_mult * (80 - norm1 / norm0) / 60)
    mupifact = 1 + 0.5*sigmas[newaxis,:] * cc_mult * (80 - energies_low[:,newaxis])/60
    mupifact *= normmupi # (Nenergies x Nsigmas)
    mupifact_1n = 1 + 0.5 * sigmas[newaxis,:] * cc_mult * (80 - energies_low_n[:,newaxis])/60
    mupifact_1n *= normmupi

    # make systematics tensors (Nen x Npdfs x Nsig)
    numbkg = 5 if use_spall else 4
    sys_shape_low = (len(energies_low), numbkg + 1, len(sigmas))
    sys_shape_low_1n = (len(energies_low_n),) + sys_shape_low[1:]

    sysmatrix_low = ones(sys_shape_low)
    sysmatrix_low_1n = ones(sys_shape_low_1n)

    sysmatrix_low[:, pdfids["mupi"], :] = mupifact

    sysmatrix_low_1n[:, pdfids["mupi"], :] = mupifact_1n

    return sysmatrix_low,sysmatrix_low_1n

def systematics_atm(energies, sknum, model, elow, ehigh, elow_1n=None,
                    energies_n=None, backgrounds=None, use_spall=False, no_nc=False):
    '''
    Compute distortion functions due to systematics (for atmospheric spectra)
    Must provide energies=[energies_mid, energies_hi] arrays if sk1/2/3.
    If sk4, energies=[lo, mid, hi], and energies_n=[lo_n, mid_n, hi_n]
    for (0 | >1) neutron region and 1 neutron region respectively.
    '''

    def pdfnorm(pdf_id, region_id, eup=ehigh, ntag=None):
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
                    return quad(pdf_en(region_id, ntag), elow_1n, eup)[0]
                else:
                    return quad(pdf_en(region_id, ntag), elow, eup)[0]
            elif region_id is not None:
                norm = quad(pdf_en(region_id, False), elow, eup)[0]
                norm += quad(pdf_en(region_id, True), elow_1n, eup)[0]
                return norm
            elif ntag is not None:
                if ntag:
                    norm = [quad(pdf_en(rid, ntag), elow_1n, eup)[0]
                            for rid in range(len(regionids))]
                else:
                    norm = [quad(pdf_en(rid, ntag), elow, eup)[0]
                            for rid in range(len(regionids))]
                return sum(norm)
            else:
                norm_1n = [quad(pdf_en(rid, True), elow_1n, eup)[0]
                           for rid in range(len(regionids))]
                norm_other = [quad(pdf_en(rid, False), elow, eup)[0]
                              for rid in range(len(regionids))]
                return sum(norm_1n) + sum(norm_other)
        else:
            return quad(pdf_en(region_id), elow, eup)[0]

    def pdfmoment(pdf_id, region_id, eup = ehigh, order = 1):
        ''' First moment of PDF for given Cherenkov region '''
        def integrand(ntag=None):
            return lambda en: en**order * pdf(en, sknum, model, elow, pdf_id,
                                        region_id, ntag, backgrounds)
        if sknum < 4:
            return quad(integrand(), elow, eup)[0]
        else:
            moment = quad(integrand(ntag=True), elow_1n, eup)[0]
            moment += quad(integrand(ntag=False), elow, eup)[0]
            return moment

    # CC distortion sigmas
    sigmas = arange(-1, 3.5, 0.5)
    numbkg = 5 if use_spall else 4
    #spacoeffs = [0.02517, -1.69954 , 43.08233, -485.46063, 2050.978 - 1]
    #spacoeffs = [0.0198, -0.4659 , 2.2375]
    spacoeffs = spacoeffs_sk[sknum - 1]
    if sknum < 4:
        assert energies_n is None
        assert len(energies) == 2
        energies_med, energies_high = energies

        # NC distortion sigmas
        sigmas2 = arange(-1, 3.5, 0.5)
        sigmas3 = arange(-2, 2.5, 0.5)

        # Normalization and correction factors for nue CC
        norm0 = pdfnorm(pdfids["nue"], regionids["medium"])
        norm1 = pdfmoment(pdfids["nue"], regionids["medium"])
        normnue = 1. / (1 + 0.5 * sigmas * cc_mult * (norm1 / norm0 - 16) / 74)
        nuefact = 1 + 0.5 * sigmas[newaxis,:] * cc_mult * (energies_med[:,newaxis] - 16)/74
        nuefact *= normnue # (Nenergies x Nsigmas)

        # Correction factors for NC
        normncmed = pdfnorm(pdfids["nc"], regionids["medium"])
        normnchigh = pdfnorm(pdfids["nc"], regionids["high"])
        ncfact_med = 1 + sigmas2 * nc_mult # (Nsigmas2)
        ncfact_high = 1 - sigmas2 * nc_mult * normncmed/normnchigh #  (Nsigmas2)

        # Correction factors for spallation
        if use_spall:
            sigmas3 = arange(-2, 2.5, 0.5)
            spafact = (1 + sigmas3[newaxis,:] * (spacoeffs[0] * energies_med[:, newaxis]**3
                                                    + spacoeffs[1] * energies_med[:, newaxis]**2 
                                                           + spacoeffs[2] * energies_med[:, newaxis] + spacoeffs[3]))
            max_pos_energies = []
            for sf in spafact.T:
                if (sf < 0).sum():
                    en = sorted(energies_med[sf < 0])[0]
                    max_pos_energies.append(sorted(energies_med[energies_med < en])[-1])
                else:
                    max_pos_energies.append(ehigh)
            max_pos_energies = array(max_pos_energies)
            #max_pos_energies = array([sorted(energies_med[energies_med < sorted(energies_med[sf < 0])[0]])[-1] for sf in spafact.T])
            norm0 = array([pdfnorm(pdfids['spall'], regionids['medium'], eup = mpe) for mpe in max_pos_energies])
            norm1 = array([pdfmoment(pdfids['spall'], regionids['medium'], eup = mpe) for mpe in max_pos_energies])
            norm2 = array([pdfmoment(pdfids['spall'], regionids['medium'], eup = mpe, order = 2) for mpe in max_pos_energies])
            norm3 = array([pdfmoment(pdfids['spall'], regionids['medium'], eup = mpe, order = 3) for mpe in max_pos_energies])
            normspa = 1./(norm0 + sigmas3 * (spacoeffs[0] * norm3 + spacoeffs[1] * norm2 + spacoeffs[2] * norm1 + spacoeffs[3] * norm0))
            spafact = where(spafact < 0, 1e-10, spafact) * normspa # Cut off pdfs when they go negative

            # make systematics tensors (Nenergies x Npdfs x Nsigmas x Nsigmas2)
            sysmatrix_med = ones((len(energies_med), numbkg + 1, len(sigmas), len(sigmas2), len(sigmas3)))
            sysmatrix_high = ones((len(energies_high), numbkg + 1, len(sigmas), len(sigmas2), len(sigmas3)))

            sysmatrix_med[:,pdfids["nue"],:,:,:] = nuefact[:, :, newaxis, newaxis]
            sysmatrix_med[:,pdfids["nc"],:,:,:] = ncfact_med[newaxis, newaxis, :, newaxis]
            sysmatrix_med[:,pdfids["spall"],:,:,:] = spafact[:, newaxis, newaxis, :]
            sysmatrix_high[:,pdfids["nc"],:,:,:] = ncfact_high[newaxis, newaxis, :, newaxis]
            return sysmatrix_med, sysmatrix_high
        else:
            # make systematics tensors (Nenergies x Npdfs x Nsigmas x Nsigmas2)
            sysmatrix_med = ones((len(energies_med), numbkg + 1, len(sigmas), len(sigmas2)))
            sysmatrix_high = ones((len(energies_high), numbkg + 1, len(sigmas), len(sigmas2)))

            sysmatrix_med[:,pdfids["nue"],:,:] = nuefact[:, :, newaxis]
            sysmatrix_med[:,pdfids["nc"],:,:] = ncfact_med[newaxis, newaxis, :]
            sysmatrix_high[:,pdfids["nc"],:,:] = ncfact_high[newaxis, newaxis, :]
            return sysmatrix_med, sysmatrix_high

    else:
        assert energies_n is not None
        assert len(energies) == len(energies_n) == 3
        energies_low, energies_med, energies_high = energies
        energies_low_n, energies_med_n, energies_high_n = energies_n

        # NC and neutron multiplicity distrtion sigmas
        sigmas2 = arange(-2, 2.5, 0.5) / 2.
        sigmas3 = arange(-2, 3.5, 0.5)

        # Normalization and correction factors for nue CC
        norm0 = pdfnorm(pdfids["nue"], regionids["medium"])
        norm1 = pdfmoment(pdfids["nue"], regionids["medium"])
        normnue = 1. / (1 + 0.5 * sigmas * cc_mult * (norm1 / norm0 - 16) / 74)
        nuefact = 1 + 0.5*sigmas[newaxis,:] * cc_mult * (energies_med[:,newaxis]-16)/74
        nuefact *= normnue # (Nenergies x Nsigmas)
        nuefact_1n = 1 + 0.5 * sigmas[newaxis,:] * cc_mult * (energies_med_n[:,newaxis]-16)/74
        nuefact_1n *= normnue

        # Correction factors for NC
        normncmed = pdfnorm(pdfids["nc"], regionids["medium"])
        normnchigh = pdfnorm(pdfids["nc"], regionids["high"])
        ncfact_med = 1 + sigmas2 * nc_mult  # (Nsigmas2)
        ncfact_high = 1 - sigmas2 * nc_mult * normncmed / normnchigh  # (Nsigmas2)
        #ncfact_high = where(ncfact_high < 0, 0, ncfact_high)

        # Neutron multiplicity correction factors
        neutnorm_1n = array([pdfnorm(pid, region_id=None, ntag=True)
                             for pid in range(len(pdfids) - 2)])[:, newaxis]
        neutnorm_other = array([pdfnorm(pid, region_id=None, ntag=False)
                                for pid in range(len(pdfids) - 2)])[:, newaxis]
        alpha_sigma3 = alpha[:, newaxis] * sigmas3[newaxis, :]
        nfact_1n = 1 + alpha_sigma3 # (Npdfs-1 x Nsigmas3)
        nfact_other = 1 - alpha_sigma3 * neutnorm_1n / neutnorm_other

        if use_spall:
            sigmas4 = arange(-2, 2.5, 0.5)
            #norm0 = pdfnorm(pdfids['spall'], regionids['medium'])
            #norm1 = pdfmoment(pdfids['spall'], regionids['medium'])
            #norm2 = pdfmoment(pdfids['spall'], regionids['medium'], order = 2)
            #normspa = 1./(1 + sigmas4 * (spacoeffs[0] * norm2 - spacoeffs[1] * norm1 + spacoeffs[2] * norm0))
            #epoly =(spacoeffs[0] * energies_med**2 - spacoeffs[1] * energies_med + spacoeffs[2]) 
            #spafact = normspa * (1 + sigmas4[newaxis,:] * epoly[:, newaxis])
            spafact = (1 + sigmas4[newaxis,:] * (spacoeffs[0] * energies_med[:, newaxis]**3
                                                    + spacoeffs[1] * energies_med[:, newaxis]**2 
                                                           + spacoeffs[2] * energies_med[:, newaxis] + spacoeffs[3]))
            max_pos_energies = []
            for sf in spafact.T:
                if (sf < 0).sum():
                    en = sorted(energies_med[sf < 0])[0]
                    max_pos_energies.append(sorted(energies_med[energies_med < en])[-1])
                else:
                    max_pos_energies.append(ehigh)
            max_pos_energies = array(max_pos_energies)
            #max_pos_energies = array([sorted(energies_med[energies_med < sorted(energies_med[sf < 0])[0]])[-1] if (sf < 0).sum()
                                      #else ehigh for sf in spafact.T])
            norm0 = array([pdfnorm(pdfids['spall'], regionids['medium'], eup = mpe) for mpe in max_pos_energies])
            norm1 = array([pdfmoment(pdfids['spall'], regionids['medium'], eup = mpe) for mpe in max_pos_energies])
            norm2 = array([pdfmoment(pdfids['spall'], regionids['medium'], eup = mpe, order = 2) for mpe in max_pos_energies])
            norm3 = array([pdfmoment(pdfids['spall'], regionids['medium'], eup = mpe, order = 3) for mpe in max_pos_energies])
            normspa = 1./(norm0 + sigmas4 * (spacoeffs[0] * norm3 + spacoeffs[1] * norm2 + spacoeffs[2] * norm1 + spacoeffs[3] * norm0))
            spafact = where(spafact < 0, 1e-10, spafact)*normspa # Cut off pdfs when they go negative

            # make systematics tensors (Nen x Npdfs x Nsig x Nsig2 x Nsig3 x Nsig4)
            sys_shape_low = (len(energies_low), numbkg + 1,
                                     len(sigmas), len(sigmas2), len(sigmas3), len(sigmas4))
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

            nuefact = nuefact[:, :, newaxis, newaxis, newaxis]
            nuefact_1n = nuefact_1n[:, :, newaxis, newaxis, newaxis]
            ncfact_med = ncfact_med[newaxis, newaxis, :, newaxis, newaxis]
            ncfact_high = ncfact_high[newaxis, newaxis, :, newaxis, newaxis]
            nfact_1n = nfact_1n[newaxis, :, newaxis, newaxis, :, newaxis]
            nfact_other = nfact_other[newaxis, :, newaxis, newaxis, :, newaxis]
            spafact = spafact[:, newaxis, newaxis, newaxis, :]

            sysmatrix_low[:, range(numbkg-1), :, :, :, :] = nfact_other

            sysmatrix_low_1n[:, range(numbkg-1), :, :, :, :] = nfact_1n

            sysmatrix_med[:, pdfids["nue"], :, :, :, :] = nuefact
            sysmatrix_med[:, pdfids["nc"], :, :, :, :] = ncfact_med
            sysmatrix_med[:, range(numbkg-1), :, :, :, :] *= nfact_other
            sysmatrix_med[:, pdfids["spall"], :, :, :, :] *= spafact

            sysmatrix_med_1n[:, pdfids["nue"], :, :, :, :] = nuefact_1n
            sysmatrix_med_1n[:, pdfids["nc"], :, :, :, :] = ncfact_med
            sysmatrix_med_1n[:, range(numbkg-1), :, :, :, :] *= nfact_1n

            sysmatrix_high[:, pdfids["nc"], :, :, :, :] = ncfact_high
            sysmatrix_high[:, range(numbkg-1), :, :, :, :] *= nfact_other

            sysmatrix_high_1n[:, pdfids["nc"], :, :, :, :] = ncfact_high
            sysmatrix_high_1n[:, range(numbkg-1), :, :, :, :] *= nfact_1n

            sysm = sysmatrix_low, sysmatrix_med, sysmatrix_high
            sysm_1n = sysmatrix_low_1n, sysmatrix_med_1n, sysmatrix_high_1n
            return sysm + sysm_1n
        elif no_nc:
            # make systematics tensors (Nen x Npdfs x Nsig x Nsig2)
            sys_shape_low = (len(energies_low), numbkg + 1,
                                     len(sigmas), len(sigmas3))
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

            nuefact = nuefact[:, :, newaxis]
            nuefact_1n = nuefact_1n[:, :, newaxis]
            nfact_1n = nfact_1n[newaxis, :, newaxis, :]
            nfact_other = nfact_other[newaxis, :, newaxis, :]

            sysmatrix_low[:, range(numbkg), :, :] = nfact_other
            sysmatrix_low_1n[:, range(numbkg), :, :] = nfact_1n

            sysmatrix_med[:, pdfids["nue"], :, :] = nuefact
            sysmatrix_med[:, range(numbkg), :, :] *= nfact_other
            sysmatrix_med_1n[:, pdfids["nue"], :, :] = nuefact_1n
            sysmatrix_med_1n[:, range(numbkg), :, :] *= nfact_1n

            sysmatrix_high[:, range(numbkg), :, :] *= nfact_other
            sysmatrix_high_1n[:, range(numbkg), :, :] *= nfact_1n

            sysm = sysmatrix_low, sysmatrix_med, sysmatrix_high
            sysm_1n = sysmatrix_low_1n, sysmatrix_med_1n, sysmatrix_high_1n
            return sysm + sysm_1n
            
        else:
            # make systematics tensors (Nen x Npdfs x Nsig x Nsig2 x Nsig3)
            sys_shape_low = (len(energies_low), numbkg + 1,
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

            sysmatrix_low[:, range(numbkg), :, :, :] = nfact_other

            sysmatrix_low_1n[:, range(numbkg), :, :, :] = nfact_1n

            sysmatrix_med[:, pdfids["nue"], :, :, :] = nuefact
            sysmatrix_med[:, pdfids["nc"], :, :, :] = ncfact_med
            sysmatrix_med[:, range(numbkg), :, :, :] *= nfact_other

            sysmatrix_med_1n[:, pdfids["nue"], :, :, :] = nuefact_1n
            sysmatrix_med_1n[:, pdfids["nc"], :, :, :] = ncfact_med
            sysmatrix_med_1n[:, range(numbkg), :, :, :] *= nfact_1n

            sysmatrix_high[:, pdfids["nc"], :, :, :] = ncfact_high
            sysmatrix_high[:, range(numbkg), :, :, :] *= nfact_other

            sysmatrix_high_1n[:, pdfids["nc"], :, :, :] = ncfact_high
            sysmatrix_high_1n[:, range(numbkg), :, :, :] *= nfact_1n

            sysm = sysmatrix_low, sysmatrix_med, sysmatrix_high
            sysm_1n = sysmatrix_low_1n, sysmatrix_med_1n, sysmatrix_high_1n
            return sysm + sysm_1n


def asym_gaussian():
    ''' Asymmetric gaussian for atm systematics weighting. '''
    return array([0.1643, 0.2517, 0.2636, 0.1888, 0.09240,
                  0.03092, 0.007076, 0.001107, 0.0001184])


def getmaxlike(nrelic, nback_ini, pdfs_low, pdfs_med, pdfs_high, sknum, sys=0, use_spall = False):
    ''' Maximum likelihood iteration '''

    def get_like_init(ncce, *nbkgs):
        ''' Likelihood without systematics, to initialize of bg numbers. '''
        ncce = ncce[0]
        ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
        nccmu = ntot - ncce - sum(nbkgs) - nrelic
        if ncce < 0 or nccmu < 0:
            return -1e10
        nevents = array([ncce, nccmu] + list(nbkgs) + [nrelic])
        totlike = log((nevents * pdfs_med).sum(axis = 1)).sum() - nevents.sum()
        return totlike

    def get_like_nosys(nbackgrounds):
        ''' Likelihood with systematics on atm spectral shapes'''
        if nbackgrounds.min() < 0:
            return -1e10
        nevents = array(list(nbackgrounds) + [nrelic])
        totlike = (log(dot(nevents,pdfs_low.T)).sum(axis = 0)
                   + log(dot(nevents,pdfs_med.T)).sum(axis = 0)
                   + log(dot(nevents,pdfs_high.T)).sum(axis = 0)
                   - nrelic - nbackgrounds.sum()) # maybe double counting?
        if isnan(totlike): raise ValueError("nan")
        return totlike

    def get_like(nbackgrounds):
        ''' Likelihood with systematics on atm spectral shapes'''
        if nbackgrounds.min() < 0:
            return -1e10
        wgauss = asym_gaussian()
        wgauss2 = asym_gaussian() if sknum < 4 else 0.20997 * exp(-arange(-2,2.5,0.5)**2/2.)
        nevents = array(list(nbackgrounds) + [nrelic])
        indexstring = "j,ijklm" if use_spall else "j,ijkl"
        totlike = (log(einsum(indexstring, nevents, pdfs_high)).sum(axis = 0)
                   + log(dot(nevents,pdfs_low.T)).sum(axis = 0)
                   + log(einsum(indexstring, nevents, pdfs_med)).sum(axis = 0)
                   - nrelic - nbackgrounds.sum()) # maybe double counting?
        totmax = totlike.max()
        if use_spall:
            likenew = log((exp(totlike - totmax) * wgauss[:,newaxis,newaxis] * wgauss2[newaxis,:,newaxis] * wgauss2[newaxis,newaxis,:]).sum()) + totmax
        else:
            likenew = log((exp(totlike - totmax) * wgauss[:,newaxis] * wgauss2[newaxis,:]).sum()) + totmax
        return likenew

    def get_like_esys(nbackgrounds):
        ''' Likelihood with systematics on energy scale and resolution'''
        if nbackgrounds.min() < 0:
            return -1e10
        nevents = array(list(nbackgrounds) + [nrelic])
        totlike = (log(einsum("j,ijk", nevents, pdfs_high)).sum(axis = 0)
                   + log(einsum("j,ijk", nevents, pdfs_med)).sum(axis = 0)
                   + log(einsum("j,ijk", nevents, pdfs_low)).sum(axis = 0)
                   - nrelic - nbackgrounds.sum()) # maybe double counting?
        totmax = totlike.max()
        gauss = exp(-arange(-4,4.5,0.5)**2/2)
        gauss /= gauss.sum()
        likenew = log((exp(totlike - totmax) * gauss).sum()) + totmax
        return likenew

    funclike = None
    if sys == 1:
        funclike = lambda nback: -get_like(nback)
        maxlike = fmin(funclike, nback_ini, full_output = True, disp = 0)
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    if sys == 2:
        funclike = lambda nback: -get_like_esys(nback)
        maxlike = fmin(funclike, nback_ini, full_output = True, disp = 0)
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    if sys == 0:
        funclike = lambda nback: -get_like_nosys(nback)
        maxlike = fmin(funclike, nback_ini, full_output = True, disp = 0)
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    if sys == -1:
        funclike = lambda nback: -get_like_init(nback, *nback_ini[1:])
        maxlike = fmin(funclike, nback_ini[0], full_output = True, disp = 0)
        ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
        nccmu = ntot - maxlike[0] - sum(nback_ini) - nrelic
        return concatenate([array([maxlike[0][0], nccmu]), nback_ini[1:],
                            array([nrelic, -maxlike[1]])])


def getmaxlike_sk4(nrelic, nback_ini, pdfs, pdfs_1n, sys=0,
                   use_spall=False, no_nc=False):
    ''' Maximum likelihood iteration '''

    def get_like_init_sk4(ncce, *nbkgs):
        '''Likelihood without systematics.
        Only used for initialization of background numbers.
        '''
        assert len(pdfs) == len(pdfs_1n) == 3
        pdfs_low, pdfs_med, pdfs_high = pdfs
        pdfs_low_1n, pdfs_med_1n, pdfs_high_1n = pdfs_1n
        ncce = ncce[0]

        ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
        ntot += len(pdfs_low_1n) + len(pdfs_high_1n) + len(pdfs_med_1n)
        nccmu = ntot - ncce - sum(nbkgs) - nrelic
        if ncce < 0 or nccmu < 0:
            return -1e10
        nevents = array([ncce, nccmu] + list(nbkgs) + [nrelic])
        totlike = log((nevents * pdfs_med).sum(axis = 1)).sum() - nevents.sum()
        return totlike

    def get_like_nosys_sk4(nbackgrounds):
        ''' Likelihood without systematics '''
        assert len(pdfs) == len(pdfs_1n) == 3

        pdfs_dist_low, pdfs_dist_med, pdfs_dist_high = [clip(p, 1e-10, None) for p in pdfs]
        pdfs_dist_low_1n, pdfs_dist_med_1n, pdfs_dist_high_1n = [clip(p, 1e-10, None) for p in pdfs_1n]

        if nbackgrounds.min() < 0:
            return -1e10
        nevents = array(list(nbackgrounds) + [nrelic])

        totlike = (log(einsum("j,ij", nevents, pdfs_dist_high)).sum(axis=0)
                   + log(einsum("j,ij", nevents, pdfs_dist_med)).sum(axis=0)
                   + log(einsum("j,ij", nevents, pdfs_dist_low)).sum(axis=0)
                   + log(einsum("j,ij", nevents, pdfs_dist_high_1n)).sum(axis=0)
                   + log(einsum("j,ij", nevents, pdfs_dist_med_1n)).sum(axis=0)
                   + log(einsum("j,ij", nevents, pdfs_dist_low_1n)).sum(axis=0)
                   - nrelic - nbackgrounds.sum())
        return totlike

    def get_like_sk4(nbackgrounds):
        ''' Likelihood with systematics '''
        assert len(pdfs) == len(pdfs_1n) == 3
        if no_nc: nbackgrounds[pdfids['nc']] = 0

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
        
        if use_spall: indexstring = "j,ijklmn"
        elif no_nc: indexstring = "j,ijkl"
        else: indexstring = "j,ijklm"

        totlike = (log(einsum(indexstring, nevents, pdfs_dist_high)).sum(axis=0)
                + log(einsum(indexstring, nevents, pdfs_dist_med)).sum(axis=0)
                + log(einsum(indexstring, nevents, pdfs_dist_low)).sum(axis=0)
                + log(einsum(indexstring, nevents, pdfs_dist_high_1n)).sum(axis=0)
                + log(einsum(indexstring, nevents, pdfs_dist_med_1n)).sum(axis=0)
                + log(einsum(indexstring, nevents, pdfs_dist_low_1n)).sum(axis=0)
                - nrelic - nbackgrounds.sum())
        totmax = totlike.max()
        if use_spall:
            likenew = log((exp(totlike - totmax)
                        * wgauss[:, newaxis, newaxis, newaxis]
                        * wgauss2[newaxis, :, newaxis, newaxis]
                        * wgauss3[newaxis, newaxis, :, newaxis]
                        * wgauss2[newaxis, newaxis, newaxis, :]).sum()) + totmax
        elif no_nc:
            likenew = log((exp(totlike - totmax)
                        * wgauss[:, newaxis]
                        * wgauss3[newaxis, :]).sum()) + totmax
        else:
            likenew = log((exp(totlike - totmax)
                        * wgauss[:, newaxis, newaxis]
                        * wgauss2[newaxis, :, newaxis]
                        * wgauss3[newaxis, newaxis, :]).sum()) + totmax
        return likenew

    def get_like_mupi_sk4(nbackgrounds):
        ''' Likelihood with systematics on atm spectral shapes'''
        assert len(pdfs) == len(pdfs_1n) == 3
        pdfs_low, pdfs_med, pdfs_high = [clip(p, 1e-10, None) for p in pdfs]
        pdfs_low_1n, pdfs_med_1n, pdfs_high_1n = [clip(p, 1e-10, None) for p in pdfs_1n]
        if nbackgrounds.min() < 0:
            return -1e10
        wgauss = asym_gaussian()
        nevents = array(list(nbackgrounds) + [nrelic])
        totlike = (log(einsum("j,ijk", nevents, pdfs_low)).sum(axis = 0)
                   + log(einsum("j,ijk", nevents, pdfs_low_1n)).sum(axis = 0)
                   + log(dot(nevents,pdfs_med.T)).sum(axis = 0)
                   + log(dot(nevents,pdfs_med_1n.T)).sum(axis = 0)
                   + log(dot(nevents,pdfs_high.T)).sum(axis = 0)
                   + log(dot(nevents,pdfs_high_1n.T)).sum(axis = 0)
                   - nrelic - nbackgrounds.sum()) # maybe double counting?
        totmax = totlike.max()
        likenew = log((exp(totlike - totmax) * wgauss).sum()) + totmax
        return likenew

    def get_like_esys_sk4(nbackgrounds):
        ''' Likelihood with systematics '''
        assert len(pdfs) == len(pdfs_1n) == 3

        pdfs_dist_low, pdfs_dist_med, pdfs_dist_high = [clip(p, 1e-10, None) for p in pdfs]
        pdfs_dist_low_1n, pdfs_dist_med_1n, pdfs_dist_high_1n = [clip(p, 1e-10, None) for p in pdfs_1n]

        if nbackgrounds.min() < 0:
            return -1e10
        nevents = array(list(nbackgrounds) + [nrelic])
        # testing = pdfs_dist_high
        # if sum(testing <= 0.0) > 0:
        #     raise ValueError("nonpositivity of distorted pdf")

        totlike = (log(einsum("j,ijk", nevents, pdfs_dist_high)).sum(axis=0)
                   + log(einsum("j,ijk", nevents, pdfs_dist_med)).sum(axis=0)
                   + log(einsum("j,ijk", nevents, pdfs_dist_low)).sum(axis=0)
                   + log(einsum("j,ijk", nevents, pdfs_dist_high_1n)).sum(axis=0)
                   + log(einsum("j,ijk", nevents, pdfs_dist_med_1n)).sum(axis=0)
                   + log(einsum("j,ijk", nevents, pdfs_dist_low_1n)).sum(axis=0)
                   - nrelic - nbackgrounds.sum())
        totmax = totlike.max()
        gauss = exp(-arange(-4,4.5,0.5)**2/2)
        gauss /= gauss.sum()
        likenew = log((exp(totlike - totmax) * gauss).sum()) + totmax
        return likenew

    if sys == 3:
        def funclike(nback):
            return -get_like_mupi_sk4(nback)
        maxlike = fmin(funclike, nback_ini, full_output=True, disp=0)
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    if sys == 2:
        def funclike(nback):
            return -get_like_esys_sk4(nback)
        maxlike = fmin(funclike, nback_ini, full_output=True, disp=0)
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    if sys == 1:
        def funclike(nback):
            return -get_like_sk4(nback)
        maxlike = fmin(funclike, nback_ini, full_output=True, disp=0)
        # print("bg, like:", maxlike[0], -maxlike[1]) # TODO: remove
        # print("") # TODO: remove
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    if sys == 0:
        def funclike(nback):
            return -get_like_nosys_sk4(nback)
        maxlike = fmin(funclike, nback_ini, full_output=True, disp=0)
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    if sys == -1:
        def funclike(nback):
            return -get_like_init_sk4(nback, *nback_ini[1:])
        maxlike = fmin(funclike, nback_ini[0], full_output=True, disp=0)
        ntot = sum([len(p) for p in pdfs])
        ntot += sum([len(p) for p in pdfs_1n])
        nccmu = ntot - maxlike[0] - sum(nback_ini) - nrelic
        return concatenate([array([maxlike[0][0], nccmu[0]]), nback_ini[1:],
                            array([nrelic, -maxlike[1]])])


def analyse(likes, final=False):
    """ Extract limits """
    lmax = likes[:, -1].max()
    bestpos = likes[:, -1].argmax()
    rel = likes[:, -2]
    best = rel[bestpos]
    # print(lmax)
    # print((likes[:, -1] - lmax).dtype)
    norm = exp(likes[:, -1] - lmax).sum()
    errminus = best - rel[searchsorted(likes[:bestpos, -1], lmax - 0.5)]
    errplus = rel[len(likes) - 1 - searchsorted(likes[bestpos:, -1][::-1], lmax - 0.5)] - best
    l90 = rel[searchsorted(exp(likes[:, -1] - lmax).cumsum(), 0.9 * norm)]
    return lmax, best, errplus, errminus, l90


def plotfit(nnue, nnumu, nnc, nmupi, model, sknum, elow, ehigh, elow_1n,
            samples, nrelic = 0, nspall = 0, samples_n=None, signal=None, background=None, use_spall = False):
    """ Plot spectral fit """
    def plotregion(region, data, ax, elow=elow, ntag=None):
        #plt.figure()
        step = 2
        en = arange(elow, ehigh, 0.1)
        nuecc = nnue * array([pdf(ee, sknum, model, elow, pdfids["nue"], region,
                                ntag=ntag, backgrounds=background) for ee in en])
        numucc = nnumu * array([pdf(ee, sknum, model, elow, pdfids["numu"], region,
                                ntag=ntag, backgrounds=background) for ee in en])
        nc = nnc * array([pdf(ee, sknum, model, elow, pdfids["nc"], region,
                            ntag=ntag, backgrounds=background) for ee in en])
        mupi = nmupi * array([pdf(ee, sknum, model, elow, pdfids["mupi"], region,
                                ntag=ntag, backgrounds=background) for ee in en])
        spall = None
        if use_spall:
            spall = nspall * array([pdf(ee, sknum, model, elow, pdfids["spall"], region,
                                ntag=ntag, backgrounds=background) for ee in en])
        relic = nrelic * array([pdf(ee, sknum, model, elow, pdfids['rel'], region,
                                    ntag=ntag,signal=signal) for ee in en])
        h = histogram(data, bins=arange(elow, ehigh, step))
        x = h[1][1:] - 0.5 * (h[1][1:] - h[1][:-1])
        ax.errorbar(x, h[0], xerr = step/2, yerr = sqrt(h[0]),
                    elinewidth=1, fmt='.', color='black')

        bglabels = [r"CC $\nu_e$", r"Decay e$^-$", r"NCQE", r"$\mu/\pi$"]
        bgs = [nuecc, numucc, nc, mupi]
        colors=["C1", "C5", "C3",  "C6"]
        bg_plots = []
        for i, bg in enumerate(bgs):
            bgplot, = ax.plot(en, step*bg, label=bglabels[i],
                              linewidth=1.5, color=colors[i])
            bg_plots += [bgplot]

        if use_spall:
            bgplot, = ax.plot(en, step*spall, label="Spallation",
                        linewidth=1.5, color="C4")
            bg_plots += [bgplot]
            bglabels += ["Spallation"]

        model_label =  "DSNB (Horiuchi+09)" if model=="horiuchi" else "DSNB"
        s_plot, = ax.plot(en, step*relic, label=model_label,
                          linewidth=3, color="C2")
        bg_total = nuecc + numucc + nc + mupi
        if use_spall: bg_total += spall
        ax.plot(en, step*(bg_total+relic), linewidth=3, color="C2")
        all_bg_plot, = ax.plot(en, step*(bg_total), label="All backgrounds",
                               linewidth=3, color="tab:blue")

        xlims = (16, 80)
        ylims = (-2, 22)
        xticks = append(arange(20, 80, 10), 16)
        yticks = arange(0, 22, 2)
        plt.sca(ax)
        plt.xlim(*xlims)
        plt.xticks(xticks)
        plt.yticks(yticks)
        plt.ylim(*ylims)
        if ax.is_first_col():
            ax.set(ylabel=f"Events / {step:d} MeV")
        else:
            plt.setp(ax.get_yticklabels(), visible=False)

        if ax.is_first_col() and ax.is_first_row():
            yshift = 0.86 if ntag is not None else 0.9
            legend1 = plt.legend([s_plot, all_bg_plot],
                        ["DSNB", "All backgrounds"],
                         loc="upper left", prop={'weight':'bold'})
            legend2 = plt.legend(bg_plots, bglabels, loc="upper left",
                       bbox_to_anchor=(0.05,yshift))
            plt.gca().add_artist(legend1)
            plt.gca().add_artist(legend2)

        if ax.is_first_col() and ntag is not None:
            nlabels = ["N$_{ntag}\\neq1$", "N$_{ntag}=1$"]
            ax.text(-0.3,0.5, nlabels[ntag], size=30,
                            transform=ax.transAxes,
                            verticalalignment='center', rotation=90)

        if ax.is_first_row() and ax.is_last_col():
            sklabels = ["SK-I\n1497 days",
                        "SK-II\n794 days",
                        "SK-III\n562 days",
                        "SK-IV\n2970 days"]
            yshift = 0.76 if ntag is not None else 0.8
            ax.text(0.88,yshift, sklabels[sknum-1], size=30,
                        transform=ax.transAxes, weight="bold",
                        horizontalalignment='right')
        

        if ax.is_first_row():
            alabels = ["$20\degree<\\theta_C<38\degree$",
               "$38\degree<\\theta_C<50\degree$",
               "$70\degree<\\theta_C<90\degree$"]

            ax.text(0.5,1.03, alabels[region], size=30,
                        transform=ax.transAxes,
                        horizontalalignment='center')
        if ax.is_first_row() and ntag is not None:
            plt.setp(ax.get_xticklabels(), visible=False)
        if ax.is_last_row():
            ax.set(xlabel="E$_\\mathrm{rec}$ [MeV]")
        plt.grid(which="minor")
    if sknum < 4:
        plt.figure(figsize=(16.0, 8.0))
        gs = gridspec.GridSpec(1, 3)
        for areg in range(3):
            ax = plt.subplot(gs[0, areg])
            plotregion(areg, samples[areg], ax)
    else:
        plt.figure(figsize=(16.0, 12.0))
        gs = gridspec.GridSpec(2, 3)
        for areg in range(3):
            ax = plt.subplot(gs[0, areg])
            plotregion(areg, samples[areg], ax, elow=elow, ntag=False)
            ax = plt.subplot(gs[1, areg])
            plotregion(areg, samples_n[areg], ax, elow=elow_1n, ntag=True)
    plt.subplots_adjust(wspace=0)
    plt.subplots_adjust(hspace=0)


def maxlike(sknum, model, elow, ehigh=90, elow_1n=16, rmin=-5, rmax=100,
            rstep=0.1, quiet=True, outdir='.', systematics = 1,
            sk4toydir=None, gd_frac=None, use_spall=False, no_nc=False,
            trig_lo=1., trig_hi=535., ineff_scale=1.0):
    '''
    Main maximum likelihood function
    sknum = 1,2,3 (SK phase)
    model = SRN model
    elow = energy threshold (here, 16MeV)
    rmin, rmax, rstep = range and step of numbers of relic events for likelihood maximization
    '''

    def load_sample(ntag=False):
        ''' Data samples for SK I-IV '''
        if sk4toydir is None:
            if ntag:
                low = loadtxt("sk{}/ntag/samplelow.txt".format(int(sknum)))[:, 1]
                med = loadtxt("sk{}/ntag/samplemed.txt".format(int(sknum)))[:, 1]
                high = loadtxt("sk{}/ntag/samplehigh.txt".format(int(sknum)))[:, 1]
                low = low[(low > elow_1n) & (low < ehigh)]
                med = med[(med > elow_1n) & (med < ehigh)]
                high = high[(high > elow_1n) & (high < ehigh)]
            else:
                #if sknum < 4:
                low = loadtxt("sk{}/samplelow.txt".format(int(sknum)))[:, 1]
                med = loadtxt("sk{}/samplemed.txt".format(int(sknum)))[:, 1]
                high = loadtxt("sk{}/samplehigh.txt".format(int(sknum)))[:, 1]
                #if sknum == 4:
                    #low = loadtxt("sk{}/tight/samplelow.txt".format(int(sknum)))[:, 1]
                    #med = loadtxt("sk{}/tight/samplemed.txt".format(int(sknum)))[:, 1]
                    #high = loadtxt("sk{}/tight/samplehigh.txt".format(int(sknum)))[:, 1]
                low = low[(low > elow) & (low < ehigh)]
                med = med[(med > elow) & (med < ehigh)]
                high = high[(high > elow) & (high < ehigh)]
        else:
            print(f"Loading toy data samples at {sk4toydir}")
            if ntag:
                low = loadtxt(f"{sk4toydir}/low_ntag.txt")
                med = loadtxt(f"{sk4toydir}/med_ntag.txt")
                high = loadtxt(f"{sk4toydir}/high_ntag.txt")
                low = low[(low > elow_1n) & (low < ehigh)]
                med = med[(med > elow_1n) & (med < ehigh)]
                high = high[(high > elow_1n) & (high < ehigh)]
            else:
                low = loadtxt(f"{sk4toydir}/low.txt")
                med = loadtxt(f"{sk4toydir}/med.txt")
                high = loadtxt(f"{sk4toydir}/high.txt")
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
        numbkg = 5 if use_spall else 4
        for e in energies:
            pdflist = [pdf(e, sknum, model, elow, i, regionids[region],
                signal=signal, backgrounds=backgrounds, ntag=ntag)
                for i in range(numbkg)]
            pdflist += [pdf(e, sknum, model, elow, pdfids['rel'], regionids[region],
                           signal=signal, backgrounds=backgrounds, ntag=ntag)]
            p += [pdflist]
        p = array(p)
        if len(p) == 0: # No events in region
            p = p.reshape((0, numbkg + 1))
        return p  # (Nen x Npdf)

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
        spall = 0

        # Maximize likelihoods over backgrounds in signal region
        if use_spall: likemax = getmaxlike(rel, array([nback/5.,nc,mupi,spall]), low, med, high, sknum, sys=-1)
        else: likemax = getmaxlike(rel, array([nback/5.,nc,mupi]), low, med, high, sknum, sys=-1)
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
        if no_nc: nc = 0
        spall = 8 if use_spall else 0

        # Maximize likelihoods over backgrounds in signal region
        if use_spall: likemax = getmaxlike_sk4(rel, array([nback/5.,nc,mupi,spall]), [low, med, high],
                                [low_1n, med_1n, high_1n], sys=-1)
        else: likemax = getmaxlike_sk4(rel, array([nback/5.,nc,mupi]), [low, med, high],
                                [low_1n, med_1n, high_1n], sys=-1)
        return likemax

    def applysys(likes, eff, rmin, rmax, rstep, sysfact = 0, maxeff=1.0):
        ''' Apply gaussian systematic efficiency error correction'''
        print(f"Signal efficiency is {eff}")
        syseff = sqrt(sys_eff[sknum - 1]**2 + sysfact**2 * sys_eff_sk4_ntag**2)
        lower = max(eff * (1 - 6*syseff), 1e-10)
        # upper = min(eff * (1 + 6*syseff), 0.999)
        upper = min(eff * (1 + 6*syseff), maxeff) # TODO: remove?
        # upper = min(eff * (1 + 6*0.001), 1.999) # TODO: remove
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
        lconv = flikes(rates[:, newaxis] * epsrange[newaxis, :] * livetimes[sknum - 1]/365.25)
        simpsoncoeff = array([step/3.] + list((1 + (arange(1,1000)%2))*2./3 * step) + [step/3.])
        ltot = (lconv * (pgaus * epsrange * simpsoncoeff)).sum(axis = 1)
        # print("negativity of ltot", any(ltot<=0))
        likenew = log(ltot) + lmax
        like_rate = interp1d(rates, likenew)
        num_rate = likes[:, -2]/(eff * livetimes[sknum - 1]) * 365.25
        return column_stack((likes[:, :-1], num_rate, like_rate(num_rate)))

    if sknum == 2:
        elow = 17.5
    
    if gd_frac is not None:
        skgd_params(gd_frac, trig_lo, trig_hi)

    samplow, sampmed, samphigh = load_sample()
    samples_n = None
    if sknum == 4:
        samplow_n, sampmed_n, samphigh_n = load_sample(ntag=True) # SK-IV ntag samples
        samples_n = [samplow_n, sampmed_n, samphigh_n]

    # Get signal and background spectra
    signal = None
    effsignal = 1.0
    signal, flux_fac, pred_rate, pred_flux = load_signal_pdf(sknum, model,
                                            elow, ehigh, elow_1n, ineff_scale)

    bgs_sk4 = None
    bg_sk4_dir = "./pdf_bg_sk4"
    maxeffsignal = signal.max_efficiency_16_90()
    effsignal = signal.overall_efficiency_16_90() # Always use calculated eff.
    ntag_sys_fact = signal._get_ntag_sysfact() if sknum >= 4 else 0
    # Load background pdfs
    # WARNING!!! 16 MeV threshold is hardcoded there!
    cut_bins_ntag, cut_effs_ntag = get_spasolbins(ntag=True)
    cut_bins, cut_effs = get_spasolbins(ntag=False)

    # ntag_rescaling = ones((4, 2))
    if gd_frac is not None:
        ntag_rescaling = ntag_rescale(ntag_ebins, ntag_effs*ntag_eff_ps,
                                            ntag_bgs, ntag_bg_ps)
        new_bins = sorted(list(set(append(cut_bins, ntag_ebins))))
        spa_new = array(cut_effs)[digitize(new_bins[:-1], cut_bins) - 1]
        rescale_new = ntag_rescaling[:,0,digitize(new_bins[:-1],ntag_ebins)-1]
        cut_bins = new_bins
        cut_effs = spa_new[newaxis,:] * rescale_new

        new_bins_n = sorted(list(set(append(cut_bins_ntag, ntag_ebins))))
        spa_new_n = array(cut_effs_ntag)[digitize(new_bins_n[:-1], cut_bins_ntag)-1]
        rescale_new_n = ntag_rescaling[:,1,digitize(new_bins_n[:-1],ntag_ebins)-1]
        cut_bins_ntag = new_bins_n
        cut_effs_ntag = spa_new_n[newaxis,:] * rescale_new_n
    else:
        cut_effs = tile(cut_effs, (4, 1))
        cut_effs_ntag = tile(cut_effs_ntag, (4, 1))

    numbkg = 4
    bgs_sk4 = [bg_sk4(i, cut_bins, cut_effs[i],
                cut_bins_ntag, cut_effs_ntag[i], bg_sk4_dir, elow,
                ehigh=ehigh, elow_n=elow_1n ) for i in range(numbkg)]
    if use_spall:
        solbins = array(list(soleff[sknum - 1][:, 0]) + [90.])
        soleffs = soleff[sknum - 1][:, 1]
        bgs_sk4 += [spall_sk(solbins, soleffs, effs_3rdred[sknum - 1], sknum, elow, ehigh=ehigh)]
    print(f"Efficiency is {effsignal}")
    if sknum == 4: print(f"Factor for ntag systematics is {ntag_sys_fact}")

    # Get pdfs. PDF matrices: (Nen x Npdf)
    pdfs_high = get_pdfmatrix(samphigh, "high", ntag=False,
                              signal=signal, backgrounds=bgs_sk4)
    pdfs_med = get_pdfmatrix(sampmed, "medium", ntag=False,
                              signal=signal, backgrounds=bgs_sk4)
    pdfs_low = get_pdfmatrix(samplow, "low", ntag=False,
                              signal=signal, backgrounds=bgs_sk4)

    if sknum < 4:
        # Set backgrounds (preliminary estimates)
        init = initialize(pdfs_low, pdfs_med, pdfs_high, rmin)
        bkglist = init[:-2]
        #nue, numu, nc, mupi, spall, _, _ = init
        #if not isscalar(numu): numu=numu[0]

        # Get systematic error matrices
        if systematics:
            sysmatrices = None
            if systematics == 1:
                sysmatrices = systematics_atm([sampmed, samphigh], sknum, model,
                                           elow, ehigh, backgrounds=bgs_sk4, use_spall = use_spall)
                # Distort pdfs: (Nen x Npdfs x Nsigma x Nsigma2 x Nsigma3)
                sysmatrix_med, sysmatrix_high = sysmatrices[:3]
                if use_spall:
                    pdfs_high = pdfs_high[...,newaxis,newaxis,newaxis] * sysmatrix_high
                    pdfs_med = pdfs_med[...,newaxis,newaxis,newaxis] * sysmatrix_med
                else:
                    pdfs_high = pdfs_high[...,newaxis,newaxis] * sysmatrix_high
                    pdfs_med = pdfs_med[...,newaxis,newaxis] * sysmatrix_med
            if systematics == 2:
                sysmatrices = systematics_escale_res([samplow, sampmed, samphigh], sknum, model,
                                           elow, ehigh, backgrounds=bgs_sk4, signal = signal, use_spall = use_spall)
                sysmatrix_low, sysmatrix_med, sysmatrix_high = sysmatrices[:3]
                # Distort pdfs: (Nen x Npdfs x Nsigma)
                pdfs_low = pdfs_low[...,newaxis] * sysmatrix_low
                pdfs_high = pdfs_high[...,newaxis] * sysmatrix_high
                pdfs_med = pdfs_med[...,newaxis] * sysmatrix_med

        # Main maximization loop
        likedata = []
        rmin = 0
        #bkglist = [nue, numu, nc, mupi, spall] if use_spall else [nue, numu, nc, mupi]
        for i,rel in enumerate(arange(rmin, rmax, rstep)):
            likeres = getmaxlike(rel, array(bkglist),
                                 pdfs_low, pdfs_med, pdfs_high,
                                 sknum, sys=systematics, use_spall = use_spall)
            likedata.append(likeres)
            # Update initial values
            bkglist = list(likedata[-1][:-2])
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

        # Set backgrounds (preliminary estimates)
        init = initialize_sk4(pdfs_low, pdfs_med, pdfs_high,
                              pdfs_low_n, pdfs_med_n, pdfs_high_n, rmin)
        bkglist = init[:-2]
        #nue, numu, nc, mupi, spall, _, _ = init

        #if sum(pdfs_high <= 0.0) > 0:
            #raise ValueError("zeros in pdf")

        # Get systematic error matrices
        if systematics:
            sysmatrices = None
            if systematics == 1:
                sysmatrices = systematics_atm([samplow, sampmed, samphigh], sknum,
                                        model, elow, ehigh, elow_1n=elow_1n,
                                        backgrounds=bgs_sk4,
                                        energies_n=[samplow_n, sampmed_n, samphigh_n],
                                        use_spall=use_spall, no_nc=no_nc)
                sysmatrix_low, sysmatrix_med, sysmatrix_high = sysmatrices[:3]
                sysmatrix_low_1n, sysmatrix_med_1n, sysmatrix_high_1n = sysmatrices[3:]
                if use_spall:
                    # Distort pdfs
                    pdfs_high = pdfs_high[...,newaxis,newaxis,newaxis,newaxis] * sysmatrix_high
                    pdfs_med = pdfs_med[...,newaxis,newaxis,newaxis,newaxis] * sysmatrix_med
                    pdfs_low = pdfs_low[...,newaxis,newaxis,newaxis,newaxis] * sysmatrix_low
                    pdfs_high_n = pdfs_high_n[...,newaxis,newaxis,newaxis,newaxis] * sysmatrix_high_1n
                    pdfs_med_n = pdfs_med_n[...,newaxis,newaxis,newaxis,newaxis] * sysmatrix_med_1n
                    pdfs_low_n = pdfs_low_n[...,newaxis,newaxis,newaxis,newaxis] * sysmatrix_low_1n
                elif no_nc:
                    # Distort pdfs
                    pdfs_high = pdfs_high[...,newaxis,newaxis] * sysmatrix_high
                    pdfs_med = pdfs_med[...,newaxis,newaxis] * sysmatrix_med
                    pdfs_low = pdfs_low[...,newaxis,newaxis] * sysmatrix_low
                    pdfs_high_n = pdfs_high_n[...,newaxis,newaxis] * sysmatrix_high_1n
                    pdfs_med_n = pdfs_med_n[...,newaxis,newaxis] * sysmatrix_med_1n
                    pdfs_low_n = pdfs_low_n[...,newaxis,newaxis] * sysmatrix_low_1n
                else:
                    # Distort pdfs
                    pdfs_high = pdfs_high[...,newaxis,newaxis,newaxis] * sysmatrix_high
                    pdfs_med = pdfs_med[...,newaxis,newaxis,newaxis] * sysmatrix_med
                    pdfs_low = pdfs_low[...,newaxis,newaxis,newaxis] * sysmatrix_low
                    pdfs_high_n = pdfs_high_n[...,newaxis,newaxis,newaxis] * sysmatrix_high_1n
                    pdfs_med_n = pdfs_med_n[...,newaxis,newaxis,newaxis] * sysmatrix_med_1n
                    pdfs_low_n = pdfs_low_n[...,newaxis,newaxis,newaxis] * sysmatrix_low_1n
            if systematics == 2:
                # Get distorsion factors
                print("Getting distorsion factors for all spectra...")
                global esys_scale_dict, esys_res_dict
                esys_scale_dict, esys_res_dict = get_distorsion_functions(elow, ehigh, elow_1n, bgs_sk4, signal, sknum)
                print("Done")
                # Get distorsion matrices
                sysmatrices = systematics_escale_res([samplow, sampmed, samphigh], sknum,
                                        model, elow, ehigh, elow_1n=elow_1n,
                                        backgrounds=bgs_sk4, signal = signal,
                                        energies_n=[samplow_n, sampmed_n, samphigh_n], use_spall = use_spall)
                sysmatrix_low, sysmatrix_med, sysmatrix_high = sysmatrices[:3]
                sysmatrix_low_1n, sysmatrix_med_1n, sysmatrix_high_1n = sysmatrices[3:]
                # Distort pdfs
                pdfs_high = pdfs_high[...,newaxis] * sysmatrix_high
                pdfs_med = pdfs_med[...,newaxis] * sysmatrix_med
                pdfs_low = pdfs_low[...,newaxis] * sysmatrix_low
                pdfs_high_n = pdfs_high_n[...,newaxis] * sysmatrix_high_1n
                pdfs_med_n = pdfs_med_n[...,newaxis] * sysmatrix_med_1n
                pdfs_low_n = pdfs_low_n[...,newaxis] * sysmatrix_low_1n
            if systematics == 3:
                sysmatrices = systematics_mupi([samplow, sampmed, samphigh], sknum,
                                        model, elow, ehigh, elow_1n=elow_1n,
                                        backgrounds=bgs_sk4,
                                        energies_n=[samplow_n, sampmed_n, samphigh_n], use_spall = use_spall)
                sysmatrix_low, sysmatrix_low_1n = sysmatrices
                # Distort pdfs
                pdfs_low = pdfs_low[...,newaxis] * sysmatrix_low
                pdfs_low_n = pdfs_low_n[...,newaxis] * sysmatrix_low_1n

        # Main maximization loop
        likedata = []
        rmin = 0
        #bkglist = [nue, numu, nc, mupi, spall] if use_spall else [nue, numu, nc, mupi]
        print(f"Livetime is {livetimes[3]}")
        for i,rel in enumerate(arange(rmin, rmax, rstep)):
            likeres = getmaxlike_sk4(rel, array(bkglist),
                            [pdfs_low, pdfs_med, pdfs_high],
                            [pdfs_low_n, pdfs_med_n, pdfs_high_n],
                            sys=systematics, use_spall=use_spall, no_nc=no_nc)
            likedata.append(likeres)
            # Update initial values
            bkglist = list(likedata[-1][:-2])
            if i % 100 == 0:
                print("Step {}/1000, like = {}".format(i, likedata[-1][-1]), flush=True)

    results = column_stack((arange(rmin, rmax, rstep), likedata))
    results = results[results[:, 0] >= 0]   # results[i] = rel, nback[0:4], rel, like

    # Systematic efficiency error correction + limits
    _, best2, errplus2, errminus2, limit2 = analyse(results)
    results_sys = applysys(results, effsignal, rmin, rmax, rstep, ntag_sys_fact, maxeffsignal)
    _, best, errplus, errminus, limit = analyse(results_sys, final=True)

    flux_best = best * flux_fac
    flux_90cl = limit * flux_fac
    lpos = results_sys[:, -1].argmax()

    # Save and display results
    savetxt(f"{outdir}/fit_sk{sknum}.txt", column_stack((results, results_sys[:, -2:])))
    print("SK-{}".format(sknum), "Best fit:")
    print("{} +{} -{} relic evts/yr".format(best, errplus, errminus))
    print("{} +{} -{} /cm^2/s > 17.3 MeV".format(flux_best, errplus*flux_fac, errminus*flux_fac))
    print("{} +{} -{} relic evts".format(best2, errplus2, errminus2))
    print("90% c.l. relic event rate: {} ev/yr {}".format(limit, limit2))
    print("90% c.l. {} /cm^2/s > 17.3 MeV".format(flux_90cl))
    print("nu-e events {}".format(results_sys[lpos, 1]))
    print("nu-mu events {}".format(results_sys[lpos, 2]))
    print("NC elastic events {}".format(results_sys[lpos, 3]))
    print("mu/pi events {}".format(results_sys[lpos, 4]))
    if use_spall: print("Spallation events {}".format(results_sys[lpos, 5]))
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
        nrelic = results_sys[lpos,6] if use_spall else results_sys[lpos,5] 
        nspall = results_sys[lpos,5] if use_spall else 0 
        plotfit(results_sys[lpos,1], results_sys[lpos,2], results_sys[lpos,3],
                results_sys[lpos,4], model, sknum,
                elow, ehigh, elow_1n, samples=[samplow, sampmed, samphigh],
                samples_n=samples_n, signal=signal, background=bgs_sk4, use_spall = use_spall,
               nrelic = nrelic, nspall = nspall)
        plt.savefig(f"{outdir}/fit_sk{sknum}.pdf")
        plt.clf()
    return limit * flux_fac, flux_fac, results_sys, pred_rate, pred_flux


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
        rels, likes = r[:, -2], r[:, -1]
        fluxes = rels * fluxfacts[i + 1]
        flike = interp1d(fluxes, likes, bounds_error=False, fill_value=1e-10)
        newlikes = flike(flux_sampling)
        liketot += newlikes - newlikes.max()
    return analyse(column_stack((flux_sampling, liketot)), final=True)


def fullike(model, elow, ehigh, elow_sk2 = 17.5, elow_sk4=None, ehigh_sk4=None, elow_sk4_1n=None,
            rmin=-5, rmax=100, rstep=0.1, quiet=False, outdir='.', systematics = 1, use_spall = False):
    """ Fit SK I-IV data """
    if elow_sk4 is None:
        elow_sk4 = elow
    if elow_sk4_1n is None:
        elow_sk4_1n = elow
    if ehigh_sk4 is None:
        ehigh_sk4 = ehigh
    like1 = maxlike(1, model, elow, ehigh, elow_sk4_1n, rmin, rmax,
                    rstep, quiet=quiet, outdir=outdir, systematics = systematics, use_spall = use_spall)
    like2 = maxlike(2, model, elow_sk2, ehigh, elow_sk4_1n, rmin, rmax,
                    rstep, quiet=quiet, outdir=outdir, systematics = systematics, use_spall = use_spall)
    like3 = maxlike(3, model, elow, ehigh, elow_sk4_1n, rmin, rmax,
                    rstep, quiet=quiet, outdir=outdir, systematics = systematics, use_spall = use_spall)
    like4 = maxlike(4, model, elow_sk4, ehigh_sk4, elow_sk4_1n, rmin, rmax,
                    rstep, quiet=quiet, outdir=outdir, systematics = systematics, use_spall = use_spall)
    fluxlims = [like1[0], like2[0], like3[0], like4[0]]
    fluxfacs = [like1[1], like2[1], like3[1], like4[1]]
    results = [like1[2], like2[2], like3[2], like4[2]]
    rate = arange(0, rmax, rstep)
    for i,r in enumerate(results):
        flike = interp1d(r[:, -2], r[:, -1], bounds_error = False, fill_value = r[:, -1].min())
        results[i] = column_stack((rate, flike(rate)))
    pred_rate = like1[3]
    pred_flux = like1[4]
    ratelims = array(fluxlims) / array(fluxfacs)
    res = combine(results)
    res_fluxes = combine_fluxes(results, fluxfacs)
    _, ratebest_comb, ratepl_comb, ratemin_comb, ratelim_comb = res
    _, fluxbest_comb, fluxpl_comb, fluxmin_comb, fluxlim_comb = res_fluxes
    print("Results: ", results[-1][results[-1][:, -1].argmax(), -2])

    if not quiet:
        plt.figure(figsize=(12.0, 8.0))
        plt.xlabel("DSNB events/year")
        plt.ylabel("Likelihood")
        x = results[0][:, -2]
        plt.plot(x, results[0][:, -1] - results[0][:,-1].max(),
                 label="SK-I", alpha=0.5)
        plt.plot(x, results[1][:, -1] - results[1][:,-1].max(),
                 '--', label="SK-II", alpha=0.5)
        plt.plot(x, results[2][:, -1] - results[2][:,-1].max(),
                 '-.', label="SK-III", alpha=0.5)
        plt.plot(x, results[3][:, -1] - results[3][:,-1].max(),
                 '--', label="SK-IV",linewidth=2)
        likesum = results[0][:, -1] + results[1][:, -1] + results[2][:, -1] + results[3][:, -1]
        likesum -= likesum.max()
        plt.plot(x, likesum, label = "Combined", color = 'black')
        plt.plot([0,20], -0.5 * ones(2), 'r')
        plt.xlim(0,20)
        plt.ylim(-2,0.2)
        # plt.grid()
        plt.legend()
        plt.savefig(outdir + "/full_like.pdf")
        plt.clf()

    print("")
    print(f"Best fit rate: {ratebest_comb} + {ratepl_comb} - {ratemin_comb}")
    print("90%% c.l. rate limits are: %f %f %f %f > 16 MeV" % tuple(ratelims))
    print(f"90%% c.l. combined rate limit: {ratelim_comb} > 16 MeV")
    print(f"Predicted rate: {pred_rate}")

    print("")
    print(f"Best fit flux: {fluxbest_comb} + {fluxpl_comb} - {fluxmin_comb}")
    print("90%% c.l. flux limits are: %f %f %f %f /cm^2/s > 17.3 MeV" % tuple(fluxlims))
    print(f"90%% c.l. combined flux limit: {fluxlim_comb} /cm^2/s > 17.3 MeV")
    print(f"Predicted flux: {pred_flux}")

def plot_results(model, elow, ehigh, elow_sk2 = 17.5, elow_sk4=None, ehigh_sk4=None, elow_sk4_1n=None,
            outdir='.', use_spall = False):
    '''
    Plot results without fitting for all SK phases
    model = SRN model
    elow = energy threshold (here, 16MeV)
    '''
    os.makedirs(f"{outdir}/replot", exist_ok=True)

    def load_sample(sknum, ntag=False):
        ''' Data samples for SK I-IV '''
        if ntag:
            low = loadtxt("sk{}/ntag/samplelow.txt".format(int(sknum)))[:, 1]
            med = loadtxt("sk{}/ntag/samplemed.txt".format(int(sknum)))[:, 1]
            high = loadtxt("sk{}/ntag/samplehigh.txt".format(int(sknum)))[:, 1]
            low = low[(low > elow_sk4_1n) & (low < ehigh_sk4)]
            med = med[(med > elow_sk4_1n) & (med < ehigh_sk4)]
            high = high[(high > elow_sk4_1n) & (high < ehigh_sk4)]
        else:
            #if sknum < 4:
            low = loadtxt("sk{}/samplelow.txt".format(int(sknum)))[:, 1]
            med = loadtxt("sk{}/samplemed.txt".format(int(sknum)))[:, 1]
            high = loadtxt("sk{}/samplehigh.txt".format(int(sknum)))[:, 1]
            print(f"sample: {(med < 18).sum()} {elow} {sknum}")
            #if sknum == 4:
                #low = loadtxt("sk{}/tight/samplelow.txt".format(int(sknum)))[:, 1]
                #med = loadtxt("sk{}/tight/samplemed.txt".format(int(sknum)))[:, 1]
                #high = loadtxt("sk{}/tight/samplehigh.txt".format(int(sknum)))[:, 1]
            low = low[(low > elow) & (low < ehigh)]
            med = med[(med > elow) & (med < ehigh)]
            high = high[(high > elow) & (high < ehigh)]
        return low,med,high

    def get_spasolbins(sknum, ntag=False):
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

    # Plot spectra for all phases of SK
    """ Fit SK I-IV data """
    if elow_sk4 is None:
        elow_sk4 = elow
    if elow_sk4_1n is None:
        elow_sk4_1n = elow
    if ehigh_sk4 is None:
        ehigh_sk4 = ehigh
    with open(f"{outdir}/fit.log") as fitfile:
        lines = fitfile.readlines()
        nnue = array([float(line.split()[-1]) for line in lines if line.startswith("nu-e")])
        nnumu = array([float(line.split()[-1]) for line in lines if line.startswith("nu-mu")])
        nnc = array([float(line.split()[-1]) for line in lines if line.startswith("NC elastic")])
        nmupi = array([float(line.split()[-1]) for line in lines if line.startswith("mu/pi")])
        nspall = array([float(line.split()[-1]) for line in lines if line.startswith("Spallation")]) if use_spall else zeros((4,))
        rates_relic = array([float(line.split()[0]) for line in lines if (line.strip()).endswith("relic evts/yr")])
        print(nnue, nnumu, nnc, nmupi, nspall, rates_relic)
        elow_sk13 = elow
        for sknum in range(1, 5):
            if sknum == 4:
                elow = elow_sk4
                ehigh = ehigh_sk4
                samplow, sampmed, samphigh = load_sample(sknum)
                samples_n = None
                samplow_n, sampmed_n, samphigh_n = load_sample(sknum, ntag=True) # SK-IV ntag samples
                samples_n = [samplow_n, sampmed_n, samphigh_n]
            elif sknum == 2:
                elow = elow_sk2
                samplow, sampmed, samphigh = load_sample(sknum)
                samples_n = None
            else:
                elow = elow_sk13
                samplow, sampmed, samphigh = load_sample(sknum)
                samples_n = None
            bg_sk4_dir = "./pdf_bg_sk4"
            cut_bins_ntag, cut_effs_ntag = get_spasolbins(sknum, ntag=True)
            cut_bins, cut_effs = get_spasolbins(sknum, ntag=False)
            cut_effs = tile(cut_effs, (4, 1))
            cut_effs_ntag = tile(cut_effs_ntag, (4, 1))

            bgs_sk4 = [bg_sk4(i, cut_bins, cut_effs[i],
                        cut_bins_ntag, cut_effs_ntag[i], bg_sk4_dir, elow,
                        ehigh=ehigh, elow_n=elow_sk4_1n ) for i in range(4)]
            if use_spall:
                solbins = array(list(soleff[sknum - 1][:, 0]) + [90.])
                soleffs = soleff[sknum - 1][:, 1]
                bgs_sk4 += [spall_sk(solbins, soleffs, effs_3rdred[sknum - 1], sknum, elow, ehigh=ehigh)]
            signal, flux_fac, pred_rate, pred_flux = load_signal_pdf(sknum, model, elow, ehigh, elow_sk4_1n)
            effsignal = signal.overall_efficiency_16_90() # Always use calculated eff.
            nrel = rates_relic[sknum - 1] * effsignal * livetimes[sknum - 1]/365.25
            print(f"sample: {(sampmed < 18).sum()} {elow} {sknum}")
            plotfit(nnue[sknum - 1], nnumu[sknum - 1], nnc[sknum - 1], 
                    nmupi[sknum - 1], model[sknum - 1], sknum,
                    elow, ehigh, elow_sk4_1n, samples=[samplow, sampmed, samphigh],
                    samples_n=samples_n, signal=signal, background=bgs_sk4, use_spall = use_spall,
                   nrelic = nrel, nspall = nspall[sknum - 1])
            plt.savefig(f"{outdir}/replot/fit_sk{sknum}.pdf")
            plt.clf()
    # Plot likelihoods
    results = [loadtxt(f"{outdir}/fit_sk{sknum}.txt") for sknum in range(1,5)]
    rate = results[0][:, 0]
    for i,r in enumerate(results):
        flike = interp1d(r[:, -2], r[:, -1], bounds_error = False, fill_value = r[:, -1].min())
        results[i] = column_stack((rate, flike(rate)))
    plt.figure(figsize=(12.0, 8.0))
    plt.xlabel("DSNB rate [events/year]")
    plt.ylabel("Likelihood")
    x = results[0][:, -2]
    plt.plot(x, results[0][:, -1] - results[0][:,-1].max(),
             label="SK-I", alpha=0.5)
    plt.plot(x, results[1][:, -1] - results[1][:,-1].max(),
             '--', label="SK-II", alpha=0.5)
    plt.plot(x, results[2][:, -1] - results[2][:,-1].max(),
             '-.', label="SK-III", alpha=0.5)
    plt.plot(x, results[3][:, -1] - results[3][:,-1].max(),
             '--', label="SK-IV",linewidth=2)
    likesum = results[0][:, -1] + results[1][:, -1] + results[2][:, -1] + results[3][:, -1]
    likesum -= likesum.max()
    plt.plot(x, likesum, label = "Combined", color = 'black')
    plt.plot([0,20], -0.5 * ones(2), 'r')
    plt.xlim(0,20)
    plt.ylim(-2,0.2)
    # plt.grid()
    plt.legend()
    plt.savefig(outdir + "/replot/full_like.pdf")
    plt.clf()

def sk4like(model, elow_sk4, ehigh_sk4, elow_sk4_1n=None, toydir=None,
            rmin=-5, rmax=100, rstep=0.1, quiet=False, outdir='.',
            systematics=1, gd_frac=None, use_spall=False, no_nc=False,
            trig_lo=1., trig_hi=535., ineff_scale=1.0, ncsys_scale=1.0):
    """ Fit SK I-IV data """
    if elow_sk4_1n is None:
        elow_sk4_1n = elow_sk4

    global nc_mult
    nc_mult *= ncsys_scale

    like4 = maxlike(4, model, elow_sk4, ehigh_sk4, elow_sk4_1n, rmin, rmax,
                    rstep, quiet=quiet, outdir=outdir, sk4toydir=toydir,
                    systematics=systematics, gd_frac=gd_frac, use_spall=use_spall,
                    no_nc=no_nc, trig_lo=trig_lo, trig_hi=trig_hi, ineff_scale=ineff_scale)
    fluxlim, fluxfac, result, = like4[0], like4[1], like4[2]
    pred_rate, pred_flux = like4[3], like4[4]
    ratelim = fluxlim / fluxfac

    print("")
    print("90%% c.l. rate limit is: %f > 16 MeV" % ratelim)
    print(f"Predicted rate: {pred_rate}")

    print("")
    print("90%% c.l. flux limits are: %f /cm^2/s > 17.3 MeV" % fluxlim)
    print(f"Predicted flux: {pred_flux}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('modelname', help="DNSB model name")
    parser.add_argument('directory', help='Fit output directory')
    parser.add_argument('--sys', help='systematics mode [-1, 0, 1, or 2]', type=int, default = 1)
    parser.add_argument('--thr', help='SK4 Energy threshold (non-IBD region)', type=float, default=16)
    parser.add_argument('--thr1n', help='SK4 Energy threshold (IBD region)', type=float, default=16)
    parser.add_argument('--toy', help='Toy dataset location (replaces data)')
    parser.add_argument('--gd', help=('Specify Gd neutron capture fraction,'
                                      ' otherwise pure water is assumed'), type=float, default=-1)
    parser.add_argument('--lt', help='Toy dataset livetime in days (default: 10yrs)', type=float, default=3652.5)
    parser.add_argument('--spall', help='Add spallation backgrounds (0/1, default: 0)', type=int, default=0)
    parser.add_argument('--nonc', help='Assume no NC background (0/1, default: 0)', type=int, default=0)
    parser.add_argument('--rmin', help="Minimum rate sampled", type=float, default=-5.)
    parser.add_argument('--rmax', help="Maximum rate sampled", type=float, default=100.)
    parser.add_argument('--rstep', help="Rate sampling step size", type=float, default=0.1)
    parser.add_argument('--triglo', help="SK-Gd trigger window start time", type=float, default=1.)
    parser.add_argument('--trighi', help="SK-Gd trigger window end time", type=float, default=535.)
    parser.add_argument('--ineff_scale', help="Signal inefficiency scaling", type=float, default=1.)
    parser.add_argument('--ncsys_scale', help="NC systematics scaling", type=float, default=1.)
    parser.add_argument('--quiet', help='Do not save plot of fit', type=int, default=0)
    parser.add_argument('--drawonly', help='Redraw plots from fit result file (0/1, default: 0)', type=int, default=0)
    args = parser.parse_args()

    modelname = args.modelname
    directory = args.directory
    systematics = args.sys
    e_thr = args.thr
    e_thr_1n = args.thr1n

    start_time = time()
    if args.toy:
        toy_data_dir = args.toy
        # quiet = False
        livetimes[3] = args.lt
        gd_fraction = args.gd
        if gd_fraction == -1:
            gd_fraction = None
        sk4like(modelname, elow_sk4=e_thr, ehigh_sk4=80, elow_sk4_1n=e_thr_1n,
            outdir=directory, toydir=toy_data_dir, systematics=systematics,
            gd_frac=gd_fraction, use_spall=args.spall, no_nc=args.nonc,
            quiet=args.quiet, rmin=args.rmin, rmax=args.rmax, rstep=args.rstep,
            trig_lo=args.triglo, trig_hi=args.trighi,
            ineff_scale=args.ineff_scale, ncsys_scale=args.ncsys_scale)
    elif args.drawonly:
        plot_results(modelname, elow=16, ehigh=90, elow_sk2=16,
            elow_sk4=e_thr, ehigh_sk4=80, elow_sk4_1n=e_thr_1n,
            outdir=directory, use_spall=args.spall)
    else:
        fullike(modelname, elow=16, ehigh=90, elow_sk2=17.5,
            elow_sk4=e_thr, ehigh_sk4=80, elow_sk4_1n=e_thr_1n,
            outdir=directory, systematics=systematics, use_spall=args.spall,
            quiet=args.quiet, rmin=args.rmin, rmax=args.rmax, rstep=args.rstep)
    
    end_time = time()
    elapsed_time = end_time - start_time
    elapsed_hr = int(elapsed_time // 3600)
    elapsed_min = int(elapsed_time % 3600 // 60)
    elapsed_sec = elapsed_time % 3600 % 60
    print("")
    print(f"Fit completed in {elapsed_hr} hr {elapsed_min} min {elapsed_sec:.2f} sec")
