from numpy import *
from fit_parameters import *

from scipy.interpolate import interp1d


def skgd_params(gd_cap_frac, srn_lo=1.0, srn_hi=535.0):
    """ Update fit parameters for fitting against SK-Gd projection.
    Efficiency correction depends on Gd capture fraction
    and neutron search time window"""
    # Variables to update
    global ntag_eff_ps, ntag_effs, ntag_bg_ps, ntag_bgs

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

    ntag_eff_ps_gdmix = (1-gd_cap_frac) * h_eff_ps + gd_cap_frac * gd_eff_ps
    ntag_effs_gdmix = gd_cap_frac * ntag_effs_gd  + (1-gd_cap_frac) * ntag_effs_h
    ntag_effs_gdmix /= ntag_eff_ps_gdmix
    print(f"Using Gd neutron capture fraction of {gd_cap_frac}%")

    # Update global variables
    ntag_eff_ps = ntag_eff_ps_gdmix
    ntag_effs = ntag_effs_gdmix
    ntag_bgs = ntag_bgs * ntag_bg_ps / ntag_bg_ps_sk6
    ntag_bg_ps = ntag_bg_ps_sk6

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