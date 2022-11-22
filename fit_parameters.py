from numpy import *
from scipy.interpolate import interp1d



######### SK6 parameters ############
##### merge with the rest later #####
#####################################
## Change MC files in get_flux_from_mc
## Update spall PDF parameters and coeffs
## Update energy resolution function (esys_pdf_distorsions.py)
#####################################

sk6_livetime = 522.2
sk6_livetime_test = 170 # SK-VI data subset before full opening of dataset
aft_eff_sk6 = 0.9
sys_fv_sk6 = 0.015
sys_eff_sk6 = 0.025
sys_eff_sk6_ntag = 0.05

# Directory with atmopsheric background PDF parametrizations
bg_sk6_dir = "./pdf_bg_sk6"

# Fractions of signal in Cherenkov angle regions
s_ch_frac_sk6 = [9.433e-04, 9.925e-01, 6.525e-03]

# Third reduction efficiencies
effs_sk6 = loadtxt("efficiencies/efficiencies_sk4.txt")
effsk6 = interp1d(effs_sk6[:,0], effs_sk6[:,1], bounds_error=False,
                  fill_value = (effs_sk6[0,1], effs_sk6[-1,1]))

# Signal efficiencies of spallation cut
spaeff_sk6_nontag = array([[12, 0.8], [16, 0.65], [18, 0.63], [20, 0.918], [24, 0.98]])
spaeff_sk6 = array([[16, 0.826], [18, 0.887], [20, 0.918], [24, 0.98]])

# Solar cut efficiencies
soleff_sk6 = array([[16, 0.731], [17, 0.822], [18, 0.883],
                    [19, 0.966], [20, 1]])

# Ntag efficiencies
ntag_ebins_sk6 = [12, 14, 16, 18, 20, 22, 24, 26, 28, 100]
bdt_cuts_sk6 = [0.997, 0.997, 0.749, 0.749, 0.712, 0.712, 0.781, 0.855, 0.937]
emin_sk6, emax_sk6 = ntag_ebins_sk6[0], ntag_ebins_sk6[-1]
bdt_roc_sk6 = genfromtxt('/disk02/usr6/giampaol/ntag-mva/models/bdt22_skg4_0.013_10M/roc_test.roc')
cuts_roc_sk6, roc_effs_sk6, roc_bg_sk6 = bdt_roc_sk6[:,0], bdt_roc_sk6[:,1], bdt_roc_sk6[:,2]
ntag_eff_sk6 = interp1d(cuts_roc_sk6, roc_effs_sk6)
ntag_bg_sk6 = interp1d(cuts_roc_sk6, roc_bg_sk6)
ntag_effs_sk6 = ntag_eff_sk6(bdt_cuts_sk6)
ntag_bgs_sk6 = ntag_bg_sk6(bdt_cuts_sk6)
ntag_eff_ps_sk6 = 0.753 # N10 > 5, 1-535 microsecs window
ntag_bg_ps_sk6 = 54 # N10 > 5, 1-535 microsecs window

#####################################
#####################################

# livetimes
livetimes = array([1497.44, 793.71, 562.04, 2970.1, sk6_livetime_test])
aft_eff = 0.94
livetimes[3] *= aft_eff
livetimes[4] *= aft_eff_sk6

# energy-independent efficiency sys
efftot = {"lma": [0.7975,0.56703,0.77969],
          "malaney": [0.7903, 0.53273, 0.766262],
          "faild": [0.7997,0.5845,0.7831] ,
          "woosley": [0.7876,0.5524,0.7764],
          "ksw" : [0.7908, 0.54527, 0.77371]}
fluxfac = {"lma": 0.535115, "malaney": 0.5845, "ksw": 0.488413,
           "faild": 0.431, "woosley": 0.47447}
# Systematics for signal: reduction efficiency, livetime, Strumia-Vissani cross-section, FV
sys_eff = array([0.0254, 0.0404, 0.0253, 0.0220, sys_eff_sk6])
sys_livetime = 0.001
sys_xsec = 0.01
sys_fv = array([0.013, 0.011, 0.010, 0.015, sys_fv_sk6])
sys_eff = sqrt(sys_eff**2 + sys_livetime**2 + sys_xsec**2 + sys_fv**2)
# Ntag cut systematics from AmBe study for SK-IV
sys_eff_sk4_ntag = 0.125

# Definitions
regionids = {"low": 0, "medium": 1, "high": 2}
pdfids = {"nue": 0, "numu": 1, "nc": 2, "mupi": 3, "spall": 4, "rel": 5}
modelids = {"lma": 0, "faild": -3, "malaney": -1, "ksw": -2, "woosley": -4}

# Directory with atmospheric background PDF parametrizations
bg_sk4_dir = "./pdf_bg_sk4"
bg_sk_dir = [bg_sk4_dir, bg_sk4_dir, bg_sk4_dir, bg_sk4_dir, bg_sk6_dir]

# Signal Cherenkov angle fractions
s_ch_frac = [9.433e-04, 9.925e-01, 6.525e-03]
s_ch_fracs = [s_ch_frac, s_ch_frac, s_ch_frac, s_ch_frac, s_ch_frac_sk6]

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
effs_3rdred = [effsk1, effsk2, effsk3, effsk4, effsk6]

# Spallation cut efficiencies
spaeff_sk1 = array([[16, 0.818], [18, 0.908], [24, 1.0]])
spaeff_sk2 = array([[17.5, 0.762], [20, 0.882], [26, 1.0]])
spaeff_sk3 = array([[16, 0.818], [18, 0.908], [24, 1.0]])
spaeff_sk4_nontag = array([[12, 0.8], [16, 0.65], [18, 0.63], [20, 0.918], [24, 0.98]])
spaeff_sk4 = array([[16, 0.826], [18, 0.887], [20, 0.918], [24, 0.98]])
spaeff = [spaeff_sk1, spaeff_sk2, spaeff_sk3, spaeff_sk4_nontag, spaeff_sk6_nontag]
spaeff_ntag = [None, None, None, spaeff_sk4, spaeff_sk6]

soleff_sk1 = array([[16, 0.738], [17, 0.821], [18, 0.878],
                    [19, 0.965], [20, 1]])
soleff_sk2 = array([[17.5, 0.738], [18.02, 0.821], [19.08, 0.878],
                    [20.14, 0.965], [21.2, 1]])
soleff_sk3 = array([[16, 0.738], [17, 0.821], [18, 0.878],
                    [19, 0.965], [20, 1]])
soleff_sk4 = array([[16, 0.731], [17, 0.822], [18, 0.883],
                    [19, 0.966], [20, 1]])
soleff = [soleff_sk1, soleff_sk2, soleff_sk3, soleff_sk4, soleff_sk6]

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

# ntag efficiencies of all SK periods
ntag_ebins_sk = [None, None, None, ntag_ebins, ntag_ebins_sk6]
bdt_cuts_sk = [None, None, None, bdt_cuts, bdt_cuts_sk6]
emin_sk = [None, None, None, emin, emin_sk6]
emax_sk = [None, None, None, emax, emax_sk6]
ntag_effs_sk = [None, None, None, ntag_effs, ntag_effs_sk6]
ntag_bgs_sk = [None, None, None, ntag_bgs, ntag_bgs_sk6]
ntag_eff_ps_sk =[None, None, None, ntag_eff_ps, ntag_eff_ps_sk6]
ntag_bg_ps_sk = [None, None, None, ntag_bg_ps, ntag_bg_ps_sk6]

# Distorsion coefficients for spallation systematics
spacoeffs_sk = array([[0.0105935649, -0.495682897, 7.93209842, -43.5523139],
                     [0.0138280665, -0.749631175, 13.6659053, -83.7150281],
                     [0.0438680847, -2.13974596, 35.1046340, -193.623584],
                     [0.0103067195, -0.475915126, 7.52829430, -40.9820538],
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
mupi_rescale_low = [1.367, 1.75, 1.34, 1.34, 1.34] # mupi from low to medium
mupi_rescale_high = [0.12777, 0.1, 0.13, 0.13, 0.13] # mupi from low to high
nc_rescale = [1.16313, 1.42, 1.14, 1.14, 1.14] # NC from high to medium

# systematics multipliers
nc_mult = 1
cc_mult = 1
n_mult = 1

# Neutron multiplicity systematcs
alpha = ones(len(pdfids) - 2) * 0.4 * n_mult
