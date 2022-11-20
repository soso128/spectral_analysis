from glob import glob
import numpy as np
import fitlike as ft
import matplotlib
from scipy.interpolate import interp1d
from ast import literal_eval
import sys
sys.path.append("spectrum_generator/")
import snspectrum as sns
matplotlib.use("pdf")
import matplotlib.pyplot as plt
plt.style.use("seaborn")

scan_dir = 'fits/'
fit_dirs = glob(scan_dir + '/fit_R*')

def combined_like(results):
    """ Combine SKI-IV likelihoods.
    Add the L - L_max for each SK era """
    liketot = results[0][:, -1] - results[0][:, -1].max() # SKI
    for r in results[1:]:
        liketot += r[:, -1] - r[:, -1].max() # SKII,III,IV
    return liketot

def analyse2d(like2d, rels, temps, final=False):
    " Calculate 2D exclusion "
    lmax = like2d.max()
    bestpos = np.argwhere(like2d == lmax)[0] # Best-fit T and rate
    fact = 0.5 if final else 1
    rels *= fact
    best_rel = rels[bestpos[1]]  # Best rate
    best_temp = temps[bestpos[0]] # Best temp
    norm = np.exp(like2d - lmax).sum()

    def likecut(cut):
        " Likelihood sum for given likelihood cut "
        lk = (like2d - lmax)
        lkcut = lk[lk > cut]
        return np.exp(lkcut).sum()

    cuts = np.arange(0.0, -3.0, -0.01) # Sample likelihood cuts
    sums = np.array([likecut(c) for c in cuts]) / norm
    cut = cuts[np.argmax(sums > 0.90)]  # Choose cut corresponding to 90% CL
    excluded = (like2d - lmax) < cut  # Boolean mask of excluded region
    return cut, excluded

def predicted_rate(Tnu, sf = 0):
    ergtoMeV = 1.6e-6
    modelname = r'{"type": 0, "tnu": ' + f'{Tnu}' + r', "lumi": ' +  f'{3e53 / ergtoMeV}' + r'}'
    imf = sns.imfs['salpeter']
    csfr = sns.csfr_fiducial
    if sf == 1:
        csfr = sns.csfr_up
    if sf == -1:
        csfr = sns.csfr_lower
    mod = literal_eval(modelname)
    en = np.arange(0.1, 100 + 0.1, 0.1)
    spec0 = np.array([sns.ibd_spectrum(ee, mod, imf, csfr) for ee in en])
    spec, en = spec0[~np.isnan(spec0)], en[~np.isnan(spec0)] # remove undefined
    en, spec_full = sns.smear_ibd_spectrum(en, np.column_stack((en, spec)), 4)
    totrate = spec_full[(en > 16) & (en < 90)].sum() * (en[1] - en[0])
    return totrate

temps, limits = [], []
like2d = []
rates = []
for dr in fit_dirs:
    tnu = float(dr[-2:])/10.0  # Temperature
    temps += [tnu]

    # Get likelihood information
    results = [np.genfromtxt(f"{dr}/fit_sk{i}.txt") for i in range(1, 5)]
    results = [r[:, [0, -1]] for r in results]
    rates = np.array(results[0][:,0])

    # Get single-model exclusion limits
    _, ratebest, ratepl, ratemin, ratelim = ft.combine(results)
    limits += [ratelim]

    # Add likelihood information to 2D array
    comblike = combined_like(results)
    like2d += [comblike]

like2d = np.array(like2d)  # 2D array of likelihoods for different temps and rates
print(like2d.shape)
np.savetxt("like2d_blackbody.txt", like2d)

# Sort temperatures
tempids = np.argsort(temps)
temps = np.array(temps)[tempids]
limits = np.array(limits)[tempids]
like2d = np.array(like2d)[tempids]

likecut, excluded2d = analyse2d(like2d, rates, temps, final=True)
limits2d = np.array([rates[np.argwhere(a == False)[-1]] for a in excluded2d])[:,0]


#### Plotting #####
matplotlib.rcParams["axes.labelsize"] = 22
matplotlib.rcParams["axes.titlesize"] = 22
matplotlib.rcParams["figure.figsize"] = (8.0, 7.0)
matplotlib.rcParams["legend.fontsize"] = 14
# matplotlib.rcParams["xtick.direction"] = 'in'
# matplotlib.rcParams["ytick.direction"] = 'in'
matplotlib.rcParams["xtick.major.size"] = 10
matplotlib.rcParams["ytick.major.size"] = 10
matplotlib.rcParams["xtick.labelsize"] = 14
matplotlib.rcParams["ytick.labelsize"] = 14
matplotlib.rcParams["xtick.major.width"] = 1.5
matplotlib.rcParams["ytick.major.width"] = 1.5
# matplotlib.rcParams["axes.axisbelow"] = 'line'
matplotlib.rcParams["axes.titleweight"] = "bold"
matplotlib.rcParams["axes.titlepad"] = 9.0
matplotlib.rcParams["axes.edgecolor"] = "black"
matplotlib.rcParams["axes.spines.top"] =   True
matplotlib.rcParams["axes.spines.right"] = True
matplotlib.rcParams["axes.linewidth"] = 1.5

plt.ylim(0,9)
plt.xlim(2,8.5)
plt.xlabel("T$_\\nu$ [MeV]")
plt.ylabel("event /22.5kton yr > 16 MeV")
plt.plot(temps, limits, linestyle="--", color='k',
         label="Individual 90% C.L. exclusion limits for each T$_\\nu$ considered separately")
plt.plot(temps, limits2d)
plt.fill_between(temps, [100]*len(temps), limits2d, color='tab:red',
                alpha=0.5, label="2D exclusion")
plt.legend()
plt.grid(which='minor')
#plt.title("Blackbody DSNB model exclusion")
plt.savefig("bexc.pdf")
plt.clf()

print("likemax: ", like2d.max(), like2d.shape)
like2d = np.where(like2d == -np.inf, like2d.max() - 100, like2d)
cs = plt.contour(temps, rates, (like2d-like2d.max()).T,
                 levels=[np.log(0.1)], colors=["tab:red"], alpha=0.0)
ymax = 9
plt.clf()
try:
    [pathlow, pathhigh] = cs.collections[0].get_paths()[0:2]
    vlow, vhigh = pathlow.vertices, pathhigh.vertices
    cxlow, cylow = vlow[:,0], vlow[:,1]
    cxhigh, cyhigh = vhigh[:,0], vhigh[:,1]
    plt.fill_between(cxlow, cylow, [-1]*len(cxlow), color="tab:red",
                     alpha=0.5)
    plt.fill_between(cxhigh, cyhigh, [25]*len(cxhigh), color="tab:red"
                     , alpha=0.5, label="2D Blackbody DSNB exclusion")
except ValueError:
    pathhigh = cs.collections[0].get_paths()[0]
    vhigh = pathhigh.vertices
    cxhigh, cyhigh = vhigh[:,0], vhigh[:,1]
    plt.fill_between(cxhigh, cyhigh, [25]*len(cxhigh), color="tab:red"
                     , alpha=0.5, label="2D Blackbody DSNB exclusion")
plt.plot(temps, limits, linestyle="--", color='k',
         label="90% C.L. upper limits for individual T$_\\nu$s ")
Tnus = np.arange(2.0, 8.6, 0.1)
#pred_low = np.array([predicted_rate(tnu, -1) for tnu in Tnus])
#pred_fid = np.array([predicted_rate(tnu, 0) for tnu in Tnus])
#pred_up = np.array([predicted_rate(tnu, 1) for tnu in Tnus])
#predictions = np.column_stack((Tnus,pred_low,pred_fid,pred_up))
#np.savetxt("predictions.txt", predictions)
predictions = np.loadtxt("predictions.txt")
pred_low, pred_fid, pred_up = predictions[:, 1], predictions[:, 2], predictions[:, 3]
plt.plot(Tnus, pred_fid, '--', color = 'tab:blue', linewidth = 1.5)
plt.fill_between(Tnus, pred_low, pred_up, alpha = 0.5, color = 'tab:blue', 
             label = "Predicted rates")
plt.ylim(0,ymax)
plt.xlim(2,8.5)
plt.xlabel("T$_\\nu$ [MeV]")
#plt.title("Blackbody DSNB model exclusion")
plt.ylabel("event /22.5kton yr > 16 MeV")
plt.legend(loc="upper left", frameon=True, fancybox=False)
plt.grid(which='minor')
plt.savefig("bexc_contour.pdf")
plt.clf()
ymax = 4e53
plt.plot(temps, limits * 5e52/pred_fid[::5], linestyle="--", color='k',
         label="90% C.L. upper limits for individual T$_\\nu$s ")
frate = interp1d(Tnus, pred_fid, bounds_error = False, fill_value = (pred_fid[0], pred_fid[-1]))
plt.fill_between(cxhigh, cyhigh * 5e52/frate(cxhigh), [2*ymax]*len(cxhigh), color="tab:red"
                 , alpha=0.5, label="2D Blackbody DSNB exclusion")
fkam = np.loadtxt("kamioka_sn1987a.csv")
fimb = np.loadtxt("imb_sn1987a.csv")
plt.fill(fkam[:, 0], fkam[:, 1]/6. * 1e53, 'tab:blue', alpha = 0.5)
plt.fill(fimb[:, 0], fimb[:, 1]/6. * 1e53, 'green', alpha = 0.5)
plt.text(3.3,1.05e53,"IMB")
plt.text(2.2,0.7e53,"Kamioka")
plt.ylim(0,ymax)
plt.xlim(2,8.5)
plt.xlabel("T$_\\nu$ [MeV]")
#plt.title("Blackbody DSNB model exclusion")
plt.ylabel("SN $\\overline{\\nu}_e$ energy in $10^{53}$erg")
plt.legend(loc="upper left", frameon=True, fancybox=False)
plt.grid(which='minor')
plt.text(5.3,3.2e53,"SK 1497+794+562+2970 Days")
plt.text(5.3, 2.9e53, "Excluded (E$_\\nu > 17.3~$MeV)")
plt.text(5.3, 2.6e53, "90% C.L.")
plt.savefig("snexc_contour.pdf")
#print("Pred rate 6 MeV: ", predicted_rate(4.0, -1))
