from glob import glob
import numpy as np
import fitlike as ft
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
plt.style.use("seaborn")

scan_dir = 'fits0/fits_Tscan'
fit_dirs = glob(scan_dir + '/*')

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
plt.xlabel("T$_\\nu$ (MeV)")
plt.ylabel("event /22.5kton yr > 16 MeV")
plt.plot(temps, limits, linestyle="--", color='k',
         label="Individual 90% C.L. exclusion limits for each T$_\\nu$ considered separately")
plt.plot(temps, limits2d)
plt.fill_between(temps, [100]*len(temps), limits2d, color='tab:red',
                alpha=0.5, label="2D exclusion")
plt.legend()
plt.grid(which='minor')
plt.title("Blackbody DSNB model exclusion")
plt.savefig("bexc.pdf")
plt.clf()

cs = plt.contour(temps, rates, (like2d-like2d.max()).T,
                 levels=[likecut], colors=["tab:red"], alpha=0.0)
[pathlow, pathhigh] = cs.collections[0].get_paths()[0:2]
vlow, vhigh = pathlow.vertices, pathhigh.vertices
cxlow, cylow = vlow[:,0], vlow[:,1]
cxhigh, cyhigh = vhigh[:,0], vhigh[:,1]
plt.plot(temps, limits, linestyle="--", color='k',
         label="90% C.L. upper limits for individual T$_\\nu$s ")
plt.fill_between(cxlow, cylow, [-1]*len(cxlow), color="tab:red",
                 alpha=0.5, label="2D Blackbody DSNB exclusion")
plt.fill_between(cxhigh, cyhigh, [10]*len(cxhigh), color="tab:red", alpha=0.5)
plt.ylim(0,9)
plt.xlim(2,8.5)
plt.xlabel("T$_\\nu$ (MeV)")
plt.title("Blackbody DSNB model exclusion")
plt.ylabel("event /22.5kton yr > 16 MeV")
plt.legend(loc="upper left", frameon=True, fancybox=False)
plt.grid(which='minor')
plt.savefig("bexc_contour.pdf")