import fitlike as ft
#from os import environ
from sys import argv

"""Simulate DSNB with blackbody emission spectrum
as in arXiv:0812:3157"""

ergtoMeV = 1.6e-6

# Get neutrino temperature
nutemp = float(argv[1][1:])/10 # Given multiplied by 10 to avoid having a dot in job name

# Make model
modelname = r'{"type": 0, "tnu": ' + f'{nutemp}' + r', "lumi": ' +  f'{3e53 * ergtoMeV}' + r'}'

# Make output folder
directory = argv[2]

# Run fit
ft.fullike(modelname, elow=16, ehigh=90, elow_sk2 = 17.5, elow_sk4=20, ehigh_sk4=80, elow_sk4_1n=16, outdir=directory)
