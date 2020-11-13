#!/bin/bash
ulimit -c 0 # Disable core dumps

model=$PJM_JOBNAME # Parse job name to get file name
outdir=fits/fit_$model
PY=~/miniconda3/envs/spectral/bin/python

spectral_dir=/disk02/usr6/giampaol/spectral/spectral_analysis
cd $spectral_dir
$PY fitlike.py $model $outdir
