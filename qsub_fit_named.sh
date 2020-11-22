#!/bin/bash
ulimit -c 0 # disable core dumps

spectral_dir=/disk02/usr6/giampaol/spectral/spectral_analysis
declare -a modelnames=("lma" "faild" "ksw" "malaney" "woosley" "horiuchi_xi2_5crit0_5" "horiuchi_xi2_5crit0_1" "nakazato_nh_min", "nakazato_ih_max" "const")

queue=ALL
maxjobs=200

cd $spectral_dir
for model in "${modelnames[@]}"; do

    jobname=$model
    outdir=fits/fit_$model
    out=$outdir/fit.log
    err=$outdir/fit.err
    mkdir -p $outdir

    jobrunning=`qstat -a $queue | grep $USER | wc -l`
    echo $jobrunning" jobs running"
    while [ "$jobrunning" -gt "$maxjobs" ]
    do
        jobrunning=`qstat -a $queue | grep $USER | wc -l`
    done
    echo "qsub -q $queue -o $out -e $err -r $jobname ./run_fit_named.sh"
    qsub -q $queue -o $out -e $err -r $jobname ./run_fit_named.sh
done
