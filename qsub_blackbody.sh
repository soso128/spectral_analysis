#!/bin/bash
ulimit -c 0 # disable core dumps

spectral_dir=/disk02/usr6/elhedri/spectral_analysis

queue=all
maxjobs=200

cd $spectral_dir
for tnu in `seq 20 5 85`; do

    jobname=$tnu
    outdir=fits/fit_R$tnu
    out=$outdir/fit.log
    err=$outdir/fit.err
    mkdir -p $outdir

    jobrunning=`qstat -a $queue | grep $USER | wc -l`
    echo $jobrunning" jobs running"
    while [ "$jobrunning" -gt "$maxjobs" ]
    do
        jobrunning=`qstat -a $queue | grep $USER | wc -l`
    done
    echo "qsub -q $queue -o $out -e $err -r $jobname ./run_blackbody.sh"
    qsub -q $queue -o $out -e $err -r $jobname ./run_blackbody.sh
done
