#!/bin/bash
ulimit -c 0 # disable core dumps

#spectral_dir=/disk02/usr6/giampaol/spectral/spectral_analysis
spectral_dir=/disk02/usr6/elhedri/spectral_analysis
#declare -a modelnames=("volpe_nh" "volpe_ih")
declare -a modelnames=("lma" "faild" "ksw" "malaney" "woosley" "horiuchi_xi2_5_crit0_5" "horiuchi_xi2_5_crit0_1" "nakazato_nh_min" "nakazato_ih_max" "const" "horiuchi" "kawasaki" "ts0" "kresse_fiducial_ih" "kresse_fiducial_nh" "kresse_high_ih" "kresse_high_nh" "kresse_low_ih" "kresse_low_nh" "galais_ih" "galais_nh" "tabrizi_nh" "horiuchi21")
#declare -a modelnames=("horiuchi_xi2_5_crit0_5" "horiuchi_xi2_5_crit0_1" "nakazato_nh_min", "nakazato_ih_max")
#declare -a modelnames=("lma" "faild" "ksw" "malaney" "woosley" "horiuchi_xi2_5_crit0_5" "horiuchi_xi2_5_crit0_1" "nakazato_nh_min" "nakazato_ih_max" "const")

queue=all
maxjobs=200

cd $spectral_dir
for model in "${modelnames[@]}"; do

    jobname=$model
    outdir=fits/nospa/fit_$model
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
