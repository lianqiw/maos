#!/bin/bash

#Sanity check maos
#NFIRAOS default
if [ -n "$1" ];then
    D=$1
    shift
    args="aper.d=$D $@"
else
    D=30
fi

case $D in
    2)
	#4/4/2017
	#REF=(90.1682 91.9371 166.001 133.074 133.682 140.137 130.599 138.509 131.078 589.324 225.003 588.795 220.674 117.774 120.16 588.795 211.419 150.461 135.644)
	#10/29/2018
	#REF=(0 88.44 88.96 143.03 140.22 140.82 155.58 137.95 155.13 137.85 589.32 226.63 588.79 218.41 117.42 119.60 588.79 209.86 149.96 138.61)
	#2/18/2020
	REF=(539.03 91.84 90.60 138.08 148.49 150.29 164.92 148.59 171.49 147.71 190.32 185.65 100.20 101.63 120.20 125.10 103.77 105.66 162.81 129.65 102.30 120.43)
	;;
    5)
	#6/7/2017
	#REF=(107.10 108.90 137.96 132.96 132.58 132.10 132.73 132.85 121.69 510.76 259.97 104.77 107.74 125.48 130.55 121.69 122.34 253.28 136.87)
	#10/15/2019
	#REF=(1078 106.37 106.75 133.98 138.84 138.46 141.10 139.85 140.25 127.97 244.87 251.29 102.67 103.99 121.63 124.67 116.55 123.19 131.01 137.78 110.64)
	#2/18/2020
	REF=(1078.03 107.77 106.77 136.44 138.98 139.46 140.80 140.25 140.17 129.38 247.20 251.37 102.50 104.05 121.59 124.61 118.64 125.63 139.27 139.40 106.61 110.21)

	;;
    10)
	#12/13/2018
	#REF=(108.80 109.09 127.83 129.54 130.46 129.81 error 129.03 115.42 353.55 287.17 104.30 104.76 139.21 143.33 195.46 165.91 217.42 128.84)
	#2/11/2019
	#REF=(109.16 109.45 132.43 129.49 129.76 129.59 129.73 129.28 114.41 365.81 277.64 106.79 107.76 143.38 145.59 193.78 166.22 215.75 126.69)
	#10/15/2019
	REF=(1208 109.13 109.46 123.31 129.72 129.97 131.27 131.55 129.88 114.93 274.21 277.10 103.01 104.65 121.70 132.94 159.73 166.28 117.76 124.77 106.46 107.93)

	;;
    30)
	#9/14/2018
	#REF=(112.55 113.02 133.92 139.07 132.86 error 116.21 365.25 364.84 109.83 109.82 145.30 146.51 361.19 361.10 154.76 125.21)
	#12/13/2018: with sim.mffocus=1 as default
	#REF=(112.68 112.93 134.55 141.49 134.51 error 118.06 375.25 370.70 109.01 109.02 146.39 147.32 343.20 368.34 155.83 143.43)
	#2/11/2019
	#REF=(112.68 112.93 134.58 141.23 133.94 146.73 117.34 375.55 370.68 108.96 109.09 146.35 147.25 343.93 367.85 156.07 138.63)
	#5/6/2019
	#REF=(112.69 112.98 134.55 141.07 134.05 140.53 117.68 374.64 370.47 109.00 109.25 146.04 146.90 343.89 367.74 154.98 136.58)
	#8/13/2019
	#REF=(112.92 113.22 134.20 140.38 136.00 139.92 119.87 368.31 373.93 109.68 111.99 128.41 143.92 346.34 367.57 113.80 141.73)
	#9/6/2019
	#REF=(112.92 113.22 133.02 138.40 136.19 140.06 138.99 120.21 368.22 373.90 110.13 110.04 127.65 134.11 346.38 367.34 113.30 142.55)
	#10/15/2019
	#REF=(1328.69 112.92 113.22 134.30 140.45 136.19 139.79 138.09 120.07 368.10 373.87 107.94 108.27 125.10 131.75 346.19 367.39 113.53 141.83 113.71 113.60)
	#5/4/2021
	REF=(1712.96 113.49 113.67 136.95 144.61 140.54 144.34 143.28 125.27 329.29 336.53 112.67 114.71 129.84 153.11 345.08 365.68 122.73 151.48 117.29 118.54 144.87 214.76 223.92 212.80 94.59)
	;;
    *)
	REF=()
esac
fnlog=maos_check_${D}.log #log of all
fntmp=maos_check_${D}.tmp #log of current simulation
fnres=maos_check_${D}.res #result summary
if [ -f $fnlog ];then
    mv $fnlog ${fnlog/.log/.bak}
fi
ii=0
printf "%-20s   Res   Ref     %%\n" "D=${D}m DM is $((D*2))x$((D*2))" | tee $fnres
function run_maos(){
	aotype=$1
	shift
    ./maos sim.end=100 $args "$*" > $fntmp
	echo $aotype $* >> $fnlog
	cat $fntmp >> $fnlog
    if [ $? == 0 ];then
		RMS[ii]=$(grep 'Mean:' $fntmp |tail -n1 |cut -d ' ' -f 2)
		a=${RMS[$ii]%.*}
    else
		RMS[ii]='error'
		a=0
		ans=1
    fi

    b=${REF[$ii]%.*}
	diff=100
    if [ "$a" != "0" -a "$b" != "error" ];then
		diff=$(((a-b)*100/a))
    fi
	printf "%-20s %5.0f %5.0f %4d%%\n" "$aotype" ${RMS[$ii]} ${REF[$ii]} "$diff" | tee -a $fnres
	ii=$((ii+1)) 
}
export MAOS_LOG_LEVEL=-1

run_maos "Openloop:        " sim.evlol=1

run_maos "Ideal fit:       " sim.idealfit=1 

run_maos "Ideal tomo:      " sim.idealtomo=1 

run_maos "LGS MCAO (inte): " recon.split=0 tomo.precond=0

run_maos "LGS MCAO (CG):   " tomo.precond=0 cn2.pair=[0 1 2 5] recon.psd=1 tomo.assemble=0 fit.assemble=1

run_maos "LGS MCAO (FDPCG):" tomo.precond=1 tomo.assemble=1 fit.assemble=0

run_maos "LGS MCAO (CBS):  " tomo.alg=0 fit.alg=0 atmr.os=[2 2 1 1 1 1 1]

if [ $D -le 10 ];then
run_maos "LGS MCAO (SVD):  " tomo.alg=2 fit.alg=2 atmr.os=[2 2 1 1 1 1 1] gpu.tomo=0
fi
run_maos "LGS MCAO (MVM):  " atmr.os=[2] tomo.precond=1 tomo.maxit=100 fit.alg=0 recon.mvm=1

run_maos "LGS MOAO:        " evl.moao=0 moao.dx=[1/2]

run_maos "LGS GLAO (inte): " dm_single.conf  recon.glao=1 recon.split=0 wfs_lgs_ttf.conf

run_maos "LGS GLAO (split):" dm_single.conf  recon.glao=1 recon.split=1 wfs_lgs_ttf.conf

run_maos "NGS SCAO (inte): " -cscao_ngs.conf recon.split=0

run_maos "NGS SCAO (split):" -cscao_ngs.conf recon.split=1

run_maos "NGS MCAO (inte): " -cmcao_ngs.conf recon.split=0

run_maos "NGS MCAO (split):" -cmcao_ngs.conf recon.split=1

run_maos "SCAO LGS (inte): " -cscao_lgs.conf recon.split=0

run_maos "SCAO LGS (split):" -cscao_lgs.conf recon.split=1

run_maos "LGS LTAO (inte): " dm_single.conf fov_oa.conf recon.split=0

run_maos "LGS LTAO (split):" dm_single.conf fov_oa.conf recon.split=1 

run_maos "NGS SCAO (lsq,inte)" -cscao_ngs.conf recon.split=0 recon.alg=1

run_maos "NGS SCAO (lsq,split)" -cscao_ngs.conf recon.split=1 recon.alg=1

run_maos "LGS MCAO PCCD:  " tomo.precond=0 cn2.pair=[0 1 2 5] recon.psd=1 powfs.radpix=[16,0,0] powfs.pixpsa=[6,0,0]

run_maos "LGS MCAO SL:    " tomo.precond=0 cn2.pair=[0 1 2 5] recon.psd=1 powfs.fnllt=['llt_SL.conf',,] powfs.pixpsa=[16,0,0]

run_maos "LGS MCAO PCCD SL:" tomo.precond=0 cn2.pair=[0 1 2 5] recon.psd=1 powfs.radpix=[18,0,0] powfs.pixpsa=[6,0,0] powfs.fnllt=['llt_SL.conf',,]
if [ $D -eq 30 ];then
run_maos "NFIRAOS LGS: " nfiraos_lgs.conf
run_maos "NFIRAOS PWFS:" nfiraos_ngs.conf
fi

echo ${RMS[*]}

exit $ans
