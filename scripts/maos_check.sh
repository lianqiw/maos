#!/bin/bash

#Sanity check maos
#NFIRAOS default
if [ -n "$1" ];then
    D=$1
    shift
    args="$@"
else
    D=30
    args=
fi
args+=" aper.d=$D"
case $D in
    2)
	#4/4/2017
	#REF=(90.1682 91.9371 166.001 133.074 133.682 140.137 130.599 138.509 131.078 589.324 225.003 588.795 220.674 117.774 120.16 588.795 211.419 150.461 135.644)
	#10/29/2018
	REF=(88.44 88.96 143.03 140.22 140.82 155.58 137.95 155.13 137.85 589.32 226.63 588.79 218.41 117.42 119.60 588.79 209.86 149.96 138.61)
	;;
    5)
	#6/7/2017
	REF=(107.10 108.90 137.96 132.96 132.58 132.10 132.73 132.85 121.69 510.76 259.97 104.77 107.74 125.48 130.55 121.69 122.34 253.28 136.87)
	;;
    30)
	#9/14/2018
	REF=(112.55 113.02 133.92 139.07 132.86 Unknown 116.21 365.25 364.84 109.83 109.82 145.30 146.51 361.19 361.10 154.76 125.21)

esac

echo > maos_check.log

function run_maos(){
    ./maos $args "$*" >> maos_check.log
    if [ $? != 0 ];then
	echo 'error'
    else
	tail -n5 maos_check.log |grep 'Mean:' |cut -d ':' -f 2
    fi
}

echo "D is ${D}m. DM order is $((D*2))."

ii=0
echo -n "Ideal fit (cpu)  "
RMS[$ii]=$(run_maos sim.idealfit=1 -g-1 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "Ideal tomo (cpu) "
RMS[$ii]=$(run_maos sim.idealtomo=1 -g-1 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (inte): "
RMS[$ii]=$(run_maos recon.split=0 tomo.precond=0 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (CG):   "
RMS[$ii]=$(run_maos tomo.precond=0 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (FDPCG):"
RMS[$ii]=$(run_maos tomo.precond=1 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (CBS):  "
RMS[$ii]=$(run_maos tomo.alg=0 fit.alg=0 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

if [ $D -le 10 ];then
echo -n "LGS MCAO (SVD):  "
RMS[$ii]=$(run_maos tomo.alg=2 fit.alg=2 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (MVM):  "
RMS[$ii]=$(run_maos atmr.os=[2] tomo.precond=1 tomo.maxit=100 fit.alg=0 recon.mvm=1 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 
fi

echo -n "LGS MOAO:        "
RMS[$ii]=$(run_maos evl.moao=0 moao.dx=[1/2] )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS GLAO (inte): "
RMS[$ii]=$(run_maos dm_single.conf  recon.glao=1 recon.split=0 wfs_lgs_ttf.conf )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS GLAO (split):"
RMS[$ii]=$(run_maos dm_single.conf  recon.glao=1 recon.split=1 wfs_lgs_ttf.conf )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "NGS SCAO (inte): "
RMS[$ii]=$(run_maos -cscao_ngs.conf recon.split=0 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "NGS SCAO (split):"
RMS[$ii]=$(run_maos -cscao_ngs.conf recon.split=1 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "NGS MCAO (inte): "
RMS[$ii]=$(run_maos -cmcao_ngs.conf recon.split=0 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "NGS MCAO (split):"
RMS[$ii]=$(run_maos -cmcao_ngs.conf recon.split=1 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "SCAO LGS (inte): "
RMS[$ii]=$(run_maos -cscao_lgs.conf recon.split=0 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "SCAO LGS (split):"
RMS[$ii]=$(run_maos -cscao_lgs.conf recon.split=1 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS LTAO (inte): "
RMS[$ii]=$(run_maos dm_single.conf fov_oa.conf recon.split=0 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS LTAO (split):"
RMS[$ii]=$(run_maos dm_single.conf fov_oa.conf recon.split=1 )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo ${RMS[*]}

