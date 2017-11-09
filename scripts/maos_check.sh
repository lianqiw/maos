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
case $D in
    2)
	#4/4/2017
	REF=(90.1682 91.9371 166.001 133.074 133.682 140.137 130.599 138.509 131.078 589.324 225.003 588.795 220.674 117.774 120.16 588.795 211.419 150.461 135.644)
	;;
    5)
	#6/7/2017
	REF=(107.10 108.90 137.96 132.96 132.58 132.10 132.73 132.85 121.69 510.76 259.97 104.77 107.74 125.48 130.55 121.69 122.34 253.28 136.87)
	;;
    30)
	#6/7/2017
	REF=(113.724 115.872 133.825 167.16 133.935 Unknown 114.086 117.495 365.903 110.422 110.985 144.098 145.37 362.133 362.086 126.993)
esac

echo > maos_check.stderr
echo "D is ${D}m. DM order is $((D*2))."
echo -n "Ideal fit (cpu)  "
ii=0
RMS[$ii]=$(./maos $args aper.d=$D sim.idealfit=1 -g-1 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "Ideal tomo (cpu) "
RMS[$ii]=$(./maos $args aper.d=$D sim.idealtomo=1 -g-1 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (inte): "
RMS[$ii]=$(./maos $args aper.d=$D recon.split=0 tomo.precond=0 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (CG):   "
RMS[$ii]=$(./maos $args aper.d=$D tomo.precond=0 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (FDPCG):"
RMS[$ii]=$(./maos $args aper.d=$D tomo.precond=1 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (CBS):  "
RMS[$ii]=$(./maos $args aper.d=$D tomo.alg=0 fit.alg=0 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

if [ $D -le 10 ];then
echo -n "LGS MCAO (SVD):  "
RMS[$ii]=$(./maos $args aper.d=$D tomo.alg=2 fit.alg=2 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (MVM):  "
RMS[$ii]=$(./maos $args aper.d=$D atmr.os=[2] tomo.precond=1 tomo.maxit=100 fit.alg=0 recon.mvm=1 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 
fi

echo -n "LGS MOAO:        "
RMS[$ii]=$(./maos $args aper.d=$D evl.moao=0 moao.dx=[1/2] 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS GLAO (inte): "
RMS[$ii]=$(./maos $args aper.d=$D dm_single.conf  recon.glao=1 recon.split=0 wfs_lgs_ttf.conf 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS GLAO (split):"
RMS[$ii]=$(./maos $args aper.d=$D dm_single.conf  recon.glao=1 recon.split=1 wfs_lgs_ttf.conf 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "NGS SCAO (inte): "
RMS[$ii]=$(./maos $args aper.d=$D  -cscao_ngs.conf recon.split=0 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "NGS SCAO (split):"
RMS[$ii]=$(./maos $args aper.d=$D  -cscao_ngs.conf recon.split=1 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "NGS MCAO (inte): "
RMS[$ii]=$(./maos $args aper.d=$D -cmcao_ngs.conf recon.split=0 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "NGS MCAO (split):"
RMS[$ii]=$(./maos $args aper.d=$D -cmcao_ngs.conf recon.split=1 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "SCAO LGS (inte): "
RMS[$ii]=$(./maos $args aper.d=$D  -cscao_lgs.conf recon.split=0 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "SCAO LGS (split):"
RMS[$ii]=$(./maos $args aper.d=$D  -cscao_lgs.conf recon.split=1 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS LTAO (inte):"
RMS[$ii]=$(./maos $args aper.d=$D dm_single.conf fov_oa.conf recon.split=0 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS LTAO (split):"
RMS[$ii]=$(./maos $args aper.d=$D dm_single.conf fov_oa.conf recon.split=1 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo ${RMS[*]}

