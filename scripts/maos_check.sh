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
	REF=(89.4582 91.262 133.036 166.016 133.689 139.809 130.262 139.625 113.996 130.983 255.921 588.795 220.674 117.581 121.81 211.419 588.795)
	;;
    5)
	#4/4/2017
	REF=(106.928 108.742 131.49 300.934 131.598 130.543 131.413 131.047 113.377 120.46 260.973 103.228 105.737 123.973 128.552 121.414 121.147)
	;;
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

echo -n "LGS MCAO (split):"
RMS[$ii]=$(./maos $args aper.d=$D  2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MCAO (inte): "
RMS[$ii]=$(./maos $args aper.d=$D recon.split=0 2>>maos_check.stderr)
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
echo -n "LGS MCAO ($((D*4))x$((D*4))):  "
RMS[$ii]=$(./maos $args aper.d=$D dm.dx=[0.25 0.25] 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS MOAO:        "
RMS[$ii]=$(./maos $args aper.d=$D evl.moao=0 moao.dx=[1/2] 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS GLAO:        "
RMS[$ii]=$(./maos $args aper.d=$D dm_single.conf  recon.glao=1 wfs_lgs_ttf.conf 2>>maos_check.stderr )
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

echo -n "SCAO LGS (split):"
RMS[$ii]=$(./maos $args aper.d=$D  -cscao_lgs.conf recon.split=1 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "SCAO LGS (inte): "
RMS[$ii]=$(./maos $args aper.d=$D  -cscao_lgs.conf recon.split=0 2>>maos_check.stderr )
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo -n "LGS LTAO (split):"
RMS[$ii]=$(./maos $args aper.d=$D dm_single.conf fov_oa.conf 2>>maos_check.stderr)
echo ${RMS[$ii]} nm, Ref: ${REF[$ii]} nm
ii=$((ii+1)) 

echo ${RMS[*]}

