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

echo > maos_check.stderr
echo "D is ${D}m. DM order is $((D*2))."
echo -n "Ideal fit (cpu)  "
echo $(./maos $args aper.d=$D sim.idealfit=1 -g-1 2>>maos_check.stderr) nm

echo -n "Ideal tomo (cpu) "
echo $(./maos $args aper.d=$D sim.idealtomo=1 -g-1 2>>maos_check.stderr) nm

echo -n "LGS MCAO (split):"
echo $(./maos $args aper.d=$D  2>>maos_check.stderr) nm

echo -n "LGS MCAO (inte): "
echo $(./maos $args aper.d=$D recon.split=0 2>>maos_check.stderr) nm

echo -n "LGS MCAO (FDPCG):"
echo $(./maos $args aper.d=$D tomo.precond=1 2>>maos_check.stderr) nm

echo -n "LGS MCAO (CBS):  "
echo $(./maos $args aper.d=$D tomo.alg=0 fit.alg=0 2>>maos_check.stderr) nm

if [ $D -le 10 ];then
echo -n "LGS MCAO (SVD):  "
echo $(./maos $args aper.d=$D tomo.alg=2 fit.alg=2 2>>maos_check.stderr) nm

echo -n "LGS MCAO (MVM):  "
echo $(./maos $args aper.d=$D atmr.os=[2 2 2 2 2 2] tomo.precond=1 tomo.maxit=100 fit.alg=0 2>>maos_check.stderr) nm
fi
echo -n "LGS MCAO ($((D*4))x$((D*4))):"
echo $(./maos $args aper.d=$D dm.dx=[0.25 0.25] 2>>maos_check.stderr ) nm

echo -n "LGS MOAO:        "
echo $(./maos $args aper.d=$D evl.moao=0 moao.order=[$D] 2>>maos_check.stderr ) nm

echo -n "LGS GLAO:        "
echo $(./maos $args aper.d=$D dm_single.conf  recon.glao=1 wfs_lgs_only.conf 2>>maos_check.stderr ) nm

echo -n "NGS SCAO (inte): "
echo $(./maos $args aper.d=$D  -cscao_ngs.conf recon.split=0 2>>maos_check.stderr ) nm

echo -n "NGS SCAO (split):"
echo $(./maos $args aper.d=$D  -cscao_ngs.conf recon.split=1 2>>maos_check.stderr ) nm

echo -n "NGS MCAO (inte): "
echo $(./maos $args aper.d=$D -cmcao_ngs.conf recon.split=0 2>>maos_check.stderr ) nm

echo -n "NGS MCAO (split):"
echo $(./maos $args aper.d=$D -cmcao_ngs.conf recon.split=1 2>>maos_check.stderr ) nm

echo -n "SCAO LGS (split):"
echo $(./maos $args aper.d=$D  -cscao_lgs.conf recon.split=1 2>>maos_check.stderr ) nm

echo -n "SCAO LGS (inte): "
echo $(./maos $args aper.d=$D  -cscao_lgs.conf recon.split=0 2>>maos_check.stderr ) nm

