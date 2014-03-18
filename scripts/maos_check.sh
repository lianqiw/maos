#!/bin/sh

#Sanity check maos
#NFIRAOS default
if [ -n "$1" ];then
    D=$1
else
    D=30
fi

echo > maos_check.stderr
echo "D is ${D}m. DM order is $((D*2))."
echo -n "LGS MCAO: "
echo $(./maos aper.d=$D dm.order=[$D*2 $D*2]  2>>maos_check.stderr) nm

echo -n "LGS MCAO $((D*4))x$((D*4)):  "
echo $(./maos aper.d=$D dm.order=[$D*4 $D*4] 2>>maos_check.stderr ) nm

echo -n "LGS MOAO:  "
echo $(./maos aper.d=$D dm.order=[$D*2 $D*2] evl.moao=0 moao.order=[$D] 2>>maos_check.stderr ) nm

echo -n "LGS GLAO:  "
echo $(./maos aper.d=$D dm_single.conf dm.order=[$D*2] recon.glao=1 wfs_lgs_only.conf 2>>maos_check.stderr ) nm

echo -n "NGS SCAO:      "
echo $(./maos aper.d=$D dm.order=[$D*2] -cscao_ngs.conf 2>>maos_check.stderr ) nm

echo -n "NGS MCAO: "
echo $(./maos aper.d=$D dm.order=[$D*2 $D*2] -cmcao_ngs.conf 2>>maos_check.stderr ) nm

echo -n "SCAO LGS: "
echo $(./maos aper.d=$D dm.order=[$D*2] -cscao_lgs.conf 2>>maos_check.stderr ) nm

echo -n "MVM mode:       "
echo $(./maos aper.d=$D dm.order=[$D*2 $D*2] atmr.os=[2 2 2 2 2 2] tomo.precond=1 tomo.maxit=100 fit.alg=0 2>>maos_check.stderr) nm
