#!/bin/bash
if [ ! -f "../bin/maos" ];then
	echo "Please run in maos build directory"
	exit
fi

unset MAOS_PYWFS_DEBUG 
export MAOS_LOG_LEVEL=-1 #reduce log output
D=30
#Sanity check maos
#NFIRAOS default
shopt -s extglob
case "$1" in
	(+([0-9.]))
	D=$1
	shift;;
esac
if [ "$D" = "30" ];then
	args="tmt.conf "
elif [ "$D" = "10" ];then
	args="keck.conf"
else
    args="aper.d=$D "
fi
args+=" $@"

ARCH=CPU
case "$@" in
  *-g-1* | *-G0*)
	;;
  *)
	if nvidia-smi -i 0 >/dev/null 2>&1 ;then
		ARCH=GPU
	fi
		;;
esac

case $D in
    5)
	REF_GPU=(1153.3 105.6 105.8 120.2 122.1 120.7 117.5 124.8 139.6 612.8 616.3 127.4 116.4 136.0 141.6 146.0 128.6 124.4 144.2 153.9 154.0 155.3 153.2 155.2 153.2 154.9 155.7 0.3 14.8) #2026-07-02 hesperus D=5
	REF_CPU=(1153.3 104.3 104.1 121.4 122.1 120.3 117.3 126.1 139.1 612.8 616.1 130.1 115.8 132.3 141.6 146.2 123.3 129.4 152.4 150.1 154.5 154.4 158.7 154.2 156.5 157.1 157.0 0.3 14.9) #2026-07-02 hesperus D=5
	;;
    10)
	REF_GPU=(1457.7 112.3 111.9 126.1 119.8 128.2 115.6 187.8 203.8 949.2 958.0 124.5 122.2 133.3 141.0 143.2 128.8 118.9 145.6 156.7 156.4 156.8 156.6 157.1 155.6 151.3 164.8 1.6 21.9) #2026-07-02 hesperus D=10
	REF_CPU=(1431.3 108.8 109.5 118.8 119.2 119.4 116.3 178.5 194.1 933.9 941.8 120.7 116.4 124.1 139.0 139.7 127.0 117.5 144.7 150.0 149.4 150.1 152.0 155.1 151.0 149.3 151.6 1.6 21.9) #2025-04-24 kepler D=10
	;;
    30)	
	REF_GPU=(1819.6 118.8 117.6 125.2 123.3 144.2 127.7 366.8 388.7 814.2 825.3 128.7 126.2 133.2 145.6 152.3 132.2 129.9 155.0 157.9 157.2 160.6 158.9 157.9 156.5 159.6 0.6 50.2 212.8 165.2) #2026-07-02 hesperus D=30
	REF_CPU=(1819.6 119.0 118.1 125.2 123.4 144.2 127.8 367.6 389.4 814.0 823.4 128.8 126.9 132.7 145.6 147.9 132.1 130.5 155.3 156.6 157.2 160.3 000.0 157.6 156.9 159.0 0.7 53.7 220 0   ) #2025-04-24 kepler D=30
	;;
    *)
	REF_GPU=()
	REF_CPU=()
	;;
esac
echo "Using $ARCH"
if [ "$ARCH" = GPU ];then
	REF=(${REF_GPU[@]})
	REF2=(${REF_CPU[@]})
else
	REF=(${REF_CPU[@]})
	REF2=(${REF_GPU[@]})
fi
fnpre=maos_check_${ARCH}_${D}
fnall=${fnpre}.all #log of all simulations
fnfail=${fnpre}.failed #log of all simulations
fnerr=${fnpre}.err #log of all failed simulations
fnlog=${fnpre}.log #log of current simulation
fnres=${fnpre}.res #result summary
fnref=maos_check.ref #all ference values
if [ -f $fnres ];then
	mv ${fnres} ${fnres}_$(date +%Y%m%d-%H%M%S)
fi
echo $(date) > $fnall
echo $(date) > $fnerr
ans=0 #result code
ii=0
s_start0=`date +%s`
printf "%-20s    Res    Ref   Ref2     %%     %%  time  total\n" "D=${D}m" | tee $fnres
function run_maos(){
	aotype=$1
	shift
	s_start=`date +%s`
    ../bin/maos sim.end=100 $* $args >$fnlog 2>$fnerr
	a=
    if [ $? -eq 0 ];then
		a=`printf %.1f $(grep 'Mean:' $fnlog |tail -n1 |cut -d ' ' -f 2)`
	else
		echo Failed to run maos sim.end=100 $* $args 
	fi
	if [ x$a = x ];then
		a=000.0
		ans=$((ans+1)) #failed to run
		echo $aotype $* >> $fnfail
		cat $fnlog >> $fnfail
		cat $fnerr >> $fnfail
	fi
	RMS[ii]=$a
    s_end=`date +%s`
	s_diff=$((s_end-s_start)) #time per simulation
	s_diff2=$((s_end-s_start0)) #cumulative time
	echo $aotype $* >> $fnall
	cat $fnlog >> $fnall
	cat $fnerr >> $fnall
	echo "" >> $fnall
	rm $fnlog
	rm $fnerr
    b=${REF[$ii]:-0}
    b2=${REF2[$ii]:-0}
	diff=$(echo "200*($a-$b)/($a+$b+10)" | bc)
	diff2=$(echo "200*($a-$b2)/($a+$b2+10)" | bc)
	if [ $diff -gt 10 -o $diff -lt -10 ];then
		ans=$((ans+1)) #mark failure
	fi
	printf "%-20s %6.1f %6.1f %6.1f %5.1f %5.1f %5d %6d\n" "$aotype" "$a" "$b" "$b2" "$diff" "$diff2" "$s_diff" "$s_diff2"| tee -a $fnres
	ii=$((ii+1)) 
}
function run_maos_gpu(){
	if [ "$ARCH" = GPU ];then
		run_maos "$@"
	else
		printf "%-20s   skipped in CPU mode.\n" "$1" | tee -a $fnres
		RMS[ii]=000.0
		ii=$((ii+1)) 
	fi
}
{
run_maos "Openloop:        " -cmcao_lgs.conf sim.evlol=1

run_maos "NGS SCAO (inte): " -cscao_ngs.conf recon.split=0

run_maos "NGS SCAO (ahst): " -cscao_ngs.conf recon.split=1

run_maos "NGS SCAO (lsq):  " -cscao_ngs.conf recon.alg=1

run_maos "NGS SCAO (modal):" -cscao_ngs.conf recon.modal=1 recon.alg=1

run_maos "NGS PWFS (zonal):" -cscao_pywfs.conf recon.modal=0 powfs.dx=[1/16] sim.end=500

run_maos "NGS PWFS (modal):" -cscao_pywfs.conf recon.modal=1 powfs.dx=[1/16] sim.end=500

run_maos "LGS SCAO (inte): " -cscao_lgs.conf recon.split=0

run_maos "LGS SCAO (ahst): " -cscao_lgs.conf recon.split=1

run_maos "LGS GLAO (inte): " -cglao.conf recon.split=0 evl.psfmean=0

run_maos "LGS GLAO (ahst): " -cglao.conf recon.split=1 evl.psfmean=0

run_maos "LGS LTAO (inte): " -cmcao_lgs.conf dm_single.conf fov_oa.conf recon.split=0

run_maos "LGS LTAO (ahst) :" -cmcao_lgs.conf dm_single.conf fov_oa.conf recon.split=1 powfs.astscale=[1 0.1 0.1] #ahst with wide NGS asterism suffers (from aliasing error?)

run_maos "LGS MOAO (ahst): " -cmcao_lgs.conf recon.split=1 evl.moao=0 moao.dx=[1/2]

run_maos "NGS MCAO (inte): " -cmcao_ngs.conf recon.split=0

run_maos "NGS MCAO (ahst): " -cmcao_ngs.conf recon.split=1 

run_maos "LGS MCAO (fit):  " -cmcao_lgs.conf sim.idealtomo=1

run_maos "LGS MCAO (tomo): " -cmcao_lgs.conf evl.tomo=1 

run_maos "LGS MCAO (inte): " -cmcao_lgs.conf recon.split=0 tomo.assemble=1 fit.assemble=1 #also test assembled matrices

run_maos "LGS MCAO (CG):   " -cmcao_lgs.conf tomo.precond=0 tomo.assemble=1 fit.assemble=1 powfs.dtrat=[1 4 2] #also test assembled matrices

run_maos "LGS MCAO (FDPCG):" -cmcao_lgs.conf tomo.precond=1 sim.fuseint=0 powfs.dtrat=[1 4 2]

run_maos "LGS MCAO (CBS):  " -cmcao_lgs.conf tomo.alg=0 fit.alg=0 atmr.os=[2 2 1 1 1 1 1] recon.split=1 #split=0 CBS does not work well in GPU

if [ ${D%.*} -le 10 ];then
run_maos "LGS MCAO (SVD):  " -cmcao_lgs.conf tomo.alg=2 fit.alg=2 atmr.os=[2 2 1 1 1 1 1]
run_maos "LGS MCAO (MVM):  " -cmcao_lgs.conf fit.alg=0 recon.mvm=1
else
run_maos_gpu "LGS MCAO (MVM):  " -cmcao_lgs.conf fit.alg=0 recon.mvm=1
fi
run_maos "LGS MCAO PCCD:  " -cmcao_lgs.conf powfs.radpix=[16,0,0] powfs.pixpsa=[6,0,0]

run_maos "SLGS MCAO (inte): " -cmcao_lgs.conf powfs.fnllt=['llt_SL.conf',,] powfs.pixpsa=[16,0,0] recon.split=0

run_maos "SLGS MCAO (ahst):  " -cmcao_lgs.conf powfs.fnllt=['llt_SL.conf',,] powfs.pixpsa=[16,0,0] recon.split=1 tomo.splitlrt=2

run_maos "LGS NCPA noatm: " -cmcao_lgs.conf sim.noatm=1 ncpa.surf=["'rms=150;L0=20;D=60;SURFEVL=0;'"] sim.wspsd= powfs.noisy=[0] atm.r0z=10 powfs0_llt.fnsurf="'rms=150;mode=5;D=0.5;dx=1/64'" #should be close to 0.

run_maos "LGS NCPA TWFS:  " -cmcao_lgs.conf powfs_shtwfs.conf ncpa.surf=["'rms=150;L0=20;D=60;SURFEVL=0;SURFWFS=[1,0,0,0]'"] powfs.dtrat=[1 1 1 10] sim.noatm=1 ncpa.calib=0 powfs.step=[0 2 2 10] powfs.noisy=[0] sim.eptwfs=0.5
if [ "$D" = 30 ];then
run_maos "NFIRAOS LGS:    " -cnfiraos_lgs.conf
run_maos_gpu "NFIRAOS PYWFS:" -cnfiraos_ngs.conf
fi

echo "REF_${ARCH}=(${RMS[*]}) #$(date +%Y-%m-%d) $HOSTNAME D=$D" | tee -a $fnres >> $fnref
if [ $ans -ne 0 ];then
	echo "$ans tests failed"
fi
};exit $ans 
