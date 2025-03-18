#The second step is to launch skyc simulations

CWD=$(pwd)
servoname=(LQG null typeI)
mrname=("" "mr")
maos_seeds=1
skyc_seed=1
skyc="skyc -d"

for servo in -1 1;do #1: integrator, -1: LQG
for multirate in 0 1;do
for prof in 25 50 75 ;do
for za in 0 30 45 50 55 60 ;do
for catscale in 1; do #scale the star density. 1 is Galactic Pole. >1 for Other Galactic Latitudes.
	
	base=" skyc.psd_ws=PSD/PSD_TMT_ws50p13mas_vib15mas_m2.bin " #maos simulation no longer included ws_vib_psd
	suffix=_v1
	
	if [ "$catscale" != 1 ];then
		base+=" skyc.catscl=$catscale"
		suffix+="_catscl$catscale"
	fi

	base+=" skyc.multirate=$multirate"
	
	cd $CWD/${prof}p_za${za}/skysim || exit	
	for servotype in $servo ;do
		fdout=${servoname[1+$servotype]}${mrname[$mr]}${suffix}
		$skyc $base maos.seeds=[$maos_seeds] skyc.servo=${servotype} skyc.seed=$skyc_seed -o $fdout
	done

done
done
done
done
done
