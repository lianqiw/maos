#!/bin/sh
test=1 #1: testing. use only 50% profile and za=30 and a single seed
test_elt=0 #1: ELT. 0: NFIRAOS
dest=202507 #result folder
atm_prefix="atm_mk13n"

zas="0 30 45 60"
profs="75 50 25"

step0_fitonly="aper.fnamp= sim.closeloop=0 sim.idealfit=1 fit.alg=0 gpu.fit=0 gpu.tomo=0 powfs_none.conf" #MCAO fit only
step0_fitonly_oa="$step0_fitonly fov_oa.conf" #On axis fit only
step0_fitonly_oa_singleDM="$step0_fitonly_oa dm_single.conf" #On axis single DM fit only
step3_geom_cl="recon.split=0 powfs.noisy=[0] powfs.phystep=[-1]" #MCAO geomtric noise free WFS
step2_geom_ol="$step3_geom_cl sim.closeloop=0 " #open loop simulation (with no servo delay)
step2_geom_ol_annular="$step2_geom_ol aper.fnamp= " #with annular aperture map

step1_scao_ngs_oa="$step2_geom_ol -c scao_ngs.conf powfs.dsa=[0.5] aper.fnamp= " #NGS AO for aliasing error

step5_phy_ny="recon.split=0 " #this is the default
step4_phy_nf_nafull="$step5_phy_ny powfs.noisy=[0]" #noise free but with elongation
step6_split="$step5_phy_ny recon.split=1"
step6_mvm="$step6_split recon.mvm=1 tomo.alg=1 tomo.precond=1 fit.alg=0" #MVM

step7_llt="$step6_mvm powfs0_llt.focus=140"  #LLT implementation
step7_hyst="$step7_llt  dm.hyst=[0.05]" #DM hystersis

steps_oa="step0_fitonly_oa step0_fitonly_oa_singleDM step1_scao_ngs_oa" #on axis only. open loop
steps_ol="step0_fitonly step2_geom_ol_annular step2_geom_ol" #open loop over fov
steps_cl="step3_geom_cl step4_phy_nf_nafull step5_phy_ny " #closed loop over fov 
steps_impl="step6_split step6_mvm step7_llt step7_hyst " 
onelayer="atmr.ht=[0] atmr.wt=[1] atmr.os=[2]"

if [ "$test_elt" = 1 ]; then #test elt
	base="-c morfeo_lgs.conf cn2.pair=[]" 
	fovs="34 52"
	if true; then
		dest+="_morfeo_elt"
		atm_prefix="atm_elt"
	else
		dest+="_morfeo_mk13n"
	fi
else
	base="-c nfiraos_lgs_noimpl.conf"
	fovs="34 "
fi

if [ "$test" = 1 ] ; then #testing
	zas="30"
	profs="50"
else
	base+=" sim.seeds=[1 10 20 30] "
fi
base+=" save.extra=0 sim.end=2000"
maos="maos -d"

for prof in $profs;do
	atm="${atm_prefix}${prof}p.conf"
	for za in $zas;do
		if test -n "$steps_oa"  ;then #on axis case.
			for step in $steps_oa;do
			$maos $base ${!step} $atm sim.zadeg=$za $onelayer -o ${dest}/${prof}p_za${za}/oa${extra_fn}/${step}  || exit
			done
		fi
		
		for fov in $fovs ;do
			for step in $steps_ol $steps_cl $steps_impl ;do
			:
			$maos $base ${!step} $atm sim.zadeg=$za evl.fov=${fov} fit.fov=${fov} -o ${dest}/${prof}p_za${za}/fov${fov}${extra_fn}/${step}  || exit
			done
		done
	done
done

