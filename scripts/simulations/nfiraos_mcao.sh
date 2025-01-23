#Launch NFIRAOS LGS MCAO simulations
save_psf="yes" #saving psf. set to no to disable

base="-c nfiraos_lgs.conf atm.size=[128 128] " #start with NFIRAOS configurations. simulated implementation errors are included.
if [ ${save_psf} = "yes" ];then
	base+=" evl.psfmean=600 sim.end=60020 " #cumulatively average and save PSF every 600 steps (1 second for 100 seconds).
	base+=" evl.psfsize=[1000]" #crop PSF to 1000x1000. the sampling is wvl/(2*D)
    base+=" evl.psfol=1" #also save open loop PSFs
fi
profs="50 25 75" #turbulence profiles
zas="0 30 45 50 55 60" #zenith angles

if true; then #testing
	profs="50"
	zas="0"
fi
for prof in $profs ;do
	for za in $zas;do
		maos -d $base atm_mk13n${prof}p.conf sim.zadeg=$za -o ${prof}p_za${za} 
	done
done
