#Launch NFIRAOS LGS MCAO simulations
test=1 #1: testing 
save_psf="2" #psf saving options
extra=""
base="-c nfiraos_lgs.conf atm.size=[256 256] " #start with NFIRAOS configurations. simulated implementation errors are included.
sys=nfiraos
if true; then
	base+=" ncpa.surf=[] sim.wspsd=" #disable surfaces and vibration
	sys="mcao"
fi
case $save_psf in
  1)
	base+=" evl.psfmean=600 sim.end=6020 " #cumulatively average and save PSF every 600 steps (1 second for 100 seconds).
	#base+=" evl.psfsize=[1000]" #crop PSF to 1000x1000. the sampling is wvl/(2*D)
    base+=" evl.psfol=1 evl_oa.conf " #also save open loop PSFs
	;;
  2) #for Ahmed (JPL/Caltech)
	base+=" evl.psfmean=600 sim.end=6020 evl.wvl=[0.75e-6] evl.psfgridsize=[9900]" #cumulatively average and save PSF every 600 steps (1 second for 100 seconds).
    base+=" evl.psfol=1 evl_oa.conf " 
	extra="_750"
	;;
esac
profs="50 25 75" #turbulence profiles
zas="0 30 45 50 55 60" #zenith angles

if [ "$test" = 1 ]; then #testing
	profs="50"
	zas="30"
fi
for prof in $profs ;do
	for za in $zas;do
		maos -d $base atm_mk13n${prof}p.conf sim.zadeg=$za -o ${sys}/${prof}p_za${za}${extra}
	done
done
