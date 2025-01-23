#script to run NFIRAOS NGSAO under varying conditions

mags=(8 9 10 11 12 13 14 15 16 ) #PWFS NGS magnitude
freq=(2000 2000 2000 1500 1000 600 300 200 50 ) #PWFS sampling freq
imag0=0
imag1=${#mags[*]}

profs="25 50 75" #turbulence profile percentile
zas="0 30 45 60" #zenith angle
maos="maos -d"
fd="ngsao"
if true; then #testing
	profs="50"
	zas="0"
	maos="maos"
	imag0=0
	imag1=1
	fnextra="_v3 sim.seeds=[1] evl.psfmean=0"
fi

for prof in $profs;do
  for za in $zas;do
	base="-c nfiraos_ngs.conf " #AO configuration: NFIRAOS NGS AO
	base+=" atm_mk13n${prof}p.conf sim.zadeg=$za " #seeing and telescope condition

	base+=" save.extra=0 evl.psfmean=1 evl.psfisim=10000" # save extra telemetry, psf time average starting at step evl.psfisim
	base+=" sim.end=20000 sim.seeds=[1,10,20,30]" #simulation length and seeds,

	for ((imag=$imag0; imag<$imag1; imag++));do
	    mag=${mags[$imag]}
	    hz=${freq[$imag]}
	    if [ $hz -gt 800 ];then
			hz=800
	    fi
	    echo mag=$mag hz=$hz
		$maos $base sim.dt=1/${hz}  powfs.mag=[$mag] -o ${fd}/${prof}p_za${za}/mag${mag}${fnextra}
	done
  done
done
