#2023-07-14: PYWFS as TTF for LGS AO
test=1 #testing
mags=(12 13 14 15 16 17 18 19 20 21 22 ) #PWFS NGS magnitude
dtratss=(1 1 1 1 1 2 4 8 16 64 80) #sampling duration of TTF PWFS over LGS WFS. determined by simulation.
imag0=0 #starting magnitude index
imag1=${#mags[*]} #ending magnitude index
lmag=7.65 #LGS magnitude
fov=1 #science FoV (diameter)
profs="25 50 75" #turbulence profile percentile
zas="0 30 45 60" #zenith angle

maos="maos -d"
fd="mcao_pyttf" #output folder name
if [ "$test" = 1 ]; then #testing with reduced parameter set
	profs="50"
	zas="0"
	maos="maos"
	imag0=0
	imag1=1
	fnextra="_v5 sim.seeds=[1] evl.psfmean=0"
fi

for prof in $profs;do
  for za in $zas;do
	base="-c nfiraos_lgs_pyttf.conf evl.fov=$fov fit.fov=$fov " #AO configuration
	base+=" atm_mk13n${prof}p.conf sim.zadeg=$za " #seeing and telescope condition
	base+=" save.extra=0 evl.psfmean=1 evl.psf=[1 0 0 0 0 0 0 0 0] evl.psfisim=2000" #save extra telemetry, on axis PSF time average from evl.psfisim step for only the on axis field.
	base+=" sim.end=5000 sim.seeds=[1,10,20,30] " #simulation length and seeds,

	for ((imag=$imag0; imag<$imag1; imag++));do
	    mag=${mags[$imag]}
		dtrat2=${dtratss[$imag]} #may loop over dtrats to determine best dtrat
		echo mag=$mag dtrat=$dtrat2
		$maos $base  powfs.dtrat=[1 $dtrat2] powfs.mag=[$lmag $mag]  -o ${fd}/${prof}p_za${za}/mag${mag}${fnextra}
	done
  done
done
