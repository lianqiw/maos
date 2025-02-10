\page skycoverage Sky Coverage

\section MAOS presimulation
Add `skyc_10.conf` to `maos` simulation to output information for sky coverage preprocessing with NGS grid spaced by 10". 
```
for prof in 50 25 75 ;do
	for za in 0 30 45 50 55 60;do
		maos -c nfiraos_lgs.conf skyc_10.conf atm_mk13n${prof}p.conf sim.zadeg=$za -o ${prof}p_za${za} sim.seeds=[1] 
	done
done
```
\section Sky coverage post processing simulation
Change directory to `${prof}p_za${za}/skysim` and run the following
```
	skyc skyc.pixpsa=[16 16] skyc.psd_ws= -d skyc.multirate=1 skyc.seed=1 -o APD_JHKs_typeImr_v2 -d  
```
The default config is `maos.conf` which is saved by `maos` presimulation. Additional parameters that can be overriden are in `${src}/config/skyc/skyc.conf`. 

\section Sky coverage results

The following file names are given for `maos seed=1` and `skyc seed=1`.

- \c Res1_1: Science field averaged RMS WFE for each for each sky of the following modes:
	- Total low order 
	- Turbulence low order
	- Turbulence tip/tilt
	- Residual windshake tip/tilt
	- Estimated error from servo analysis
- \c Res1_1_oa: On axis RMS WFE in the same format as Res1_1
- \c Res1_1_aster: An array of all aterisms that have been evaluated. Each sky field can have up to skyc.maxaster number of entries. Each column is for an asterism. Each row contains the properties of this astierms:
	- RMS WFE in physical optics simulations
	- RMS WFE in servo optimization estimate
	- Sky field index
	- Number of active WFS
	- Multiple numbers per WFS
		- sampling rate (0 if no star is available)
		- signal to noise ratio at this sampling rate
		- x coordinate
		- y coordinate
		- magnitude in first wavelength
		- magnitude in second wavelength (if applicable)
		- magnitude in third wavelength (if applicable)
	
