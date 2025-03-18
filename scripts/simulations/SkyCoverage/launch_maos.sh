#!/bin/sh
#First launch MAOS simulations with as much implementations as possible
#- LGS spot size
#- 200 nm additional off axis
# windshake/vibration can be done in post processing 

base="-c nfiraos_lgs.conf sim.wspsd= " #do not include windshake and vibration psd to have more flexibility in sky coverage

maos='maos -d'
for prof in 50 25 75 ;do
	for za in 0 30 45 50 55 60;do
		$maos $base skyc_10.conf atm_mk13n${prof}p.conf sim.zadeg=$za -o ${prof}p_za${za} sim.seeds=[1] 
	done
done
