{
	for D in 5 10 30 ;do
		for dev in '' '-g-1' ;do
			maos_check.sh $D $dev
		done
	done
};exit 
