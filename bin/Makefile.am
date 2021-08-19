SUBDIRS=../maos ../skyc ../tools 

TESTS=maos.test
maos.test:
	@echo ./maos aper.d=5 sim.end=10 > maos.test && chmod u+x maos.test

CLEANFILES=maos.test
	
MAKEOVERRIDES=--no-print-directory
MAKEFLAGS=-j8 #enable parallel build
