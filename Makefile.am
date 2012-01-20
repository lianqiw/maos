ACLOCAL_AMFLAGS=-I m4 --install
SUBDIRS=
if COMP_CHOL
SUBDIRS+=external/cholmod
endif
SUBDIRS+=lib
if USE_CUDA
SUBDIRS+=cuda
endif
SUBDIRS+=maos skyc tools
if USE_MEX
SUBDIRS+=mex
endif
include $(top_srcdir)/aminclude.am

dist-hook:
	rm -rf $(distdir)/external $(distdir)/test

all-local:$(srcdir)/lib/tmp.c
$(srcdir)/tmp.c:$(srcdir)/lib/mat.c $(srcdir)/lib/cell.c $(srcdir)/lib/matbin.c \
	$(srcdir)/lib/spbin.c $(srcdir)/lib/sp.c
	@echo "#define USE_COMPLEX" > $(srcdir)/lib/tmp.c
	@cat $(srcdir)/lib/mat.c $(srcdir)/lib/cell.c $(srcdir)/lib/matbin.c $(srcdir)/lib/spbin.c $(srcdir)/lib/sp.c >> $(srcdir)/lib/tmp.c

EXTRA_DIST = doxygen.cfg doxygen.xml doxygen.css docs 

#Note: in EXTRA_DIST, don't specify a folder. If you do, everything in that
#folder, including compiled files (in the compiled folder), svn directory, will
#be included.

EXTRA_DIST += scripts/readbin.pro scripts/readdata.pro scripts/readbin.m scripts/writebin.m\
scripts/save_surf.m scripts/save_tsurf.m
nobase_dist_pkgdata_DATA=\
config/maos/bin/TMT_CurDesign.bin\
config/maos/bin/KECK.bin\
config/maos/bin/NapMean.bin\
config/maos/bin/PSD_TMT_WindShake20mas_rad2.bin\
config/maos/atm_mk13n25p.conf\
config/maos/atm_mk13n50p.conf\
config/maos/atm_mk13n75p.conf\
config/maos/atm_srd.conf\
config/maos/atm_single.conf\
config/maos/dbg.conf\
config/maos/dm_dual.conf\
config/maos/dm_none.conf\
config/maos/dm_single.conf\
config/maos/dm_triple.conf\
config/maos/fov_iris.conf\
config/maos/fov_wirc.conf\
config/maos/fov_oa.conf\
config/maos/fov_pmos.conf\
config/maos/llt_CL.conf\
config/maos/llt_SL.conf\
config/maos/nfiraos.conf\
config/maos/recon.conf\
config/maos/scao_base.conf\
config/maos/scao_lgs.conf\
config/maos/scao_ngs.conf\
config/maos/sim.conf\
config/maos/skyc_10.conf\
config/maos/skyc_20.conf\
config/maos/wfs_lgs_tt.conf\
config/maos/wfs_lgs_ttf_tt.conf\
config/maos/wfs_ngs.conf\
config/maos/default.conf\
config/maos/change.log\
config/skyc/maos_test.conf\
config/skyc/skyc.conf\
config/skyc/skyc_za0.conf\
config/skyc/skyc_za15.conf\
config/skyc/skyc_za30.conf\
config/skyc/skyc_za45.conf\
config/skyc/skyc_za55.conf\
config/skyc/skyc_za60.conf\
config/skyc/PSD/PSD_NGS.bin\
config/skyc/PSD/PSD_PS.bin\
config/skyc/PSD/PSD_TMT_WindShake20mas_m2.bin\
config/skyc/PSD/PSD_TT.bin\
config/skyc/besancon/JH_5sqdeg_lat90_lon200_besancon.bin
