ACLOCAL_AMFLAGS=-I m4 --install
SUBDIRS=
SUBDIRS+=sys math lib
if USE_CUDA
SUBDIRS+=cuda
endif
SUBDIRS+=maos skyc
if USE_MEX
SUBDIRS+=mex
endif
SUBDIRS+=tools test 
include $(top_srcdir)/aminclude.am

#dist-hook:
#	rm -rf $(distdir)/configure.ac 

all-local:$(srcdir)/math/tmp.c
$(srcdir)/math/tmp.c:$(srcdir)/math/mat.c $(srcdir)/math/cell.c $(srcdir)/math/matbin.c \
	$(srcdir)/math/spbin.c $(srcdir)/math/sp.c
	@echo "#define USE_COMPLEX" > $(srcdir)/math/tmp.c
	@cat $(srcdir)/math/mat.c $(srcdir)/math/cell.c $(srcdir)/math/matbin.c $(srcdir)/math/spbin.c $(srcdir)/math/sp.c >> $(srcdir)/math/tmp.c

EXTRA_DIST = doxygen.cfg docs 

#Note: in EXTRA_DIST, don't specify a folder. If you do, everything in that
#folder, including compiled files (in the compiled folder), svn directory, will
#be included.

EXTRA_DIST += scripts/
nobase_dist_pkgdata_DATA=\
config/maos/bin/TMT_30cm3.bin\
config/maos/bin/TMT_23cm6.bin\
config/maos/bin/KECK.bin\
config/maos/bin/NapMean.bin\
config/maos/bin/PSD_TMT_WindShake20mas_rad2.bin\
config/maos/atm_mk13n25p.conf\
config/maos/atm_mk13n50p.conf\
config/maos/atm_mk13n75p.conf\
config/maos/atm_single.conf\
config/maos/dbg.conf\
config/maos/default.conf\
config/maos/dm_none.conf\
config/maos/dm_dual.conf\
config/maos/dm_single.conf\
config/maos/dm_triple.conf\
config/maos/fov_iris.conf\
config/maos/fov_irms.conf\
config/maos/fov_oa.conf\
config/maos/fov_pmos.conf\
config/maos/fov_wirc.conf\
config/maos/llt_CL.conf\
config/maos/llt_SL.conf\
config/maos/mcao_ngs.conf\
config/maos/nfiraos.conf\
config/maos/recon.conf\
config/maos/scao_base.conf\
config/maos/scao_lgs.conf\
config/maos/scao_ngs.conf\
config/maos/sim.conf\
config/maos/skyc_10.conf\
config/maos/skyc_20.conf\
config/maos/powfs_common.conf\
config/maos/powfs_lgs.conf\
config/maos/powfs_ngs.conf\
config/maos/powfs_none.conf\
config/maos/powfs_tt.conf\
config/maos/powfs_ttf.conf\
config/maos/powfs_twfs.conf\
config/maos/wfs_lgs_hitt.conf\
config/maos/wfs_lgs_ngs.conf\
config/maos/wfs_lgs_only.conf\
config/maos/wfs_lgs_ttf_tt.conf\
config/maos/wfs_lgs_ttf_tt_tt.conf\
config/maos/wfs_lgs_ttf_tt_twfs.conf\
config/maos/wfs_lgs_tt_tt_tt.conf\
config/maos/wfs_ngs.conf\
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
config/skyc/besancon/JH_5sqdeg_lat90_lon0_besancon.bin\
config/skyc/besancon/JHK_5sqdeg_lat90_lon0_besancon.bin
