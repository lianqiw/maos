ACLOCAL_AMFLAGS=-I m4 --install
#SUBDIRS=
#if COMP_ARPACK
#SUBDIRS+=external/arpack/
#endif
#if COMP_CHOL
#SUBDIRS+=external/cholmod/
#endif
SUBDIRS=external/arpack/ external/cholmod/ lib/ maos/ skyc/ 

include $(top_srcdir)/aminclude.am

EXTRA_DIST = doxygen.cfg doxygen_layout.xml doxygen.css index.c

#Note: in EXTRA_DIST, don't specify a folder. If you do, everything in that
#folder, including compiled files (in the compiled folder), svn directory, will
#be included.

EXTRA_DIST += mex/Makefile mex/write.c mex/read.c mex/io.c mex/io.h

nobase_dist_pkgdata_DATA=\
config/maos/bin/TMT_CurDesign.bin.gz\
config/maos/bin/KECK.bin.gz\
config/maos/bin/NapMean.bin.gz\
config/maos/aper_tmt.conf\
config/maos/atm_mk13n25p.conf\
config/maos/atm_mk13n50p.conf\
config/maos/atm_scao.conf\
config/maos/atm_srd.conf\
config/maos/dbg.conf\
config/maos/dm_dual.conf\
config/maos/dm_single.conf\
config/maos/fov_oa.conf\
config/maos/fov_iris.conf\
config/maos/llt_CL.conf\
config/maos/mcao_ngs.conf\
config/maos/nfiraos.conf\
config/maos/nfiraos45.conf\
config/maos/recon.conf\
config/maos/scao_base.conf\
config/maos/scao_lgs.conf\
config/maos/scao_ngs.conf\
config/maos/sim.conf\
config/maos/wfs_multi_lgs.conf\
config/maos/wfs_single_lgs.conf\
config/maos/wfs_single_ngs.conf\
config/maos/skyc_10.conf\
config/maos/skyc_20.conf\
config/maos/keck_lgs.conf\
config/skyc/maos.conf\
config/skyc/skyc.conf\
config/skyc/skyc_za0.conf\
config/skyc/skyc_za15.conf\
config/skyc/skyc_za30.conf\
config/skyc/skyc_za45.conf\
config/skyc/skyc_za60.conf\
config/skyc/PSD_NGS.bin.gz\
config/skyc/PSD_PS.bin.gz\
config/skyc/PSD_TMT_WindShake20mas.bin.gz\
config/skyc/PSD_TT.bin.gz\
config/skyc/besancon/JH_5sqdeg_lat90_lon200_besancon.bin.gz
