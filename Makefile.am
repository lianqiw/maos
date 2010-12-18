ACLOCAL_AMFLAGS=-I m4 --install
SUBDIRS=external/cholmod/ lib/ maos/ skyc/ tools/

include $(top_srcdir)/aminclude.am

EXTRA_DIST = doxygen.cfg doxygen.xml doxygen.css docs

#Note: in EXTRA_DIST, don't specify a folder. If you do, everything in that
#folder, including compiled files (in the compiled folder), svn directory, will
#be included.

EXTRA_DIST += mex/Makefile mex/write.c mex/read.c mex/io.c mex/io.h
EXTRA_DIST += scripts/readbin.pro scripts/readdata.pro scripts/readbin.m scripts/writebin.m\
scripts/save_surf.m scripts/save_tsurf.m
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
config/maos/dm_none.conf\
config/maos/dm_single.conf\
config/maos/dm_triple.conf\
config/maos/fov_iris.conf\
config/maos/fov_oa.conf\
config/maos/fov_pmos.conf\
config/maos/keck_lgs.conf\
config/maos/llt_CL.conf\
config/maos/mcao_ngs.conf\
config/maos/nfiraos45.conf\
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
config/skyc/maos.conf\
config/skyc/skyc.conf\
config/skyc/skyc_za0.conf\
config/skyc/skyc_za15.conf\
config/skyc/skyc_za30.conf\
config/skyc/skyc_za45.conf\
config/skyc/skyc_za55.conf\
config/skyc/skyc_za60.conf\
config/skyc/PSD_NGS.bin.gz\
config/skyc/PSD_PS.bin.gz\
config/skyc/PSD_TMT_WindShake20mas.bin.gz\
config/skyc/PSD_TT.bin.gz\
config/skyc/besancon/JH_5sqdeg_lat90_lon200_besancon.bin.gz
