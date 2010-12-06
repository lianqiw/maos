EXTRA_DIST=\
./Include\
./Supernodal/t_cholmod_super_numeric.c\
./Supernodal/t_cholmod_super_solve.c\
./Core/t_cholmod_triplet.c\
./Core/t_cholmod_transpose.c\
./Core/t_cholmod_dense.c\
./Core/t_cholmod_change_factor.c\
./MatrixOps/t_cholmod_sdmult.c\
./Cholesky/t_cholmod_rowfac.c\
./Cholesky/t_cholmod_lsolve.c\
./Cholesky/t_cholmod_solve.c\
./Cholesky/t_cholmod_ltsolve.c\
./Modify/t_cholmod_updown.c\
./Modify/t_cholmod_updown_numkr.c\
./Supernodal/gpl.txt\
./Supernodal/License.txt\
./Include/README.txt\
./Include/License.txt\
./Core/lesser.txt\
./Core/License.txt\
./MatrixOps/gpl.txt\
./MatrixOps/License.txt\
./Check/lesser.txt\
./Check/License.txt\
./Cholesky/lesser.txt\
./Cholesky/License.txt\
./Modify/gpl.txt\
./Modify/License.txt

#remove _DNPARTITION is METIS is installed 
AM_CFLAGS=-DNDEBUG=1 -O3 -DUSE_MEM=0 -DHAVE_SSE2=1 -DNPARTITION 
if USE_LONG
AM_CFLAGS+=-DDLONG
endif
if COMP_CHOL
lib_LTLIBRARIES=libcholmod.la
libcholmod_la_CFLAGS=$(AM_CFLAGS) -I$(srcdir)/Include
libcholmod_la_SOURCES=$(CHOL_SRC)
libcholmod_la_LIBADD=-lm
else
EXTRA_DIST+=$(CHOL_SRC)
endif
CHOL_SRC=\
./Supernodal/cholmod_super_solve.c \
./Supernodal/cholmod_super_numeric.c \
./Supernodal/cholmod_super_symbolic.c \
./COLAMD/colamd_global.c \
./COLAMD/colamd.c \
./Core/cholmod_factor.c \
./Core/cholmod_complex.c \
./Core/cholmod_band.c \
./Core/cholmod_common.c \
./Core/cholmod_aat.c \
./Core/cholmod_copy.c \
./Core/cholmod_dense.c \
./Core/cholmod_change_factor.c \
./Core/cholmod_error.c \
./Core/cholmod_transpose.c \
./Core/cholmod_memory.c \
./Core/cholmod_add.c \
./Core/cholmod_sparse.c \
./Core/cholmod_triplet.c \
./MatrixOps/cholmod_scale.c \
./MatrixOps/cholmod_submatrix.c \
./MatrixOps/cholmod_symmetry.c \
./MatrixOps/cholmod_sdmult.c \
./MatrixOps/cholmod_horzcat.c \
./MatrixOps/cholmod_norm.c \
./MatrixOps/cholmod_ssmult.c \
./MatrixOps/cholmod_vertcat.c \
./MatrixOps/cholmod_drop.c \
./Check/cholmod_read.c \
./Check/cholmod_check.c \
./Check/cholmod_write.c \
./CCOLAMD/ccolamd.c \
./CCOLAMD/ccolamd_global.c \
./CAMD/Source/amd_info.c \
./CAMD/Source/amd_order.c \
./CAMD/Source/amd_postorder.c \
./CAMD/Source/amd_preprocess.c \
./CAMD/Source/amd_dump.c \
./CAMD/Source/amd_post_tree.c \
./CAMD/Source/amd_1.c \
./CAMD/Source/amd_aat.c \
./CAMD/Source/amd_control.c \
./CAMD/Source/amd_defaults.c \
./CAMD/Source/amd_valid.c \
./CAMD/Source/amd_2.c \
./CAMD/Source/amd_global.c \
./CAMD/camd_1.c \
./CAMD/camd_info.c \
./CAMD/camd_preprocess.c \
./CAMD/camd_defaults.c \
./CAMD/camd_postorder.c \
./CAMD/camd_dump.c \
./CAMD/camd_control.c \
./CAMD/camd_order.c \
./CAMD/camd_2.c \
./CAMD/camd_global.c \
./CAMD/camd_valid.c \
./CAMD/camd_aat.c \
./Cholesky/cholmod_postorder.c \
./Cholesky/cholmod_solve.c \
./Cholesky/cholmod_rowfac.c \
./Cholesky/cholmod_spsolve.c \
./Cholesky/cholmod_analyze.c \
./Cholesky/cholmod_resymbol.c \
./Cholesky/cholmod_etree.c \
./Cholesky/cholmod_factorize.c \
./Cholesky/cholmod_rowcolcounts.c \
./Cholesky/cholmod_colamd.c \
./Cholesky/cholmod_amd.c \
./Cholesky/cholmod_rcond.c \
./Modify/cholmod_updown.c \
./Modify/cholmod_rowdel.c \
./Modify/cholmod_rowadd.c 

#CC=$(CCsave)