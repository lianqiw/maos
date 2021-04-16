\page page40_results
\section geometry Geometry Data

When \c save.setup=1, MAOS will save the geometry data created before
entering simulation to files in folder \c setup. When \c save.recon=1, the
tomography and DM fitting matrices are also saved, which takes space and are
not enabled by \c save.setup=1.

The following explains each file. Not all of them may exist, depending on set
of options. The suffix of the file names are removed.

File name | Description
---------|--------------
\c actslave  |The slaving operator
\c ahst_GM   |ad hoc split tomography, model to gradient operator
\c ahst_MCC  |ahst, model cross coupling matrix.
\c ahst_Mdm  |ahst, the NGS mode defined on DM grid
\c ahst_Pngs |ahst, ngs mode removal operator from DM commands
\c ahst_Ptt  |ahst, tip/tilt removal operator from DM commands
\c ahst_Rngs |ahst, ngs mode reconstructor
\c ahst_Wa   |ahst, the weighting using science field.
\c aloc      |DM actuator grid.
\c ploc      |The coarse sampled grid on aperture for reconstruction
\c xloc      |The tomography grid.
\c TT        |The global tip/tilt modes for wfs
\c PTT       |The global tip/tilt removal operator for wfs
\c saneai    |The inverse of nea^2 used in tomography (rad^2)
\c W0        |The W0 weighting defined on ploc. 
\c W1        |The W1 weighting defined on ploc.
\c aper_locs |The fine sampled grid on telescope aperture.
\c aper_amp  |Telescope aperture amplitude defined on aper_locs
\c aper_mcc  |modal cross-coupling matrix for modes defined on aper_locs
\c FLM       |Fit left hand side operator, sparse matrix
\c FLU       |Fit left hand side operator, Low rank U matrix
\c FLV       |Fit left hand side operator, Low rank V matrix
\c FRM       |Fit right hand side operator, sparse matrix
\c FRU       |Fit right hand side operator, Low rank U matrix
\c FRV       |Fit right hand side operator, Low rank V matrix
\c GA        |Gradient operator from aloc.
\c GP        |Gradient operator from ploc.
\c GX        |Gradient operator from xloc.
\c HA        |Ray tracing from aloc to ploc along fit directions
\c HXF       |Ray tracing from xloc to ploc along fit directions
\c HXW       |Ray tracing from xloc to ploc along wfs directions
\c L2        |Laplacian regularization on xloc
\c NW        |Low rank terms in fitting
\c NW2       |Adjusted low rank terms by slaving.
\c powfs0_area         |The subaperture area
\c powfs0_dtf0_nominal |The nominal of DTF
\c powfs0_dtf0_si      |The si of DTF
\c powfs0_etfprep0_2d  |The elongation transfer function
\c powfs0_GP           |The gradient operator from ploc.
\c powfs0_GS0          |The gradient operator from aper_locs.
\c powfs0_i0           |The subaperture time averaged image for matched filter
\c powfs0_gx           |The pixel by pixel x gradient of i0
\c powfs0_gy           |The pixel by pixel y gradient of i0
\c powfs0_imcc         |The inverse of model cross-coupling matrix of piston/tip/tilt modes
\c powfs0_loc          |The grid for all subapertures (grouped by subapertures)
\c powfs0_amp          |The amplitude defined on powfs0_loc
\c powfs0_llt_loc      |The aperture grid of the uplink laser launch telescope (LLT)
\c powfs0_llt_amp      |The aperture amptidue of LLT defined on powfs0_llt_loc
\c powfs0_llt_imcc     |The inverse of model cross-coupling matrix of p/t/t modes for LLT
\c powfs0_srot         |The orientation of each subaperture wrt LLT
\c powfs0_srsa         |The distance of each subaperture from the LLT
\c powfs0_mtche        |The matched filter gradient estimator
\c powfs0_pts          |The lower left grid point of each subaperture.
\c powfs0_saloc        |The lower left corner of each subaperture
\c powfs0_sanea        |The subaperture noise equivalent angle(nea) in rad^2
\c powfs0_saneaxy      |The nea along x/y directions in rad^2
\c powfs0_saneaixy     |The inverse of nea along x/y directions in rad^-2
\c powfs0_saneara      |The subaperture nea along r/a directions in rad^-2
\c powfs0_sepsf        |The subaperture short exposure psf (tip/tilt removed)
\c powfs0_sodium       |The sodium layer profile
\c powfs0_sprint       |which subapertures to print
\c powfs0_GS0          |The averaging gradient operator from pts.
\c powfs1_ZS0          |The ztilt best fit gradient operator from pts.

</table>

