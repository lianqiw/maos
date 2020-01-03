/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \page page41 Geometry Data

   When \c save.setup=1, MAOS will save the geometry data created before
   entering simulation to files in folder \c setup. When \c save.recon=1, the
   tomography and DM fitting matrices are also saved, which takes space and are
   not enabled by \c save.setup=1.

   The following explains each file. Not all of them may exist, depending on set
   of options. 

   The suffix of the file names are removed.

<table border="0">
             <tr><td>\c actslave  </td><td>The slaving operator
   </td></tr><tr><td>\c ahst_GM   </td><td>ad hoc split tomography, model to gradient operator
   </td></tr><tr><td>\c ahst_MCC  </td><td>ahst, model cross coupling matrix.
   </td></tr><tr><td>\c ahst_Mdm  </td><td>ahst, the NGS mode defined on DM grid
   </td></tr><tr><td>\c ahst_Pngs </td><td>ahst, ngs mode removal operator from DM commands
   </td></tr><tr><td>\c ahst_Ptt  </td><td>ahst, tip/tilt removal operator from DM commands
   </td></tr><tr><td>\c ahst_Rngs </td><td>ahst, ngs mode reconstructor
   </td></tr><tr><td>\c ahst_Wa   </td><td>ahst, the weighting using science field.
   </td></tr><tr><td>\c aloc      </td><td>DM actuator grid.
   </td></tr><tr><td>\c ploc      </td><td>The coarse sampled grid on aperture for reconstruction
   </td></tr><tr><td>\c xloc      </td><td>The tomography grid.
   </td></tr><tr><td>\c TT        </td><td>The global tip/tilt modes for wfs
   </td></tr><tr><td>\c PTT       </td><td>The global tip/tilt removal operator for wfs
   </td></tr><tr><td>\c saneai    </td><td>The inverse of nea^2 used in tomography (rad^2)
   </td></tr><tr><td>\c W0        </td><td>The W0 weighting defined on ploc. 
   </td></tr><tr><td>\c W1        </td><td>The W1 weighting defined on ploc.
   </td></tr><tr><td>\c aper_locs </td><td>The fine sampled grid on telescope aperture.
   </td></tr><tr><td>\c aper_amp  </td><td>Telescope aperture amplitude defined on aper_locs
   </td></tr><tr><td>\c aper_mcc  </td><td>modal cross-coupling matrix for modes defined on aper_locs
   </td></tr><tr><td>\c FLM       </td><td>Fit left hand side operator, sparse matrix
   </td></tr><tr><td>\c FLU       </td><td>Fit left hand side operator, Low rank U matrix
   </td></tr><tr><td>\c FLV       </td><td>Fit left hand side operator, Low rank V matrix
   </td></tr><tr><td>\c FRM       </td><td>Fit right hand side operator, sparse matrix
   </td></tr><tr><td>\c FRU       </td><td>Fit right hand side operator, Low rank U matrix
   </td></tr><tr><td>\c FRV       </td><td>Fit right hand side operator, Low rank V matrix
   </td></tr><tr><td>\c GA        </td><td>Gradient operator from aloc.
   </td></tr><tr><td>\c GP        </td><td>Gradient operator from ploc.
   </td></tr><tr><td>\c GX        </td><td>Gradient operator from xloc.
   </td></tr><tr><td>\c HA        </td><td>Ray tracing from aloc to ploc along fit directions
   </td></tr><tr><td>\c HXF       </td><td>Ray tracing from xloc to ploc along fit directions
   </td></tr><tr><td>\c HXW       </td><td>Ray tracing from xloc to ploc along wfs directions
   </td></tr><tr><td>\c L2        </td><td>Laplacian regularization on xloc
   </td></tr><tr><td>\c NW        </td><td>Low rank terms in fitting
   </td></tr><tr><td>\c NW2       </td><td>Adjusted low rank terms by slaving.
   </td></tr><tr><td>\c powfs0_area         </td><td>The subaperture area
   </td></tr><tr><td>\c powfs0_dtf0_nominal </td><td>The nominal of DTF
   </td></tr><tr><td>\c powfs0_dtf0_si      </td><td>The si of DTF
   </td></tr><tr><td>\c powfs0_etfprep0_2d  </td><td>The elongation transfer function
   </td></tr><tr><td>\c powfs0_GP           </td><td>The gradient operator from ploc.
   </td></tr><tr><td>\c powfs0_GS0          </td><td>The gradient operator from aper_locs.
   </td></tr><tr><td>\c powfs0_i0           </td><td>The subaperture time averaged image for matched filter
   </td></tr><tr><td>\c powfs0_gx           </td><td>The pixel by pixel x gradient of i0
   </td></tr><tr><td>\c powfs0_gy           </td><td>The pixel by pixel y gradient of i0
   </td></tr><tr><td>\c powfs0_imcc         </td><td>The inverse of model cross-coupling matrix of piston/tip/tilt modes
   </td></tr><tr><td>\c powfs0_loc          </td><td>The grid for all subapertures (grouped by subapertures)
   </td></tr><tr><td>\c powfs0_amp          </td><td>The amplitude defined on powfs0_loc
   </td></tr><tr><td>\c powfs0_llt_loc      </td><td>The aperture grid of the uplink laser launch telescope (LLT)
   </td></tr><tr><td>\c powfs0_llt_amp      </td><td>The aperture amptidue of LLT defined on powfs0_llt_loc
   </td></tr><tr><td>\c powfs0_llt_imcc     </td><td>The inverse of model cross-coupling matrix of p/t/t modes for LLT
   </td></tr><tr><td>\c powfs0_srot         </td><td>The orientation of each subaperture wrt LLT
   </td></tr><tr><td>\c powfs0_srsa         </td><td>The distance of each subaperture from the LLT
   </td></tr><tr><td>\c powfs0_mtche        </td><td>The matched filter gradient estimator
   </td></tr><tr><td>\c powfs0_pts          </td><td>The lower left grid point of each subaperture.
   </td></tr><tr><td>\c powfs0_saloc        </td><td>The lower left corner of each subaperture
   </td></tr><tr><td>\c powfs0_sanea        </td><td>The subaperture noise equivalent angle(nea) in rad^2
   </td></tr><tr><td>\c powfs0_saneaxy      </td><td>The nea along x/y directions in rad^2
   </td></tr><tr><td>\c powfs0_saneaixy     </td><td>The inverse of nea along x/y directions in rad^-2
   </td></tr><tr><td>\c powfs0_saneara      </td><td>The subaperture nea along r/a directions in rad^-2
   </td></tr><tr><td>\c powfs0_sepsf        </td><td>The subaperture short exposure psf (tip/tilt removed)
   </td></tr><tr><td>\c powfs0_sodium       </td><td>The sodium layer profile
   </td></tr><tr><td>\c powfs0_sprint       </td><td>which subapertures to print
   </td></tr><tr><td>\c powfs0_GS0          </td><td>The averaging gradient operator from pts.
   </td></tr><tr><td>\c powfs1_ZS0          </td><td>The ztilt best fit gradient operator from pts.
   </td></tr>
</table>
*/
