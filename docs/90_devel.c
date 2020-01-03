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
   \page devel Development

   The following shows how various components are implemented.
 

   \section Fundementals

   \section Auxiliary
   
   \subsection Sodium range variation

   The sodium range variation is specified through the LLT configuration key
   llt.fnrange in llt_CL.song or llt_SL.conf. When supplied by user, use
   powfs0_llt.fnrange. The file should contain additional change (wrt fnprof) of
   sodium height in meters with dimension nx1 or nxnLGS where n can be less
   than, equal to or greater than the length of simulation and nLGS is the
   number of LGS in this powfs (6 for NFIRAOS).
   

*/
