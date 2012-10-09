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
