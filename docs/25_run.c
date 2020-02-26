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
   \page page25 Running simulations
   
   \section sect-run Introduction
   
   We assume that maos is already in the PATH so we can type "maos" to launch it.
   
   Use the following command to get the help message:
   \verbatim
   maos -h  #will print the help message and exit.
   \endverbatim

   If maos is launched without any parameters, it will use default.conf which
   contains the baseline NFIRAOS configuration nfiraos.conf. Users can change
   the content of default.conf to include another file to make that the default
   option. The valid command line arguments are listed as follows

   \verbatim
   maos   # Run the default configuration of NFIRAOS: nfiaros.conf as the baseline.
   maos -c scao_ngs.conf -s 2 -n 2 -d -o scao_ngs override.conf chol.conf
   # Run a single conjugate natural guide star case, with seed 2, 2 threads
   # detach from the terminal and output results to folder scao_ngs
   # and read in overriding parameters stored in override.conf and chol.conf

   Options: 
   -h, --help        to print out this message
   -d, --detach      to detach from terminal and run in background
   -f, --force       force starting simulation without scheduler
   -n, --nthread=N   Use N threads, default is 1
   -s, --seed=N      Use seed N instead of the numbers specified in .conf files
   -o, --output=DIR  output the results to DIR.
   -c, --conf=FILE.conf
   Use FILE.conf as the baseline config instead of nfiraos.conf
   -p, --path=dir    Add dir to the internal PATH
   -P, --pause       paulse simulation in the end of every time step
   -g, --gpu=i       Use the i'th gpu. 0 for the first. -1 to disable. default: automatic
   -G, --ngpu=N'     Use a total of N gpus.
   -r, --run=host    Run the job in another host.

   \endverbatim

   Beside this, key=value pairs, as found in the configuration files, can be
   supplied in command line to override default values. Please be aware that
   such options that appear later will override those that appear earlier if the
   key is the same.
   
   \section sect-exe AO Configurations

   With a 10 meter outer diameter, 1 meter inner diameter annular telescope, the
   default subaperture size is dm.dx=0.5, implying order 20x20 AO system. the
   following shows how to run each AO mode:
   
   - LGS MCAO 
   
   \verbatim
   maos aper.d=[10 1] -o lgsmcao20 #default is dual conjugate
   maos aper.d=[10 1] dm_triple.conf -o lgsmcao20 #triple conjugate
   \endverbatim

   - NGS MCAO
   
   \verbatim
   maos -cmcao_ngs.conf aper.d=[10 1] -o ngsmcao20
   \endverbatim

   - NGS SCAO
    
   \verbatim
   maos -cscao_ngs.conf aper.d=[10 1] -o ngsscao20 
   \endverbatim

   - LGS SCAO
   
   \verbatim
   maos -cscao_lgs.conf aper.d=[10 1] -o lgsscao20 
   \endverbatim
   
   - LGS GLAO
   
   \verbatim
   maos dm_single.conf aper.d=[10 1] wfs_lgs_only.conf recon.glao=1 -o lgsglao20
   \endverbatim
   
   \section sect-details Change Parameters

   This section describes how to change simulation detailed parameters. The
   configurations can be combined to make multiple changes in a single run.

   - Change turbulence

   Override keys listed in `config/maos/atm_mk13n50p.conf`

   E.g., change r0 at zenith to 0.1m 
   \verbatim
   maos atm.r0z=0.1
   \endverbatim

   - Change zenith angle to 30 degrees

   \verbatim
   maos sim.zadeg=30
   \endverbatim

   - Use noise free wfs

   \verbatim
   maos powfs.noisy=[0]
   \endverbatim

   - Use geometric wfs (sense wavefront directly) instead of physical optics wfs

   \verbatim 
   maos powfs.phystep=[-1]
   \endverbatim

   Adjust the number of elements depending on how many powfs is in use. A powfs
   is a type of WFS. Each type can have multiple WFS.

   - Change gain of the main integrator to 0.3
   
   \verbatim
   maos sim.ephi=0.3
   \endverbatim

   \section sect-perfevl Performance Evaluation

   Maos only computes RMS averaged wavefront errors by default, which are saved
   to Res_seed.bin. When desired, PSFs computing can be enabled for some or all
   of the science evaluation directions. Modify \c evl.thetax , \c evl.thetay , and \c
   evl.wt to modify evaluation directions and relative weighting for field
   averaging. Use the following parameters to enable PSF computation.
   \verbatim
   evl.psfisim=20 #time step to start accumulating PSFs.
   evl.psfmean=1 #1: output time averaged PSF in the end. n>1: save cumulative average every n steps.
   evl.psfol=1 #1:also compute time averaged open loop PSF. 
   evl.psf=[1 0] #a flag for each direction
   evl.wvl=[] #specify a list of wavelength in meter
   evl.psfsize=[] #size of the psf to output (center crop)
   \endverbatim
   See \c sim.conf for details.


   See \ref simulation for further explanation.

   See \ref page40 for result interpretation.

*/
