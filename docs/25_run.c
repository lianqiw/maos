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
   
   \section sect-exe Examples
   With a 10 meter outer diameter, 1 meter inner diameter annular telescope, the defaule subaperture size is dm.dx=0.5, implying order 20x20 AO system. the following shows how to run each AO mode:
   
   - LGS MCAO 

   \verbatim
   maos aper.d=[10 1] -o lgsmcao20
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
   
   \section sect-details Details

   This section describes how to change simulation detailed parameters. The
   configurations can be combined to make multiple changes in a single run.

   - Change turbulence

   Override keys listed in `config/maos/atm_mk13n50p.conf`

   E.g., change r0 at zenith to 0.1m 
   
   `maos atm.r0z=0.1 `

   - Change zenith angle to 30 degrees

   `maos sim.zadeg=30`

   - Use noise free wfs

   `maos powfs.noisy=[0]`

   - Use geometric wfs (sense wavefront directly) instead of physical optics wfs

   `maos powfs.phystep=[-1]` 

   Adjust the number of elements depending on how many powfs is in use. A powfs
   is a type of WFS. Each type can have multiple WFS.

   - Change gain of the main integrator to 0.3

   `maos sim.epdm=0.3`

   See \ref page30 for further explanation.

   See \ref page40 for result interpretation.

*/
