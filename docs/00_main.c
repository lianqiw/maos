/** 
    \mainpage MAOS Documentation

    Multi-Thread Adaptive Optics Simulator (MAOS) is a end-to-end adaptive
    optics system simulator. It has the capability to simulation many different
    type of adaptive optics systems, including conventional single conjugate AO
    (SCAO), multi-conjugate AO (MCAO), laser tomography AO (LTAO), multi-object
    AO (MOAO), and ground layer AO (GLAO).

    Some benchmarking results using the TMT NFIRAOS system:
    Dual Intel Xeon W5590 at 3.33 Ghz: 1.3 s per time step with 8 threads.
    Intel Core i7-2600 at 3.4 Ghz: 1.9s per time step with 4 threads.
    Intel Core i7-2600 at 3.4 Ghz: 1.67s per time step with 8 hyper threads (4 physical cores).
    <p>

    - \ref page1 

    - \ref page2 

    - \ref page3 

    - \ref page4 
    
    - \ref page41 
    
    - \ref page42 

    - \ref page5 

    - \ref page6 

    - \ref page7 

    \section sect-structure Code Structure    

    - main() is the entrance into the program. It calls
        - setup_parms(): to setup the parameters of types PARMS_T and check for errors.
        - maos(): calls the following to do simulation.
            - setup_aper():  to setup the aperture (of type APER_T)
            - setup_powfs(): to setup the wfs type information (of type POWFS_T)
            - setup_recon(): to setup the wavefront reconstructor and dm fit structs (of type RECON_T)
            - sim(): to start the simulation. It then calls
                - save_skyc(): (optional) to save information for sky coverage postproc
                - FOR EACH SEED
                    - init_simu(): to initialize the run time structs
                    - genscreen(): (in closed loop) to generate atmosphere turbulence
                    - FOR EACH TIME STEP
                         - sim_update_etf(): (optional) to update the sodium profile
                         - genscreen(): (in open loop) to generate atmosphere turbulence
                         - perfevl(): (in closed loop) to evaluate the performance ins cience field
                         - wfsgrad(): to compute WFS gradients
                         - tomofit(): to do tomography and DM fit
                         - filter(): to do temporal filtering of DM commands.
                         - moao_recon(): (optional) to compute MOAO DM commands.
                         - perfevl(): (in open loop) to performance evaluation
                         - save_simu(): to save simulation telemetry data.
                         - print_progress(): to display progress data.
                    - free_simu(): to delete run time structs and close files


    \author Lianqi Wang <lianqiw@gmail.com> at TMT Corporation www.tmt.org 
    
    The source code can be obtained in github.com/lianqiw/maos 

 */
