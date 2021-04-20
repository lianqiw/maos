\page page93_maos Architecture

\section sect-structure Simulation Flow

- main() is the entrance into the program. It calls
    - setup_parms(): to setup the parameters of types parms_t and check for errors.
    - setup_parms_gpu(): to setup the GPU usage for various tasks
    - maos_setup(): to setup the AO system geometry and data structs. It calls
        - setup_aper():  to setup the aperture grid and amplitude map (of type aper_t)
        - setup_powfs_init(): to setup the WFS subapertures geometry
        - setup_recon_prep(): to setup reconstruction grid (ploc, aloc, floc, GP, GX, etc.)
        - setup_surf(): to setup NCPA
        - setup_powfs_phy(): to setup physical optics (dtf, etf, i0, mtch, cog)
        - setup_powfs_neasim(): to setup grad noise model if not using photon noise.
        - setup_powfs_calib(): NCPA calibration
        - setup_recon_prep2(): to setup GA, GF, GR, ngsmod geometry, t/t removal
        - setup_recon(): to setup sanea, t/t removal, read MVM, ngsmod recon
        - setup_recon_fit(): to setup DM fitting algorithm
        - setup_powfs_fit(): to setup fitting to WFS subaperture algorithm
        - setup_recon_moao(): to setup MOAO
        - setup_recon_mvm(): to calculate MVM control matrix 
        - gpu_setup_recon_mvm(): #to setup recon in all GPUs and calculate MVM control matrix 
        - gpu_setup_recon(): to setup recon in the preferred GPU
        - gpu_perfevl_init(): to initialize GPU performance evaluation parameters
        - gpu_wfssurf2gpu(): to copy NCPA OPDs to gpu
        - setup_recon_post(): to compute PSD, etc.
        - plot_setup(): to plot the geometry (ploc, floc, aloc, amplitude map, etc.).
    - maos_sim(): to start the simulation. It then calls
        - maos_iseed(): to initialize runtime data for each seed including turbulence
        - In event driven simulation mode (PARALLEL=2). It calls the following in parallel each in a loop
            - perfevl() & print_progress(): performance simulation
            - wfsgrad() and shift_grad(): compute gradients and shift to gradlast
            - reconstruct() & filter_dm: wfs reconstruction and servo filtering
        - In parallel simulation, calls maos_isim() for each step. It then calls
            - sim_update_etf(): update sodium profile if needed
            - genatm(): update turbulence if not using frozen-flow
            - In parallel mode (PARALLEL=1)
                - perfevl_pre: wfsgrad_prep, reconstruct() in parallel and then wait
                - perfevl(): wfsgrad(), in parallel and then wait
                - shift_grad(): copy from grad to gradlast for reconstruct()
                - filter_dm(): servo filtering (dmreal)
            - In serial mode (PARALLEL=0), calls the following in sequence
                - perfevl(): (in closed loop) to evaluate the performance in science field
                - wfsgrad(): to compute WFS gradients
                - reconstruct():
                - shift_grad():
                - filter_dm():
                - perfevl():
            - setup_recon_tomo_update(): Update tomography parameter when needed
            - setup_recon_update(): Update recon parameters when needed.
            - print_progress(): show progress
        - free_simu(): to delete run time structs and close files


