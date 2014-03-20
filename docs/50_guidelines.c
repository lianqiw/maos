/**
   \page page50 Programming Guidelines

   \section sect-guide Guidelines

   When modifying the code, please adhere to the following guidelines as
   strictly as possible.

   Memory allocation:

   - When declaring a struct, always initialize it to zero by calling calloc or
   memset, unless you have a good reason not to.

   - Use calloc instead of malloc to initialize the memory unless for large
   arrays that you will initialize immediately.

   Functions:
      
      
   - The input parameters that is also the output should be grouped in the
   beginning.

   - In any function, there should usually be at most 1 return statement unless
   NULL or error code are returned before matural.

   - One utility function should be handle only one major mask to maximize
   reusability.

   - Input arguments to utility functions (in lib folder) should be basic types
   to maximize resuability. Input arguments to simulation functions can be
   wrapped (like simu) to maximum readability.
   
   Others:
   
   - Avoid non constant static variables/pointers unless you have a good reason
   to do so. static variables are bad for multi-threading and hard to free
   after usage.

   - whenever modify something temporarily for debugging, make a warning. let
   it easily identifiable. not make a hidden bug.
   
   - Do not hard code adjustable parameters.
   
   - Always declare unchanged variables to constants. declare pointers as
   restrict if possible.

   - Do not include system headers in header files unless necessary. include
   thoese headers in .c file. 

   - Avoid function casting. It will hide data type check and hide bugs.

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


 */
