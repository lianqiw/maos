/**
   \page page5 Programming Guidelines

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

 */
