#include "../lib/aos.h"
#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif
#if __STDC_VERSION__ < 199901L /*C99*/
#define inline __inline
#endif
/**
   Register the signal handler so that error() in the library does not cause matlab to crash.
 */
static void mex_signal_handler(int sig){
    if(sig){
	mexErrMsgTxt("Signal caught.\n");
    }else{
	info("signal 0 caught\n");
    }
}
static void(*default_handler)(int)=NULL;
static __attribute__((constructor)) void init(){
    if(!default_handler){
	default_handler=signal(SIGTERM, mex_signal_handler);
    }
}
static __attribute__((destructor)) void deinit(){
    if(default_handler){
	signal(SIGTERM, default_handler);
    }else{
	signal(SIGTERM, SIG_DFL);
    }
}
