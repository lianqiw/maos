#if defined (__CYGWIN__)
#include "proc.h"
//Largely not implemented.
const char *get_job_progname(void){
    return "maos";
}
int get_job_mem(void){
    return 0;
}
double get_job_launchtime(int pid){
    (void) pid;
    return 0;
}
int get_usage_running(void){
    return 0;
}
double get_usage_load(void){
    return 0;
}
double get_usage_mem(void){
    return 0;
}
double read_self_cpu(void){
    return 0;
}
int read_usage_cpu(long *user, long *tot){
    *user=0;
    *tot=0;
    return 0;
}
int get_ncpu(void){
    return 1;
}
#endif
