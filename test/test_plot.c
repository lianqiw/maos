



#include "../lib/aos.h"
#define USE_SHM 0
int main(){
#if USE_SHM==1
    key_t key=12345;
    int m=300;
    shm_data *data=shm_getdata_rw(&key, m, m, m*m, SHM_CHAR);
    unsigned char *p=(unsigned char*)data->x;
    for(int i=0; i<m*m; i++){
	p[i]=(unsigned char)(256.*drand48());
    }
    fprintf(stderr, "Finish preparing data\n");
    shm_detach_rw(key);
    fprintf(stderr, "Detached\n");
    if(!shm_attach_ro(key,4)){
	fprintf(stderr, "Unable to reattach");
    }else{
	fprintf(stderr, "Reattached read only\n");
    }
    p=(unsigned char*)data->x;
    fprintf(stderr, "p=%p",p);
    /*for(int i=0; i<m*m; i++){ */
    /*	fprintf(stderr,"%d",p[i]); */
    /*} */
   
#endif
}
