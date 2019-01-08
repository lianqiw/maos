/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
