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

rand_t rstat;/*a struct for each sequence*/
rand_t rstat2;


static void test_number(){
    rand_t a,b;
    seed_rand(&a,1);
    dbg("First number is %g\n",randn(&a));
    seed_rand(&a,1);
    unsigned int i=lrand(&a);
    dbg("First number is %u\n",i);
    i=lrand(&a);
    dbg("First number is %u\n",i);
    seed_rand(&b,i);
    dbg("First number is %g\n",randn(&b));
}
static void test_speed(){
    double *p;
    int N,i;
    clock_t ck;
    int seed=2;
    N=1000000;
    p=mymalloc(N,double);
    memset(p, 0, sizeof(double)*N);
    seed_rand(&rstat, seed);
    srand(seed);
    ck=clock();
    seed_rand(&rstat, seed);
    ck=clock();
    for(i=0; i<N; i++){
	p[i]=randu(&rstat);
    }
    printf("randu elapsed %f seconds\n", (double)(clock()-ck)/CLOCKS_PER_SEC);
    writedbl(p, N, 1,"randu.bin");
    srand(1);
    ck=clock();
    for(i=0; i<N; i++){
	p[i]=rand();
    }
    printf("Rand elapsed %f seconds\n", (double)(clock()-ck)/CLOCKS_PER_SEC);
    /*ck=clock();
    for(i=0; i<N; i++){
	p[i]=slow_randn(&rstat);
    }
    printf("Slow Randn elapsed %f seconds\n", (double)(clock()-ck)/CLOCKS_PER_SEC);
    writedbl("slow_randn.bin", p, N, 1);
    */
    ck=clock();
    for(i=0; i<N; i++){
	p[i]=randn(&rstat)*10;
    }
    printf("Fast Randn elapsed %f seconds\n", (double)(clock()-ck)/CLOCKS_PER_SEC);
    writedbl(p, N, 1, "randn.bin");

    seed_rand(&rstat2,seed);
    ck=clock();
    for(i=0; i<N; i++){
	p[i]=(double)randp(&rstat2,30);
    }
    printf("Randp elapsed %f seconds\n", (double)(clock()-ck)/CLOCKS_PER_SEC);
    writedbl(p, N, 1,"randp.bin");
    free(p);
}
int main(){
    test_number();
    test_speed();
    return 0;
}
