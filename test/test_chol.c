/*
  Copyright 2009-2026 Lianqi Wang
  
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
/*
int main(){
    dspcell *FLM=dspcellread("FLM.bin");
    dsp *F=P(FLM,0);
    spchol *Fchol=chol_factorize(F);
 
    rand_t rrand;
    seed_rand(&rrand,1);
    dmat *y=dnew(F->m,10);
    drandu(y,1,&rrand);
    dmat *x=NULL;
    chol_solve(&x,Fchol,y);
    writebin(y,"y");
    writebin(x,"x");
    chol_save(Fchol,"Chol");
    chol_free(Fchol);
    dfree(x);
    dfree(y);
    dspcellfree(FLM);
}
*/
TIC;
int main(int argc, char* argv[]){
    /*dsp *RLMc1=dspread("RLMc_old.bin"); */
    if(argc!=2){
	info("Need 1 argument for RLM file.\n");
	return 0;
    }
    dspcell *RLM=dspcellread("%s",argv[1]);
    dsp *RLMc=dspcell2sp(RLM);
    tic;info("chol ...");
    spchol *R1=chol_factorize(RLMc);
    toc("done");
    chol_save(R1, "%s_chol.bin", argv[1]);
    rand_t rstat;
    seed_rand(&rstat,1);
    dmat *y=dnew(RLMc->nx, 1);
    drandn(y, 1, &rstat);
    dmat *x=NULL, *x2=NULL, *x3=NULL;
    tic;
    chol_solve(&x, R1, y);
    toc("cholmod");tic;
    chol_solve(&x, R1, y);
    toc("cholmod");tic;
    chol_convert(R1, 0);//convert to Cl
    
    chol_solve(&x2, R1, y);
    toc("lower");tic;
    chol_solve(&x2, R1, y);
    toc("lower");tic;

    R1->Cu=dsptrans(R1->Cl);
    dspfree(R1->Cl);
    chol_solve(&x3, R1, y);
    toc("upper");tic;
    chol_solve(&x3, R1, y);
    toc("upper");tic;
    
    writebin(y,"y");
    writebin(x,"x");
    writebin(x2,"x2");
    writebin(x3,"x3");
    chol_free(R1);
    dspfree(RLMc);
    dspcellfree(RLM);
}
