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

#include "hexpupil.h"
/**
    Generate center of segments */
dmat *hexcent(const dmat *segment){
    if(!segment || P(segment, 0)!=1 || PN(segment)<4 || P(segment, 4)<=0){
        warning("segment configuration is invalid.\n");
        dshow(segment, "segment");
        return NULL;
    }
    const real dseg=P(segment, 1);//vertex to vertex diameter
    const real dgap=P(segment, 2);
    const real dout=P(segment, 3);
    const real din=(PN(segment)>4 && P(segment, 4)<dout)?P(segment, 4):0;
    const real roc=2*dout*(PN(segment)>5?P(segment, 5):1); //radius of curvature of the mirror. used to scale dgap. default is 2*dout
    const real drot=PN(segment)>6?P(segment, 6)*M_PI/180.:0;//pupil rotation CCW
    
    const real dseg1=dseg+dgap*2/sqrt(3.);//vertex to vertex segment separation (along y)
    const real dseg2=dseg1*sqrt(3)*0.5;//edge to edge segment separation (along x)
    const real dseg3=0.75*dseg1;//height increase of stacked hexagons. (vertex on top)
    const real yr=dseg*0.5;//segment y radius, without gap
    const real xr=yr*sqrt(3)*0.5;//segment x radius, without gap
    const real routsq=dout*dout*0.25;
    const int rskip=(din<=0)?0:round((din+dseg1)*0.5/dseg3);
    const int rtot=round((dout-dseg1)*0.5/dseg3)+1;
    const int ntot=3*rtot*rtot-3*rtot+1;
    const int nskip=rskip==0?0:3*rskip*rskip-3*rskip+1;
    dmat *cents=dnew(ntot-nskip, 2);
    const real cs=cos(M_PI/3);
    const real ss=sin(M_PI/3);
    int ic=0;//segment counter
    real r=rskip*dseg2;
    if(rskip==0) ic+=1; //first ring has 1 segment
    dmat *edge=dnew(rtot, 2);
    //info("din=%g, dout=%g, rskip=%d, rtot=%d, nskip=%d, ntot=%d\n", din, dout, rskip, rtot, nskip, ntot);
    for(int ring=rskip+1; ring<rtot+1; ring++){//rings
        real x0=r*cs;
        real y0=r*ss;
        real dseg2r;
        if(roc>0 && isfinite(roc)){
            real dgap2=dgap*sqrt(1-pow(r/roc, 2));
            dseg2r=(dseg+dgap2*2/sqrt(3))*sqrt(3)*0.5;
        }else{
            dseg2r=dseg2;
        }
        r+=dseg2r;
        const int nxsave=edge->nx;
        edge->nx=ring-1;
        //create an edge
        for(int xsub=0; xsub<ring-1; xsub++){
            P(edge, xsub, 0)=x0-xsub*dseg2r;
            P(edge, xsub, 1)=y0;
        }
        for(int arot=0; arot<6; arot++){
            if(ic+ring-1>ntot-nskip){
                warning("hexcent overflows. ic=%d, ring-1=%d, ntot-nskip=%d\n", ic, ring-1, ntot-nskip);
                break;
            }
            for(int xsub=0; xsub<ring-1; xsub++){
                real xc=fabs(P(edge, xsub, 0));
                real yc=fabs(P(edge, xsub, 1));
                real r1sq=(xc+xr)*(xc+xr)+yc*yc;
                real r2sq=xc*xc+(yc+yr)*(yc+yr);
                if(r1sq<routsq && r2sq<routsq){
                    P(cents, ic, 0)=P(edge, xsub, 0);
                    P(cents, ic, 1)=P(edge, xsub, 1);
                    ic++;
                }
            }
            //rotate the edge to the right angle
            if(arot<6) drotvec(edge, M_PI/3);
        }
        edge->nx=nxsave;
    
    }
    dresize(cents, ic, 2);
    dfree(edge);
    if(drot) drotvec(cents, drot);
    return cents;
}
/**
    Fill a hexagon to A
*/
void hexfill(map_t *A, real cx, real cy, real dseg, real value){
    const real ox=A->ox;
    const real oy=A->oy;
    const real dx=A->dx;
    //convert to unit spacing
    cx=(cx-ox)/dx;
    cy=(cy-oy)/dx;
    dseg=dseg/dx;
    const real dsegx=dseg*sqrt(3)/2;//edge to edge
    const real xthres=dsegx/2;//along x
    const real ythres=dseg/2;//along y
    int ix0=MAX(0, ceil(cx-xthres));
    int ix1=MIN(A->nx, floor(cx+xthres)+1);
    int iy0=MAX(0, ceil(cy-ythres));
    int iy1=MIN(A->ny, floor(cy+ythres)+1);
    const real slope=1./sqrt(3.);
    for(int ix=ix0; ix<ix1; ix++){
        real ixc=fabs(ix-cx);//x distance to center
        real ythres2=ythres/2+fabs(xthres-ixc)*slope;
        for (int iy=iy0; iy<iy1; iy++){
            real iyc=fabs(iy-cy);//y distance to center
            if(iyc<ythres2){
                P(A, ix, iy)=value;
            }
        }
    }
}
/**
    Fill a hexagon to A using gray pixel approach
*/
void hexfill_gray(map_t *A, real cx, real cy, real dseg, real value){
    const real ox=A->ox;
    const real oy=A->oy;
    const real dx=A->dx;
    //convert to unit spacing
    cx=(cx-ox)/dx;
    cy=(cy-oy)/dx;
    dseg=dseg/dx;
    const real dsegx=dseg*sqrt(3)/2;//edge to edge
    const real xthres=dsegx/2;//along x
    const real ythres=dseg/2;//along y
    int ix0=MAX(0, floor(cx-xthres-0.5));
    int ix1=MIN(A->nx, ceil(cx+xthres+0.5)+1);
    int iy0=MAX(0, floor(cy-ythres-0.5));
    int iy1=MIN(A->ny, ceil(cy+ythres+0.5)+1);
    const real slope=1./sqrt(3.);
    const real xedge=0.5;
    const real yedge=slope;//project distance to slanted edge to y axis 
    for(int ix=ix0; ix<ix1; ix++){
        real ixc=fabs(ix-cx);//x distance to center
        real xwt=1;
        if(ixc+xedge<=xthres){
            xwt=1;
        }else if(ixc>=xthres+xedge){
            xwt=0; continue;
        }else{
            xwt=1-fabs(ixc-xthres)/xedge;
        }
        real ythres2=ythres/2+fabs(xthres-ixc)*slope;
        for (int iy=iy0; iy<iy1; iy++){
            real iyc=fabs(iy-cy);//y distance to center
            real ywt=1;
            if(iyc+yedge<ythres2){
                ywt=1;
            }else if(iyc>ythres2+yedge){
                ywt=0; continue;
            }else{
                ywt=1-fabs(iyc-ythres2)/yedge;
            }
            P(A, ix, iy)=value*xwt*ywt;
        }
    }
}
/**
    Generate amplitude map of segments */
map_t *hexamp(const dmat *segment, real dx, int nx, int ny, int gray){
    if(!segment || P(segment, 0)!=1 || PN(segment)<4 || P(segment, 4)<=0){
        warning("segment configuration is invalid.\n");
        return NULL;
    }
    dmat *cents=hexcent(segment);
    const real dseg=P(segment, 1);
    const real dout=P(segment, 3);
    const real drot=PN(segment)>6?P(segment, 6)*M_PI/180.:0;
    if(drot){
        drotvec(cents, -drot);//undo rotation
    }
    if(!nx || !ny){
        nx=nextfftsize(dout/dx);
        ny=nextfftsize(dout/dx);
    }
    map_t *opd=mapnew(nx, ny, dx, dx);
    for(int iseg=0; iseg<NX(cents); iseg++){
        if(gray){
            hexfill_gray(opd, P(cents, iseg, 0), P(cents, iseg, 1), dseg, 1);
        }else{
            hexfill(opd, P(cents, iseg, 0), P(cents, iseg, 1), dseg, 1);
        }
    }
    if(drot){
        dmaprot(opd, drot);//pupil rotation
    }
    return opd;
}
/**
    Generate p/t/t mode of segments in transposed mode. faster but not accuracy when there is rotation due to interpolation.*/
dsp *hexptt_trans(const dmat *segment, real dx, int nx, int ny){
    if(!segment || P(segment, 0)!=1 || PN(segment)<4 || P(segment, 4)<=0){
        warning("segment configuration is invalid.\n");
        return NULL;
    }
    dmat *cents=hexcent(segment);
    const real dseg=P(segment, 1);
    const real dout=P(segment, 3);
    const real drot=PN(segment)>6?P(segment, 6)*M_PI/180.:0;
    dmat *cents2=NULL;//unrotated cents
    if(drot){
        cents2=ddup(cents);
        drotvec(cents2, -drot);//undo rotation
    }else{
        cents2=dref(cents);
    }
    if(!nx || !ny){
        nx=nextfftsize(dout/dx);
        ny=nextfftsize(dout/dx);
    }
    map_t *opd=mapnew(nx, ny, dx, dx);

    int iseg=0;
    hexfill(opd, P(cents2, iseg, 0), P(cents2, iseg, 1), dseg, iseg+1);
    const int nzseg=floor(dsum(DMAT(opd))*1.1);//number of non-zeros per segment
    const int nseg=NX(cents);
    for(iseg=1; iseg<nseg; iseg++){
        hexfill(opd, P(cents2, iseg, 0), P(cents2, iseg, 1), dseg, iseg+1);
    }
    if(drot){
        dmaprot(opd, drot);
    }
    dsp *ptt=dspnew(3*nseg, nx*ny, nzseg*nseg*3);//transposed version first
    spint *pp=ptt->pp;
    spint *pi=ptt->pi;
    real *px=ptt->px;
    long count=0;    
    const real ox=opd->ox;
    const real oy=opd->oy;
    const real ss=sin(drot);
    const real cs=cos(drot);
    for(int iy=0; iy<ny; iy++){
        for(int ix=0; ix<nx; ix++){
            real v=(int)round(P(opd, ix, iy));
            if(fabs(P(opd, ix, iy)-v)<1e-6){
                P(opd, ix, iy)=v;
            }else{
                P(opd, ix, iy)=0;
            }
        }
    }

    for(int iy=0; iy<ny; iy++){
        for(int ix=0; ix<nx; ix++){
            pp[ix+iy*nx]=count;
            iseg=(int)P(opd, ix, iy)-1;
            if(iseg>=0 && iseg<nseg){
                real cx=(ix*dx)+ox-P(cents, iseg, 0);
                real cy=(iy*dx)+oy-P(cents, iseg, 1);
                pi[count]=iseg;//piston
                px[count]=1;
                count++;
                pi[count]=iseg+nseg;//tip
                px[count]=cx*cs+cy*ss;
                count++;
                pi[count]=iseg+nseg*2;//tilt
                px[count]=cx*(-ss)+cy*cs;
                count++;
            }
        }
    }
    pp[nx*ny]=count;
    if(count>ptt->nzmax){
        error("overflow\n");
    }
    dspsetnzmax(ptt, count);
    dfree(cents);
    dfree(cents2);
    mapfree(opd);
    return ptt;
}
/**
    Generate p/t/t mode of segments */
dsp *hexptt(const dmat *segment, real dx, int nx, int ny){
    if(!segment || P(segment, 0)!=1 || PN(segment)<4 || P(segment, 4)<=0){
        warning("segment configuration is invalid.\n");
        return NULL;
    }
    const real drot=PN(segment)>6?P(segment, 6)*M_PI/180.:0;
    if(fmod(drot, (M_PI/2))<1e-5){//faster
        dsp *ptt_t=hexptt_trans(segment, dx, nx, ny);
        dsp *ptt=dsptrans(ptt_t);
        dspfree(ptt_t);
        return ptt;
    }
    dmat *cents=hexcent(segment);
    const real dseg=P(segment, 1);
    const real dout=P(segment, 3);
    dmat *cents2=NULL;//unrotated cents
    if(drot){
        cents2=ddup(cents);
        drotvec(cents2, -drot);//undo rotation
    }else{
        cents2=dref(cents);
    }
    if(!nx || !ny){
        nx=nextfftsize(dout/dx);
        ny=nextfftsize(dout/dx);
    }
    map_t *opd=mapnew(nx, ny, dx, dx);

    int iseg=0;
    hexfill(opd, P(cents2, iseg, 0), P(cents2, iseg, 1), dseg, 1);
    const int nzseg=floor(dsum(DMAT(opd))*1.1);//number of non-zeros per segment
    const int nseg=NX(cents);
    if(drot){
        dmaprot(opd, drot);
    }
    dsp *ps=dspnew(nx*ny, nseg, nzseg*nseg);//transposed version first
    dsp *tx=dspnew(nx*ny, nseg, nzseg*nseg);//transposed version first
    dsp *ty=dspnew(nx*ny, nseg, nzseg*nseg);//transposed version first
    long count=0;    
    const real ox=opd->ox;
    const real oy=opd->oy;
    const real ss=sin(drot);
    const real cs=cos(drot);

    for(iseg=0; iseg<nseg; iseg++){
        ps->pp[iseg]=count;
        tx->pp[iseg]=count;
        ty->pp[iseg]=count;
        dzero(DMAT(opd));
        hexfill(opd, P(cents2, iseg, 0), P(cents2, iseg, 1), dseg, 1);
        if(drot){
            dmaprot(opd, drot);
        }
        for(int iy=0; iy<ny; iy++){
            for(int ix=0; ix<nx; ix++){
                if(P(opd, ix, iy)>0.5){//don't use >= to avoid double counting
                    real cx=(ix*dx)+ox-P(cents, iseg, 0);
                    real cy=(iy*dx)+oy-P(cents, iseg, 1);
                    ps->pi[count]=ix+nx*iy;
                    tx->pi[count]=ix+nx*iy;
                    ty->pi[count]=ix+nx*iy;

                    ps->px[count]=1;
                    tx->px[count]=   cs*cx+ss*cy;
                    ty->px[count]=(-ss)*cx+cs*cy;
                    count++;
                }
            }
        }
    }
    ps->pp[nseg]=count;
    tx->pp[nseg]=count;
    ty->pp[nseg]=count;
    if(count>ps->nzmax){
        error("overflow\n");
    }
    dspsetnzmax(ps, count);
    dspsetnzmax(tx, count);
    dspsetnzmax(ty, count);
    dfree(cents);
    dfree(cents2);
    mapfree(opd);
    dsp* tmp=dspcat(ps, tx, 1);
    dspfree(ps);
    dspfree(tx);
    dsp* ptt=dspcat(tmp, ty, 1);
    dspfree(tmp);
    dspfree(ty);
    return ptt;
}
