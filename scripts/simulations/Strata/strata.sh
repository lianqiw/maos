#!/bin/bash
lgs4ngs0="powfs.nwfs=[4 0] wfs.thetax=[-1 -1 1 1]/2 wfs.thetay=[-1 1 -1 1]/2 powfs0_llt.ox=[-1 -1 1 1]*6.5  powfs0_llt.oy=[-1 1 -1 1]*1"
lgs4ngs1="powfs.nwfs=[4 1] wfs.thetax=[-1 -1 1 1 1]/2 wfs.thetay=[-1 1 -1 1 0]/2 powfs0_llt.ox=[-1 -1 1 1]*6.5  powfs0_llt.oy=[-1 1 -1 1]*1"
lgs6ngs1="powfs.nwfs=[6 1] wfs.thetax=[1 0.5 -0.5 -1 -0.5 0.5 0]/2 wfs.thetay=[0 0.87 0.87 0 -0.87 -0.87 1]/2 \
    powfs0_llt.ox=[1 0.5 -0.5 -1 -0.5 0.5]*6.5  powfs0_llt.oy=[0 0.87 0.87 0 -0.87 -0.87]*6.5"
lgs8ngs3="wfs_lgs_ttf_tt.conf powfs.nwfs=[8 1 2] wfs.thetax=[1 0.71 0 -0.71 -1 -0.71 -0 0.71 0 -0.87 0.87]/2 wfs.thetay=[0 0.71 1 0.71 0 -0.71 -1 -0.71 1 -0.5 -0.5]/2 \
    powfs0_llt.ox=[1 0.71 0 -0.71 -1 -0.71 -0 0.71]*6.5  powfs0_llt.oy=[0 0.71 1 0.71 0 -0.71 -1 -0.71 ]*6.5"
declare -A fit #field of view points for DM fitting
declare -A evl #field of view points for evaluation
fit[iris2rec]="fit.thetax = [-0.5 -0.5 -0.5 -0.5 -0.5 -0.25 -0.25 -0.25 -0.25 -0.25 0 0 0 0 0 0.25 0.25 0.25 0.25 0.25 0.5 0.5 0.5 0.5 0.5]*600 \
            fit.thetay = [-0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5]*300 \
            fit.wt     = [1 1 1 1 1 1 1 1 1 1 1 1 0.5 1 1 1 1 1 1 1 1 1 1 1 1]  \
            fit.fov    = 1"
#Only use a quadrant
evl[iris2rec]="evl.thetax = [0 0.25 0.5 0 0.25 0.5 0 0.25 0.5]*600 \
            evl.thetay = [0 0 0 0.25 0.25 0.25 0.5 0.5 0.5]*300 \
            evl.wt     = [1 1 1 1 1 1 1 1 1] \
            evl.fov    = 1" 

fit[iris2ifu]="fit.thetax = [0 0.5 0 -0.5 0 0.5 -0.5 -0.5 0.5]*20 \
            fit.thetay = [0 0 0.5 0 -0.5 0.5 0.5 -0.5 -0.5]*7.2 \
            fit.wt     = [0.5 1 1 1 1 1 1 1 1]  \
            fit.fov    = 1"

evl[iris2ifu]="evl.thetax = [0 0.25 0.5 0 0.25 0.5 0 0.25 0.5]*20 \
            evl.thetay = [0 0 0 0.25 0.25 0.25 0.5 0.5 0.5]*7.2 \
            evl.wt     = [1 1 1 1 1 1 1 1 1] \
            evl.fov    = 1" 
declare -A config
config[lris2rec]="$lgs6ngs1 powfs.astscale=[804 960] ${fit[iris2rec]} ${evl[iris2rec]} " #LRIS2 with rectangular FoV 10'x5'
config[lris2cir]="$lgs6ngs1 powfs.astscale=[804 960] fit_cir60.conf fit.fov=600 evl_x.conf evl.fov=600 " #LRIS2 with circular FoV D=10'
config[lris2ifu]="$lgs6ngs1 powfs.astscale=[804 960] ${fit[iris2ifu]} ${evl[iris2ifu]} " #LRIS2 IFU with 20x7.2 arcsec
config[mosfire]="$lgs6ngs1 powfs.astscale=[816 720] fit_cir60.conf fit.fov=367 evl_x.conf evl.fov=367 " #MOSFIRE with circular FoV D=6.12' 
config[zimager]="$lgs4ngs0 powfs.astscale=[36  240] fit_cir60.conf fit.fov=180 evl_x.conf evl.fov=180 " #Zimager with circular FoV D=3'
config[kola]="$lgs8ngs3 powfs.astscale=[60 30 30] fit_cir60.conf fit.fov=60 evl_x.conf evl.fov=60 dm_triple.conf dm.dx=[0.185 0.185 0.185] dm.ht=[-80 6000 10000] sim.dt=1/1500 " #Example 
declare -A algs
algs[idealfit]="sim.idealfit=1" #best performance by fitting turbulence directly to DM
algs[glao]="recon.glao=1"
algs[tomo]="recon.glao=0"
#atms="25pGL_25pFA 25pGL_50pFA 25pGL_75pFA 50pGL_25pFA 50pGL_50pFA 50pGL_75pFA 75pGL_25pFA 75pGL_50pFA 75pGL_75pFA "
atms="50pGL_50pFA"
#zas="30 45 55 60"
zas="30"
systems="lris2rec lris2cir lris2ifu mosfire zimager kola"
for sys in $systems;do
    if [ "$sys" = "kola" ];then
        alg0="idealfit tomo "
    else
        alg0="idealfit tomo glao"
    fi
    for za in $zas;do
        for atm in $atms;do
            for alg in $alg0; do
                maos -c strata.conf ${config[${sys}]} atm/atm_mk${atm}.conf ${algs[${alg}]} -o strata/${sys}/za${za}/${atm}/${alg} -d
                #append powfs.step=[0 1e5] to disable NGS WFS
            done
        done
    done
done
