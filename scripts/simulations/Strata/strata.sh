#!/bin/bash
lgs4ngs1="powfs.nwfs=[4 1] wfs.thetax=[-1 -1 1 1 1]/2 wfs.thetay=[-1 1 -1 1 0]/2 "
lgs6ngs1="powfs.nwfs=[6 1] wfs.thetax=[1 0.5 -0.5 -1 -0.5 0.5 0]/2 wfs.thetay=[0 0.87 0.87 0 -0.87 -0.87 1]/2 \
    powfs0_llt.ox=[1 0.5 -0.5 -1 -0.5 0.5]*6.5  powfs0_llt.oy=[0 0.87 0.87 0 -0.87 -0.87]*6.5 "
lgs8ngs3="wfs_lgs_ttf_tt.conf powfs.nwfs=[8 1 2] wfs.thetax=[1 0.71 0 -0.71 -1 -0.71 -0 0.71 0 -0.87 0.87]/2 wfs.thetay=[0 0.71 1 0.71 0 -0.71 -1 -0.71 1 -0.5 -0.5]/2 \
    powfs0_llt.ox=[1 0.71 0 -0.71 -1 -0.71 -0 0.71]*6.5  powfs0_llt.oy=[0 0.71 1 0.71 0 -0.71 -1 -0.71 ]*6.5 "
declare -A fit #field of view points for DM fitting
declare -A evl #field of view points for evaluation

fit[lris2csu]="fit.thetax = [-0.5 -0.5 -0.5 -0.5 -0.5 -0.25 -0.25 -0.25 -0.25 -0.25 0 0 0 0 0 0.25 0.25 0.25 0.25 0.25 0.5 0.5 0.5 0.5 0.5]*600 \
            fit.thetay = [-0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5]*300 \
            fit.wt     = [1 1 1 1 1 1 1 1 1 1 1 1 0.5 1 1 1 1 1 1 1 1 1 1 1 1]  \
            fit.fov    = 1"

#Only use a quadrant
evl[lris2csu]="evl.thetax = [0 0.25 0.5 0 0.25 0.5 0 0.25 0.5]*600 \
            evl.thetay = [0 0 0 0.25 0.25 0.25 0.5 0.5 0.5]*300 \
            evl.wt     = [1 1 1 1 1 1 1 1 1] \
            evl.fov    = 1" 

fit[lris2ifu]="fit.thetax = [0 0.5 0 -0.5 0 0.5 -0.5 -0.5 0.5]*20 \
            fit.thetay = [0 0 0.5 0 -0.5 0.5 0.5 -0.5 -0.5]*7.2 \
            fit.wt     = [0.5 1 1 1 1 1 1 1 1]  \
            fit.fov    = 1"

evl[lris2ifu]="evl.thetax = [0 0.25 0.5 0 0.25 0.5 0 0.25 0.5]*20 \
            evl.thetay = [0 0 0 0.25 0.25 0.25 0.5 0.5 0.5]*7.2 \
            evl.wt     = [1 1 1 1 1 1 1 1 1] \
            evl.fov    = 1" 

fit[mosfire]="fit.thetax = [-0.5 -0.5 -0.5 -0.5 -0.5 -0.25 -0.25 -0.25 -0.25 -0.25 0 0 0 0 0 0.25 0.25 0.25 0.25 0.25 0.5 0.5 0.5 0.5 0.5]*360 \
            fit.thetay = [-0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5 -0.5 -0.25 0 0.25 0.5]*180 \
            fit.wt     = [1 1 1 1 1 1 1 1 1 1 1 1 0.5 1 1 1 1 1 1 1 1 1 1 1 1]  \
            fit.fov    = 1"

evl[mosfire]="evl.thetax = [0 0.25 0.5 0 0.25 0.5 0 0.25 0.5]*360\
            evl.thetay = [0 0 0 0.25 0.25 0.25 0.5 0.5 0.5]*180 \
            evl.wt     = [1 1 1 1 1 1 1 1 1] \
            evl.fov    = 1" 
declare -A config
#For GLAO, set powfs.dtrat to 600 for NGS to make it slow guider only 
config[lris2csu]=" powfs.astscale=[600 540] ${fit[lris2csu]} ${evl[lris2csu]} powfs.dtrat=[1 600] " #LRIS2 FoV 10'x5'. LGS is 10'x10'. Guider is at 4.5' off axis
config[lris2ifu]=" powfs.astscale=[150 540] ${fit[lris2ifu]} ${evl[lris2ifu]} powfs.dtrat=[1 600] " #LRIS2 IFU FoV 20x7.2 arcsec. LGS is 2.5'x2.5'
config[mosfire]=" powfs.astscale=[540 800] ${fit[mosfire]} ${evl[mosfire]} powfs.dtrat=[1 600] " #MOSFIRE spectrograph 6'x3' fov. Guider is 6.7' off axis.
config[zimager]=" powfs.astscale=[36  240] fit_cir60.conf fit.fov=180 evl_x.conf evl.fov=180 powfs.dtrat=[1 600] " #Zshooter imager with circular FoV D=3'
config[zspec]=" powfs.astscale=[15.2  240] fit_oa.conf fit.fov=0 evl_x.conf evl.fov=180 powfs.dtrat=[1 600] " #Zshooter spectrograph in LTAO mode
config[kapa]=" powfs.astscale=[15.2   0] powfs.wvl=[0.589 1.25]  powfs.pixtheta=[1 0.05] fit_oa.conf fit.fov=0 evl_x.conf evl.fov=20 " #Kapa LTAO
#config[kola]="$lgs8ngs3 powfs.astscale=[60 30 30] powfs.wvl=[0.589 1.25 1.25] powfs.nwvl=[1 1 1] powfs.pixtheta=[1 0.05 0.05] fit_cir60.conf fit.fov=60 evl_x.conf evl.fov=60 dm.offset+=[1/3 2/3] dm.ht+=[6000 10000] sim.dt=1/1500 " #Example 
declare -A algs
algs[idealfit]="sim.idealfit=1" #best performance by fitting turbulence directly to DM
algs[glao]="recon.alg=1" #instead of averaging gradients which cannot handle misregistration or differential measurement noise, we use LSQ method.
algs[tomo]="" #default
atms="25pGL_25pFA 25pGL_50pFA 25pGL_75pFA 50pGL_25pFA 50pGL_50pFA 50pGL_75pFA 75pGL_25pFA 75pGL_50pFA 75pGL_75pFA "
atms="50pGL_50pFA"
#zas="30 45 55 60"
zas="30"
extra="sim.seeds=[1 10 20 30]"
systems="zspec zimager lris2csu lris2ifu mosfire kapa "
#systems="zspec"
alg0="idealfit tomo"
rne=0.1
for sys in $systems;do   
    for za in $zas;do
        for atm in $atms;do
            for alg in $alg0 ;do
                if [ $alg = idealfit ];then
                    lmags=0
                    nlgss=0
                else
                    lmags="6 7.5 8"
                    nlgss="4 "
                fi
                for nlgs in $nlgss ;do
                    for lmag in $lmags ;do
                        ast="lgs${nlgs}ngs1"
                        if [ $alg = idealfit ];then
                            fd=${alg}
                        else
                            fd=${alg}_rne${rne}_${nlgs}_${lmag}
                        fi
                        maos -c strata.conf ${config[${sys}]} ${!ast} powfs.mag=[${lmag} 8.5] atm/atm_mk${atm}.conf powfs.rne=[$rne 0.5] ${algs[${alg}]} $extra -o strata/${sys}/za${za}/${atm}/${fd} -d
                    done
                done
            done
        done
    done
done
