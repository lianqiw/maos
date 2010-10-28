addpath('../mex');
addpath('~/work/laos/LAOS_SVN_SKY_TEST/recons');

xloc=read('xloc.dblc.gz');
parms.r0=0.15;
dxts=[1/4,1/4,1/2,1/2,1/2,1/2];
wts=[0.6523 0.1723 0.0551 0.0248 0.0736 0.0219];
for ips=1:length(xloc)
    parms.atmr(ips).dxt=dxts(ips);
    parms.atmr(ips).wt=wts(ips);
    trp{ips}=mktrp(parms, xloc{ips},ips, 1);
end
L2=read('L2.spc.gz');
