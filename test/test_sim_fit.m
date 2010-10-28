addpath('../mex');
opdr=read('opdr.bin.gz');
HX=read('HX.bin.gz');
HA=read('HA.bin.gz');
ploc=read('ploc.bin.gz');
aloc0=read('dm0_loc.bin.gz');
aloc1=read('dm1_loc.bin.gz');

opd=cell_mm(HX, opdr);
W0=read('W0.bin.gz');
W1=read('W1.bin.gz');
HAT=celltransp(HA);
wtevl=[4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
for ievl=1:9
    opd{ievl}=opd{ievl}*wtevl(ievl);
    W0opd{ievl}=W0*opd{ievl};
    W1opd{ievl}=-W1*(W1'*opd{ievl});
    Wopd{ievl}=W0opd{ievl}+W1opd{ievl};
    opd{ievl}'*Wopd{ievl}
end



adm=cell_mm(HAT,opd);

fit_rhs=read('fit_rhs.bin.gz');

FRM=read('FRM.bin.gz');
FRU=read('FRU.bin.gz');
FRV=read('FRV.bin.gz');
adm0=cell_mm(FRM,opdr);
FRVT=celltransp(FRV);
tmp=cell_mm(FRVT,opdr);
adm01=cell_mm(FRU,tmp);

