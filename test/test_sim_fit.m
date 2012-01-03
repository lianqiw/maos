addpath('../mex');
opdr=read('opdr.bin');
HX=read('HX.bin');
HA=read('HA.bin');
ploc=read('ploc.bin');
aloc0=read('dm0_loc.bin');
aloc1=read('dm1_loc.bin');

opd=cell_mm(HX, opdr);
W0=read('W0.bin');
W1=read('W1.bin');
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

fit_rhs=read('fit_rhs.bin');

FRM=read('FRM.bin');
FRU=read('FRU.bin');
FRV=read('FRV.bin');
adm0=cell_mm(FRM,opdr);
FRVT=celltransp(FRV);
tmp=cell_mm(FRVT,opdr);
adm01=cell_mm(FRU,tmp);

