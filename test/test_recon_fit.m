
HX=read('HX.bin');
HA=read('HA.bin');

W0=read('W0.bin'); 
W1=read('W1.bin');

wtevl=[4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
%wtevl=[1 zeros(1,8)];
FRM=read('FRM.bin');
FRM2=cell(2,6);
for ips=1:6
    for idm=1:2
        FRM2{idm,ips}=0;
        for ievl=1:9
            FRM2{idm,ips}=sparse(HA{ievl,idm}'*W0*HX{ievl,ips}*wtevl(ievl)+FRM2{idm,ips});
        end
    end
end
diff=full(cell2mat(FRM))-full(cell2mat(FRM2));
norm(diff(:))
FRU2=cell(2,1);
for idm=1:2
    for ievl=1:9
        FRU2{idm}(:,ievl)=HA{ievl,idm}'*W1*sqrt(wtevl(ievl));
    end
end
FRU=read('FRU.bin');
norm(cell2mat(FRU)-cell2mat(FRU2))


FRV2=cell(6,1);
for ips=1:6
    for ievl=1:9
        FRV2{ips}(:,ievl)=HX{ievl,ips}'*W1*sqrt(wtevl(ievl));
    end
end
FRV=read('FRV.bin');
norm(cell2mat(FRV)-cell2mat(FRV2))


%%Order is different!
FLM2=cell(2,2);
for idm=1:2
    for jdm=1:2
        FLM2{jdm,idm}=0;
        for ievl=1:9
            FLM2{jdm, idm}=sparse(FLM2{jdm, idm}+ HA{ievl,jdm}'*W0*HA{ievl, ...
                                idm}*wtevl(ievl));
        end
    end
end
FLM=read('FLM.bin');
flm=full(cell2mat(FLM));
flm2=full(cell2mat(FLM2));
norm(flm(:)-flm2(:))
