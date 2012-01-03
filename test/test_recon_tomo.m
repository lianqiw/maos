%test the setting up of tomography in setup_recon
split=1
%right hand side
pts=read('powfs0_loc.bin');
amp=read('powfs0_amp.bin');
orig=read('powfs0_orig.bin');
nsa=size(orig,1);
%TT=zeros(nsa*2,2);
%TT(1:end/2,1)=1;
%TT(end/2+1:end,2)=1;
%sanea=ones(nsa*2,1)*(20/206265000)^2;
%saneai=(1./sanea);
%TTw=zeros(size(TT));
%TTw(:,1)=TT(:,1).*saneai;
%TTw(:,2)=TT(:,2).*saneai;
%PTT=inv(TTw'*TT)*TTw';
TT=read('TT.bin');
PTT=read('PTT.bin');

G0=read('G0.bin');
[nwfs,nps]=size(G0);
RRM=cell(nps,nwfs);
RRM2=read('RRM.bin');

RRV=cell(nwfs,nwfs);
PTTt=PTT';
RRV2 = read('RRV.bin');

RRU=cell(nps,nwfs);
RRU2 = read('RRU.bin');
saneai=read('saneai.bin');
disp('RHS');
if split==1
    skipwfs=[zeros(6,1); ones(3,1);];
else
    skipwfs=zeros(9,1);
end
for iwfs=1:nwfs
    if skipwfs(iwfs)
        continue
    end
    nsa2=size(G0(iwfs,1),1);
    for ips=1:nps
        RRM{ips,iwfs}=G0{iwfs,ips}'*saneai{iwfs,iwfs};
        diff=norm(nonzeros(RRM{ips,iwfs}-RRM2{ips,iwfs}));
        fprintf('iwfs=%d,ips=%d,RRM diff=%g\n',iwfs,ips,diff);
    end
    if iwfs<7
        RRV{iwfs,iwfs}=PTT{iwfs}';
        diff=norm(RRV{iwfs,iwfs}-RRV2{iwfs,iwfs});
        fprintf('iwfs=%d,RRV diff=%g\n',iwfs,diff);
        for ips=1:nps
            RRU{ips,iwfs}=RRM{ips,iwfs}*TT{iwfs};
            diff=norm(RRU{ips,iwfs}-RRU2{ips,iwfs});
            fprintf('iwfs=%d,ips=%d,RRU diff=%g\n',iwfs,ips,diff);
        end
    end
   
end
disp('LHS');
%left hand side
RLU=cell(nps,nwfs);
RLV=cell(nps,nwfs);
RLM2=read('RLM.bin');
%RLM20=read('RLM_0.bin');
RLU2=read('RLU.bin');
RLV2=read('RLV.bin');
L2=read('L2.bin');
xloc=read('xloc.bin');
Z=cell(nps,1);
%when compare sparse matrices, becareful about the ordering.
%for ips=1:nps
%    for jps=1:nps
%        RLM2{ips,jps}=sparse(full(RLM2{ips,jps}));
%        %RLM20{ips,jps}=sparse(full(RLM20{ips,jps}));
%    end
%end
ZZT=read('ZZT.bin');
RLM=cell(nps,nps);
mt0=read('/home/lianqiw/work/laos/runs/nominal/test_cmp_aos_geom_noisefree_noreg/mt_0.bin');
mt1=read('/home/lianqiw/work/laos/runs/nominal/test_cmp_aos_geom_noisefree_noreg/mt_1.bin');
mt2=read('/home/lianqiw/work/laos/runs/nominal/test_cmp_aos_geom_noisefree_noreg/mt_2.bin');
for ips=1:nps
    disp(ips)
    for jps=1:nps
        for iwfs=1:6
            if iwfs==1
                RLM{ips, jps}=RRM{ips,iwfs}*G0{iwfs,jps};
            else
                RLM{ips, jps}=RLM{ips,jps}+RRM{ips,iwfs}*G0{iwfs,jps};
            end
        end
    end
    keyboard
    if L2
        RLM{ips,ips}=RLM{ips,ips}+L2{ips,ips}'*L2{ips,ips};
    if ZZT
        RLM{ips,ips}=RLM{ips,ips}+ZZT{ips,ips};
end
for ips=1:nps
    for jps=1:nps
        v=rand(size(RLM2{ips,jps},2),1);
        diff=norm(RLM{ips,jps}*v-RLM2{ips,jps}*v);
        fprintf('jps=%d,ips=%d, RLM diff=%g\n',jps,ips,diff);
    end 
end

for ips=1:nps
    for iwfs=1:9
        if skipwfs(iwfs)
            continue
        end
        if iwfs<7 
            RLU{ips,iwfs}=RRU{ips,iwfs};
            RLV{ips,iwfs}=G0{iwfs,ips}'*RRV{iwfs,iwfs};
        else
            RLU{ips,iwfs}=-RRM{ips,iwfs};
            RLV{ips,iwfs}=G0{iwfs,ips}';
        end
        diff1=norm(RLU{ips,iwfs}-RLU2{ips,iwfs});
        diff2=norm(RLV{ips,iwfs}-RLV2{ips,iwfs});
        fprintf('ips=%d,iwfs=%d,diff=%g, %g\n',ips,iwfs,diff1,diff2);
    end
end

RLMc=cell2mat(RLM);
X=rand(size(RLMc,2),1);
RLUc=cell2mat(RLU);
RLVc=cell2mat(RLV);
Y=RLMc*X-RLVc*(RLUc'*X);
norm(Y)
