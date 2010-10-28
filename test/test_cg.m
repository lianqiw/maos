
%addpath('../mex/');
RR.M=cell2mat(read('RRM.bin.gz'));
RR.U=cell2mat(read('RRU.bin.gz'));
RR.V=cell2mat(read('RRV.bin.gz'));

RL.M=cell2mat(read('RLM.bin'));
RL.U=cell2mat(read('RLU.bin.gz'));
RL.V=cell2mat(read('RLV.bin.gz'));

gradcl=read('gradcl_1.bin.gz');
grad=cell2mat(gradcl(1:6));

rhs=RR.M*grad-RR.U*(RR.V'*grad);
%rhs2=cell2mat(read('rhs.dblc'));
%norm(rhs-rhs2)/norm(rhs)
disp('CG');
tol=0;
maxit=30;
[x,flag,resres,iter,resvec]=pcg(@MUV,rhs,tol,maxit,[],[],zeros(size(rhs)),RL);
resvec
keyboard
return
x0=zeros(size(rhs));
r0=rhs-(RL.M*x0-RL.U*(RL.V'*x0));
p0=r0;
rn0=r0'*r0;
for k=1:300
    ap=(RL.M*p0-RL.U*(RL.V'*p0));
    ak=rn0/(p0'*ap);
    x0=x0+ak*p0;
    r0=r0-ak*ap;
    rn1=r0'*r0;
    disp(rn1);
    bk=rn1/rn0;
    rn0=rn1;
    p0=r0+bk*p0;
end

