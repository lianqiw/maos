
if 0
denom=read('denom');
ind=denom<0.5;
b=read('atm_cov_100'); 
b(ind)=0;
c=real(b);
cov=(c(end/2+1:end,end/2+1));
d=2*(cov(1)-cov);

r0=0.20;
r0i=1/r0;
lambda=0.5e-6;
N=2;
dx=1./64;
r=(0:N-1)'*dx;
D=6.88*(2*pi/lambda)^(-2)*(r*r0i).^(5./3);



dx=1/64.;
r=dx;
diam=dx*2;
r0=0.20;
r0i=1/r0;
lambda=0.5e-6;
sigma2=0.5*6.88*(2*pi/lambda)^(-2)*(diam*sqrt(2)*r0i).^(5./3);

c0=sigma2;
c1=sigma2-0.5*6.88*(2*pi/lambda)^(-2)*(r*r0i).^(5./3);
c2=sigma2-0.5*6.88*(2*pi/lambda)^(-2)*(r*sqrt(2)*r0i).^(5./3);

C=[c0 c1 c2 c1;
   c1 c0 c1 c2;
   c2 c1 c0 c1;
   c1 c2 c1 c0];

ind=[1 4 2 3];
C=C(ind,ind);
a=read('cov');
a
C
D0=6.88*(2*pi/lambda)^(-2)*(dx*r0i).^(5./3);

end
nframe=10000;

%Calculate structure function
istep=1;
atm=read(sprintf('atm_%d.bin',istep-1));
atm=atm(1:end-1,1:end-1);
N=size(atm,1);
nsep=N;
d=zeros(nframe,nsep);
dx=1/64;
amp=ones(size(atm));
r=[1:nsep]*dx;
D2=0;
for istep=1:nframe
    istep
    atm=read(sprintf('atm_%d.bin',istep-1));
    atm=atm(1:end-1,1:end-1);
    D2=D2+struct_function(atm, amp, dx, r);
end
D2=D2/nframe;
Dr=azimuthal_average(D2, 1:nsep);
save D2 D2;

r0=0.20;
r0i=1/r0;
lambda=0.5e-6;
N=nsep;
dx=1./64;
r=(1:N)'*dx;
D=6.88*(2*pi/lambda)^(-2)*(r*r0i).^(5./3);



%Calculate structure function
istep=1;
atm=read(sprintf('atmfft_%d.bin',istep-1));
N=size(atm,1);
nsep=N;
d=zeros(nframe,nsep);
dx=1/64;
amp=ones(size(atm));
r=[1:nsep]*dx;
D3=0;
for istep=1:nframe
    istep
    atm=read(sprintf('atmfft_%d.bin',istep-1));
    D3=D3+struct_function(atm, amp, dx, r);
end
D3=D3/nframe;
Dr3=azimuthal_average(D3, 1:nsep);
save D3 D3





lambda=0.5d-6;
r0=0.2;
r0i=1./r0;
D=16;
sigma2=6.88*(2*pi/lambda)^(-2)*(D*sqrt(2)*r0i).^(5./3)/2
c0=sigma2;
c1=sigma2-6.88*(2*pi/lambda)^(-2)*(D*r0i).^(5./3)/2
c2=sigma2-6.88*(2*pi/lambda)^(-2)*(sqrt(2)*D*r0i).^(5./3)/2

C=[c0 c1 c2 c1; c1 c0 c1 c2; c2 c1 c0 c1; c1 c2 c1 c0];
a=sqrt(c0+2*c1+c2);
b=sqrt(c0-2*c1+c2);
c=sqrt(2*(c0-c2));
K=0.5*[a -b -c 0; a b 0 -c; a -b c 0; a b 0 c];

nframe=10000;
cov=zeros(16,nframe);

for iframe=1:nframe
    v=randn(4,1);
    atm=K*v;
    tmp=atm*atm';
    cov(:,iframe)=tmp(:);
end
covm=zeros(nframe,16);
stdm=zeros(16,1);
for i=1:16
    covm(:,i)=cummean(cov(i,:));
    stdm(i)=std(cov(i,:));
end
K2=chol(C)';

cov2=zeros(16,nframe);

for iframe=1:nframe
    v=randn(4,1);
    atm=K2*v;
    tmp=atm*atm';
    cov2(:,iframe)=tmp(:);
end 

covm2=zeros(nframe,16);
stdm2=zeros(16,1);

for i=1:16
    covm2(:,i)=cummean(cov2(i,:));
    stdm2(i)=std(cov2(i,:));
end

[u,s,v]=svd(C);
K=u*sqrt(s);

cov3=zeros(16,nframe);

for iframe=1:nframe
    v=randn(4,1);
    atm=K2*v;
    tmp=atm*atm';
    cov3(:,iframe)=tmp(:);
end 

covm3=zeros(nframe,16);
stdm3=zeros(16,1);

for i=1:16
    covm3(:,i)=cummean(cov3(i,:));
    stdm3(i)=std(cov3(i,:));
end

C=read('vkcov_C');
K=read('vkcov_K');
nframe=5000;
cov4=zeros(81*81,nframe);

for iframe=1:nframe
    v=randn(81,1);
    atm=K*v;
    tmp=atm*atm';
    cov4(:,iframe)=tmp(:);
end 

covm4=zeros(nframe,81*81);
stdm4=zeros(81*81,1);

for i=1:81*81
    covm4(:,i)=cummean(cov4(i,:));
    stdm4(i)=std(cov4(i,:));
end

covmm4=mean(cov4,2);
C4=reshape(covmm4, 81, 81);