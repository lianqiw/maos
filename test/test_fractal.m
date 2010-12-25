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