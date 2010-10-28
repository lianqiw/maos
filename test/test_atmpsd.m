%function atm=test_atmpsd(nx)
n2=64;
nx=16;

dx=1/64.;
dx2=dx*dx;
dk=1/(nx*dx);
dk2=1/(nx^2*dx2);
r0=0.15;
L0=30;
L0eff=L0;
L0eff=(L0.^-2+dk^2/2)^-0.5;%effective L0
coeff=0.0229*(2*pi/0.5e-6)^-2*r0^(-5/3.)*dk2;
[KY,KX]=meshgrid([-nx/2:nx/2-1]*dk);
K2=KX.*KX+KY.*KY;

L0isq=1/L0^2;
PSD=coeff*(K2+L0isq).^-(11/6.);
PSD(end/2+1,end/2+1)=0;
PSD(1,:)=0;
PSD(:,1)=0;
%COV=fftshift(ifft2(fftshift(PSD)));
%cov=COV((nx-n2)/2+1:(nx+n2)/2,(nx-n2)/2+1:(nx+n2)/2);
%cov=cov-cov(end/2+1,end/2+1);
%atm=real(fft2(fftshift(sqrt(PSD)).*(randn(nx,nx)+1i*randn(nx,nx))));

%integral of PSD is 0.0229k^-2*(L0/r0)^5/3*pi*6/5
PINT=0.0229*(2*pi/0.5e-6)^-2*(L0eff/r0)^(5/3)*pi*6/5
sum(PSD(:))
