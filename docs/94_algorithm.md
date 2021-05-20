
\page algorithm Algorithm Details
The following shows how various algorithms are implemented.

<!--
pandoc -t pdf -V geometry:margin=2cm -o ~/Sync/fresnel.pdf thisfile.md
-->

# Physical Optics Beam Propagation
## Maxwell Equation

The Maxwell equation predicts that 

$$\nabla^2 \bold{E} = \mu_0 \epsilon_0 \frac{\partial^2 \bold{E}}{\partial t^2} $$
where $\mu_0 \epsilon_0 = 1/c_0^2$

It has a general solution of the form

$$\bold{E}=\bold{E_0}f(\hat{\bold{k}}\cdot{\bold{x}} - c_0 t)$$ where
$\hat{\bold{k}}$ is the wave propagation direction. We have
$\bold{E}\cdot\hat{\bold{k}}=0$ and $c_0\bold{B}=\hat{\bold{k}}\times \bold{E}$. 

## Fresnel diffraction integral

The electric field at a point $(x,y,z)$ following an aperture at $z=0$ is given by

$$ E(x,y,z)=\frac{1}{i\lambda}\iint E(x^\prime, y^\prime, 0)\frac{e^{ikr}}{r} \frac{z}{r} dx^\prime dy^\prime $$

where $k=2\pi/\lambda$ is the wave number and
$r=\sqrt{(x-x^\prime)^2+(y-y^\prime)^2+z^2}$. This is the exact solution, but
lacks analytical or convenient numerical solution.

## Fresnel approximation

Introduce $\rho^2=(x-x^\prime)^2+(y-y^\prime)^2$, and substitute $r$ with its
Taylor expansion:
$$r=\sqrt{\rho^2+z^2}=z\sqrt{1+\frac{\rho^2}{z^2}}=z+\frac{\rho^2}{2z}-\frac{\rho^4}{8z^3}+...$$
The equation can be simplified if we can ignore the third term, which requires
that the phase $k\frac{\rho^4}{8z^3}\ll2\pi$ or
$$z^3\gg\frac{\rho^4}{8\lambda}.$$ For the $r$ outside of the phase, we make
another approximate by keeping only the first term,  which is valid if the
aperture is significantly smaller than the propagation distance: $\rho\ll z$.
The electric field can then be simplified as

$$E(x,y,z)=\frac{e^{ikz}}{i\lambda z}\iint E(x^\prime, y^\prime,0){e^{\frac{i\pi}{z\lambda}[(x-x^\prime)^2+(y-y^\prime)^2]}} dx^\prime dy^\prime$$ 

This is the Fresnel diffraction integral. It means that, if the Fresnel
approximation is valid, the propagating field is a spherical wave, originating
at the aperture and moving along z. The term before the integral represents
constant phase evolution and energy normalization. 

### Angular Spectrum

The integral can be viewed as a convolution between $E(x,y,0)$ and 
$$h(x,y)=e^{\frac{i\pi}{z\lambda}(x^2+y^2)}.$$ 

The convolution can be done with Fourier transform (there is no need for embedding): $$E=\mathcal{F}^{-1}[\mathcal{F}[E(x,y,0)]\mathcal{F}[h(x,y)]].$$ 

The fourier transform of $h(x,y)$ can be computed analytically. Recall  the Fourier transform $$\mathcal{F}[\exp[ia(x^2+y^2)]]=\frac{i\pi}{a}\exp[\frac{-i\pi^2}{a}(f_x^2+f_y^2)]$$ for $k>0$. We have 
$$\mathcal{F}[h(x,y)]=iz\lambda\exp[-iz\pi\lambda(f_x^2+f_y^2)].$$ Using this requires upsampling since FFT convolution assumes the function is periodic. We do not use this method.

This convolution view is commonly referred to as the angular spectrum method.
The original field $E(x,y,0)$ is only FFT'ed once. For each propagation
distance, the FFT of $h()$ is multiplied to $\hat{E}$ and then inverse FFT'ed to
obtain the field. This approach is mostly useful in very near field propagations
when $x$ and $x^\prime$ has the same sampling. 


### Single FFT

The integral can also be rewritten as a single Fourier transform: 

$$E(x,y,z)=\frac{e^{ikz}}{i\lambda z} h(x,y)\iint E(x^\prime, y^\prime, 0){e^{\frac{i\pi}{z\lambda}[{x^\prime}^2+{y^\prime}^2]}}{e^{-i2\pi(\frac{x}{z\lambda}x^\prime+\frac{y}{z\lambda}y^\prime)}} dx^\prime dy^\prime$$

Define a new function $G(x^\prime, y^\prime)$ as 
$$G(x^\prime, y^\prime)=E(x^\prime, y^\prime, 0){e^{\frac{i\pi}{z\lambda}({x^\prime}^2+{y^\prime}^2)}}$$
Its fourier transform is
$$\hat{G}(p,q)=\mathcal{F}\{G(x^\prime,y^\prime)\}.$$ We have $p=\frac{x}{z\lambda}$ and $q=\frac{y}{z\lambda}$ and 
$$E(x,y,z)=\frac{e^{ikz}}{i\lambda z} h(x,y) \hat{G}(\frac{x}{z\lambda},\frac{y}{z\lambda}).$$ 

The FFT view is more suitable when the propagation is near far field as the FFT
of $G()$ gives coordinate in the angular space. Replace $x,y$ by the angular
coordinate $p,q$ and renormalize the function to angular coordinate, we have 

$$ E(p,q,z)=h(p z \lambda, q z \lambda)\hat{G}(p,q)$$. The scaling of
$\frac{1}{z\lambda}$ is dropped due to change of variable and preservation of
total energy. In numerical calculation, the sampling of $p,q$ and $x,y$ follows
$dp=1/(n*dx),dq=1/(n*dy)$ where $n$ is the array dimension along $x$ and $y$.

## Fraunhofer approximation

In a special case when the propagation is at far field (from pupil to focus or
vice versa), so that $z\gg\frac{\rho^2}{2\lambda}$ the extra phase term also be dropped:

$$E(x,y,z)=\frac{1}{i\lambda z}e^{ikz}\iint E(x^\prime, y^\prime,
0){e^{-i2\pi(\frac{x}{z\lambda}x^\prime+\frac{y}{z\lambda}y^\prime)}} dx^\prime
dy^\prime$$

Replace $x,y$ by angular coordinate $p=\frac{x}{z\lambda},q=\frac{y}{z\lambda}$, ignore the constant phase, and renormalize the function to angular coordinate, we have
$$E^\prime(p,q)=\hat{E}(p,q)$$ and the energy is preserved via
$$\iint{||E^\prime(p,q)||^2 dp dq}=\iint{||E(x,y,0)||^2dxdy}.$$

## Sphere to sphere propagation

When the wavefront at the starting aperture is highly curved, regular numerical
solution may not be able to sample the wavefront sufficiently, resuling in poor
numerical accuracy. 
