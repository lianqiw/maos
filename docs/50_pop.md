Physical Optics Beam Propagation {#pop}
======

<!--
Use the following to convert to pdf
sed 's/\\f\[/$$/g' 50_pop.md |sed 's/\\f\]/$$/g' | sed 's/\\f\$/$/g' | pandoc -t pdf -V geometry:margin=2cm -o ~/Downloads/fresnel.pdf
-->

# Maxwell Equation {#pop-maxwell}

The Maxwell equation predicts that 

\f[
	\nabla^2 \mathbf{E} = \mu_0 \epsilon_0 \frac{\partial^2 \mathbf{E}}{\partial t^2} \f]

where \f${\mu}_0{\epsilon}_0=1/{c}_0^2\f$

It has a general solution of the form

\f[\mathbf{E}=\mathbf{E_0}f(\mathbf{\hat{k}}\cdot{\mathbf{x}} - c_0 t)\f] where
\f$\mathbf{\hat{k}}\f$ is the wave propagation direction. We have
\f$\mathbf{E}\cdot\mathbf{\hat{k}}=0\f$ and \f$c_0\mathbf{B}=\mathbf{\hat{k}}\times \mathbf{E}\f$. 

# Fresnel diffraction integral {#pop-fresnel-integral}

The electric field at a point \f$(x,y,z)\f$ following an aperture at \f$z=0\f$ is given by

\f[ E(x,y,z)=\frac{1}{i\lambda}\iint E(x^\prime, y^\prime, 0)\frac{e^{ikr}}{r} \frac{z}{r} dx^\prime dy^\prime \f]

where \f$k=2\pi/\lambda\f$ is the wave number and
\f$r=\sqrt{(x-x^\prime)^2+(y-y^\prime)^2+z^2}\f$. This is the exact solution, but
lacks analytical or convenient numerical solution.

# Fresnel approximation {#pop-fresnel}

Introduce \f$\rho^2=(x-x^\prime)^2+(y-y^\prime)^2\f$, and substitute \f$r\f$ with its
Taylor expansion:
\f[r=\sqrt{\rho^2+z^2}=z\sqrt{1+\frac{\rho^2}{z^2}}=z+\frac{\rho^2}{2z}-\frac{\rho^4}{8z^3}+...\f]
The equation can be simplified if we can ignore the third term, which requires
that the phase \f$k\frac{\rho^4}{8z^3}\ll2\pi\f$, or
\f[z^3\gg\frac{\rho^4}{8\lambda}.\f] For the \f$r\f$ in the denominator (not within the exponential), we make
another approximation by keeping only the first term, which is valid if the
aperture is significantly smaller than the propagation distance: \f$\rho\ll z\f$.
The electric field can then be simplified as

\f[E(x,y,z)=\frac{e^{ikz}}{i\lambda z}\iint E(x^\prime, y^\prime,0){e^{\frac{i\pi}{z\lambda}[(x-x^\prime)^2+(y-y^\prime)^2]}} dx^\prime dy^\prime\f] 

This is the Fresnel diffraction integral. It means that, if the Fresnel
approximation is valid, the field is a superposition of spherical waves
(approximated by parabolic waves), originating at the aperture and moving
along z. The term before the integral represents
constant phase evolution and energy normalization. 

## Angular Spectrum {#pop-angular}

The integral can be viewed as a convolution between \f$E(x,y,0)\f$ and 
\f[h(x,y)=e^{\frac{i\pi}{z\lambda}(x^2+y^2)}.\f] 

The convolution can be done with Fourier transforms: \f[E=\mathscr{F}^{-1}[\mathscr{F}[E(x,y,0)]\mathscr{F}[h(x,y)]].\f] 

The Fourier transform of \f$h(x,y)\f$ can be computed analytically. Recall the Fourier transform \f[\mathscr{F}[\exp[ia(x^2+y^2)]]=\frac{i\pi}{a}\exp[\frac{-i\pi^2}{a}(f_x^2+f_y^2)]\f] for \f$a>0\f$. With \f$a=\frac{\pi}{z\lambda}\f$, we have
\f[\mathscr{F}[h(x,y)]=iz\lambda\exp[-iz\pi\lambda(f_x^2+f_y^2)].\f] The transfer function is applied directly in the frequency domain, instead
of sampling \f$h(x,y)\f$ in the spatial domain, which would be aliased. Since
the FFT treats the field as periodic, the output is periodic and the field
should be padded (embedded in a larger array) to avoid wrapping around.

This convolution view is commonly referred to as the angular spectrum method.
The original field \f$E(x,y,0)\f$ is only FFT'ed once. For each propagation
distance, the FFT of \f$h()\f$ is multiplied to \f$\hat{E}\f$ and then inverse FFT'ed to
obtain the field. The output sampling equals the input sampling. This approach is mostly useful
in near field propagation, where the output sampling of the single FFT method
below, \f$z\lambda/(n\,dx)\f$, would be finer than the input sampling.


## Single FFT {#pop-fft}

The integral can also be rewritten as a single Fourier transform: 

\f[E(x,y,z)=\frac{e^{ikz}}{i\lambda z} h(x,y)\iint E(x^\prime, y^\prime, 0){e^{\frac{i\pi}{z\lambda}[{x^\prime}^2+{y^\prime}^2]}}{e^{-i2\pi(\frac{x}{z\lambda}x^\prime+\frac{y}{z\lambda}y^\prime)}} dx^\prime dy^\prime\f]

Define a new function \f$G(x^\prime, y^\prime)\f$ as 
\f[G(x^\prime, y^\prime)=E(x^\prime, y^\prime, 0){e^{\frac{i\pi}{z\lambda}({x^\prime}^2+{y^\prime}^2)}}\f]
Its fourier transform is
\f[\hat{G}(p,q)=\mathscr{F}\{G(x^\prime,y^\prime)\}.\f] We have \f$p=\frac{x}{z\lambda}\f$ and \f$q=\frac{y}{z\lambda}\f$ and 
\f[E(x,y,z)=\frac{e^{ikz}}{i\lambda z} h(x,y) \hat{G}(\frac{x}{z\lambda},\frac{y}{z\lambda}).\f] 

The FFT view is more suitable when the propagation is near far field as the FFT
of \f$G()\f$ gives coordinate in the angular space. Replace \f$x,y\f$ by the angular
coordinate \f$p,q\f$ and renormalize the function to angular coordinate, we have 

\f[ E(p,q,z)=h(p z \lambda, q z \lambda)\hat{G}(p,q)\f]. The scaling of
\f$\frac{1}{z\lambda}\f$ is dropped due to change of variable and preservation of
total energy. In numerical calculation, the sampling of \f$p,q\f$ and \f$x^\prime,y^\prime\f$ follows
\f$dp=1/(n\,dx^\prime),dq=1/(n\,dy^\prime)\f$ where \f$n\f$ is the array dimension along \f$x\f$ and \f$y\f$.
The output sampling in the spatial domain is therefore \f$z\lambda/(n\,dx^\prime)\f$.

# Fraunhofer approximation {#pop-fraunhofer}

In a special case when the propagation is at far field (from pupil to focus or
vice versa), so that \f$z\gg\frac{\rho^2}{2\lambda}\f$, the quadratic phase term inside the integral can also be dropped:

\f[E(x,y,z)=\frac{1}{i\lambda z}e^{ikz}\iint E(x^\prime, y^\prime,
0){e^{-i2\pi(\frac{x}{z\lambda}x^\prime+\frac{y}{z\lambda}y^\prime)}} dx^\prime
dy^\prime\f]

Replace \f$x,y\f$ by angular coordinate \f$p=\frac{x}{z\lambda},q=\frac{y}{z\lambda}\f$, ignore the constant phase, and renormalize the function to angular coordinate, we have
\f[E^\prime(p,q)=\hat{E}(p,q)\f] and the energy is preserved via
\f[\iint{||E^\prime(p,q)||^2 dp dq}=\iint{||E(x,y,0)||^2dxdy}.\f]

# Sphere to sphere propagation {#pop-sphere}

When the wavefront at the starting aperture is highly curved, regular numerical
solution may not be able to sample the wavefront sufficiently, resulting in poor
numerical accuracy. This case is not handled by fresnel_prop(). 

# Implementation {#pop-impl}

fresnel_prop() in \c lib/locfft.c implements the Fresnel propagation of a
complex field (\c in, sampled at \c dxin) by a distance \c z at wavelength
\c wvl. The argument \c method selects the algorithm:

- 0: direct (brute force) summation of the Fresnel diffraction integral
  without the Fresnel approximation. It is very slow and mostly used for
  verification.
- 1: the single FFT method. The quadratic phase \f$h\f$ is applied before and
  after the FFT. The output sampling is \f$z\lambda/(n\,dx)\f$.
- 2: the angular spectrum method. The transfer function
  \f$\exp[-i\pi\lambda z(f_x^2+f_y^2)]\f$ is applied in the frequency domain.
  The output sampling equals the input sampling.
- -1: automatic. The angular spectrum method is used if the output sampling of
  the single FFT method would be smaller than the input sampling, otherwise
  the single FFT method.

For methods 1 and 2, if the Fresnel condition
\f$z^3\gg\frac{\rho^4}{8\lambda}\f$ is not satisfied by a large margin (with
\f$\rho\f$ taken as the array width \f$n\,dx\f$), the propagation is split
into two steps of \f$z+z_{min}\f$ and \f$-z_{min}\f$ that each satisfy it.
The output is multiplied by \c scale; a value of 0 or less computes it from the
input. fresnel_prop() is also available in Python through
\c scripts/libaos.py.

The far field (Fraunhofer) propagation from the pupil to the focal plane is
done with FFTs, for example by locfft_psf() (\c lib/locfft.c) to compute the
science PSFs.
