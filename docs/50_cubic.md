Algorithms {#algorithm}
======

This section describes some of the models and numerical algorithms used in
MAOS.

- \ref sect-dm-actuator (below)
- \subpage hysteresis
- \subpage pop

# DM Actuator Influence Function {#sect-dm-actuator}
## Linear influence function

The linear influence function is defined as
\f[
h(x;x_{i};\delta)=\begin{cases}
1-|x-x_{i}|/\delta & \textrm{if }|x-x_{i}|\leq\delta\\
0 & \textrm{else}
\end{cases}
\f]
where \f$\delta\f$ is usually equal to the grid spacing. A bilinear influence
function is simply the product of two linear influence functions with
variables \f$x_{1}\f$ and \f$x_{2}\f$.

## Cubic influence function

A cubic influence function that can reproduce piston/tip/tilt was proposed
by Ellerbroek to model the piezo stack DM actuator. The influence function
has the form
\f[
h(x;x_{i};\delta)=h_{0}((x-x_{i})/\delta)
\f]
where \f$\delta\f$ is the same as the grid spacing, and \f$h_{0}\f$ is the
influence function defined in the normalized coordinates
\f[
h_{0}(x)=\frac{1}{1+2c}\begin{cases}
1+(4c-\frac{5}{2})|x|^{2}+(\frac{3}{2}-3c)|x|^{3} & |x|\leq1\\
(2c-\frac{1}{2})(2-|x|)^{2}+(\frac{1}{2}-c)(2-|x|)^{3} & 1<|x|\leq2\\
0 & |x|>2
\end{cases}
\f]
where \f$c\f$ is the inter-actuator coupling coefficient: the ratio
\f$h_0(1)/h_0(0)\f$ of the response at the nearest neighbor actuator to that
at the actuator itself. The leading coefficient normalizes the influence
function so that the sum over the actuator grid points,
\f$h_0(0)+2h_0(1)\f$, equals 1. The 2-d influence function is the product of
the influence functions along \f$x\f$ and \f$y\f$.

## Implementation

The influence function is selected per grid by the \c iac member of #loc_t
and #map_t: 0 uses the bilinear influence function, and a positive value uses
the cubic influence function with \f$c=\f$\c iac.

- Ray tracing from a grid or point set with \c iac set uses the cubic
  interpolation in \c lib/accphi.c (\c prop_grid_cubic(),
  \c prop_nongrid_cubic(), and related functions, called by prop()). Transpose ray
  tracing of a #map_t is not available with the cubic influence function.
- The sparse ray tracing operators built by mkh() and mkht() (\c lib/mkh.c)
  use mkht_cubic() automatically when the input #loc_t has \c iac set.
- The coefficients \f$c_0\f$ to \f$c_4\f$ of the two polynomial pieces are
  computed in the macro \c PREP_CUBIC_PARAM (\c lib/accphi.c) and in
  mkht_cubic().

The coupling is set by the following parameters:

- \c dm.iac: the DM actuators (default 0.223 in \c sim.conf; 0 uses the
  bilinear influence function). It is assigned to \c recon->aloc.
- \c moao.iac: the MOAO DM actuators.
- \c tomo.iac: the tomography grid \c xloc (default 0, for testing).
- \c dbg.wfs_iac: the grid used to fit turbulence onto the WFS subapertures
  (default 0). It is a real number like the others.
