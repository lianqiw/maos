Algorithms {#algorithm}
======

# DM Actuator Influence Function {#sect-dm-actuator}
##	Linear influence function

The linear influence function is defined as
\f[
h(x;x_{i};\delta)=\begin{cases}
1-|x-x_{i}|/\sigma & \textrm{if }|x-x_{i}|\leq\delta\\
0 & \textrm{else}
\end{cases}
\f]
where \f$\delta\f$ is usually equal to the grid spacing. A bilinear influence
function is simply the product of two linear influence functions with
variables \f$x_{1}\f$and \f$x_{2}\f$.

## Cubic influence function

A cubic influence function that can reproduce piston/tip/tilt is coined
by Ellerbroek to model the piezo stack DM actuator. The influence
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
where c is the nearest neighbor coupling frequency. The leading coefficient
is to normalize the influence function so that it sums to 1.

