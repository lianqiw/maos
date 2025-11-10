Algorithms {#algorithm}
======

# DM Hysteresis {#hysteresis}

The DM hysteresis modeling is based on "Modeling the hysteresis of a scanning probe microscope"
J. Vac. Sci. Technol. B 18 (2), Mar/Apr 2000

Formula (2) is used with x replaced by actual DM position y and V replaced by command x, 
\f[
	\frac{dy}{dx}=\alpha sign(dx) (\beta x-y)+(\beta-u)
\f]


For closed hyteresis loop within the stroke limit \f$S\f$, the linear ordinary differential equation can be solved analytically:
\f{eqnarray*}{
	y_{t/r}(x)=\beta x \mp \frac{u}{\alpha}\left(1-\frac{2\exp(\mp\alpha x)}{\exp(-\alpha S)+\exp(\alpha S)}\right)
\f}

where t and r represents trace and retrace directions. It is desirable to have y be
within the same range as x to simulate calibrated command:

\f$y_{t/r}(S)=S\f$. Let the hysteresis ratio be \f$h\f$, we have

\f{eqnarray*}{
	\frac{y_r(0)}{S}&=&\frac{u}{\alpha S}\left(1-\frac{2}{\exp(-\alpha S)+exp(\alpha S)}\right)=h \\
	\frac{y_{t/r}(S)}{S}&=&\beta-\frac{u}{\alpha S}\left(1-\frac{2\exp(-\alpha S)}{\exp(-\alpha S)+\exp(\alpha S)}\right) = 1
\f}

The equation can be solved uniquely when \f$\alpha S\f$ is fixed.
The hysteresis curve is a weak function of \f$\alpha S\f$ and having \f$\alpha S=2\f$ produces a good looking curve.
Smaller \f$\alpha S\f$ has slower convergence of the hysteresis curve.

The hysteresis model is simulated using simple integration. For each step with a new x, the position y is updated as
\f{eqnarray*}{
	dx &=& x-x_{last} \\
	dy &=& \alpha*abs(dx)*(\beta*x_{last}-y_{last})+dx*(\beta-u) \\
	y_{last}&\mathrel{+}=&dy \\
	x_{last}&=&x \\
	y&=&y_{last}
\f}

