DM Hysteresis {#hysteresis}
======

The DM hysteresis model is based on "Modeling the hysteresis of a scanning
probe microscope", J. Vac. Sci. Technol. B 18 (2), Mar/Apr 2000.

Formula (2) of the paper is used with \f$x\f$ replaced by the actual DM
position \f$y\f$ and \f$V\f$ replaced by the command \f$x\f$:
\f[
	\frac{dy}{dx}=\alpha\, \mathrm{sign}(dx) (\beta x-y)+(\beta-u)
\f]

For a closed hysteresis loop within the stroke limit \f$S\f$, the linear
ordinary differential equation can be solved analytically:
\f{eqnarray*}{
	y_{t/r}(x)=\beta x \mp \frac{u}{\alpha}\left(1-\frac{2\exp(\mp\alpha x)}{\exp(-\alpha S)+\exp(\alpha S)}\right)
\f}

where \f$t\f$ and \f$r\f$ represent the trace (increasing \f$x\f$, upper sign)
and retrace (decreasing \f$x\f$, lower sign) directions. It is desirable to
have \f$y\f$ within the same range as \f$x\f$ to simulate a calibrated
command: \f$y_{t/r}(S)=S\f$. Let the hysteresis ratio be \f$h\f$, defined as
the maximum difference between the trace and retrace curves, which occurs at
\f$x=0\f$, divided by \f$2S\f$. We have

\f{eqnarray*}{
	\frac{y_r(0)}{S}&=&\frac{u}{\alpha S}\left(1-\frac{2}{\exp(-\alpha S)+\exp(\alpha S)}\right)=h \\
	\frac{y_{t/r}(S)}{S}&=&\beta-\frac{u}{\alpha S}\left(1-\frac{2\exp(-\alpha S)}{\exp(-\alpha S)+\exp(\alpha S)}\right) = 1
\f}

The equations can be solved uniquely for \f$u\f$ and \f$\beta\f$ when
\f$\alpha S\f$ is fixed. The hysteresis curve is a weak function of
\f$\alpha S\f$, and \f$\alpha S=2\f$ produces a good looking curve. A smaller
\f$\alpha S\f$ has slower convergence of the hysteresis curve.

The hysteresis model is simulated using simple integration. For each step with
a new command \f$x\f$, the position \f$y\f$ is updated as
\f{eqnarray*}{
	dx &=& x-x_{last} \\
	dy &=& \alpha |dx| (\beta x_{last}-y_{last})+dx(\beta-u) \\
	y_{last}&\mathrel{+}=&dy \\
	x_{last}&=&x \\
	y&=&y_{last}
\f}

# Implementation {#hyst-impl}

The model is implemented in \c lib/hyst.c:

- hyst_new() creates the model for one DM. It takes the hysteresis ratio
  \f$h\f$, the dimensionless \f$\alpha S\f$ (a negative value uses the default
  of 2), the stroke \f$S\f$ (infinity is replaced by 10 µm), and the number of
  actuators, and computes \f$u\f$, \f$\beta\f$, and \f$\alpha\f$ from the
  equations above.
- hyst_dmat() and hyst_dcell() apply the integration above to each actuator
  and keep \f$x_{last}\f$ and \f$y_{last}\f$ as the state. hyst_reset()
  clears the state.

In \c maos, the model is enabled for a DM when \c dm.hyst is nonzero, using
the parameters

- \c dm.hyst: the hysteresis ratio \f$h\f$, between 0 and 1 (default 0, i.e.,
  disabled),
- \c dm.hyst_alpha: \f$\alpha S\f$ (default 2),
- \c dm.hyst_stroke: the stroke \f$S\f$ in meters at which the hysteresis is
  measured (default 10e-6),

with defaults in \c config/maos/sim.conf. The models are created in
\c maos/sim_utils.c, and filter_dm() (\c maos/filter.c) applies them to
convert the DM commands \c simu->dmcmd to the actual DM shape
\c simu->dmreal.
