/*
  Copyright 2009-2026 Lianqi Wang
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <unistd.h>
#include "bessel_mode.h"


typedef struct {
    int m;
    double alpha;//spatial frequency
    double scale;//normalization constant for unit energy
} bmode_t;

// Target function: returns Jm(x) or Jm'(x)
double target_func(int m, double x, int type) {
    if (type == 0) return jn(m, x); // Dirichlet
    if (x < 1e-9) return (m == 1) ? -0.5 : 0.0; // Neumann
    return 0.5 * (jn(m - 1, x) - jn(m + 1, x));
}

double find_next_root(int m, double start, double k_cutoff, int type) {
    // Neumann (type 1) has a root at exactly 0 for m=0
    if (type == 1 && m == 0 && start < 0) return 0.0;
    
    // Dirichlet (type 0) or m > 0 Neumann: roots are > 0
    double a = (start < 0) ? 0.1 : start + 1e-5;
    double b = a + 0.3; 
    double fa = target_func(m, a, type);

    while (fa * target_func(m, b, type) > 0) {
        a = b;
        fa = target_func(m, a, type);
        b += 0.3;
        if (b > k_cutoff + 5.0) return 1e9;
    }

    for (int i = 0; i < 60; i++) {
        double mid = (a + b) / 2.0;
        if (fa * target_func(m, mid, type) < 0) b = mid; else a = mid;
    }
    return (a + b) / 2.0;
}
int cmp_alpha(const void *a, const void *b) {
    double d = ((bmode_t*)a)->alpha - ((bmode_t*)b)->alpha;
    return (d > 0) - (d < 0);
}
/**
 * Computes the Disk Harmonic Modes, aka Fourier Bessel modes, for a given set of coordinates and a cutoff frequency.
 * The modes are ordered by their spatial frequency (alpha) and returned in a column-major format
 * Inputs:
 *   - loc: loc_t structure containing the coordinates of the points where modes are evaluated
 *   - k_cutoff: maximum spatial frequency (alpha) to include in the output modes
 *   - type: 0 for Dirichlet (Jm=0 at boundary), 1 for Neumann (Jm'=0 at boundary)
 * Outputs: OPD
 */

dmat* bessel_modes(loc_t *loc, double D, double k_cutoff, int type) {
    const int nloc = loc->nloc;
    //double **out_phi, int *out_n_modes
    // 1. Pre-calculate polar coordinates
    double *r = (double*)malloc(nloc * sizeof(double));
    double *theta = (double*)malloc(nloc * sizeof(double));
	if(D==0){
		D=loc_diam(loc);
	}
    for (int i = 0; i < nloc; i++) {
        r[i] = sqrt(loc->locx[i] * loc->locx[i] + loc->locy[i] * loc->locy[i])*2./D;
        theta[i] = atan2(loc->locy[i], loc->locx[i]);
    }

    // 2. Build mode list and calculate normalization constants
	int mode_max=k_cutoff*k_cutoff;
    bmode_t *modes = (bmode_t*)malloc(mode_max * sizeof(bmode_t));
    int mode_count = 0;
    int total_cols = 0;

    for (int m = 0; m <= (int)(k_cutoff * 2); m++) {
        int found_in_m = 0;
		double alpha = -1;
        while (1) {
            alpha = find_next_root(m, alpha, k_cutoff, type);
            if (alpha > k_cutoff) break;
			if(mode_count==mode_max){
				mode_max*=2;
				bmode_t* modes_new=(bmode_t*)realloc(modes, mode_max*sizeof(bmode_t));
				if(modes_new) modes=modes_new; else error("Memory allocation for bmode_t failed\n");
			}
            modes[mode_count].m = m;
            modes[mode_count].alpha = alpha;
            
            // Normalization integral: 0.5 * (1 - m^2/alpha^2) * Jm(alpha)^2

			double J_alpha = jn(m, alpha);
			double r_int;
			if (type == 0) { // Dirichlet: 0.5 * J_{m+1}(alpha)^2
				r_int = 0.5 * pow(jn(m + 1, alpha), 2);
			} else { // Neumann: 0.5 * (1 - m^2/alpha^2) * Jm(alpha)^2
				r_int = (alpha < 1e-9) ? 0.5 : 0.5 * (1.0 - (double)m*m/(alpha*alpha))*J_alpha*J_alpha;
			}
			// Radial normalization + Angular normalization (1/pi or 1/2pi)
			double a_int = (m == 0) ? 2. : 1.;
			modes[mode_count].scale = sqrt(1./(r_int * a_int));

            mode_count++;
            total_cols += (m == 0) ? 1 : 2;
            found_in_m = 1;
        }
        if (!found_in_m && m > 0) break;
    }
	qsort(modes, mode_count, sizeof(bmode_t), cmp_alpha);
    // 3. Evaluate modes in Column-Major format: Phi[col * nloc + row]
	dmat *phi=dnew(nloc, total_cols);
    int current_col = 0;

    for (int i = 0; i < mode_count; i++) {
        int m = modes[i].m;
        double alpha = modes[i].alpha;
        double scale = modes[i].scale;
        if (m == 0) {
            for (int j = 0; j < nloc; j++) {
                P(phi, j, current_col) = scale * jn(0, alpha * r[j]);
            }
            current_col++;
        } else {
            // Cosine column
            for (int j = 0; j < nloc; j++) {
				double jnj = scale*jn(m, alpha * r[j]);
                P(phi, j, current_col  ) = jnj * cos(m * theta[j]);
				P(phi, j, current_col+1) = jnj * sin(m * theta[j]);
            }
            current_col+=2;
        }
    }

    free(r); free(theta); free(modes);
	return phi;
}
