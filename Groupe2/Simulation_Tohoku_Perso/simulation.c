// qcc -autolink -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DDUMBGL solver.c -o simulation -l

#include "grid/quadtree.h"
#include "paraview2d.h"
#include "run.h"

#include "terrain.h"

//==============================Variables==============================

scalar h[]; // Height of the water
scalar zb[]; // Bathymetrie
vector u[]; // Speed in x and y direction

double g = 9.81;

double TMAX = 1.;

double dry = 1.e-10;
int level = 7;

double Radius = 6371220.; // Earth radius in metres

//==============================Main function==============================

int main(int argc, char * argv[]) {
    // resolution
    if (argc > 1) level = atoi(argv[1]);


    L0 = sqrt(73.); // the domain is 73 degrees squared
    // centered on 142,38 longitude,latitude
    X0 = 142. - L0/2.;
    Y0 = 38. - L0/2.;

    init_grid(1 << level);
    CFL = 0.4;
    run();
}

//==============================boundary conditions==============================
u.n[left] = neumann(0.); u.n[right] = neumann(0.);
u.n[bottom] = neumann(0.); u.n[top] = neumann(0.);

//==============================Functions for the solver==============================

void flux (double h, double u, double *flux1, double *flux2) {
	*flux1 = h*u;
	*flux2 = h*u*u + g*h*h/2.;
}

// Roe's state
void WRoe (double hL, double uL, double hR, double uR, // inputs
           double * hRoe, double * uRoe ) {
    // Roe's state
    *hRoe = (hL + hR)/2.;
    *uRoe = (uL*sqrt(hL) + uR*sqrt(hR))/(sqrt(hL) + sqrt(hR));
}

void HLL (double hL, double uL, double hR, double uR, double *F1, double *F2, double *vmax) {
    // Roe's state
    double hRoe, uRoe;
    WRoe (hL, uL, hR, uR, &hRoe, &uRoe);

	// wave speeds
	double sL = min(uL - sqrt(g*hL), uRoe - sqrt(g*hRoe));
	double sR = max(uR + sqrt(g*hR), uRoe + sqrt(g*hRoe));
	*vmax = max(fabs(sL), fabs(sR));

	// fluxes
	double f1L, f2L;
	flux(hL, uL, &f1L, &f2L);
	double f1R, f2R;
	flux(hR, uR, &f1R, &f2R);

	if (0. < sL) {
		*F1 = f1L;
		*F2 = f2L;
	}
	else if (sR < 0.) {
		*F1 = f1R;
		*F2 = f2R;
	}
	else {
	    *F1 = (sR*f1L - sL*f1R)/(sR-sL) + sL*sR*(hR - hL)/(sR-sL);
        *F2 = (sR*f2L - sL*f2R)/(sR-sL) + sL*sR*(hR*uR - hL*uL)/(sR-sL);
	}
}

//==============================Events to solve==============================

event init (i = 0)
{
    terrain(zb, "~/Documents/terrain/etopo2", "~/Documents/terrain/srtm_japan", NULL);

    if (restore (file = "restart"))
        conserve_elevation();
    else {
        conserve_elevation();

        foreach()
        h[] = max(0., - zb[]);
        boundary ({h});
        #include "tohoku/faults.h"
    }
}

event plot (t <= TMAX; t += 0.05) {
	fprintf (stdout, "# t = %g\n", t); // time
	// paraview
	scalar zs[];
	foreach() zs[] = h[] + zb[];
	output_paraview (slist = {h, zs, zb}, vlist = {u});
}


event adapt (i++) {
    astats s = adapt_wavelet ({h}, (double[]){1.e-3}, maxlevel=9, minlevel=7);
    fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}

event integration (i++) {
	// Compute the solution for this iteration

	// Neumann boundary condition
	boundary ({h, u, zb});

	// Fluxes at interfaces
	face vector Fh[];	// Fh.x,	Fh.y
	tensor Fhu[];   // FhuL.x.x, Fhu.x.y
                    // Fhu.y.x,  FhuR.y.y
    face vector Fhv[]; // Fhv = FhuR.x.x, FhuR.y.y

	// Fluxes
	double dtmax = 1.E+10;
	foreach_face(reduction (min:dtmax)) {
		// Face in x AND y directions
		double hL = h[-1], zL = zb[-1], uL = u.x[-1], vL = u.y[-1]; // left value
		double hR = h[], zR = zb[], uR = u.x[], vR = u.y[];   // right value

		// Audusse2004 --> M = Minus ; P = Plus
		double zi = max(zL, zR); // eq 2.9
		double hM = max(0., hL + zL - zi), uM = uL, vM = vL; // eq 2013
		double hP = max(0., hR + zR -zi), uP = uR, vP = vR;
		// source term
		double sM = g/2. * (sq(hM) - sq(hL)); // eq 2.14
		double sP = g/2. * (sq(hR) - sq(hP));
		if (hM > dry || hP > dry) {
			// F(W_L,W_R) = (F1, F2)
			double F1, F2, F3, vmax;
			// numerical flux
			HLL (hM, uM, hP, uP, &F1, &F2, &vmax);
			// get fluxes
            F3 = (F1 > 0. ? F1*vM : F1*vP);  // HLLC
			Fh.x[] = fm.x[] * F1;
			Fhu.x.x[] = fm.x[] * (F2 - sM); // FhuL
            Fhu.y.x[] = fm.x[] * F3;
			Fhv.x[] = fm.x[] * (F2 + sP); // FhuR

			// Time step
			double dxloc = Delta*cm[]/fm.x[];
			dtmax = min(dtmax, CFL*dxloc/vmax);
		} else {
			Fh.x[] = Fhu.x.x[] = Fhu.y.x[] = Fhv.x[] = 0.;
		}
	}

	boundary_flux ({Fh, Fhu, Fhv});

	// next timestep
	double dt = dtnext(dtmax);

	// Update
	foreach() {
	    double dxloc = cm[]*Delta; // size of the cell
		// Fh.x,    Fh.y
		// Fhu.x.x, Fhu.x.y
		// Fhu.y.x, Fhu.y.y

		// update h
		double h_old = h[]; // backup h at t = t^n 
		// note that: FL = F[], FR = F[1]
		double FhL = Fh.x[]  + Fh.y[]; // left + bottom
		double FhR = Fh.x[1,0] + Fh.y[0,1]; // right + top 
		h[] = h_old - (dt/dxloc)*(FhR - FhL); // h at t = t^{n+1}

		// update u,v
		if (h[] > dry) {
			foreach_dimension() { // u, v
				double FhuL = Fhu.x.y[]  + Fhv.x[]; // left + bottom
				double FhuR = Fhu.x.x[1,0] + Fhu.x.y[0,1]; // right + top
				double hu = h_old*u.x[] - (dt/dxloc)*(FhuR - FhuL);
				u.x[] = hu/h[];
			}
		} else { // dry ==> set u = v = 0
			u.x[] = u.y[] = 0.;
		}
	}
}

/* FOR MESH ADAPTATION (see http://basilisk.fr/src/saint-venant.h) */
event defaults (i=0) {
    // default initial condition
    foreach () {
        for (scalar s in {h,u,zb}) {
            s[] = 0.; // initial condition
        }
    }

    // for refine/coarsen
    for (scalar s in {h,u,zb}) {
        // see http://basilisk.fr/src/grid/multigrid-common.h
        s.refine = s.prolongation = refine_linear; // refine
        s.restriction = restriction_volume_average; // coarsen
    }
}
/* END MESH ADAPTATION */
