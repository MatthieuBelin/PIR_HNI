//#include "grid/cartesian.h"  // cartesian 2d
#include "grid/quadtree.h"
#include "paraview2d.h"
#include "run.h"


// variables
scalar h[];  // water height
vector u[];  // velocity (u.x, u.y)
scalar zb[]; // topography

// parameters
double TMAX = 1.; // final time
double G = 9.81; // gravity
double dry = 1.e-10;
int level = 7; // resolution 2^8 = 256

// boundary conditions
//h[left] = neumann(0.); h[right] = neumann(0.);
u.n[left] = neumann(0.); u.n[right] = neumann(0.);
u.n[bottom] = neumann(0.); u.n[top] = neumann(0.);

int main (int argc, char * argv[]) {
  // resolution
  if (argc > 1) level = atoi(argv[1]);

  // domain
  //G = 1.; 
  //L0 = 3.; X0 = -L0/2.; Y0 = -L0/2.;
  L0 = 4.; X0 = 0.; Y0 = 0.;

  init_grid (1 << level);
  CFL = 0.4;
  run();
}

//bid obstacle; // obstacle

event init (i = 0) {

  /* dam-break */
  /* foreach() { // for each cell x_i */
  /*   double dist = sqrt(sq(x) + sq(y)); */
  /*   h[] = dist < 1. ? 4. : 1.; */
  /*   u.x[] = u.y[] = 0.; */
  /* } */
  
  /* Thacker without friction */
  const double h0 = 0.1, a = 1., r0 = 0.8;
  double A = (a*a - r0*r0)/(a*a + r0*r0);
  double omega = sqrt(8.*G*h0)/a;
  TMAX = 3.*2.*M_PI/omega;
  foreach() {
    double r2 = (x-L0/2.)*(x-L0/2.) + (y-L0/2.)*(y-L0/2.);
    double eta = h0*(sqrt(1.-A*A)/(1.-A) - 1. - r2/(a*a)
		     * ((1.-A*A)/((1.-A)*(1.-A)) - 1.));
    zb[] = -h0 * (1.-(r2/(a*a)));
    h[] = max(0., eta - zb[]);
    u.x[] = u.y[] = 0.;
  }
  
}

event result (t <= TMAX; t += 0.1) {
  fprintf (stdout, "# t = %g\n", t); // time
  // paraview
  scalar eta[];
  foreach() {
    eta[] = h[] + zb[];
  }
  output_paraview (slist = {h,eta,zb}, vlist = {u});
}

event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1.e-3}, maxlevel=9, minlevel=7);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}

// F(W) = (f1,f2)
void flux (double h, double u, // inputs
	   double * f1, double * f2){
  *f1 = h*u;
  *f2 = h*u*u + G*h*h/2.;
}

// Roe's state
void WRoe (double hL, double uL, double hR, double uR, // inputs
	   double * hRoe, double * uRoe ) {
  // Roe's state
  *hRoe = (hL + hR)/2.;
  *uRoe = (uL*sqrt(hL) + uR*sqrt(hR))/(sqrt(hL) + sqrt(hR));
}

void HLL (double hL, double uL, double hR, double uR, // inputs
	  double * F1, double * F2, double * vmax) {
  // Roe's state
  double hRoe, uRoe;
  WRoe (hL, uL, hR, uR, &hRoe, &uRoe);
  
  // wave speeds
  double sL, sR;
  sL = min(uL - sqrt(G*hL), uRoe - sqrt(G*hRoe));
  sR = max(uR + sqrt(G*hR), uRoe + sqrt(G*hRoe));
  *vmax = max(fabs(sL), fabs(sR));
  
  // fluxes
  double f1L, f2L; // F(WL)
  flux(hL, uL, &f1L, &f2L);
  double f1R, f2R; // F(WR)
  flux(hR, uR, &f1R, &f2R);
  
  if (0. < sL) { // F(WL)
    *F1 = f1L;
    *F2 = f2L;
  } else if (sR < 0.) { // F(WR)
    *F1 = f1R;
    *F2 = f2R;
  } else { // FHLL
    *F1 = (sR*f1L - sL*f1R)/(sR-sL) + sL*sR*(hR - hL)/(sR-sL);
    *F2 = (sR*f2L - sL*f2R)/(sR-sL) + sL*sR*(hR*uR - hL*uL)/(sR-sL);
  }
}

event integration (i++) {
  // boundary condition 
  boundary({h, u, zb});

  // |----------|----------|
  //    (i-1)  [i]   (i)       
  
  // fluxes at interfaces 
  face vector Fh[]; // Fh.x,     Fh.y
  tensor Fhu[];     // FhuL.x.x, Fhu.x.y
                    // Fhu.y.x,  FhuR.y.y
  face vector S[];  // S = FhuR.x.x, FhuR.y.y
  double dtmax = 1.E+10;   
  // for each interface x_{i+1/2}, y_{j+1.2}
  foreach_face(reduction (min:dtmax)) { // face in x and y directions
    double hL = h[-1], uL = u.x[-1], vL = u.y[-1], zL = zb[-1]; // left value
    double hR = h[],   uR = u.x[],   vR = u.y[],   zR = zb[];   // right value
    // hydrostatic reconstruction - Audusse2004: 
    double zi = max(zL, zR); // eq 2.9 - M = minus, P = plus
    double hM = max(0., hL + zL - zi), uM = uL, vM = vL; // eq 2.13
    double hP = max(0., hR + zR - zi), uP = uR, vP = vR;
    // topographic source term 
    double sM = G/2.*(sq(hM) - sq(hL)); // eq 2.14
    double sP = G/2.*(sq(hR) - sq(hP));
    if (hM > dry || hP > dry) {
      // F(W_L,W_R) = (F1, F2)
      double F1, F2, F3, vmax;
      // homogeneous numerical flux applied on (hM,uM) and (hP,uP)
      HLL (hM, uM, hP, uP, &F1, &F2, &vmax);
      F3 = (F1 > 0. ? F1*vM : F1*vP); // HLLC
      // get flux at interface
      Fh.x[]    = fm.x[]*F1;
      Fhu.x.x[] = fm.x[]*(F2 - sM); // FhuL - eq 2.16 
      S.x[]     = fm.x[]*(F2 + sP); // FhuR 
      Fhu.y.x[] = fm.x[]*F3;
      // time step
      double dxloc = Delta*cm[]/fm.x[]; // size of face
      dtmax = min(dtmax, CFL*dxloc/vmax);
    } else {
      Fh.x[] = Fhu.x.x[] = Fhu.y.x[] = S.x[] = 0.;
    }
  }
  boundary_flux ({Fh, Fhu, S});
    
  // next timestep
  double dt = dtnext(dtmax);
  
  // update
  foreach() {
    // note that: FL = F[], FR = F[1]
    double dxloc = cm[]*Delta; // size of cell

    // update h
    double h_old = h[]; // backup h at t = t^n
    double FhL = Fh.x[] + Fh.y[]; // entering fluxes: left + bottom
    double FhR = Fh.x[1,0] + Fh.y[0,1]; // leaving fluxes: right + top
    h[] = h_old - (dt/dxloc)*(FhR - FhL);
    
    // update u,v
    if (h[] > dry) { // wet ==> update by finite volume scheme
      foreach_dimension() { // recall: FhuL = Fhu.x.x, S = FhuR.x.x
	double FhuL = S.x[]  + Fhu.x.y[]; // entering fluxes: left + bottom
	double FhuR = Fhu.x.x[1,0] + Fhu.x.y[0,1]; // leaving fluxes: right + top
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
