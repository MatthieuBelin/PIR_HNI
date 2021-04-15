
#include "saint-venant.h"

#define LEVEL 6

event init (t=0) {
  // Initial conditions
	foreach()
		h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
}

event graphs (i++) {
  // To use with gnuplot
  // It shows the minimal and maximal height of the water for each iteration
	stats s = statsf (h);
	fprintf (stderr, "%g %g %g\n", t, s.min, s.max);
}

event images (t += 8./300.) {
  // To use with animate
  
  // It shows the height of water each 8./300. seconds
  // The linear option implies an linar interpolation between the meshes
	output_ppm (h, linear = true);

  // It allows us to see the grid
	scalar l[];
	foreach()
		l[] = level;
	static FILE * fp = fopen("grid.ppm", "w");
	output_ppm(l, fp, min = 0, max = LEVEL);
}

event end (t = 8) {
  // It ends the simulation
	printf ("i = %d t = %g\n", i, t);
}

event adapt (i++) {
  // Adapt the "quatree" mesh
	adapt_wavelet ({h}, (double []){4e-3}, maxlevel = LEVEL);
}

int main() {
  // The solution is calculated on a symetric field
	origin(-3., -3.);
	init_grid(1 << LEVEL);
  // Run the code
	run();
}
