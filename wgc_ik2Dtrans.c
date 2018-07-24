/**
# Gravito-Capillary Waves with a Moving Pressure Source in 3D

A similar test to the [capillary wave](capwave.c) but for a pure
gravity wave, using the [reduced gravity](/src/reduced.h) approach.

We use a constant-resolution grid, the Navier--Stokes solver for
two-phase flows and reduced gravity. */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "tension.h"
#include "iforce.h"
#include "curvature.h"

/**
We make sure that the boundary conditions for the face-centered
velocity field are consistent with the centered velocity field (this
affects the advection term). */

uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;
uf.n[front] = 0.;
uf.n[back] = 0.;

/* uf.n[left] =  dirichlet(0); */
/* uf.n[right]  =  dirichlet(0); */
/* //uf.n[top] = dirichlet(0); */
/* uf.n[bottom]  = dirichlet(0); */
/* uf.n[front] = dirichlet(0); */
/* uf.n[back]  =  dirichlet(0); */


/* uf.t[left]   = 0.; */
/* uf.t[right]  = 0.; */
/* uf.t[top]    = 0.; */
/* uf.t[bottom] = 0.; */

/**
We will store the accumulated error in *se* and the number of samples
in *ne*. */

int main() {
  /**
  The domain is 2x2 to minimise finite-size effects. The viscosity is
  constant. The acceleration of gravity is 50. */
  L0 = 50.;
  X0 = -L0/2;
  Y0 = -L0/2;
  Z0 = -L0/10; 
  G.z = 50;       //'G' has an attribute 'z' to specify direction
  f.sigma = 1;  //'f' interface has an attribute 'sigma' for surface tension
  rho1 = 1, rho2 = 0.1225;
  mu1 = 0.018; mu2 = 0.0001; 
  TOLERANCE = 1e-6;
  N = 64;
  init_grid(N);
  //for(f.sigma = 0; f.sigma <=  5; f.sigma = f.sigma + 0.1)
    // {  
  run();
  // }
}

event init (t = 0) {
  //fraction (f, y - 0.1*cos(2.*pi*x));
  //double ainit = 0.2; 
  /* double b = 1; */
  /* double c = 1; */
  fraction(f, z);// + 0.1*exp(-(x-2)*(x-2)*(y-2)*(y-2)/(2*ainit*ainit)));
}

/**
By default tracers are defined at $t-\Delta t/2$. We use the *first*
keyword to move VOF advection before the *amplitude* output i.e. at
$t+\Delta/2$. This improves the results. */

event vof (i++, first);
/**
We output the amplitude at times matching exactly those in the
reference file. */

event amplitude (t += 0.005; t <= 0.5) {
  scalar pos[];
  vector h[];

  position (f, pos);
  heights(f,h);

  char nameposfac[100];
  sprintf (nameposfac, "Kakade/wgc_2d/data2307/R3S2/datafacqcc-%g.txt",t);
  FILE * fac = fopen (nameposfac, "w");
  output_facets (f, fac);

  //double max = statsf(pos).max;
  output_ppm(pos,linear = true);

  /* char namepos[100]; */
  /* sprintf (namepos, "Kakade/wgc_2d/data1607/2D/10/dataallqcc-%g.txt",t); */
  /* FILE * fpos = fopen (namepos, "w"); */
  /* output_field(all, fpos, N); //Does not work in 3D */
  /* fflush (fpos); */
  //fprintf (stderr, h[]);

  /* /\** */
  /* We output the corresponding evolution in a file indexed with the */
  /* number of grid points N and the current numerical surface tension. * */
  /* char name[100]; */
  /* double st; */
  /* st = f.sigma; */
  /* sprintf (name, "Kakade/wgc_2d/data0207/wave-%d-%g.txt", N, st);    //CHANGE DATE BEFORE RUNNING AND CREATE FOLDER */
  /* static FILE * fp = fopen (name, "w"); */
  /* fprintf (fp, "%g %g\n", t*16.032448313657, max); */
  /* fflush (fp); */
  fprintf (stderr, "It. No.: %d\n",i);
}

event display (t = end)
{
  //double cap_length = sqrt(f.sigma/(rho1*G.y));
  fprintf (stderr, "Number of Grid Points per Wavelength: %g Numeric Surface Tension: %g\nNumeric Capillary Length:  Numeric Liquid Viscosity: %g\n",N/L0, f.sigma, mu1);
}


event acceleration (i++)
{
  scalar phi = f.phi;
  double a = 0.1;     //Width of Gaussian Distributed Force
  double Po = -150; //Magnitude of Force Applied
  double s = 2;     //Angular Velocity
  double b = 5;    //Radius of Rotation
  coord G1;
  G1.z = -G.z;      //Direction of Applied
 if (phi.i)
   positionbumptrans2D(f, phi, a, Po, b, t, s, G1, Z, add = true);
     else {
     phi = new scalar;
     positionbumptrans2D(f, phi, a, Po, b, t, s, G1, Z, add = false);
     f.phi = phi;
    }
}

event adapt(i++)
{
  adapt_wavelet ({f}, (double []){1e-4}, maxlevel = 10, minlevel = 6);
  
  scalar l[];
  foreach()
    l[] = level;
  static FILE * fp = fopen ("gridrefine12.ppm", "w");
  output_ppm (l, fp, min = 0, max = 8);
}

/* #if 1 */
/* event gfsview (i += 1) { */
/*   static FILE * fp = popen ("gfsview2D -s Kakade/wgc_2d/data1607/2D/1/gravity.gfv", "w"); */
/*   output_gfs (fp); */
/* } */
/* #endif */
