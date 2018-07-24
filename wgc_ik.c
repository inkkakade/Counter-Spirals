/**
# Gravito-Capillary Waves with a Moving Pressure Source

A similar test to the [capillary wave](capwave.c) but for a pure
gravity wave, using the [reduced gravity](/src/reduced.h) approach.

We use a constant-resolution grid, the Navier--Stokes solver for
two-phase flows and reduced gravity. */

#include "grid/multigrid.h"
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
/* uf.n[front]    = 0.; */
/* uf.n[back] = 0.; */
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
  
  L0 = 20.;
  Y0 = -L0/2.;
  G.y = 50.; //'G' has an attribute 'y' to specify direction
  //f.sigma = 1; //'f' interface has an attribute 'sigma' for surface tension
  rho1 = 1, rho2 = 0.1225;
  mu1 = 30; mu2 = 0.0001; 
  TOLERANCE = 1e-6;
  N = 512;
  // for(G.y = 10; G.y <= 50; G.y = G.y +10)
  // {
  for(f.sigma = 5.1; f.sigma <=  6.5; f.sigma = f.sigma + 0.1)
      {  
       run();
      }
	  //  }
  
}

/**
The initial condition is a small amplitude plane wave of wavelength
unity. */

event init (t = 0) {
  //fraction (f, y - 0.1*cos(4.*pi*x));
  //double ainit = 0.5;
  fraction(f, y);//+ 0.5*exp(-(x-10)*(x-10)/(2*ainit*ainit)));
}

/**
By default tracers are defined at $t-\Delta t/2$. We use the *first*
keyword to move VOF advection before the *amplitude* output i.e. at
$t+\Delta/2$. This improves the results. */

event vof (i++, first);

/**
clearWe output the amplitude at times matching exactly those in the
reference file. */

event amplitude (t += 0.005; t <= 20)
{
  scalar pos[];
  position (f, pos, {0,1});
  double max = statsf(pos).max;
  //output_ppm (pos,linear = true);
  /**
  We output the corresponding evolution in a file indexed with the
  number of grid points N and the current numerical surface tension. */
  char name[100];
  double st;
  st = f.sigma;
  sprintf (name, "Kakade/wgc_2d/data0907/BondHigh/wave-%d-%g.txt", N, st); //CHANGE DATE BEFORE RUNNING AND CREATE FOLDER
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t*16.032448313657, max);
  fflush (fp);
  fprintf (stderr, "Iteration Number %d\n",i);

  char namepos[150];
  sprintf (namepos, "Kakade/wgc_2d/data0907/BondHigh/wavepos-%g.txt", t*1000);
  FILE * fpos = fopen (namepos, "w");
  output_facets (f, fpos);
  fflush(fpos);
}

event display (t = end)
{
  double cap_length = sqrt(f.sigma/(rho1*G.y));
  fprintf (stderr, "Number of Grid Points per Wavelength: %g Numeric Surface Tension: %g\nNumeric Capillary Length: %g Numeric Liquid Viscosity: %g\n",N/L0,f.sigma,cap_length, mu1);
}

/* #if 0 */
/* event gfsview (i += 1) { */
/*   static FILE * fp = popen ("gfsview2D -s gravity.gfv", "w"); */
/*   output_gfs (fp); */
/* } */
/* #endif */

/* event images (t += 0.05) { */
/*   scalar pos[]; */
/*   position (f, pos, {0,1}); */
/*   output_ppm (pos); */
/* } */

event acceleration (i++)
{
  scalar phi = f.phi;
  double a = 0.5;   //Width of Gaussian Distributed Force
  double Po = -100; //Magnitude of Force Applied
  //double s = 2;     //Speed of Movement of Source
  double b = 10;  //Position of Applied Force
  coord G1;

  G1.y = -G.y;      //Direction of Applied Force
   if (phi.i)
     positionbumpstat(f, phi, a, Po, b, G1, Z, add = true);
     else {
     phi = new scalar;
     positionbumpstat(f, phi, a, Po, b, G1, Z, add = false);
     f.phi = phi;
    }
}

/* event adapt (i++) { */
/*   adapt({f}, (double[]){5e-3}, LEVEL + 1); */
/* } */
