// Author: Radu Cimpeanu
// Date: 24/10/2024

#include "axi.h"                     // axisymmetric geometry
#include "navier-stokes/centered.h"  // solve NS equations
#define FILTERED // Smear density and viscosity jumps
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "tension.h"                 // include surf tension between phases
#include "vof.h"                     // solve using VoF method
#include "fractions.h"               // initially define fraction of fluid phases
#include "view.h"                    // need to make the animations
#include "draw.h"                    // visualisation helper
#include "tag.h"                     // helps track droplet properties

// dimensional quantities (fluids at 20 degrees C)
double rhoLiquid; // liquid phase density (kg/m^3)
double rhoGas;    // gas phase density (kg/m^3)

double muLiquid;  // liquid dynamic viscosity (kg/ms)
double muGas;     // gas dynamic viscosity(kg/ms)

double sig;       // surface tension (N/m)

double g_accel;   // gravitational acceleration (m/s^2)

double dRadius;   // drop radius (m)

double v_init;    // initial drop velocity (m/s)

// dimensionless specifications (key groupings defined in main below)
#define rho_ratio   (rhoGas/rhoLiquid) // density ratio
#define mu_ratio    (muGas/muLiquid)   // viscosity ratio 

#define poolHeight 0.0             // Pool height (in radii)
#define domainSize 8.0               // Computational box size (in radii)

face vector av[];

FILE * fp_stats;
FILE * fp_vol;
FILE * fp_droplets;

double ND_Weber;
double ND_Reynolds;
double ND_Froude;
double ND_Bond;
double ND_Ohnesorge;

double filmHeight;

int minLevel = 4;
int maxLevel;// = 10;

double tEnd;

// Solid BCs for when using Gravity

// Bottom of domain = LEFT of domain
// No slip, impermeability, droplet not touching
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
f[left] = dirichlet(0.);

// Laterally
// Free slip and impermeability
u.n[top] = dirichlet(0.); // Impermeability
u.t[top] = neumann(0.); // Slip

// Top of domain
// Outflow
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

// Default for bottom is symmetry

int main(int argc, char * argv[]) {

  rhoLiquid = atof(argv[1]);  // prescribed liquid density 
  rhoGas = atof(argv[2]);     // prescribed gas density
  muLiquid = atof(argv[3]);   // prescribed liquid (dynamic) viscosity
  muGas = atof(argv[4]);      // prescribed gas (dynamic) viscosity
  sig = atof(argv[5]);      // prescribed surface tension coefficient
  g_accel = atof(argv[6]);    // prescribed gravitational acceleration
  dRadius = atof(argv[7]);    // prescribed drop radius
  v_init = atof(argv[8]);     // prescribed initial velocity
  tEnd = atof(argv[9]);       // prescribed simulation end time
  maxLevel = atof(argv[10]);  // prescribed maximum resolution level
  
  ND_Weber = (rhoLiquid*pow(v_init,2.0)*dRadius)/sig;
  ND_Reynolds = (rhoLiquid*v_init*dRadius)/muLiquid;
  ND_Froude = v_init/pow(dRadius*g_accel,0.5);
  ND_Bond = rhoLiquid*g_accel*pow(dRadius,2.0)/sig;
  ND_Ohnesorge = muLiquid/pow(rhoLiquid*sig*dRadius,0.5);
    
  init_grid(1 << 6);
  
  size(domainSize);                     
  origin(-0.5*domainSize, 0.0);

  //make folders
  mkdir("Slices", 0700);
  mkdir("Animations", 0700);
  mkdir("Interfaces", 0700);
  //mkdir("Data", 0700);

  //print dimensionless numbers
  fprintf(stdout, "Reynolds number = %0.6f \n", ND_Reynolds); fflush(stdout);
  fprintf(stdout, "Weber number = %0.6f \n", ND_Weber); fflush(stdout);
  fprintf(stdout, "Froude number = %0.6f \n", ND_Froude); fflush(stdout);
  fprintf(stdout, "Bond number = %0.6f \n", ND_Bond); fflush(stdout);
  fprintf(stdout, "Ohnesorge number = %0.6f \n", ND_Ohnesorge); fflush(stdout);

  rho1 = 1.;
  rho2 = rho_ratio;
  
  mu1 = 1./ND_Reynolds;
  mu2 = mu_ratio*mu1;
  
  f.sigma = 1./ND_Weber;

  a = av;

  // Pointer of the file to save stats
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");
  }

  // Pointer of the file to droplet count info
  {
    char name[200];
    sprintf(name, "logdroplets.dat");
    fp_droplets = fopen(name, "w");
  }

  DT = 1e-3;
  NITERMIN = 1; // default 1
  NITERMAX = 300; // default 100
  TOLERANCE = 1e-4; // default 1e-3
  
  run();

  fclose(fp_stats);
  //fclose(fp_vol);
  fclose(fp_droplets);
}

event acceleration (i++) {
  foreach_face(x)  
    av.x[] -= 1./pow(ND_Froude,2.0);
  foreach_face(y)  
    av.y[] += 0.0;
}

scalar omega[], viewingfield[], mylevel[];

event init (t = 0.0) {

  filmHeight = -domainSize/2. + poolHeight;

  // strong refinement around the interfacial regions
  refine (((sq(x - (filmHeight + 1.5 + 0.5)) + sq(y) < sq(1.0*1.05) && sq(x - (filmHeight + 1.5 + 0.5)) + sq(y) > sq(1.0*0.95)) || fabs(x - filmHeight) <= 0.005) && level < maxLevel);
  
  // creative active liquid phase as union between drop and film
  fraction (f, sq(1.0) - sq(x - (filmHeight + 1.5 + 0.5)) - sq(y));
  
  // initialise uniform velocity field inside droplet
  foreach()
  {
  	u.x[] = -1.0*f[];
        u.y[] = 0.0;
        p[] = 0.0;
	omega[] = 0.0;
  }
}

event adapt (i++) {

  // refine only with respect to interfacial shape location and velocity component magnitude
  adapt_wavelet ((scalar *){f, u}, (double[]){1e-6, 4e-3, 4e-3}, maxLevel, minLevel);

}

event gfsview (t = 0.0; t += 1.0; t <= tEnd) {
    char name_gfs[200];
    sprintf(name_gfs,"Slices/DropImpact-%0.1f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}

event saveInterfaces (t += 0.01) {

    char nameInterfaces1[200];

    sprintf(nameInterfaces1,"Interfaces/interfaceDrop-%0.2f.dat",t);

    FILE * fp1 = fopen(nameInterfaces1, "w");
    output_facets (f, fp1);	
    fclose(fp1);
}

event small_droplet_removal (i++) {
// Removes any small droplets that have formed, that are smaller than a specific size
    remove_droplets(f, 8);       // Removes droplets of diameter 8 cells or less
    remove_droplets(f, 8, true); // Removes bubbles of diameter 8 cells or less
}

event droplets (t += 0.01)
{
  scalar m[];
  foreach()
    m[] = f[] > 1e-2;
  int n = tag (m);

  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      coord p = {x,y,z};
      foreach_dimension()
	b[j].x += dv()*f[]*p.x;
    }

  #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  for (int j = 0; j < n; j++)
    fprintf (fp_droplets, "%d %g %d %g %g %g\n", i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j]);
  fflush (fp_droplets);
}

event movies (t += 0.01){

  char timestring[100];
  
  foreach(){
  	viewingfield[] = 1.0 - f[];
  	mylevel[] = level;
  }
  
  view(width=1200, height=800, fov=30.0, ty = 0.0, quat = { 0, 0, -0.707, 0.707 });
  clear();
  
  draw_vof("f", lw=2);
  squares("viewingfield", map = cool_warm, min = -0.5, max = 2.5);
  mirror({0,1}) {
  	draw_vof("f", lw=2);	
	cells(lw=0.5);
	squares("mylevel", map = cool_warm, min = minLevel, max = maxLevel);
  } 

  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/ImpactSummary.mp4");

}

event logstats (t += 0.01; t <= tEnd) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}
