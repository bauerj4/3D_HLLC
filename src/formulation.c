#include "../include/allvars.h"
#include "../include/hllc_defs.h"
#include "../include/proto.h"
#include <math.h>

int toConservative(double *** density, double *** vx, double *** vy, double *** vz, \
		   double *** pressure, double *** internal_energy, double *** C1, \
		   double *** C2, double *** C3, double *** C4, double *** C5)
{
  double v_dot_v;
  int i,j,k;
  for (i=0; i < (int)(PARENT_RESOLUTION + 2*ghost_cell_offset); i++)
    for (j=0; j < (int)(PARENT_RESOLUTION + 2*ghost_cell_offset); j++)
      for (k=0; k < (int)(PARENT_RESOLUTION + 2*ghost_cell_offset); k++)
	{
	  v_dot_v = vx[i][j][k] * vx[i][j][k] +  vy[i][j][k] * vy[i][j][k] +\
	    vz[i][j][k] * vz[i][j][k];
	  C1[i][j][k] = density[i][j][k];
	  C2[i][j][k] = density[i][j][k] * vx[i][j][k]; 
          C3[i][j][k] = density[i][j][k] * vy[i][j][k];
          C4[i][j][k] = density[i][j][k] * vz[i][j][k];
	  C5[i][j][k] = density[i][j][k] * (0.5 * v_dot_v +\
					    internal_energy[i][j][k]);
	}

  
  return 0;
}
int toPrimitive(double *** density, double *** vx, double *** vy, double *** vz, \
		double *** pressure, double *** internal_energy, double *** C1,\
		double *** C2, double *** C3, double *** C4, double *** C5)
{
  double v_dot_v, p;
  int i,j,k;
  for (i=0; i < (int)(PARENT_RESOLUTION + 2*ghost_cell_offset); i++)
    for (j=0; j < (int)(PARENT_RESOLUTION + 2*ghost_cell_offset); j++)
      for (k=0; k < (int)(PARENT_RESOLUTION + 2*ghost_cell_offset); k++)
        {
	  density[i][j][k] = C1[i][j][k];
	  vx[i][j][k] = C2[i][j][k] / C1[i][j][k];
          vy[i][j][k] = C3[i][j][k] / C1[i][j][k];
          vz[i][j][k] =C2[i][j][k] / C1[i][j][k];
	  v_dot_v = vx[i][j][k] * vx[i][j][k] +  vy[i][j][k] * vy[i][j][k] + \
	    vz[i][j][k] * vz[i][j][k];
	  internal_energy[i][j][k] = C5[i][j][k] / density[i][j][k] - 0.5 * v_dot_v;

	  // Get pressure from EoS
#ifdef IDEAL_GAS
	  p = pressure_ideal_gas(internal_energy[i][j][k], density[i][j][k], GAMMA);
#endif // IDEAL_GAS
	  pressure[i][j][k] = p;
        }
  return 0;
}
