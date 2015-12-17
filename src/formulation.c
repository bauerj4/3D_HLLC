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
	  //printf("vx.vy.vz = %f, %f, %f\n", vx[i][j][k], vy[i][j][k], vz[i][j][k]);
	  C5[i][j][k] = density[i][j][k] * (0.5 * v_dot_v +\
					    internal_energy[i][j][k]);

	  //printf("C1, C2, C3, C4, C5 = [%f, %f, %f, %f, %f]\n", C1[i][j][k], \
	    //	 C2[i][j][k], C3[i][j][k], C4[i][j][k], C5[i][j][k]);
	  if(C1[0] == 0)
	    {
	      printf("Warning: 0 density!\n");
	    }
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
          vz[i][j][k] =C4[i][j][k] / C1[i][j][k];
	  v_dot_v = vx[i][j][k] * vx[i][j][k] +  vy[i][j][k] * vy[i][j][k] + \
	    vz[i][j][k] * vz[i][j][k];
	  internal_energy[i][j][k] = C5[i][j][k] / density[i][j][k] - 0.5 * v_dot_v;
	  //printf("rho, vx, vy, vz, p, e = [%f, %f, %f, %f, %f, %f]\n", density[i][j][k], \
	    //      vx[i][j][k], vy[i][j][k], vz[i][j][k], pressure[i][j][k], internal_energy[i][j][k]);

	  // Get pressure from EoS
#ifdef IDEAL_GAS
	  p = pressure_ideal_gas(internal_energy[i][j][k], density[i][j][k], GAMMA);
#endif // IDEAL_GAS
	  pressure[i][j][k] = p;
        }
  return 0;
}
