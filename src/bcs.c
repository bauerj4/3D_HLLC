#include "../include/proto.h"
#include "../include/allvars.h"
#include "../include/hllc_defs.h"
#include <math.h>
#include <stdlib.h>

int set_HLLC_1_transmissive(double *** density, double *** vx, double *** vy, double *** vz,\
			    double *** pressure, double *** internal_energy)
{
  int i,j,k;
  for (i=0; i < (int)(PARENT_RESOLUTION + 2*ghost_cell_offset); i++)
    for (j=0; j < (int)(PARENT_RESOLUTION + 2*ghost_cell_offset); j++)
      for (k=0; k < (int)(PARENT_RESOLUTION + 2*ghost_cell_offset); k++)
	{
	  // X boundary
	  
	  if (i == 0)
	    {
	      density[i][j][k] = density[i+1][j][k];
	      vx[i][j][k] = vx[i+1][j][k];
	      vy[i][j][k] = vy[i+1][j][k];
	      vz[i][j][k] = vz[i+1][j][k];
	      pressure[i][j][k] = pressure[i+1][j][k];
	      internal_energy[i][j][k] = internal_energy[i+1][j][k];
	    }

          if (i == PARENT_RESOLUTION + ghost_cell_offset)
            {
              density[i][j][k] = density[i-1][j][k];
              vx[i][j][k] = vx[i-1][j][k];
              vy[i][j][k] = vy[i-1][j][k];
              vz[i][j][k] = vz[i-1][j][k];
              pressure[i][j][k] = pressure[i-1][j][k];
              internal_energy[i][j][k] = internal_energy[i-1][j][k];
            }

	  // Y boundary

          if (j == 0)
            {
              density[i][j][k] = density[i][j+1][k];
              vx[i][j][k] = vx[i][j+1][k];
              vy[i][j][k] = vy[i][j+1][k];
              vz[i][j][k] = vz[i][j+1][k];
              pressure[i][j][k] = pressure[i][j+1][k];
              internal_energy[i][j][k] = internal_energy[i][j+1][k];
            }

          if (j == PARENT_RESOLUTION + ghost_cell_offset)
            {
              density[i][j][k] = density[i][j-1][k];
              vx[i][j][k] = vx[i][j-1][k];
              vy[i][j][k] = vy[i][j-1][k];
              vz[i][j][k] = vz[i][j-1][k];
              pressure[i][j][k] = pressure[i][j-1][k];
              internal_energy[i][j][k] = internal_energy[i][j-1][k];
            }

	  // Z boundary
          if (k == 0)
            {
              density[i][j][k] = density[i][j][k+1];
              vx[i][j][k] = vx[i][j][k+1];
              vy[i][j][k] = vy[i][j][k+1];
              vz[i][j][k] = vz[i][j][k+1];
              pressure[i][j][k] = pressure[i][j][k+1];
              internal_energy[i][j][k] = internal_energy[i][j][k+1];
            }

          if (k == PARENT_RESOLUTION + ghost_cell_offset)
            {
              density[i][j][k] = density[i][j][k-1];
              vx[i][j][k] = vx[i][j][k-1];
              vy[i][j][k] = vy[i][j][k-1];
              vz[i][j][k] = vz[i][j][k-1];
              pressure[i][j][k] = pressure[i][j][k-1];
              internal_energy[i][j][k] = internal_energy[i][j][k-1];
            }



	}

  return 0;
}

int set_HLLC_1_periodic(double *** density, double *** vx, double *** vy, double *** vz, \
			double *** pressure, double *** internal_energy)
{
  return 0;
}

