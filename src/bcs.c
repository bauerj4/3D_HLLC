#include "../include/proto.h"
#include "../include/allvars.h"
#include "../include/hllc_defs.h"
#include <math.h>
#include <stdlib.h>

int set_HLLC_1_transmissive(double *** density, double *** vx, double *** vy, double *** vz,\
			    double *** pressure, double *** internal_energy, double *** C1, \
			    double *** C2, double *** C3, double *** C4, double *** C5)
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

	      C1[i][j][k] = C1[i+1][j][k];
	      C2[i][j][k] = C2[i+1][j][k];
	      C3[i][j][k] = C3[i+1][j][k];
	      C4[i][j][k] = C4[i+1][j][k];
	      C5[i][j][k] = C5[i+1][j][k];

	    }

          if (i == PARENT_RESOLUTION + ghost_cell_offset)
            {
              density[i][j][k] = density[i-1][j][k];
              vx[i][j][k] = vx[i-1][j][k];
              vy[i][j][k] = vy[i-1][j][k];
              vz[i][j][k] = vz[i-1][j][k];
              pressure[i][j][k] = pressure[i-1][j][k];
              internal_energy[i][j][k] = internal_energy[i-1][j][k];

	      C1[i][j][k] = C1[i-1][j][k];
              C2[i][j][k] = C2[i-1][j][k];
              C3[i][j][k] = C3[i-1][j][k];
              C4[i][j][k] = C4[i-1][j][k];
              C5[i][j][k] = C5[i-1][j][k];

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

	      C1[i][j][k] = C1[i][j+1][k];
              C2[i][j][k] = C2[i][j+1][k];
              C3[i][j][k] = C3[i][j+1][k];
              C4[i][j][k] = C4[i][j+1][k];
              C5[i][j][k] = C5[i][j+1][k];

            }

          if (j == PARENT_RESOLUTION + ghost_cell_offset)
            {
              density[i][j][k] = density[i][j-1][k];
              vx[i][j][k] = vx[i][j-1][k];
              vy[i][j][k] = vy[i][j-1][k];
              vz[i][j][k] = vz[i][j-1][k];
              pressure[i][j][k] = pressure[i][j-1][k];
              internal_energy[i][j][k] = internal_energy[i][j-1][k];

	      C1[i][j][k] = C1[i][j-1][k];
              C2[i][j][k] = C2[i][j-1][k];
              C3[i][j][k] = C3[i][j-1][k];
              C4[i][j][k] = C4[i][j-1][k];
              C5[i][j][k] = C5[i][j-1][k];

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

	      C1[i][j][k] = C1[i][j][k+1];
              C2[i][j][k] = C2[i][j][k+1];
              C3[i][j][k] = C3[i][j][k+1];
              C4[i][j][k] = C4[i][j][k+1];
              C5[i][j][k] = C5[i][j][k+1];

            }

          if (k == PARENT_RESOLUTION + ghost_cell_offset)
            {
              density[i][j][k] = density[i][j][k-1];
              vx[i][j][k] = vx[i][j][k-1];
              vy[i][j][k] = vy[i][j][k-1];
              vz[i][j][k] = vz[i][j][k-1];
              pressure[i][j][k] = pressure[i][j][k-1];
              internal_energy[i][j][k] = internal_energy[i][j][k-1];

	      C1[i][j][k] = C1[i][j][k-1];
              C2[i][j][k] = C2[i][j][k-1];
              C3[i][j][k] = C3[i][j][k-1];
              C4[i][j][k] = C4[i][j][k-1];
              C5[i][j][k] = C5[i][j][k-1];

            }

          if (j==0 && k==0)
            {
              density[i][j][k] = density[i][j+1][k+1];
              vx[i][j][k] = vx[i][j+1][k+1];
              vy[i][j][k] = vy[i][j+1][k+1];
              vz[i][j][k] = vz[i][j+1][k+1];
              pressure[i][j][k] = pressure[i][j+1][k+1];
              internal_energy[i][j][k] = internal_energy[i][j+1][k+1];

              C1[i][j][k] = C1[i][j+1][k+1];
              C2[i][j][k] = C2[i][j+1][k+1];
              C3[i][j][k] = C3[i][j+1][k+1];
              C4[i][j][k] = C4[i][j+1][k+1];
              C5[i][j][k] = C5[i][j+1][k+1];

            }

          if (i==0 && j==0)
            {
              density[i][j][k] = density[i+1][j+1][k];
              vx[i][j][k] = vx[i+1][j+1][k];
              vy[i][j][k] = vy[i+1][j+1][k];
              vz[i][j][k] = vz[i+1][j+1][k];
              pressure[i][j][k] = pressure[i+1][j+1][k];
              internal_energy[i][j][k] = internal_energy[i+1][j+1][k];

              C1[i][j][k] = C1[i+1][j+1][k];
              C2[i][j][k] = C2[i+1][j+1][k];
              C3[i][j][k] = C3[i+1][j+1][k];
              C4[i][j][k] = C4[i+1][j+1][k];
              C5[i][j][k] = C5[i+1][j+1][k];

            }

          if (i==0 && k==0)
            {
              density[i][j][k] = density[i+1][j][k+1];
              vx[i][j][k] = vx[i+1][j][k+1];
              vy[i][j][k] = vy[i+1][j][k+1];
              vz[i][j][k] = vz[i+1][j][k+1];
              pressure[i][j][k] = pressure[i+1][j][k+1];
              internal_energy[i][j][k] = internal_energy[i+1][j][k+1];

              C1[i][j][k] = C1[i+1][j][k+1];
              C2[i][j][k] = C2[i+1][j][k+1];
              C3[i][j][k] = C3[i+1][j][k+1];
              C4[i][j][k] = C4[i+1][j][k+1];
              C5[i][j][k] = C5[i+1][j][k+1];

            }



	  if (i==0 && j==0 && k==0)
	    {
	      density[i][j][k] = density[i+1][j+1][k+1];
              vx[i][j][k] = vx[i+1][j+1][k+1];
              vy[i][j][k] = vy[i+1][j+1][k+1];
              vz[i][j][k] = vz[i+1][j+1][k+1];
              pressure[i][j][k] = pressure[i+1][j+1][k+1];
              internal_energy[i][j][k] = internal_energy[i+1][j+1][k+1];

              C1[i][j][k] = C1[i+1][j+1][k+1];
              C2[i][j][k] = C2[i+1][j+1][k+1];
              C3[i][j][k] = C3[i+1][j+1][k+1];
              C4[i][j][k] = C4[i+1][j+1][k+1];
              C5[i][j][k] = C5[i+1][j+1][k+1];

	    }

          if (i==PARENT_RESOLUTION && j==PARENT_RESOLUTION)
            {
              density[i][j][k] = density[i-1][j-1][k];
              vx[i][j][k] = vx[i-1][j-1][k];
              vy[i][j][k] = vy[i-1][j-1][k];
              vz[i][j][k] = vz[i-1][j-1][k];
              pressure[i][j][k] = pressure[i-1][j-1][k];
              internal_energy[i][j][k] = internal_energy[i-1][j-1][k];

              C1[i][j][k] = C1[i-1][j-1][k];
              C2[i][j][k] = C2[i-1][j-1][k];
	      C3[i][j][k] = C3[i-1][j-1][k];
              C4[i][j][k] = C4[i-1][j-1][k];
              C5[i][j][k] = C5[i-1][j-1][k];

            }

          if (i==PARENT_RESOLUTION && k==PARENT_RESOLUTION)
            {
              density[i][j][k] = density[i-1][j][k-1];
              vx[i][j][k] = vx[i-1][j][k-1];
              vy[i][j][k] = vy[i-1][j][k-1];
              vz[i][j][k] = vz[i-1][j][k-1];
              pressure[i][j][k] = pressure[i-1][j][k-1];
              internal_energy[i][j][k] = internal_energy[i-1][j][k-1];

              C1[i][j][k] = C1[i-1][j][k-1];
              C2[i][j][k] = C2[i-1][j][k-1];
	      C3[i][j][k] = C3[i-1][j][k-1];
              C4[i][j][k] = C4[i-1][j][k-1];
              C5[i][j][k] = C5[i-1][j][k-1];

            }

          if (j==PARENT_RESOLUTION && k==PARENT_RESOLUTION)
            {
              density[i][j][k] = density[i][j-1][k-1];
              vx[i][j][k] = vx[i][j-1][k-1];
              vy[i][j][k] = vy[i][j-1][k-1];
              vz[i][j][k] = vz[i][j-1][k-1];
              pressure[i][j][k] = pressure[i][j-1][k-1];
              internal_energy[i][j][k] = internal_energy[i][j-1][k-1];

              C1[i][j][k] = C1[i][j-1][k-1];
              C2[i][j][k] = C2[i][j-1][k-1];
	      C3[i][j][k] = C3[i][j-1][k-1];
              C4[i][j][k] = C4[i][j-1][k-1];
              C5[i][j][k] = C5[i][j-1][k-1];

            }



          if (i==PARENT_RESOLUTION && j==PARENT_RESOLUTION && k==PARENT_RESOLUTION)
            {
              density[i][j][k] = density[i-1][j-1][k-1];
              vx[i][j][k] = vx[i-1][j-1][k-1];
              vy[i][j][k] = vy[i-1][j-1][k-1];
              vz[i][j][k] = vz[i-1][j-1][k-1];
              pressure[i][j][k] = pressure[i-1][j-1][k-1];
              internal_energy[i][j][k] = internal_energy[i-1][j-1][k-1];

              C1[i][j][k] = C1[i-1][j-1][k-1];
              C2[i][j][k] = C2[i-1][j-1][k-1];
              C3[i][j][k] = C3[i-1][j-1][k-1];
              C4[i][j][k] = C4[i-1][j-1][k-1];
              C5[i][j][k] = C5[i-1][j-1][k-1];

            }


	}

  return 0;
}

int set_HLLC_1_periodic(double *** density, double *** vx, double *** vy, double *** vz, \
			double *** pressure, double *** internal_energy)
{
  return 0;
}

