#include "../include/allvars.h"
#include "../include/proto.h"
#include "../include/hllc_defs.h"
#include <stdio.h>

int read_parent_mesh()
{
  parent_mesh = fopen("parent_mesh.dat","r");
  int i,j,k;
  printf("Allocating memory...\n");
  printf("The top level block is %d bytes\n",(int)PARENT_RESOLUTION*sizeof(double**));

  // Allocate top level memory and set ghost cell size.
#ifdef HLLC_1
  pressure = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  vx = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  vy = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  vz = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  density = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  internal_energy = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);

  C1 = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  C2 = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  C3 = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  C4 = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  C5 = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  //C6 = (double ***)malloc((int)PARENT_RESOLUTION*sizeof(double**) + 2);
  ghost_cell_offset = 1;
#endif

  // Allocate intermediate level memory
  for (i=0; i < (int)(PARENT_RESOLUTION + 2*ghost_cell_offset); i++)
    {
    pressure[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
    vx[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
    vy[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
    vz[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
    density[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
    internal_energy[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));

    C1[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
    C2[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
    C3[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
    C4[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
    C5[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
    //C6[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));

    for (j=0; j < (int)(PARENT_RESOLUTION + 2 * ghost_cell_offset); j++)
      {
	pressure[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
        vx[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
        vy[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
        vz[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
        density[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
	internal_energy[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));

        C1[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
        C2[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
        C3[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
        C4[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
        C5[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
        //C6[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));

      }
    }
  byte_size += PARENT_RESOLUTION*PARENT_RESOLUTION*PARENT_RESOLUTION*sizeof(double) * 12;
  byte_size += ghost_cell_offset * PARENT_RESOLUTION * PARENT_RESOLUTION * (12 * 4);
  for (i=ghost_cell_offset; i < PARENT_RESOLUTION + ghost_cell_offset; i++)
    {
      for (j=ghost_cell_offset; j < PARENT_RESOLUTION + ghost_cell_offset; j++)
	{
	  for (k=ghost_cell_offset; k < PARENT_RESOLUTION + ghost_cell_offset; k++)
	    {
	      fread(&density[i][j][k], sizeof(double), 1,  parent_mesh);
              fread(&vx[i][j][k], sizeof(double), 1,  parent_mesh);
              fread(&vy[i][j][k], sizeof(double), 1,  parent_mesh);
              fread(&vz[i][j][k], sizeof(double), 1,  parent_mesh);
              fread(&pressure[i][j][k], sizeof(double), 1,  parent_mesh);
#ifdef IDEAL_GAS
	      internal_energy[i][j][k] = internal_energy_ideal_gas(pressure[i][j][k], \
								  density[i][j][k], GAMMA);
#endif
	    }
	}
    }

  printf("Allocated %10.1f MB of mesh data.\n", (int)byte_size/1e6);
  return 0;
}

