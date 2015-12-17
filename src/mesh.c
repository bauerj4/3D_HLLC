#include "../include/allvars.h"
#include "../include/proto.h"
#include "../include/hllc_defs.h"
#include <stdio.h>

int write_mesh(char * FILE_NAME)
{
  int i,j,k;
  snapshot = fopen(FILE_NAME, "w");
  for (i=0; i < PARENT_RESOLUTION + 2*ghost_cell_offset; i++)
    for (j=0; j < PARENT_RESOLUTION + 2*ghost_cell_offset; j++)
      for (k=0; k < PARENT_RESOLUTION + 2*ghost_cell_offset; k++)
	{
	  /*
	  fwrite(&density[i][j][k], sizeof(double), 1,  snapshot);
	  fwrite(&vx[i][j][k], sizeof(double), 1,  snapshot);
	  fwrite(&vy[i][j][k], sizeof(double), 1,  snapshot);
	  fwrite(&vz[i][j][k], sizeof(double), 1,  snapshot);
	  fwrite(&pressure[i][j][k], sizeof(double), 1,  snapshot);
	  fwrite(&internal_energy[i][j][k], sizeof(double), 1, snapshot);
	  */
	  
          fwrite(&C1[i][j][k], sizeof(double), 1,  snapshot);
          fwrite(&C2[i][j][k], sizeof(double), 1,  snapshot);
          fwrite(&C3[i][j][k], sizeof(double), 1,  snapshot);
          fwrite(&C4[i][j][k], sizeof(double), 1,  snapshot);
          fwrite(&C5[i][j][k], sizeof(double), 1,  snapshot);
	  
          //fwrite(&internal_energy[i][j][k], sizeof(double), 1, snapshot);

	}
  return 0;
}

int HLLC1()
{
  write_mesh("ics.dat");
  while (CurrentTime < SIMULATION_TIME)
    {
      advanceTimestep_HLLC1();
      toPrimitive(density, vx,vy,vz, pressure, internal_energy,	\
			  C1, C2, C3, C4, C5);
      toConservative(density, vx,vy,vz, pressure, internal_energy, \
		     C1, C2, C3, C4, C5);

      set_HLLC_1_transmissive(density, vx,vy,vz, pressure, internal_energy, C1, C2,C3, C4,C5);
      //toConservative(density, vx,vy,vz, pressure, internal_energy,	\
	//	     C1, C2, C3, C4, C5);

    }
  printf("\n\nSimulation time %f reached.  Simulation ends.\n", SIMULATION_TIME);
  write_mesh("final_snapshot.dat");
  return 0;
}

int advanceTimestep_HLLC1()
{
  compute_X_HLLC_Fluxes(density, vx,vy,vz,pressure,internal_energy,\
			C1,C2,C3,C4,C5,F_R, F_RSTAR, F_L,F_LSTAR);
  dx = X_MAX / ((double) PARENT_RESOLUTION);
  dy = Y_MAX / ((double) PARENT_RESOLUTION);
  dz = Z_MAX / ((double) PARENT_RESOLUTION);

  dt = CFL * dx / S_max;

  // if dt y < dt x -> dt = dty 

  CurrentTime += dt;
  printf("\n\nMax Signal Speed: %10.10f\n", S_max);
  printf("Time: %10.5f \t dt = %5.10f\n", CurrentTime, dt);
  
  int i,j,k;
  //printf("%f\n",F_TOTAL[1][1][1][0]);
  for (i=ghost_cell_offset; i < PARENT_RESOLUTION + ghost_cell_offset; i++)
    for (j=ghost_cell_offset; j < PARENT_RESOLUTION + ghost_cell_offset; j++)
      for (k=ghost_cell_offset; k < PARENT_RESOLUTION + ghost_cell_offset; k++)
	{
	  //printf("k = %d < %d\n", k, (int)PARENT_RESOLUTION + ghost_cell_offset);
	  //printf("[%f, %f, %f, %f, %f]\n", F_TOTAL[i][j][k][0], F_TOTAL[i][j][k][1], \
	    // F_TOTAL[i][j][k][2],  F_TOTAL[i][j][k][3],  F_TOTAL[i][j][k][4]);
	  C1[i][j][k] -= dt/dx * (F_TOTAL[i][j][k][0] - F_TOTAL[i-1][j][k][0]);
	  C2[i][j][k] -= dt/dx * (F_TOTAL[i][j][k][1] - F_TOTAL[i-1][j][k][1]);
          C3[i][j][k] -= dt/dx * (F_TOTAL[i][j][k][2] - F_TOTAL[i-1][j][k][2]);
          C4[i][j][k] -= dt/dx * (F_TOTAL[i][j][k][3] - F_TOTAL[i-1][j][k][3]);
          C5[i][j][k] -= dt/dx * (F_TOTAL[i][j][k][4] - F_TOTAL[i-1][j][k][4]);	  
	}

  return 0;
}

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
  double v_dot_v;
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
	      C1[i][j][k] = density[i][j][k];
	      C2[i][j][k] = density[i][j][k] * vx[i][j][k];
	      C3[i][j][k] = density[i][j][k] * vy[i][j][k];
	      C4[i][j][k] = density[i][j][k] * vz[i][j][k];
	      v_dot_v = vx[i][j][k] * vx[i][j][k] + vy[i][j][k] * vy[i][j][k] +\
		vz[i][j][k] * vz[i][j][k];

	      C5[i][j][k] = density[i][j][k] * (0.5 * v_dot_v + internal_energy[i][j][k]);

	      if (density[i][j][k] == 0)
		{
		  printf("WARNING: %f DENSITY\n",density[i][j][k]);
		}

              if (pressure[i][j][k] == 0)
                {
                  printf("WARNING: 0 PRESSURE\n");
                }

              if (internal_energy[i][j][k] == 0)
                {
                  printf("WARNING: 0 ENERGY\n");
                }


	    }
	}
    }

  printf("Allocated %10.1f MB of mesh data.\n", (int)byte_size/1e6);
  return 0;
}

