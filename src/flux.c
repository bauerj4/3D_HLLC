#include "../include/hllc_defs.h"
#include "../include/hllc_defs.h"
#include "../include/allvars.h"
#include "../include/proto.h"
#include <stdlib.h>
#include <math.h>

int compute_USTAR_XL(double * U_STAR,int i, int j, int k)
{
  double coeff = density[i][j][k] * (S_L - vx[i][j][k])/(S_L - S_star);
  U_STAR[0] = coeff;
  U_STAR[1] = coeff * S_star;
  U_STAR[2] = coeff * vy[i][j][k];
  U_STAR[3] = coeff * vz[i][j][k];
  double fn = C5[i][j][k] / density[i][j][k] + (S_star - vx[i][j][k] *\
						 (S_star + pressure[i][j][k]/\
						  (density[i][j][k] * (S_L - vx[i][j][k]))));
  U_STAR[4] = coeff * fn;
  return 0;
}

int compute_USTAR_XR(double * U_STAR,int i, int j, int k)
{
  double coeff = density[i+1][j][k] * (S_L - vx[i+1][j][k])/(S_R - S_star);
  U_STAR[0] = coeff;
  U_STAR[1] = coeff * S_star;
  U_STAR[2] = coeff * vy[i+1][j][k];
  U_STAR[3] = coeff * vz[i+1][j][k];
  double fn = C5[i+1][j][k] / density[i+1][j][k] + (S_star - vx[i+1][j][k] *\
						    (S_star + pressure[i+1][j][k]/ \
						     (density[i+1][j][k] * (S_R - vx[i+1][j][k]))));
  U_STAR[4] = coeff * fn;
  return 0;
}


int compute_X_HLLC_Fluxes(double *** density, double *** vx, double *** vy, double *** vz,\
                     double *** pressure, double *** internal_energy, double *** C1,\
		     double *** C2, double *** C3, double *** C4, double *** C5,\
                     double *** F_R, double *** F_RSTAR, double *** F_L, double *** F_LSTAR)
{
  int i,j,k;
  for (i=0; i < (int)(PARENT_RESOLUTION + ghost_cell_offset); i++)
    for (j=0; j < (int)(PARENT_RESOLUTION + ghost_cell_offset); j++)
      for (k=0; k < (int)(PARENT_RESOLUTION + ghost_cell_offset); k++)
	{
	  
	}

  return 0;
}

int compute_Y_HLLC_Fluxes(double *** density, double *** vx, double *** vy, double *** vz,\
                     double *** pressure, double *** internal_energy, double *** C1,\
		     double *** C2, double *** C3, double *** C4, double *** C5,\
                     double *** G_R, double *** G_RSTAR, double *** G_L, double *** G_LSTAR)
{
  return 0;
}

int compute_Z_HLLC_Fluxes(double *** density, double *** vx, double *** vy, double *** vz,\
                     double *** pressure, double *** internal_energy, double *** C1,\
		     double *** C2, double *** C3, double *** C4, double *** C5,\
                     double *** H_R, double *** H_RSTAR, double *** H_L, double *** H_LSTAR)
{
  return 0;
}




int initialize_flux_arrays()
{
  printf("Allocating flux arrays...\n");
  int i,j,k;
  U_STAR =(double *) malloc(5 * sizeof(double));

  F_R = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));
  F_RSTAR = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));
  F_LSTAR = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));
  F_L = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));

  G_R = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));
  G_RSTAR = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));
  G_LSTAR = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));
  G_L = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));


  H_R = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));
  H_RSTAR = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));
  H_LSTAR = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));
  H_L = (double ***) malloc((int) PARENT_RESOLUTION*sizeof(double **));


  for (i=0; i < (int)(PARENT_RESOLUTION + ghost_cell_offset); i++)
    {
      F_R[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
      F_RSTAR[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
      F_LSTAR[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
      F_L[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));

      G_R[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
      G_RSTAR[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
      G_LSTAR[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
      G_L[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));

      H_R[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
      H_RSTAR[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
      H_LSTAR[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
      H_L[i] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));

      for (j=0; j < (int)(PARENT_RESOLUTION +  ghost_cell_offset); j++)
	{
	  F_R[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
	  F_RSTAR[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
	  F_LSTAR[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
	  F_L[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));

          G_R[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
          G_RSTAR[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
          G_LSTAR[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
          G_L[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));

          H_R[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
          H_RSTAR[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
          H_LSTAR[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
          H_L[i][j] = (double*)malloc((int) PARENT_RESOLUTION * sizeof(double));
	}
    }
  byte_size += PARENT_RESOLUTION*PARENT_RESOLUTION*PARENT_RESOLUTION*sizeof(double) * 4 * 3;
  printf("Currently allocated %10.1f MB for simulation data.\n",byte_size/1e6);

  return 0;
}
