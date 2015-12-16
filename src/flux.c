#include "../include/hllc_defs.h"
#include "../include/hllc_defs.h"
#include "../include/allvars.h"
#include "../include/proto.h"
#include <stdlib.h>
#include <math.h>

int conservedToXFlux(double * U, double * F)
{
  double rho = U[0];
  double u = U[1]/rho;
  double v = U[2]/rho;
  double w = U[3]/rho;
  double u_dot_u = u*u + v*v + w*w;
  double e = U[4]/rho - 0.5 * u_dot_u;
#ifdef IDEAL_GAS
  double p = pressure_ideal_gas(e,rho,GAMMA);
#endif// IDEAL_GAS
  F[0] = rho * u;
  F[1] = rho *u * u + p;
  F[2] = rho  * u * v;
  F[3] = rho * u * w;
  F[4] = u * (U[4] + p);

  return 0;
}

int compute_X_STAR_state_vars(int i, int j, int k)
{
  rho_bar = 0.5 * (density[i][j][k] + density[i+1][j][k]);
  a_L = sound_speed_ideal_gas(pressure[i][j][k], density[i][j][k], GAMMA);
  a_R =sound_speed_ideal_gas(pressure[i+1][j][k], density[i+1][j][k], GAMMA);
  //  printf("rho_bar, a_L, a_R = %f, %f, %f\n", rho_bar, a_L, a_R);
  a_bar = 0.5 *(a_L + a_R);
  p_bar = 0.5 * (pressure[i][j][k] + pressure[i+1][j][k]);
  pvrs = p_bar - 0.5 * (vx[i+1][j][k] - vx[i][j][k])*rho_bar*a_bar;
  if (0 > pvrs)
    {
      p_star = 0;
    }
  else
    {
      p_star = pvrs;
    }

#ifdef WAVESPEED_ENTHALPY
  QL = 1;
  QR = 1;
  
  if (p_star > pressure[i][j][k])
    {
      QL = 1 + (GAMMA + 1)/(2 * GAMMA) *\
	(p_star / pressure[i][j][k] - 1);
      QL = sqrt(QL);
    }
  if (p_star > pressure[i+1][j][k])
    {
      QR = 1 + (GAMMA +1)/(2 * GAMMA) *		\
	(p_star / pressure[i+1][j][k] - 1);
      QR = sqrt(QL);
    }
  
  S_L = vx[i][j][k] - a_L * QL;
  S_R = vx[i+1][j][k] + a_R * QR;

#endif

#ifdef WAVESPEED_SIMPLE
  S_L = vx[i][j][k] - a_L;
  if (vx[i+1][j][k] - a_R < S_L)
    {
      S_L = vx[i+1][j][k] - a_R;
    }

  S_R =vx[i][j][k] + a_L;
  if (vx[i+1][j][k] - a_R > S_R)
    {
      S_R = vx[i+1][j][k] + a_R;
    }

#endif // WAVESPEED_SIMPLE

  S_star = pressure[i+1][j][k] - pressure[i][j][k] + density[i][j][k] * vx[i][j][k] *\
    (S_L - vx[i][j][k]) - density[i+1][j][k] * vx[i+1][j][k] * (S_R - vx[i+1][j][k]);
  S_star /= (density[i][j][k]*(S_L - vx[i][j][k]) - density[i+1][j][k] * (S_R - vx[i+1][j][k]));
  
  //printf("S_L. S_R, S_star = %f, %f, %f\n", S_L, S_R, S_star);
  return 0;
}

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
                     double **** F_R, double **** F_RSTAR, double **** F_L, double **** F_LSTAR)
{
  int i,j,k;
  S_max = 0;
  for (i=0; i < (int)(PARENT_RESOLUTION + ghost_cell_offset); i++)
    for (j=0; j < (int)(PARENT_RESOLUTION + ghost_cell_offset); j++)
      for (k=0; k < (int)(PARENT_RESOLUTION + ghost_cell_offset); k++)
	{
	  compute_X_STAR_state_vars(i,j,k);
	  //printf("S_star = %10.10f\n",S_star);

	  if (fabs(S_star) > S_max)
	    {
	      //printf("S_star = %10.10f\n",S_star);
	      S_max = S_star;
	    }
	  U_L[0] = C1[i][j][k];
	  U_L[1] = C2[i][j][k];
	  U_L[2] = C3[i][j][k];
	  U_L[3] = C4[i][j][k];
	  U_L[4] = C5[i][j][k];
	  
	  U_R[0] = C1[i+1][j][k];
	  U_R[1] = C2[i+1][j][k];
	  U_R[2] = C3[i+1][j][k];
	  U_R[3] = C4[i+1][j][k];
	  U_R[4] = C5[i+1][j][k];
	  
	  compute_USTAR_XL(U_LSTAR,i, j,k);
	  compute_USTAR_XR(U_RSTAR,i, j,k);
	  
	  conservedToXFlux(U_R, F_R[i][j][k]);
	  conservedToXFlux(U_L, F_L[i][j][k]);
	  conservedToXFlux(U_LSTAR, F_LSTAR[i][j][k]);
	  conservedToXFlux(U_RSTAR, F_RSTAR[i][j][k]);
#ifdef HLLC_1
	  if (S_L >= 0)
	    {
	      F_TOTAL[i][j][k] = F_L[i][j][k];
	    }
	  else if(S_L <= 0 && S_star >= 0)
	    {
	      F_TOTAL[i][j][k] = F_LSTAR[i][j][k];
	    }
	  else if (S_star <= 0 && S_R >= 0)
	    {
	      F_TOTAL[i][j][k] = F_RSTAR[i][j][k];
	    }
	  else if (S_R <= 0)
	    {
	      F_TOTAL[i][j][k] = F_R[i][j][k];
	    }
#endif // HLLC_1
	    
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
  U_LSTAR =(double *) malloc(5 * sizeof(double));
  U_RSTAR =(double *) malloc(5 * sizeof(double));
  U_L =(double *) malloc(5 * sizeof(double));
  U_R =(double *) malloc(5 * sizeof(double));


  F_R = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  F_RSTAR = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  F_LSTAR = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  F_L = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  F_TOTAL = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));


  G_R = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  G_RSTAR = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  G_LSTAR = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  G_L = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  G_TOTAL = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));


  H_R = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  H_RSTAR = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  H_LSTAR = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  H_L = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));
  H_TOTAL = (double ****) malloc((int) PARENT_RESOLUTION*sizeof(double ***));


  for (i=0; i < (int)(PARENT_RESOLUTION + ghost_cell_offset); i++)
    {
      F_R[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      F_RSTAR[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      F_LSTAR[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      F_L[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      F_TOTAL[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));

      G_R[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      G_RSTAR[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      G_LSTAR[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      G_L[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      G_TOTAL[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));

      H_R[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      H_RSTAR[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      H_LSTAR[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      H_L[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));
      H_TOTAL[i] = (double***)malloc((int) PARENT_RESOLUTION * sizeof(double**));

      for (j=0; j < (int)(PARENT_RESOLUTION +  ghost_cell_offset); j++)
	{
	  F_R[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
	  F_RSTAR[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
	  F_LSTAR[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
	  F_L[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
          F_TOTAL[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));

          G_R[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
          G_RSTAR[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
          G_LSTAR[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
          G_L[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
          G_TOTAL[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));

          H_R[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
          H_RSTAR[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
          H_LSTAR[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
          H_L[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));
          H_TOTAL[i][j] = (double**)malloc((int) PARENT_RESOLUTION * sizeof(double*));

	
	  for (k=0; k < (int)(PARENT_RESOLUTION +  ghost_cell_offset); k++)
	    {
	      //printf("k = %d\n",k);
	      F_R[i][j][k] = (double*)malloc(5 * sizeof(double));
	      F_L[i][j][k] = (double*)malloc(5 * sizeof(double));
	      F_LSTAR[i][j][k] = (double*)malloc(5 * sizeof(double));
	      F_RSTAR[i][j][k] = (double*)malloc(5 * sizeof(double));
	      F_TOTAL[i][j][k] = (double*)malloc(5 * sizeof(double));
	      
	      //printf("F allocated.\n");
	      
	      G_R[i][j][k] = (double*)malloc(5 * sizeof(double));
	      G_L[i][j][k] = (double*)malloc(5 * sizeof(double));
	      G_LSTAR[i][j][k] = (double*)malloc(5 * sizeof(double));
	      G_RSTAR[i][j][k] = (double*)malloc(5 * sizeof(double));
	      G_TOTAL[i][j][k] = (double*)malloc(5 * sizeof(double));
	      
	      H_R[i][j][k] = (double*)malloc(5 * sizeof(double));
	      H_L[i][j][k] = (double*)malloc(5 * sizeof(double));
	      H_LSTAR[i][j][k] = (double*)malloc(5 * sizeof(double));
	      H_RSTAR[i][j][k] = (double*)malloc(5 * sizeof(double));
	      H_TOTAL[i][j][k] = (double*)malloc(5 * sizeof(double));
	    }
	}

    }
  byte_size += PARENT_RESOLUTION*PARENT_RESOLUTION*PARENT_RESOLUTION*sizeof(double) * 5 * 3;
  printf("Currently allocated %10.1f MB for simulation data.\n",byte_size/1e6);

  return 0;
}
