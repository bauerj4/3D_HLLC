#include "../include/hllc_defs.h"
#include "../include/hllc_defs.h"
#include "../include/allvars.h"
#include "../include/proto.h"
#include <stdlib.h>
#include <math.h>

int compute_X_STAR_state_vars(int i, int j, int k)
{
  P_R = pressure[i+1][j+1][k+1];
  P_L = pressure[i][j+1][k+1];
  rho_R =density[i+1][j+1][k+1];
  rho_L =density[i][j+1][k+1];
  ux_R =vx[i+1][j+1][k+1];
  ux_L =vx[i][j+1][k+1];


  rho_bar = 0.5 * (rho_L + rho_R);
  a_L = sound_speed_ideal_gas(P_L, rho_L, GAMMA);
  a_R =sound_speed_ideal_gas(P_R, rho_R, GAMMA);
  //printf("rho_bar, a_L, a_R = %f, %f, %f\n", rho_bar, a_L, a_R);
  a_bar = 0.5 *(a_L + a_R);
  p_bar = 0.5 * (P_L + P_R);
  pvrs = p_bar - 0.5 * (ux_R - ux_L)*rho_bar*a_bar;
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
  
  if (p_star * p_star > P_L * P_L)
    {
      QL = 1 + (GAMMA + 1)/(2 * GAMMA) *\
	(p_star/P_L - 1);
      QL = sqrt(QL);
    }
  if (p_star*p_star > P_R * P_R)
    {
      QR = 1 + (GAMMA +1)/(2 * GAMMA) *		\
	(p_star/P_R - 1);
      QR = sqrt(QR);
    }
  
  S_L = ux_L - (a_L * QL);
  S_R = ux_R + (a_R * QR);

#endif

#ifdef WAVESPEED_SIMPLE
  S_L = ux_L - a_L;
  if (ux_R - a_R < S_L)
    {
      S_L = ux_R - a_R;
    }

  S_R =ux_L + a_L;
  if (ux_R - a_R > S_R)
    {
      S_R = ux_R + a_R;
    }

#endif // WAVESPEED_SIMPLE
  
  
  double numerator = P_R - P_L + (rho_L * ux_L * (S_L - ux_L)) - (rho_R * ux_R * (S_R - ux_R));
  double denominator = (rho_L * (S_L - ux_L)) - (rho_R * (S_R - ux_R));

  S_star = numerator/denominator;
  /*
  if (denominator != denominator)
    {
      printf("rho_L, rho_R, ux_L, rho_R, S_R, ux_R = %f, %f, %f, %f, %f, %f\n", rho_L, rho_R, ux_L, rho_R, S_R, ux_R);
      }*/
  /*
  S_star = pressure[i+1][j+1][k+1] - pressure[i][j+1][k+1] + density[i][j+1][k+1] * vx[i][j+1][k+1] *\
    (S_L - vx[i][j+1][k+1]) - density[i+1][j+1][k+1] * vx[i+1][j+1][k+1] * (S_R - vx[i+1][j+1][k+1]);
  S_star /= (density[i][j+1][k+1]*(S_L - vx[i][j+1][k+1]) - density[i+1][j+1][k+1] * (S_R - vx[i+1][j+1][k+1]));
  */
  //printf("S_L. S_R, S_star = %f, %f, %f\n", S_L, S_R, S_star);
  return 0;
}

int compute_USTAR_XL(double * U_STAR,int i, int j, int k)
{
  P_R = pressure[i+1][j+1][k+1];
  P_L = pressure[i][j+1][k+1];
  rho_R =density[i+1][j+1][k+1];
  rho_L =density[i][j+1][k+1];
  ux_R =vx[i+1][j+1][k+1];
  ux_L =vx[i][j+1][k+1];
  uy_R =vy[i+1][j+1][k+1];
  uy_L =vy[i][j+1][k+1];
  uz_R =vz[i+1][j+1][k+1];
  uz_L =vz[i][j+1][k+1];

  double uL_dot_uL = ux_L*ux_L + uy_L*uy_L + uz_L * uz_L;
  double uR_dot_uR = ux_R*ux_R + uy_R*uy_R + uz_R * uz_R;

  E_L = rho_L *(P_L / (rho_L*(GAMMA - 1)) + 0.5 *uL_dot_uL);
  E_R = rho_R *(P_R / (rho_R*(GAMMA - 1)) + 0.5 * uR_dot_uR);

  double coeff = rho_L * (S_L - ux_L)/(S_L - S_star);
  U_STAR[0] = coeff;
  U_STAR[1] = coeff * S_star;
  U_STAR[2] = coeff * vy[i][j+1][k+1];
  U_STAR[3] = coeff * vz[i][j+1][k+1];
  double fn = ((E_R / rho_R) + (S_star - ux_R) * (S_star + (P_R)/(rho_R*(S_R - ux_R))));
    //(C5[i][j+1][k+1] / rho_L) + (S_star - ux_L) * (S_star + P_L/(rho_L * (S_L - ux_L)));
  U_STAR[4] = coeff * fn;

  //printf("S_L, S_star = %f,%f\n", S_L, S_star);
  return 0;
}

int compute_USTAR_XR(double * U_STAR,int i, int j, int k)
{

  P_R = pressure[i+1][j+1][k+1];
  P_L = pressure[i][j+1][k+1];
  rho_R =density[i+1][j+1][k+1];
  rho_L =density[i][j+1][k+1];
  ux_R =vx[i+1][j+1][k+1];
  ux_L =vx[i][j+1][k+1];
  uy_R =vy[i+1][j+1][k+1];
  uy_L =vy[i][j+1][k+1];
  uz_R =vz[i+1][j+1][k+1];
  uz_L =vz[i][j+1][k+1];


  double uL_dot_uL = ux_L*ux_L + uy_L*uy_L + uz_L * uz_L;
  double uR_dot_uR = ux_R*ux_R + uy_R*uy_R + uz_R * uz_R;

  E_L = rho_L *(P_L / (rho_L*(GAMMA - 1)) + 0.5 *uL_dot_uL);
  E_R = rho_R *(P_R / (rho_R*(GAMMA - 1)) + 0.5 * uR_dot_uR);

  double coeff = rho_R * (S_R - ux_R)/(S_R - S_star);
  U_STAR[0] = coeff;
  U_STAR[1] = coeff * S_star;
  U_STAR[2] = coeff * vy[i+1][j+1][k+1];
  U_STAR[3] = coeff * vz[i+1][j+1][k+1];
  double fn = ((E_R / rho_R) + (S_star - ux_R) * (S_star + (P_R)/(rho_R*(S_R - ux_R))));
    //(C5[i+1][j+1][k+1] / rho_R) + (S_star - ux_R) * (S_star + P_R / (rho_R * (S_R - ux_R)));
  U_STAR[4] = coeff * fn;
  return 0;
}


int compute_X_HLLC_Fluxes(double *** density, double *** vx, double *** vy, double *** vz,\
                     double *** pressure, double *** internal_energy, double *** C1,\
		     double *** C2, double *** C3, double *** C4, double *** C5,\
                     double **** F_R, double **** F_RSTAR, double **** F_L, double **** F_LSTAR)
{
  int i,j,k;
  /*
  P_R = pressure[i+1][j+1][k+1];
  P_L = pressure[i][j+1][k+1];
  rho_R =density[i+1][j+1][k+1];
  rho_L =density[i][j+1][k+1];
  ux_R =vx[i+1][j+1][k+1];
  ux_L =vx[i][j+1][k+1];
  */
  S_max = 0;
  for (i=0; i < (int)(PARENT_RESOLUTION + ghost_cell_offset); i++)
    for (j=0; j < (int)(PARENT_RESOLUTION + ghost_cell_offset); j++)
      for (k=0; k < (int)(PARENT_RESOLUTION + ghost_cell_offset); k++)
	{
	  compute_X_STAR_state_vars(i,j,k);
	  P_R = pressure[i+1][j+1][k+1];
	  P_L = pressure[i][j+1][k+1];
	  rho_R =density[i+1][j+1][k+1];
	  rho_L =density[i][j+1][k+1];
	  ux_R =vx[i+1][j+1][k+1];
	  ux_L =vx[i][j+1][k+1];
          uy_R =vy[i+1][j+1][k+1];
          uy_L =vy[i][j+1][k+1];
          uz_R =vz[i+1][j+1][k+1];
          uz_L =vz[i][j+1][k+1];
	  
	  double uL_dot_uL = ux_L*ux_L + uy_L*uy_L + uz_L * uz_L;
          double uR_dot_uR = ux_R*ux_R + uy_R*uy_R + uz_R * uz_R;

	  E_L = rho_L *(P_L / (rho_L*(GAMMA - 1)) + 0.5 *uL_dot_uL);
	  E_R = rho_R *(P_R / (rho_R*(GAMMA - 1)) + 0.5 * uR_dot_uR);

	  //printf("S_star = %10.10f\n",S_star);

	  if (fabs(S_star) > S_max)
	    {
	      //printf("S_star = %10.10f\n",S_star);
	      S_max = S_star;
	    }
	  U_L[0] = C1[i][j+1][k+1];
	  U_L[1] = C2[i][j+1][k+1];
	  U_L[2] = C3[i][j+1][k+1];
	  U_L[3] = C4[i][j+1][k+1];
	  U_L[4] = C5[i][j+1][k+1];
	  
	  U_R[0] = C1[i+1][j+1][k+1];
	  U_R[1] = C2[i+1][j+1][k+1];
	  U_R[2] = C3[i+1][j+1][k+1];
	  U_R[3] = C4[i+1][j+1][k+1];
	  U_R[4] = C5[i+1][j+1][k+1];
	  
	  //if (U_L[2] != 0 || U_L[3] != 0)
	    {
	      /*
	      printf("U_L = [%f, %f, %f, %f, %f]\n", U_L[0], U_L[1], U_L[2], U_L[3], U_L[4]);
	      printf("U_R = [%f, %f, %f, %f, %f]\n", U_R[0], U_R[1], U_R[2], U_R[3], U_R[4]);
	      printf("U_RSTAR = [%f, %f, %f, %f, %f]\n", U_RSTAR[0], U_RSTAR[1], U_RSTAR[2], U_RSTAR[3], U_RSTAR[4]); 
	      printf("U_LSTAR = [%f, %f, %f, %f, %f]\n", U_LSTAR[0], U_LSTAR[1], U_LSTAR[2], U_LSTAR[3], U_LSTAR[4]); 
	      */
	    }
	  compute_USTAR_XL(U_LSTAR,i, j,k);
	  compute_USTAR_XR(U_RSTAR,i, j,k);
	  
	  if (C1[i][j+1][k+1] == 0)
	    {
	      printf("Warning: %f density\n", C1[i][j][k]);
	    }

	  if (density[i][j+1][k+1] == 0)
	    {
	      printf("WARNING: 0 DENSITY at %d %d %d\n",i,j,k);
	    }
	  if (pressure[i][j+1][k+1] == 0)
	    {
	      printf("WARNING: 0 PRESSURE\n");
	    }
	  if(internal_energy[i][j+1][k+1] == 0)
	    {
	      printf("WARNING: 0 ENERGY\n");
	    }


	  F_R[i][j][k][0] = rho_R * ux_R;
	  F_R[i][j][k][1] = rho_R * ux_R*ux_R + P_R;
	  F_R[i][j][k][2] = rho_R * ux_R * uy_R;
	  F_R[i][j][k][3] = rho_R * ux_R * uz_R;
	  F_R[i][j][k][4] = ux_R * (E_R + P_R);

          F_L[i][j][k][0] = rho_L * ux_L;
          F_L[i][j][k][1] = rho_L * ux_L*ux_L + P_L;
          F_L[i][j][k][2] = rho_L * ux_L * uy_L;
          F_L[i][j][k][3] = rho_L * ux_L * uz_L;
          F_L[i][j][k][4] = ux_L * (E_L + P_L);

	  
	  F_LSTAR[i][j][k][0] = F_L[i][j][k][0] + S_L * (U_LSTAR[0] - U_L[0]);
          F_LSTAR[i][j][k][1] = F_L[i][j][k][1] + S_L * (U_LSTAR[1] - U_L[1]);
          F_LSTAR[i][j][k][2] = F_L[i][j][k][2] + S_L * (U_LSTAR[2] - U_L[2]);
          F_LSTAR[i][j][k][3] = F_L[i][j][k][3] + S_L * (U_LSTAR[3] - U_L[3]);
          F_LSTAR[i][j][k][4] = F_L[i][j][k][4] + S_L * (U_LSTAR[4] - U_L[4]);

          F_RSTAR[i][j][k][0] = F_R[i][j][k][0] + S_R * (U_RSTAR[0] - U_R[0]);
          F_RSTAR[i][j][k][1] = F_R[i][j][k][1] + S_R * (U_RSTAR[1] - U_R[1]);
          F_RSTAR[i][j][k][2] = F_R[i][j][k][2] + S_R * (U_RSTAR[2] - U_R[2]);
          F_RSTAR[i][j][k][3] = F_R[i][j][k][3] + S_R * (U_RSTAR[3] - U_R[3]);
          F_RSTAR[i][j][k][4] = F_R[i][j][k][4] + S_R * (U_RSTAR[4] - U_R[4]);
	  
#ifdef HLLC_1
	  if (S_L >= 0)
	    {
	      //printf("Condition 1.\n");
	      F_TOTAL[i][j][k][0] = F_L[i][j][k][0];
	      F_TOTAL[i][j][k][1] = F_L[i][j][k][1];
              F_TOTAL[i][j][k][2] = F_L[i][j][k][2];
              F_TOTAL[i][j][k][3] = F_L[i][j][k][3];
              F_TOTAL[i][j][k][4] = F_L[i][j][k][4];
	    }
	  else if(S_L <= 0 && S_star >= 0)
	    {
	      //printf("Condition 2.\n");
	      F_TOTAL[i][j][k][0] = F_LSTAR[i][j][k][0];
	      F_TOTAL[i][j][k][1] = F_LSTAR[i][j][k][1];
              F_TOTAL[i][j][k][2] = F_LSTAR[i][j][k][2];
              F_TOTAL[i][j][k][3] = F_LSTAR[i][j][k][3];
              F_TOTAL[i][j][k][4] = F_LSTAR[i][j][k][4];
	    }
	  else if (S_star <= 0 && S_R >= 0)
	    {
	      //printf("Condition 3.\n");
	      F_TOTAL[i][j][k][0] = F_RSTAR[i][j][k][0];
	      F_TOTAL[i][j][k][1] = F_RSTAR[i][j][k][1];
	      F_TOTAL[i][j][k][2] = F_RSTAR[i][j][k][2];
	      F_TOTAL[i][j][k][3] = F_RSTAR[i][j][k][3];
	      F_TOTAL[i][j][k][4] = F_RSTAR[i][j][k][4];
	    }
	  else if (S_R <= 0)
	    {
	      //printf("Condition 4.\n");
	      F_TOTAL[i][j][k][0] = F_R[i][j][k][0];
              F_TOTAL[i][j][k][1] = F_R[i][j][k][1];
              F_TOTAL[i][j][k][2] = F_R[i][j][k][2];
              F_TOTAL[i][j][k][3] = F_R[i][j][k][3];
              F_TOTAL[i][j][k][4] = F_R[i][j][k][4];
	    }
	  
	  /*
	  printf("F_L = [%f, %f, %f, %f, %f]\n", F_L[i][j][k][0], F_L[i][j][k][1], F_L[i][j][k][2], F_L[i][j][k][3], F_L[i][j][k][4]);
	  printf("F_R = [%f, %f, %f, %f, %f]\n", F_R[i][j][k][0], F_R[i][j][k][1], F_R[i][j][k][2], F_R[i][j][k][3], F_R[i][j][k][4]);
	  printf("F_RSTAR = [%f, %f, %f, %f, %f]\n", F_RSTAR[i][j][k][0], F_RSTAR[i][j][k][1], F_RSTAR[i][j][k][2],\
		 F_RSTAR[i][j][k][3], F_RSTAR[i][j][k][4]);
	  printf("F_LSTAR = [%f, %f, %f, %f, %f]\n", F_LSTAR[i][j][k][0], F_LSTAR[i][j][k][1], F_LSTAR[i][j][k][2],\
	  F_LSTAR[i][j][k][3], F_LSTAR[i][j][k][4]);*/

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
