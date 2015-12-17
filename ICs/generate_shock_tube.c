#include <stdio.h>
#include <math.h>
/*
int LEVEL = 4;
int PARENT_RESOLUTION = (int)pow((pow(2,LEVEL)), 3);
double rho_l = 1.0;
double ux_l = 0.75;
double uy_l = 0.;
double uz_l = 0.;
double p_l = 1.0;
double rho_r = 0.125;
double ux_r = 0.;
double uy_r = 0.;
double uz_r = 0.;
double p_r = 0.1;
double barrier = (int)(pow(2,LEVEL))/4;

FILE * parent_mesh = open("parent_mesh.dat",'w');
*/
int main()
{
  int LEVEL = 7;
  int PARENT_RESOLUTION = (int)pow((pow(2,LEVEL)),1);
  double rho_l = 1.0;
  double ux_l = 0.75;
  double uy_l = 0.;
  double uz_l = 0.;
  double p_l = 1.0;
  double rho_r = 0.125;
  double ux_r = 0.;
  double uy_r = 0.;
  double uz_r = 0.;
  double p_r = 0.1;
  double barrier = (int)(pow(2,LEVEL))/4;

  FILE * parent_mesh = fopen("parent_mesh.dat","wb");


  int i,j,k;
  printf("Writing IC file on level %d...\n", LEVEL);
  for (i=0; i<PARENT_RESOLUTION; i++)
    for (j=0; j<PARENT_RESOLUTION; j++)
      for (k=0; k<PARENT_RESOLUTION; k++)
	{
	  if (i < barrier)
	    {	      
	      fwrite(&rho_l, sizeof(double), 1, parent_mesh);
	      fwrite(&ux_l, sizeof(double), 1,  parent_mesh);
	      fwrite(&uy_l, sizeof(double), 1,  parent_mesh);
	      fwrite(&uz_l, sizeof(double), 1,  parent_mesh);
	      fwrite(&p_l, sizeof(double), 1,  parent_mesh);

	    }
	  else
	    {
              fwrite(&rho_r, sizeof(double), 1, parent_mesh);
              fwrite(&ux_r, sizeof(double), 1,  parent_mesh);
              fwrite(&uy_r, sizeof(double), 1,  parent_mesh);
              fwrite(&uz_r, sizeof(double), 1,  parent_mesh);
              fwrite(&p_r, sizeof(double), 1,  parent_mesh);

	    }
	}
  printf("ICs written.\n");
  return 0;
}
