
/*

  Function prototypes.
 
*/

// mesh.c
int read_parent_mesh();
int advanceTimestep_HLLC1();

// flux.c
int initialize_flux_arrays();

// formulation.c
int toConservative(double***, double ***, double ***, double ***, double ***, double ***,\
		   double***, double ***, double ***, double ***, double ***);
int toPrimitive(double***, double ***, double ***, double ***, double ***, double ***,\
		double***, double ***, double ***, double ***, double ***);

// eos.c
double internal_energy_ideal_gas(double,double,double);
double internal_energy_covolume_gas(double, double, double, double);
double pressure_ideal_gas(double,double,double);
double pressure_covolume_gas(double, double, double, double);
double sound_speed_ideal_gas(double, double, double);

// bcs.c

int set_HLLC_1_transmissive(double***, double ***, double ***, double ***, double ***, double ***);
int set_HLLC_1_periodic(double***, double ***, double ***, double ***, double ***, double ***);
