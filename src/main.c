#include <math.h>
#include <stdio.h>
#include "../include/hllc_defs.h"
#include "../include/allvars.h"
#include "../include/proto.h"

int main()
{

#ifdef MPICH
#endif // MPICH

#ifdef OPENCL
#endif // OPENCL

#if !defined(OPENCL) && !defined(MPICH)
  // Serial mesh initialization
  byte_size = 0;
  CurrentTime = 0.;
  printf("HLLC3 is initializing mesh...\n");
  printf("Parent mesh is (%d X %d X %d).\n", (int)PARENT_RESOLUTION, (int)PARENT_RESOLUTION, (int)PARENT_RESOLUTION);
  
  read_parent_mesh();

  initialize_flux_arrays();

#ifdef HLLC_1

  printf("Method chosen: 1st order HLLC (TVD)\n");

#ifdef TRANSMISSIVE
  printf("Transmissive boundaries selected.\n");
  set_HLLC_1_transmissive(density, vx,vy,vz, pressure,internal_energy,C1,C2,C3,C4,C5);
#endif // TRANSMISSIVE

  // Convert to conservative variables for calculation                                                                                                               
  toConservative(density, vx,vy,vz,pressure,internal_energy,
		 C1, C2, C3, C4, C5);

#ifdef TRANSMISSIVE
  set_HLLC_1_transmissive(density, vx,vy,vz, pressure,internal_energy,C1,C2,C3,C4,C5);
#endif // TRANSMISSIVE
  printf("Now working in conservative form.\n");
  HLLC1();

#endif // HLLC_1
#endif //!defined(OPENCL) && !defined(MPICH)

  return 0;
}
