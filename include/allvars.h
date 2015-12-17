#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hllc_defs.h"


#ifndef PARENT_RESOLUTION
#define PARENT_RESOLUTION  pow(2, PARENT_MESH_LEVEL)
#endif

#ifndef parent_mesh
FILE * parent_mesh;
FILE * snapshot;
int snapshot_count;

// Primitive
double *** pressure;
double *** density;
double *** vx;
double *** vy;
double *** vz;
double *** internal_energy;

// Conserved
double *** C1;
double *** C2;
double *** C3;
double *** C4;
double *** C5;
//double *** C6;

double byte_size;
double CurrentTime;
double rho_L, ux_L, uy_L, uz_L, P_L, E_L;
double rho_R, ux_R, uy_R, uz_R,P_R, E_R;
double dt;
double dx;
double dy;
double dz;
#endif // parent_mesh

#ifdef HLLC_1
int ghost_cell_offset;
double S_max;
double S_star;
double S_L;
double S_R;
double pvrs;
double a_L;
double a_R;
double a_bar;
double rho_bar;
double p_bar;
double p_star;
double QL;
double QR;
double * U_RSTAR;
double * U_LSTAR;
double * U_R;
double * U_L;

// X-dir
double **** F_R; // right fan flux state
double **** F_RSTAR; // right star fan flux state
double **** F_LSTAR; // left star fan flux state
double **** F_L; // left fan flux state
double **** F_TOTAL;

// Y-dir
double **** G_R; // right fan flux state                                                                                                                            
double **** G_RSTAR; // right star fan flux state                                                                                                                   
double **** G_LSTAR; // left star fan flux state                                                                                                                    
double **** G_L; // left fan flux state                
double **** G_TOTAL;

// Z-dir
double **** H_R; // right fan flux state                                                                                                                            
double **** H_RSTAR; // right star fan flux state                                                                                                                   
double **** H_LSTAR; // left star fan flux state                                                                                                                    
double **** H_L; // left fan flux state                
double **** H_TOTAL;
#endif // HLLC

