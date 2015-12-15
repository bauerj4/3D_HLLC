#include "../include/allvars.h"
#include "../include/hllc_defs.h"
#include "../include/proto.h"
#include <stdlib.h>
#include <math.h>

double internal_energy_ideal_gas(double p, double rho, double gamma)
{
  double e = p / ((gamma - 1)  * rho);
  return e;
}

double internal_energy_covolume_gas(double p, double rho, double b, double gamma)
{
  double e = p * (1 - b * rho) / ((gamma - 1) * rho);
  return e;
}

double pressure_ideal_gas(double e, double rho, double gamma)
{
  double p = e * ((gamma - 1) * rho);
  return p;
}

double pressure_covolume_gas(double e, double rho, double b, double gamma)
{
  double p = e * ((gamma - 1) * rho)/(1 - b*rho);
  return p;
}
