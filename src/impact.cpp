#include "./impact.h"

using namespace std;

namespace Impact {

 double Get_impact_parameter(gsl_rng* random, double min_b, double max_b)
 {
   double impact_parameter = gsl_rng_uniform(random)*(max_b*max_b - min_b*min_b);
   impact_parameter +=  min_b*min_b;
   impact_parameter  = sqrt(impact_parameter);

   return impact_parameter;
 }

 void Shift_nuclei_transverse(double impact_parameter_b, nucleon* Target, nucleon* Projectile)
 {
   for (int inucleon = 0; inucleon < 197; inucleon++)
   {
     Target[inucleon].x     += impact_parameter_b/2.;
     Projectile[inucleon].x -= impact_parameter_b/2.;
   }
   return;
 }

} 

