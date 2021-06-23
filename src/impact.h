#ifndef _SRC_IMPACT_H_
#define _SRC_IMPACT_H_

#include<string>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<math.h>
#include<sstream>
#include<cmath>
#include<vector>
#include<memory>
#include<gsl/gsl_rng.h>
#include "./globals.h"
#include "./struct.h"

using namespace std;

namespace Impact {

 double Get_impact_parameter(gsl_rng* random, double min_b, double max_b);
 void Shift_nuclei_transverse(double impact_parameter_b, nucleon* Target, nucleon* Projectile, int TargetA, int ProjectileA);
}


#endif
