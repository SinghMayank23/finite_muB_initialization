#ifndef _SRC_NUCLEI_H_
#define _SRC_NUCLEI_H_

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
#include<algorithm>
#include<array>
#include<gsl/gsl_rng.h>
#include "./globals.h"
#include "./struct.h"

namespace Nuclei {

 void Initialize_projectile_nucleus(std::ifstream& nucleusfile, int nuclei_number, nucleon* Projectile, double sqrtsNN);
 void Initialize_target_nucleus(std::ifstream& nucleusfile, int nuclei_number, nucleon* Target, double sqrtsNN);
 void Apply_random_rotation(nucleon*Target, nucleon* Projectile, gsl_rng* random);
 void Apply_Lorentz_Contraction(nucleon* Target, nucleon* Projectile, double sqrtNN);
 void Shift_and_sort_nuclei(nucleon* Target, nucleon* Projectile);
 void sortincreasing(nucleon* A);
 void sortdecreasing(nucleon* A);
}


#endif
