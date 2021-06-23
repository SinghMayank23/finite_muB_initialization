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

using namespace std;

namespace Nuclei {

 void Initialize_T_P_properties(string nucleusname, nucleus_properties &tpproperties);
 void Initialize_projectile_nucleus(std::ifstream& nucleusfile, int nuclei_number, nucleon* Projectile, double sqrtsNN);
 void Initialize_target_nucleus(std::ifstream& nucleusfile, int nuclei_number, nucleon* Target, double sqrtsNN);
 void Apply_random_rotation(nucleon*Target, nucleon* Projectile, gsl_rng* random, int TargetA, int ProjectileA);
 void Apply_Lorentz_Contraction(nucleon* Target, nucleon* Projectile, double sqrtNN, int TargetA, int ProjectileA);
 void Shift_and_sort_nuclei(nucleon* Target, nucleon* Projectile, int TargetA, int ProjectileA);
 void sortincreasing(nucleon* A, int nucleonnum);
 void sortdecreasing(nucleon* A, int nucleonnum);
 void Sample_WS_nucleus(nucleon* Nucleus, double sqrtsNN, gsl_rng* random, nucleus_properties &properties, bool is_nucleus_going_right);
 bool Nucleon_is_close(nucleon* Nucleus, int inucleon, double x, double y, double z);
}


#endif
