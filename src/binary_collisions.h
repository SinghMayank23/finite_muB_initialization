#ifndef _SRC_BINARY_COLLISIONS_H_
#define _SRC_BINARY_COLLISIONS_H_

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
#include<gsl/gsl_rng.h>
#include "./struct.h"
#include "./globals.h"
#include "./propagation.h"

namespace Binary_Collisions {

 void Determine_strings(std::vector<std::shared_ptr<binary_coll>>& binary_list, nucleon* Target, nucleon* Projectile,
                                 double sqrtsNN, double sigmaNN, gsl_rng* random);
 bool sort_cond(std::shared_ptr<binary_coll> a, std::shared_ptr<binary_coll> b);
 void Determine_all_collisions(std::vector<std::shared_ptr<binary_coll>>& binary_list, nucleon* Target, nucleon* Projectile,
                                 double sqrtsNN, double sigmaNN, gsl_rng* random);
} 


#endif
