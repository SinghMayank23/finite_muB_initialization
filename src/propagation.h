#ifndef _SRC_PROPAGATION_H_
#define _SRC_PROPAGATION_H_

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
#include "./struct.h"
#include "./globals.h"

namespace Propagation {

 void Propagate_strings(std::vector<std::shared_ptr<binary_coll>>& string_list, std::vector<std::shared_ptr<string_final>>& string_final_list,
                                 double sqrtsNN, gsl_rng* random);
 double Sample_y_loss(double yinT, double yinP, gsl_rng* random);
 double Calculate_total_final_energy(std::vector<std::shared_ptr<string_final>>& string_final_list);
 void Propagate_gaussians(std::vector<std::shared_ptr<binary_coll>>& binary_list, std::vector<std::shared_ptr<gaussians>>& gaussian_list,double sqrtsNN);
 void Propagate_remnants(nucleon* Target, nucleon* Projectile, std::vector<std::shared_ptr<remnant>>& remnant_list, double sqrtsNN);
} 


#endif
