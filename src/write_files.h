#ifndef _SRC_WRITE_FILES_H_
#define _SRC_WRITE_FILES_H_

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

namespace Write_files {

 void Write_nucleon_positions_and_binary_collisions(nucleon* Target, nucleon* Projectile, std::vector<std::shared_ptr<binary_coll>>& string_list, FILE *fileout1, 
                                                    int TargetA, int ProjectileA);
 void Write_final_terms(std::vector<std::shared_ptr<string_final>>& string_final_list, FILE *fileout2, double total_energy);
 void Write_gaussians(std::vector<std::shared_ptr<gaussians>>& gaussian_list, FILE *fileout1);
 void Write_remnants(std::vector<std::shared_ptr<remnant>>& remnant_list, FILE *fileout2);

} 


#endif
