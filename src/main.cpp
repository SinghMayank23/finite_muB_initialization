//© Mayank Singh
// Wounded nucleon string initial state model based on work by Schenke and Shen (arxiv:1710.00881)
#include<string>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<math.h>
#include<sstream>
#include<cmath>
#include<gsl/gsl_rng.h>
#include "./struct.h"
#include "./globals.h"
#include "./binary_collisions.h"
#include "./nuclei.h"
#include "./propagation.h"

using namespace std;

int main(int argc, char *argv[]) {

  std::string sseed;
  if (argc > 1)
  {
    sseed = argv[1];
  } else {
    cout << "itau not defined" << endl;
    exit(1);
  }


  int nevents = 1;
  double sqrtsNN = 200.;
  double sigmaNN = 1.;
  
  double rnum = time(0) + (double)(atoi(sseed.c_str()));

  gsl_rng* random = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(random,static_cast <unsigned long int> (rnum));

  int nuclei_number = 1;//set 0 for start reading from current line
  std::ifstream nucleusfile("../Nuclei_DATA/pb208-01.dat");

  for (int ievents = 0; ievents < nevents; ievents++)
  {

    nucleon Target[197];
    nucleon Projectile[197];
    std::vector<std::shared_ptr<string_initial>> string_list     ;
    std::vector<std::shared_ptr<string_final>> string_final_list     ;

    Nuclei::Initialize_projectile_nucleus(nucleusfile, nuclei_number, Projectile);
    Nuclei::Initialize_target_nucleus(nucleusfile, nuclei_number, Target);

    Nuclei::Apply_random_rotation(Target, Projectile, random);
    Nuclei::Apply_Lorentz_Contraction(Target, Projectile, sqrtsNN);
    Nuclei::Shift_nuclei(Target, Projectile);

//    Impact::Get_impact_parameter();
//    Impact::Shift_nuclei_transverse();

    Binary_Collisions::Determine_strings(string_list, Target, Projectile, sqrtsNN, sigmaNN, random);

    Propagation::Propagate_strings(string_list, string_final_list, sqrtsNN, random);
  }

  gsl_rng_free(random);

}
