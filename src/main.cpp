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
#include "./impact.h"
#include "./propagation.h"

using namespace std;

int main(int argc, char *argv[]) {

  std::string sseed;
  if (argc > 1)
  {
    sseed = argv[1];
  } else 
  {
    sseed = "0";
  }

  bool fixed_impact_parameter = true;
  int nevents = 1;
  double sqrtsNN = 200.;
  double sigmaNN = 1.;
  double impact_parameter_b = 0.;
  double min_b = 0.;
  double max_b = 7.5;
  
  double rnum = time(0) + (double)(atoi(sseed.c_str()));

  gsl_rng* random = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(random,static_cast <unsigned long int> (rnum));

  int nuclei_number = 1;//set 0 for start reading from current line
  std::ifstream nucleusfile("../au197-sw-full_3Bchains-conf1820.dat");

  for (int ievents = 0; ievents < nevents; ievents++)
  {

    nucleon Target[197];
    nucleon Projectile[197];
    std::vector<std::shared_ptr<string_initial>> string_list     ;
    std::vector<std::shared_ptr<string_final>> string_final_list     ;

    Nuclei::Initialize_projectile_nucleus(nucleusfile, nuclei_number, Projectile);
    Nuclei::Initialize_target_nucleus(nucleusfile, nuclei_number, Target);

    FILE *fileout1, *fileout2;
    string outfile1 = "file1.dat";
    string outfile2 = "file2.dat";
    fileout1 = fopen(outfile1.c_str(), "w");
    fileout2 = fopen(outfile2.c_str(), "w");

    Nuclei::Apply_random_rotation(Target, Projectile, random);
    Nuclei::Apply_Lorentz_Contraction(Target, Projectile, sqrtsNN);
    for (int inucleon = 0; inucleon < 197; inucleon++)
    {
//      fprintf(fileout1, "%e %e %e \n", Target[inucleon].x, Target[inucleon].y, Target[inucleon].z);
      fprintf(fileout1, "%e %e %e \n", Projectile[inucleon].x, Projectile[inucleon].y, Projectile[inucleon].z);
    }
    Nuclei::Shift_nuclei(Target, Projectile);
    for (int inucleon = 0; inucleon < 197; inucleon++)
    {
//      fprintf(fileout2, "%e %e %e \n", Target[inucleon].x, Target[inucleon].y, Target[inucleon].z);
      fprintf(fileout2, "%e %e %e \n", Projectile[inucleon].x, Projectile[inucleon].y, Projectile[inucleon].z);
    }

    fflush(fileout1);
    fflush(fileout2);
    fclose(fileout1);
    fclose(fileout2);

    exit(1);


    if (! fixed_impact_parameter)
    impact_parameter_b = Impact::Get_impact_parameter(random, min_b, max_b);

    Impact::Shift_nuclei_transverse(impact_parameter_b, Target, Projectile);

    Binary_Collisions::Determine_strings(string_list, Target, Projectile, sqrtsNN, sigmaNN, random);

    Propagation::Propagate_strings(string_list, string_final_list, sqrtsNN, random);
  }

  gsl_rng_free(random);

}
