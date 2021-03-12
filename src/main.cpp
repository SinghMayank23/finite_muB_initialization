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

  bool fixed_impact_parameter = false;
  bool save_nucleon_pos_and_bin_collisions = true;
  int nevents = 1;
  double sqrtsNN = 200.;
  double sigmaNN = 0.5;//4.2 fm^2 (for 200 GeV)
  double impact_parameter_b = 0.;
  double min_b = 5.;
  double max_b = 7.5;
  
  double rnum = time(0) + (double)(atoi(sseed.c_str()));

  gsl_rng* random = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(random,static_cast <unsigned long int> (rnum));

  int nuclei_number = 1;//set 0 for start reading from current line
  std::ifstream nucleusfile("../au197-sw-full_3Bchains-conf1820.dat");

   FILE *fileout1, *fileout2;

  for (int ievents = 0; ievents < nevents; ievents++)
  {

    string event_num = to_string(ievent);
    string outfile1 = "n_positions_and_binary_collisions_" + event_num.c_str() + ".dat";
    string outfile2 = "initial_" + event_num.c_str() + ".dat";
    if (save_nucleon_pos_and_bin_collisions)
    fileout1 = fopen(outfile1.c_str(), "w");
    fileout2 = fopen(outfile2.c_str(), "w");

    nucleon Target[197];
    nucleon Projectile[197];
    auto string_list = std::vector<std::shared_ptr<string_initial>>();
    auto string_final_list = std::vector<std::shared_ptr<string_final>>();

    Nuclei::Initialize_projectile_nucleus(nucleusfile, nuclei_number, Projectile);
    Nuclei::Initialize_target_nucleus(nucleusfile, nuclei_number, Target);

    Nuclei::Apply_random_rotation(Target, Projectile, random);
    Nuclei::Apply_Lorentz_Contraction(Target, Projectile, sqrtsNN);
    Nuclei::Shift_nuclei(Target, Projectile);

    if (! fixed_impact_parameter)
    impact_parameter_b = Impact::Get_impact_parameter(random, min_b, max_b);

    Impact::Shift_nuclei_transverse(impact_parameter_b, Target, Projectile);

    cout << "Impact_parameter = " << impact_parameter_b << endl;

    Binary_Collisions::Determine_strings(string_list, Target, Projectile, sqrtsNN, sigmaNN, random);

    if (save_nucleon_pos_and_bin_collisions) Write_files::Write_nucleon_positions_and_binary_collisions(Target, Projectile, string_list, fileout1);

    Propagation::Propagate_strings(string_list, string_final_list, sqrtsNN, random);
    double total_energy = Propagation::Calculate_total_final_energy(string_final_list);

    if (save_nucleon_pos_and_bin_collisions)
    {
      fflush(fileout1);
      fclose(fileout1);
    }
    fflush(fileout2);
    fclose(fileout2);
  }

  gsl_rng_free(random);

}
