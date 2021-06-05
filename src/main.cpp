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
#include "./write_files.h"

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

  int model_choice = 1;
  bool fixed_impact_parameter = false;
  bool save_nucleon_pos_and_bin_collisions = true;
  int nevents = 10;
  double sqrtsNN = 200.;
  double sigmaNN = 4.2;//4.2 fm^2 (for 200 GeV)
  sigmaNN *= 0.6;
  double impact_parameter_b = 7.;
  double min_b = 5.;
  double max_b = 7.5;
  
  double rnum = time(0) + (double)(atoi(sseed.c_str()));
  rnum = 20.;

  gsl_rng* random = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(random,static_cast <unsigned long int> (rnum));

  int nuclei_number = 1;//set 0 for start reading from current line
  std::ifstream nucleusfile("../au197-sw-full_3Bchains-conf1820.dat");

   FILE *fileout1, *fileout2;

  for (int ievents = 0; ievents < nevents; ievents++)
  {


    nucleon Target[197];
    nucleon Projectile[197];

    Nuclei::Initialize_projectile_nucleus(nucleusfile, nuclei_number, Projectile, sqrtsNN);
    Nuclei::Initialize_target_nucleus(nucleusfile, nuclei_number, Target, sqrtsNN);

    Nuclei::Apply_random_rotation(Target, Projectile, random);
    Nuclei::Apply_Lorentz_Contraction(Target, Projectile, sqrtsNN);
    Nuclei::Shift_and_sort_nuclei(Target, Projectile);

    if (! fixed_impact_parameter)
    impact_parameter_b = Impact::Get_impact_parameter(random, min_b, max_b);

    Impact::Shift_nuclei_transverse(impact_parameter_b, Target, Projectile);

    cout << "Impact_parameter = " << impact_parameter_b << endl;

    if (model_choice == 1)
    {
      string event_num = to_string(ievents);
      event_num = event_num.c_str();
      string outfile1 = "gaussians_" + event_num + ".dat";
      string outfile2 = "remnants_" + event_num + ".dat";
      fileout1 = fopen(outfile1.c_str(), "w");
      fileout2 = fopen(outfile2.c_str(), "w");

      auto binary_list = std::vector<std::shared_ptr<binary_coll>>();
      auto gaussian_list = std::vector<std::shared_ptr<gaussians>>();
      auto remnant_list  = std::vector<std::shared_ptr<remnant>>();
      Binary_Collisions::Determine_all_collisions(binary_list, Target, Projectile, sqrtsNN, sigmaNN, random);
      Propagation::Propagate_gaussians(binary_list, gaussian_list);
      Propagation::Propagate_remnants(Target, Projectile, remnant_list, sqrtsNN);
      Write_files::Write_gaussians(gaussian_list, fileout1);
      Write_files::Write_remnants(remnant_list, fileout2);

      fflush(fileout1);
      fflush(fileout2);
      fclose(fileout1);
      fclose(fileout2);
    }
    else if (model_choice == 2)
    {
      string event_num = to_string(ievents);
      event_num = event_num.c_str();
      string outfile1 = "n_positions_and_binary_collisions_" + event_num + ".dat";
      string outfile2 = "initial_" + event_num + ".dat";
      if (save_nucleon_pos_and_bin_collisions)
      fileout1 = fopen(outfile1.c_str(), "w");
      fileout2 = fopen(outfile2.c_str(), "w");

      auto binary_list = std::vector<std::shared_ptr<binary_coll>>();
      auto string_final_list = std::vector<std::shared_ptr<string_final>>();
      Binary_Collisions::Determine_strings(binary_list, Target, Projectile, sqrtsNN, sigmaNN, random);
  
      if (save_nucleon_pos_and_bin_collisions) Write_files::Write_nucleon_positions_and_binary_collisions(Target, Projectile, binary_list, fileout1);
  
      Propagation::Propagate_strings(binary_list, string_final_list, sqrtsNN, random);
      double total_energy = Propagation::Calculate_total_final_energy(string_final_list);
      if (save_nucleon_pos_and_bin_collisions)
      {
        fflush(fileout1);
        fclose(fileout1);
      }
      fflush(fileout2);
      fclose(fileout2);
    }

  }

  gsl_rng_free(random);

}
