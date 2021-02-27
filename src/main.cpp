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

int main(int argc, char *argv[]) {

  string sseed;
  if (argc > 1)
  {
    sseed = argv[1];
  } else {
    cout << "itau not defined" << endl;
    exit(1);
  }

  gsl_rng* random = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(random,static_cast <unsigned long int> (rum));

  int nuclei_number = 1;//set 0 for start reading from current line
  ifstream nucleusfile("../Nuclei_DATA/pb208-01.dat")

  for (ievents = 0; ievents < nevents; ievents++)
  {
    Impact::Get_impact_parameter();

    nucleon Target[197];
    nucleon Projectile[197];
    Nuclei::Get_projectile_nuleus(nucleusfile, nuclei_number, Projectile);
    Nuclei::Get_target_nucleus(nucleusfile, nuclei_number, Target);

    Nuclei::Apply_random_rotation(Target, Projectile, random);
    Nuclei::Apply_Lorentz_Contraction(Target, Projectile, sNN);
  }

  gsl_rng_free(random);

}
