#include<string>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<math.h>
#include<sstream>
#include<cmath>

using namespace std;

void sortdecreasing(double** A, int eventnum);

int main() {
  int number_of_events = 10000;
  double **energies = new double* [2];
  for (int idx = 0; idx < 2; idx++)
  {
    energies[idx] = new double [number_of_events];
  }

  for (int ievent = 0; ievent < number_of_events; ievent++)
  {
    string event_num = to_string(ievent);
    event_num = event_num.c_str();
    string filename1 = "gaussians_" + event_num + ".dat";
    string filename2 = "remnants_" + event_num + ".dat";
    ifstream gaussianfile(filename1.c_str());
    ifstream remnantsfile(filename2.c_str());
  
    string line;
    double dummy, energy;
    double total_energy = 0;

    while (getline(gaussianfile, line))
    {
      istringstream iss(line);
      iss >> dummy >> dummy >> dummy >> dummy >> energy >> dummy >> dummy;
      total_energy += energy;
    }
    gaussianfile.close();

    while (getline(remnantsfile, line))
    {
      istringstream iss(line);
      iss >> dummy >> dummy >> dummy >> dummy >> energy >> dummy;
      total_energy += energy;
    }
    remnantsfile.close();

    energies[0][ievent] = (double)ievent;
    energies[1][ievent] = total_energy;
  }

  sortdecreasing(energies, number_of_events);

  FILE* outfile  = fopen("Ordered_energies.dat", "w");
  FILE* outfile2 = fopen("change_ordering.sh", "w");

  for (int ievent = 0; ievent < number_of_events; ievent++)
  {
    fprintf(outfile , "%d %e \n", (int)energies[0][ievent], energies[1][ievent]);
    fprintf(outfile2, "%s%d%s%d%s \n","mv gaussians_", (int)energies[0][ievent], ".dat gaussians_", ievent, ".dat");
    fprintf(outfile2, "%s%d%s%d%s \n","mv remnants_", (int)energies[0][ievent], ".dat remnants_", ievent, ".dat");
  } 
  fclose(outfile);

  for(int idx = 0; idx < 2; idx++) delete [] energies[idx];
  delete [] energies;
}

void sortdecreasing(double** A, int eventnum)
{
  double temp1, temp2;
  for (int ievent = 0; ievent < eventnum; ievent++)
  {
    for (int jevent = ievent; jevent < eventnum; jevent++)
    {
      if (A[1][ievent] < A[1][jevent])
      {
        temp1 = A[0][jevent];
        temp2 = A[1][jevent];
        A[0][jevent] = A[0][ievent];
        A[1][jevent] = A[1][ievent];
        A[0][ievent] = temp1;
        A[1][ievent] = temp2;
      }
    } 
  }
  return;
}
