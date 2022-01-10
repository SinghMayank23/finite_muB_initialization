#include<string>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<math.h>
#include<sstream>
#include<cmath>

using namespace std;

int main() {
  int number_of_events = 10000;
  double eta_min = -0.5;
  double eta_max =  0.5;

  FILE* outfile  = fopen("Epsilon_n.dat", "w");

  for (int ievent = 0; ievent < number_of_events; ievent++)
  {
    string line;
    string event_num = to_string(ievent);
    event_num = event_num.c_str();
    string filename1 = "gaussians_" + event_num + ".dat";
    string filename2 = "remnants_" + event_num + ".dat";

    double dummy, energy, x_c, y_c, eta_c;
    double total_energy = 0;
    double centre_x = 0.0;
    double centre_y = 0.0;

    ifstream gaussianfile(filename1.c_str());
    while (getline(gaussianfile, line))
    {
      istringstream iss(line);
      iss >> dummy >> x_c >> y_c >> eta_c >> dummy >> energy >> dummy;
      if (eta_c >= eta_min && eta_c <= eta_max)
      {
        centre_x += energy*x_c;
        centre_y += energy*y_c;
        total_energy += energy;
      }
    }
    gaussianfile.close();

    ifstream remnantsfile(filename2.c_str());
    while (getline(remnantsfile, line))
    {
      istringstream iss(line);
      iss >> dummy >> x_c >> y_c >> eta_c >> energy >> dummy;
      if (eta_c >= eta_min && eta_c <= eta_max)
      {
        centre_x += energy*x_c;
        centre_y += energy*y_c;
        total_energy += energy;
      }
    }
    remnantsfile.close();

    centre_x /= total_energy;
    centre_y /= total_energy;

    double epsilon_2_cos = 0.0;
    double epsilon_3_cos = 0.0;
    double epsilon_4_cos = 0.0;
    double epsilon_5_cos = 0.0;
    double epsilon_2_sin = 0.0;
    double epsilon_3_sin = 0.0;
    double epsilon_4_sin = 0.0;
    double epsilon_5_sin = 0.0;
    double epsilon_2_den = 0.0;
    double epsilon_3_den = 0.0;
    double epsilon_4_den = 0.0;
    double epsilon_5_den = 0.0;

    ifstream gaussianfile_again(filename1.c_str());
    while (getline(gaussianfile_again, line))
    {
      istringstream iss(line);
      iss >> dummy >> x_c >> y_c >> eta_c >> dummy >> energy >> dummy;
      if (eta_c >= eta_min && eta_c <= eta_max)
      {
        double dist = sqrt((x_c - centre_x)*(x_c - centre_x) + (y_c - centre_y)*(y_c - centre_y));
        double phi  = atan((y_c - centre_y)/(x_c - centre_x + 1e-10));

        epsilon_2_cos += pow(dist, 2.)*energy*cos(2.*phi);
        epsilon_3_cos += pow(dist, 3.)*energy*cos(3.*phi);
        epsilon_4_cos += pow(dist, 4.)*energy*cos(4.*phi);
        epsilon_5_cos += pow(dist, 5.)*energy*cos(5.*phi);

        epsilon_2_sin += pow(dist, 2.)*energy*sin(2.*phi);
        epsilon_3_sin += pow(dist, 3.)*energy*sin(3.*phi);
        epsilon_4_sin += pow(dist, 4.)*energy*sin(4.*phi);
        epsilon_5_sin += pow(dist, 5.)*energy*sin(5.*phi);

        epsilon_2_den += pow(dist, 2.)*energy;
        epsilon_3_den += pow(dist, 3.)*energy;
        epsilon_4_den += pow(dist, 4.)*energy;
        epsilon_5_den += pow(dist, 5.)*energy;
      }
    }
    gaussianfile_again.close();

    ifstream remnantsfile_again(filename2.c_str());
    while (getline(remnantsfile_again, line))
    {
      istringstream iss(line);
      iss >> dummy >> x_c >> y_c >> eta_c >> energy >> dummy;
      if (eta_c >= eta_min && eta_c <= eta_max)
      {
        double dist = sqrt((x_c - centre_x)*(x_c - centre_x) + (y_c - centre_y)*(y_c - centre_y));
        double phi  = atan((y_c - centre_y)/(x_c - centre_x + 1e-10));

        epsilon_2_cos += pow(dist, 2.)*energy*cos(2.*phi);
        epsilon_3_cos += pow(dist, 3.)*energy*cos(3.*phi);
        epsilon_4_cos += pow(dist, 4.)*energy*cos(4.*phi);
        epsilon_5_cos += pow(dist, 5.)*energy*cos(5.*phi);

        epsilon_2_sin += pow(dist, 2.)*energy*sin(2.*phi);
        epsilon_3_sin += pow(dist, 3.)*energy*sin(3.*phi);
        epsilon_4_sin += pow(dist, 4.)*energy*sin(4.*phi);
        epsilon_5_sin += pow(dist, 5.)*energy*sin(5.*phi);

        epsilon_2_den += pow(dist, 2.)*energy;
        epsilon_3_den += pow(dist, 3.)*energy;
        epsilon_4_den += pow(dist, 4.)*energy;
        epsilon_5_den += pow(dist, 5.)*energy;
      }
    }
    remnantsfile_again.close();

    epsilon_2_cos /= epsilon_2_den;
    epsilon_3_cos /= epsilon_3_den;
    epsilon_4_cos /= epsilon_4_den;
    epsilon_5_cos /= epsilon_5_den;
    epsilon_2_sin /= epsilon_2_den;
    epsilon_3_sin /= epsilon_3_den;
    epsilon_4_sin /= epsilon_4_den;
    epsilon_5_sin /= epsilon_5_den;

    fprintf(outfile, "%d %e %e %e %e %e %e %e %e %e %e %e \n", ievent, total_energy, centre_x, centre_y,
           epsilon_2_cos, epsilon_2_sin, epsilon_3_cos, epsilon_3_sin, epsilon_4_cos, epsilon_4_sin, epsilon_5_cos, epsilon_5_sin); 
  }

  fclose(outfile);
}
