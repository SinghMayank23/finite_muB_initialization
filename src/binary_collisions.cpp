#include "./binary_collisions.h"

using namespace std;

namespace Binary_Collisions {

 void Determine_strings(std::vector<std::shared_ptr<string_initial>>& string_list, nucleon* Target, nucleon* Projectile,
                                 double sqrtsNN, double sigmaNN, gsl_rng* random)
 {
   int Nconn_T[197], Nconn_P[197];
   for (int inucleon = 0; inucleon < 197; inucleon++)
   {
     Nconn_T[inucleon] = 0;
     Nconn_P[inucleon] = 0;
   }

   int sizeofstringlist; 

   for (int inucleon_T = 0; inucleon_T < 197; inucleon_T++)
   {
   for(int inucleon_P = 0; inucleon_P < 197; inucleon_P++)
   {
     double distance = sqrt(pow((Target[inucleon_T].x - Projectile[inucleon_P].x), 2.)
                          + pow((Target[inucleon_T].y - Projectile[inucleon_P].y), 2.));
     if (distance < 2.*sqrt(sigmaNN/M_PI))
     {
       double num = gsl_rng_uniform(random);
       int Nconn = Nconn_T[inucleon_T] + Nconn_P[inucleon_P];
       double probability = exp(-0.5*(double)Nconn);

       if (num < probability)
       {
         std::shared_ptr<string_initial> new_string(new string_initial);
         new_string->x = (Target[inucleon_T].x + Projectile[inucleon_P].x)/2.;
         new_string->y = (Target[inucleon_T].y + Projectile[inucleon_P].y)/2.;

         double delta_t = (Target[inucleon_T].z - Projectile[inucleon_P].z)/(2.*tanh(acosh(sqrtsNN/2.*nucleon_mass)));
         new_string->time = delta_t;

         new_string->z = Projectile[inucleon_P].z + delta_t*tanh(acosh(sqrtsNN/2.*nucleon_mass));

         new_string->itarget     = inucleon_T;
         new_string->iprojectile = inucleon_P;

         Nconn_T[inucleon_T]++;
         Nconn_P[inucleon_P]++;

         string_list.push_back(new_string);
         sizeofstringlist = string_list.size();
       }
     }
   }
   }

   std::sort(string_list.begin(), string_list.end(), sort_cond);
   return;
 }

 bool sort_cond(std::shared_ptr<string_initial> a, std::shared_ptr<string_initial> b)
 {
   return a->time < b->time;
 }
} 

