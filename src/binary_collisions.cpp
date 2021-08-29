#include "./binary_collisions.h"

using namespace std;

namespace Binary_Collisions {

 void Determine_strings(std::vector<std::shared_ptr<binary_coll>>& binary_list, nucleon* Target, nucleon* Projectile,
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
         std::shared_ptr<binary_coll> new_string(new binary_coll);
         new_string->x = (Target[inucleon_T].x + Projectile[inucleon_P].x)/2.;
         new_string->y = (Target[inucleon_T].y + Projectile[inucleon_P].y)/2.;

         double delta_t = (Target[inucleon_T].z - Projectile[inucleon_P].z)/(2.*tanh(acosh(sqrtsNN/2.*nucleon_mass)));
         new_string->time = delta_t;

         new_string->z = Projectile[inucleon_P].z + delta_t*tanh(acosh(sqrtsNN/2.*nucleon_mass));

         new_string->itarget     = inucleon_T;
         new_string->iprojectile = inucleon_P;

         Nconn_T[inucleon_T]++;
         Nconn_P[inucleon_P]++;

         binary_list.push_back(new_string);
         sizeofstringlist = binary_list.size();
       }
     }
   }
   }

   std::sort(binary_list.begin(), binary_list.end(), sort_cond);
   return;
 }

 void Determine_all_collisions(std::vector<std::shared_ptr<binary_coll>>& binary_list, nucleon* Target, nucleon* Projectile,
      double sqrtsNN, double sigmaNN, gsl_rng* random, int TargetA, int ProjectileA, double inelastic_fraction)
 {
   int sizeofbinarylist; 

   for (int inucleon_T = 0; inucleon_T < TargetA; inucleon_T++)
   {
   for(int inucleon_P = 0; inucleon_P < ProjectileA; inucleon_P++)
   {
     double distance = sqrt(pow((Target[inucleon_T].x - Projectile[inucleon_P].x), 2.)
                          + pow((Target[inucleon_T].y - Projectile[inucleon_P].y), 2.));
     if (distance < sqrt(sigmaNN/M_PI))
     {
       double multiply_rapidities = Target[inucleon_T].rapidity*Projectile[inucleon_P].rapidity;
       double rand  = gsl_rng_uniform(random);
       if (rand < inelastic_fraction && multiply_rapidities < 0.)
       {
         double target_velocity     = tanh(Target[inucleon_T].rapidity)    ;
         double projectile_velocity = tanh(Projectile[inucleon_P].rapidity);

         double target_z, projectile_z, time_shift;
         if (Target[inucleon_T].time > Projectile[inucleon_P].time)
         {
           target_z = Target[inucleon_T].z;
           projectile_z = Projectile[inucleon_P].z + (Target[inucleon_T].time - Projectile[inucleon_P].time)*projectile_velocity;
           time_shift = Target[inucleon_T].time;
         } else
         {
           projectile_z = Projectile[inucleon_P].z;
           target_z = Target[inucleon_T].z + (Projectile[inucleon_P].time - Target[inucleon_T].time)*target_velocity;
           time_shift = Projectile[inucleon_P].time;
         }
         double delta_t = (target_z - projectile_z)/(projectile_velocity-target_velocity);
//         cout << Target[inucleon_T].z << " " << Projectile[inucleon_P].z << " " << target_z << "  " <<  projectile_z << endl;

         double mid_z =  projectile_z + delta_t*projectile_velocity;

         double y_loss = Propagation::Sample_y_loss(Target[inucleon_T].rapidity, Projectile[inucleon_P].rapidity, random);

         std::shared_ptr<binary_coll> new_coll(new binary_coll);
         new_coll->time = time_shift + delta_t;
         new_coll->x = (Target[inucleon_T].x + Projectile[inucleon_P].x)/2.;
         new_coll->y = (Target[inucleon_T].y + Projectile[inucleon_P].y)/2.;
         new_coll->z = mid_z;
         new_coll->itarget     = inucleon_T;
         new_coll->iprojectile = inucleon_P;

         new_coll->rapidity = y_loss;
         new_coll->y_com = (Target[inucleon_T].rapidity + Projectile[inucleon_P].rapidity)/2.;

         double y_T = Target[inucleon_T].rapidity;
         double y_P = Projectile[inucleon_P].rapidity;

         new_coll->sqrtsNN = 2.*nucleon_mass*cosh((y_P - y_T)/2.);

         if (Target[inucleon_T].rapidity < Projectile[inucleon_P].rapidity)
         {
           Target[inucleon_T].rapidity     += y_loss/2.;
           Projectile[inucleon_P].rapidity -= y_loss/2.;
         } else {
           Target[inucleon_T].rapidity     -= y_loss/2.;
           Projectile[inucleon_P].rapidity += y_loss/2.;
         }
//         Target[inucleon_T].rapidity     -= y_loss*y_T/(y_P - y_T);
//         Projectile[inucleon_P].rapidity -= y_loss*y_P/(y_P - y_T);
         Target[inucleon_T].z     = mid_z;
         Projectile[inucleon_P].z = mid_z;
         Target[inucleon_T].time     = time_shift + delta_t;
         Projectile[inucleon_P].time = time_shift + delta_t;

//         y_loss += y_com;

         new_coll->energy = nucleon_mass*(cosh(y_T) + cosh(y_P)
                           - cosh(Target[inucleon_T].rapidity) - cosh(Projectile[inucleon_P].rapidity));
         new_coll->p_z    = nucleon_mass*(sinh(y_T) + sinh(y_P) 
                           - sinh(Target[inucleon_T].rapidity) - sinh(Projectile[inucleon_P].rapidity));

//         if (new_coll->energy < 0.)
//         cout << new_coll->energy << "  " << y_loss << " " << (y_loss - y_com) << endl;

         binary_list.push_back(new_coll);
       }
     }
   }
   }

   std::sort(binary_list.begin(), binary_list.end(), sort_cond);
   return;
 }

 bool sort_cond(std::shared_ptr<binary_coll> a, std::shared_ptr<binary_coll> b)
 {
   return a->time < b->time;
 }
} 

