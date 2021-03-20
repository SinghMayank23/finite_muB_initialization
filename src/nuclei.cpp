#include "./nuclei.h"

namespace Nuclei {

 void Initialize_projectile_nucleus(std::ifstream& nucleusfile, int nuclei_number, nucleon* Projectile, double sqrtsNN)
 {
   std::string dummy;
   for (int iline = 0; iline < (nuclei_number-1)*197; iline++)//Works only for Au nuclei
   {
     std::getline(nucleusfile, dummy);
   }
   
   for (int inucleon = 0; inucleon < 197; inucleon++)//Works only for Au-197
   {
     double posx, posy, posz, spin, isospin;//check is last 2 columns are spin and isospin
     nucleusfile >> posx >> posy >> posz >> spin >> isospin;

     Projectile[inucleon].x = posx;
     Projectile[inucleon].y = posy;
     Projectile[inucleon].z = posz;
     Projectile[inucleon].time = 0.;
     Projectile[inucleon].rapidity = acosh(sqrtsNN/(2.*nucleon_mass));
   }
   return; 
 }

 void Initialize_target_nucleus(std::ifstream& nucleusfile, int nuclei_number, nucleon* Target, double sqrtsNN)
 {
   std::string dummy;
   for (int iline = 0; iline < (nuclei_number-1)*197; iline++)//Works only for Au nuclei
   {
     std::getline(nucleusfile, dummy);
   }
   
   for (int inucleon = 0; inucleon < 197; inucleon++)//Works only for Au-197
   {
     double posx, posy, posz, spin, isospin;//check is last 2 columns are spin and isospin
     nucleusfile >> posx >> posy >> posz >> spin >> isospin;

     Target[inucleon].x = posx;
     Target[inucleon].y = posy;
     Target[inucleon].z = posz; 
     Target[inucleon].time = 0.; 
     Target[inucleon].rapidity = -1.*acosh(sqrtsNN/(2.*nucleon_mass));
   }
   return; 
 }

 void Apply_random_rotation(nucleon* Target, nucleon* Projectile, gsl_rng* random)
 {
   double Talpha = 2.*M_PI*gsl_rng_uniform(random);//Check Euler angles
   double Tbeta  =    M_PI*gsl_rng_uniform(random);
   double Tgamma = 2.*M_PI*gsl_rng_uniform(random);

   double Tcosalpha  = cos(Talpha);
   double Tsinalpha  = sin(Talpha);
   double Tcosbeta   = cos(Tbeta );
   double Tsinbeta   = sin(Tbeta );
   double Tcosgamma  = cos(Tgamma);
   double Tsingamma  = sin(Tgamma);

   double Palpha = 2.*M_PI*gsl_rng_uniform(random);//Check Euler angles
   double Pbeta  =    M_PI*gsl_rng_uniform(random);
   double Pgamma = 2.*M_PI*gsl_rng_uniform(random);

   double Pcosalpha  = cos(Palpha);
   double Psinalpha  = sin(Palpha);
   double Pcosbeta   = cos(Pbeta );
   double Psinbeta   = sin(Pbeta );
   double Pcosgamma  = cos(Pgamma);
   double Psingamma  = sin(Pgamma);

   for (int inucleon = 0; inucleon < 197; inucleon++)
   {
     double Tx, Ty, Tz, Px, Py, Pz;

     Tx = Target[inucleon].x;
     Ty = Target[inucleon].y;
     Tz = Target[inucleon].z;

     Px = Projectile[inucleon].x;
     Py = Projectile[inucleon].y;
     Pz = Projectile[inucleon].z;

     Target[inucleon].x =    Tcosalpha*Tcosbeta*Tx 
                          + (Tcosalpha*Tsinbeta*Tsingamma - Tsinalpha*Tcosgamma)*Ty
                          + (Tcosalpha*Tsinbeta*Tcosgamma + Tsinalpha*Tsingamma)*Tz;
     Target[inucleon].y =    Tsinalpha*Tcosbeta*Tx
                          + (Tsinalpha*Tsinbeta*Tsingamma + Tcosalpha*Tcosgamma)*Ty
                          + (Tsinalpha*Tsinbeta*Tcosgamma - Tcosalpha*Tsingamma)*Tz;
     Target[inucleon].z = -1.* Tsinbeta*Tx
                          + Tcosbeta*Tsingamma*Ty
                          + Tcosbeta*Tcosgamma*Tz;

     Projectile[inucleon].x =    Pcosalpha*Pcosbeta*Px 
                              + (Pcosalpha*Psinbeta*Psingamma - Psinalpha*Pcosgamma)*Py
                              + (Pcosalpha*Psinbeta*Pcosgamma + Psinalpha*Psingamma)*Pz;
     Projectile[inucleon].y =    Psinalpha*Pcosbeta*Px
                              + (Psinalpha*Psinbeta*Psingamma + Pcosalpha*Pcosgamma)*Py
                              + (Psinalpha*Psinbeta*Pcosgamma - Pcosalpha*Psingamma)*Pz;
     Projectile[inucleon].z = -1.* Psinbeta*Px
                              + Pcosbeta*Psingamma*Py
                              + Pcosbeta*Pcosgamma*Pz;
   }
   return;
 }

 void Apply_Lorentz_Contraction(nucleon* Target, nucleon* Projectile, double sqrtNN)
 {
   double alpha = sinh(acosh(sqrtNN/(2.*nucleon_mass)));
   double Lorentz_gamma = sqrt(1.+ alpha*alpha);

   for (int inucleon = 0; inucleon < 197; inucleon++)
   {
     Target[inucleon].z     /= Lorentz_gamma;
     Projectile[inucleon].z /= Lorentz_gamma;
   }
   return;
 }

 void Shift_and_sort_nuclei(nucleon* Target, nucleon* Projectile)
 {
   double zshift_T = 0.;
   double zshift_P = 0.;
   for (int inucleon = 0; inucleon < 197; inucleon++)
   {
     if (zshift_T > Target[inucleon].z    ) zshift_T = Target[inucleon].z    ;
     if (zshift_P < Projectile[inucleon].z) zshift_P = Projectile[inucleon].z;
   }

   for (int inucleon = 0; inucleon < 197; inucleon++)
   {
     Target[inucleon].z     -= zshift_T;
     Projectile[inucleon].z -= zshift_P;
   }

   sortincreasing(Target);
   sortdecreasing(Projectile);
   return;
 }

 void sortincreasing(nucleon* A)
 {
   nucleon temporary;
   for (int inucleon = 0; inucleon < 197; inucleon++)
   {
     for (int jnucleon = inucleon; jnucleon < 197; jnucleon++)
     {
       if (A[inucleon].z > A[jnucleon].z)
       {
         temporary = A[jnucleon];
         A[jnucleon] = A[inucleon];
         A[inucleon] = temporary;
       }
     } 
   }
   return;
 }

 void sortdecreasing(nucleon* A)
 {
   nucleon temporary;
   for (int inucleon = 0; inucleon < 197; inucleon++)
   {
     for (int jnucleon = inucleon; jnucleon < 197; jnucleon++)
     {
       if (A[inucleon].z < A[jnucleon].z)
       {
         temporary = A[jnucleon];
         A[jnucleon] = A[inucleon];
         A[inucleon] = temporary;
       }
     } 
   }
   return;
 }
}

