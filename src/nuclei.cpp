#include "./grid.h"

namespace Nuclei {

   double epsabs = 1.e-15;
   double epsrel = 1.e-10;
//   struct cell
//   {
//     double   ;
//     double heavymomentum;
//     int    int_index    ;
//   };

 void Initialize_projectile_nucleus(ifstream& nucleusfile, int nuclei_number, nucleon Projectile)
 {
   std::string dummy;
   for (int iline = 0; iline < (nuclei_number-1)*197; iline++)//Works only for Au nuclei
   {
     std::getline(nucleus, dummy);
   }
   
   for (int inucleon = 0; inucleon < 197; inucleon++)//Works only for Au-197
   {
     double posx, posy, posz, spin, isospin;//check is last 2 columns are spin and isospin
     nucleusfile >> posx >> posy >> posz >> spin >> isospin;

     Projectile[inucleon].x = posx;
     Projectile[inucleon].y = posy;
     Projectile[inucleon].z = posz;
   }
   return; 
 }

 void Initialize_target_nucleus(ifstream& nucleusfile, int nuclei_number, nucleon Target)
 {
   std::string dummy;
   for (int iline = 0; iline < (nuclei_number-1)*197; iline++)//Works only for Au nuclei
   {
     std::getline(nucleus, dummy);
   }
   
   for (int inucleon = 0; inucleon < 197; inucleon++)//Works only for Au-197
   {
     double posx, posy, posz, spin, isospin;//check is last 2 columns are spin and isospin
     nucleusfile >> posx >> posy >> posz >> spin >> isospin;

     Target[inucleon].x = posx;
     Target[inucleon].y = posy;
     Target[inucleon].z = posz;
   }
   return; 
 }

 void Apply_random_rotation(nucleon Target, nucleon Projectile, gsl_rng* random)
 {
   double Talpha = 2.*M_PI*gsl_rng_uniform(random);//Check Euler angles
   double Tbeta  =    M_PI*gsl_rng_uniform(random);
   double Tgamma = 2.*M_PI*gsl_rng_uniform(random);

   double Tcosalpha  = cos(alpha)
   double Tsinalpha  = sin(alpha);
   double Tcosbeta   = cos(beta )
   double Tsinbeta   = sin(beta );
   double Tcosgamma  = cos(gamma)
   double Tsingamma  = sin(gamma);

   double Palpha = 2.*M_PI*gsl_rng_uniform(random);//Check Euler angles
   double Pbeta  =    M_PI*gsl_rng_uniform(random);
   double Pgamma = 2.*M_PI*gsl_rng_uniform(random);

   double Pcosalpha  = cos(alpha)
   double Psinalpha  = sin(alpha);
   double Pcosbeta   = cos(beta )
   double Psinbeta   = sin(beta );
   double Pcosgamma  = cos(gamma)
   double Psingamma  = sin(gamma);

   for (int inuleon = 0; inucleon < 197; inucleon++)
   {
     double Tx, Ty, Tz, Px, Py, Px;

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

 void Apply_Lorentz_Contraction(nucleon Target, nucleon Projectile, double sqrtNN)
 {
   double nucleon_mass = 0.939;
   double alpha = sinh(acosh(sqrtNN/(2.*nucleon_mass)));
   double Lorentz_gamma = sqrt(1.+ alpha*alpha);

   for (int inuleon = 0; inucleon < 197; inucleon++)
   {
     Target[inucleon].z     /= Lorentz_gamma;
     Projectile[inucleon].z /= Lorentz_gamma;
   }
   return;
 }
}

