#include "./nuclei.h"

using namespace std;

namespace Nuclei {

 void Initialize_T_P_properties(string nucleusname, nucleus_properties tpproperties)
 {
   if (nucleusname.compare("p") == 0)
   {
     tpproperties.Z = 1;
     tpproperties.A = 1;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 1.0;
     tpproperties.w_WS = 0.0;
     tpproperties.a_WS = 1.0;
   } else if (nucleusname.compare("C") == 0)
   {
     tpproperties.Z = 6;
     tpproperties.A = 12;
     tpproperties.density_func = 3;
     tpproperties.R_WS = 2.44;
     tpproperties.w_WS = 1.403;
     tpproperties.a_WS = 1.635;
   } else if (nucleusname.compare("O") == 0)
   {
     tpproperties.Z = 8;
     tpproperties.A = 16;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 2.608;
     tpproperties.w_WS = -0.051;
     tpproperties.a_WS = 0.513;
   } else if (nucleusname.compare("S") == 0)
   {
     tpproperties.Z = 16;
     tpproperties.A = 32;
     tpproperties.density_func = 2;
     tpproperties.R_WS = 2.54;
     tpproperties.w_WS = 0.16;
     tpproperties.a_WS = 2.191;
   } else if (nucleusname.compare("W") == 0)
   {
     tpproperties.Z = 74;
     tpproperties.A = 184;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 6.51;
     tpproperties.w_WS = 0.0;
     tpproperties.a_WS = 0.535;
   } else if (nucleusname.compare("Al") == 0)
   {
     tpproperties.Z = 13;
     tpproperties.A = 27;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 3.07;
     tpproperties.w_WS = 0.0;
     tpproperties.a_WS = 0.519;
   } else if (nucleusname.compare("Au") == 0)
   {
     tpproperties.Z = 79;
     tpproperties.A = 197;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 6.38;
     tpproperties.w_WS = 0.0;
     tpproperties.a_WS = 0.535;
   } else if (nucleusname.compare("Ca") == 0)
   {
     tpproperties.Z = 20;
     tpproperties.A = 40;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 3.766;
     tpproperties.w_WS = -0.161;
     tpproperties.a_WS = 0.586;
   } else if (nucleusname.compare("Cu") == 0)
   {
     tpproperties.Z = 29;
     tpproperties.A = 63;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 4.163;
     tpproperties.w_WS = 0.0;
     tpproperties.a_WS = 0.606;
   } else if (nucleusname.compare("Fe") == 0)
   {
     tpproperties.Z = 26;
     tpproperties.A = 56;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 4.106;
     tpproperties.w_WS = 0.0;
     tpproperties.a_WS = 0.519;
   } else if (nucleusname.compare("Pb") == 0)
   {
     tpproperties.Z = 82;
     tpproperties.A = 208;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 6.62;
     tpproperties.w_WS = 0.0;
     tpproperties.a_WS = 0.546;
   } else if (nucleusname.compare("Pt") == 0)
   {
     tpproperties.Z = 78;
     tpproperties.A = 195;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 6.78;
     tpproperties.w_WS = 0.0;
     tpproperties.a_WS = 0.54;
   } else if (nucleusname.compare("Pt") == 0)
   {
     tpproperties.Z = 92;
     tpproperties.A = 238;
     tpproperties.density_func = 1;
     tpproperties.R_WS = 6.874;
     tpproperties.w_WS = 0.0;
     tpproperties.a_WS = 0.556;
   } else
   {
     cout << "Nucleus " << nucleusname << " is not known" << endl;
     exit(1);
   }
 }

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

 void Apply_random_rotation(nucleon* Target, nucleon* Projectile, gsl_rng* random, int TargetA, int ProjectileA)
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

   for (int inucleon = 0; inucleon < TargetA; inucleon++)
   {
     double Tx, Ty, Tz;

     Tx = Target[inucleon].x;
     Ty = Target[inucleon].y;
     Tz = Target[inucleon].z;

     Target[inucleon].x =    Tcosalpha*Tcosbeta*Tx 
                          + (Tcosalpha*Tsinbeta*Tsingamma - Tsinalpha*Tcosgamma)*Ty
                          + (Tcosalpha*Tsinbeta*Tcosgamma + Tsinalpha*Tsingamma)*Tz;
     Target[inucleon].y =    Tsinalpha*Tcosbeta*Tx
                          + (Tsinalpha*Tsinbeta*Tsingamma + Tcosalpha*Tcosgamma)*Ty
                          + (Tsinalpha*Tsinbeta*Tcosgamma - Tcosalpha*Tsingamma)*Tz;
     Target[inucleon].z = -1.* Tsinbeta*Tx
                          + Tcosbeta*Tsingamma*Ty
                          + Tcosbeta*Tcosgamma*Tz;
   }
   for (int inucleon = 0; inucleon < ProjectileA; inucleon++)
   {
     double Px, Py, Pz;

     Px = Projectile[inucleon].x;
     Py = Projectile[inucleon].y;
     Pz = Projectile[inucleon].z;

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

 void Apply_Lorentz_Contraction(nucleon* Target, nucleon* Projectile, double sqrtNN, int TargetA, int ProjectileA)
 {
   double alpha = sinh(acosh(sqrtNN/(2.*nucleon_mass)));
   double Lorentz_gamma = sqrt(1.+ alpha*alpha);

   for (int inucleon = 0; inucleon < TargetA    ; inucleon++) Target[inucleon].z     /= Lorentz_gamma;
   for (int inucleon = 0; inucleon < ProjectileA; inucleon++) Projectile[inucleon].z /= Lorentz_gamma;
   return;
 }

 void Shift_and_sort_nuclei(nucleon* Target, nucleon* Projectile, int TargetA, int ProjectileA)
 {
   double zshift_T = 0.;
   double zshift_P = 0.;
   for (int inucleon = 0; inucleon < TargetA    ; inucleon++) if (zshift_T > Target[inucleon].z    ) zshift_T = Target[inucleon].z    ;
   for (int inucleon = 0; inucleon < ProjectileA; inucleon++) if (zshift_P > Projectile[inucleon].z) zshift_P = Projectile[inucleon].z;

   for (int inucleon = 0; inucleon < TargetA    ; inucleon++) Target[inucleon].z     -= zshift_T;
   for (int inucleon = 0; inucleon < ProjectileA; inucleon++) Projectile[inucleon].z -= zshift_P;

   sortincreasing(Target, TargetA);
   sortdecreasing(Projectile, ProjectileA);
   return;
 }

 void sortincreasing(nucleon* A, int nucleonnum)
 {
   nucleon temporary;
   for (int inucleon = 0; inucleon < nucleonnum; inucleon++)
   {
     for (int jnucleon = inucleon; jnucleon < nucleonnum; jnucleon++)
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

 void sortdecreasing(nucleon* A, int nucleonnum)
 {
   nucleon temporary;
   for (int inucleon = 0; inucleon < nucleonnum; inucleon++)
   {
     for (int jnucleon = inucleon; jnucleon < nucleonnum; jnucleon++)
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

 void Sample_WS_nucleus(nucleon* Nucleus, double sqrtsNN, gsl_rng* random, nucleus_properties properties)
 {
   int nucleon_num = properties.A;
   double Radii[nucleon_num];
   double radius;

   //Sample nucleon radii. Taking cue from Trento, we sample radii in the range [0,R+10a]
   for (int inucleon = 0; inucleon < nucleon_num; inucleon++)
   {
     do {
        radius = (properties.R_WS + 10.*properties.a_WS)*pow(gsl_rng_uniform(random),1./3.);
     } while (gsl_rng_uniform(random) > (1. + properties.w_WS*pow(radius/properties.R_WS,2.))/(1.+exp((radius - properties.R_WS)/properties.a_WS)));
     Radii[inucleon] = radius;
   }

   int num = sizeof(Radii)/sizeof(Radii[0]);
   sort(Radii, Radii + num);
  
   double theta, phi, posx, posy, posz; 
   for (int inucleon = 0; inucleon < nucleon_num; inucleon++)
   {

     do {
        theta = acos(1. - 2.*gsl_rng_uniform(random));
        phi   = 2.*M_PI*gsl_rng_uniform(random);

        posx = Radii[inucleon]*cos(theta)*cos(phi);
        posy = Radii[inucleon]*cos(theta)*sin(phi);
        posz = Radii[inucleon]*sin(theta);

     } while (Nucleon_is_close(Nucleus, inucleon, posx, posy, posz));

     Nucleus[inucleon].x = posx;
     Nucleus[inucleon].y = posy;
     Nucleus[inucleon].z = posz;
     Nucleus[inucleon].time = 0.;
     Nucleus[inucleon].rapidity = acosh(sqrtsNN/(2.*nucleon_mass));
   }
   return; 
 }

 bool Nucleon_is_close(nucleon* Nucleus, int inucleon, double x, double y, double z)
 {
   bool check = false;
   double dmin_sq = dmin*dmin;
   for (int idx = 0; idx < inucleon; idx++)
   {
     double dx = Nucleus[inucleon].x - x;
     double dy = Nucleus[inucleon].y - y;
     double dz = Nucleus[inucleon].z - z;

     double distsq = dx*dx + dy*dy + dz*dz;

     if (distsq < dmin_sq) 
     {
       check = true;
       break;
     }
   }
   return check;
 }

}

