#include "./propagation.h"

using namespace std;

namespace Propagation {

 void Propagate_strings(std::vector<std::shared_ptr<binary_coll>>& string_list, std::vector<std::shared_ptr<string_final>>& string_final_list,
                                 double sqrtsNN, gsl_rng* random)
 {
   double yin_T[197], yin_P[197];
   for (int inucleon = 0; inucleon < 197; inucleon++)
   {
     yin_T[inucleon] = -1.*acosh(sqrtsNN/(2.*nucleon_mass));
     yin_P[inucleon] = acosh(sqrtsNN/(2.*nucleon_mass));
   }

   for (int istring = 0; istring < string_list.size(); istring++)
   {
     std::shared_ptr<string_final> new_string(new string_final);
     new_string->x_c = string_list[istring]->x;
     new_string->y_c = string_list[istring]->y;
     int itarget     = string_list[istring]->itarget    ;
     int iprojectile = string_list[istring]->iprojectile;

     double y_loss = Sample_y_loss(yin_T[itarget], yin_P[iprojectile], random);

     double delta_tau_f = (nucleon_mass/tension)*sqrt(2.*(cosh(y_loss)-1.));
     new_string->delta_tau_f = delta_tau_f;

     double Y_com_string = (yin_T[itarget] + yin_P[iprojectile])/2.;

     double y_deceleration = acosh((pow(delta_tau_f*tension/nucleon_mass,2.)/2.) + 1.);
     double y_l = yin_T[itarget]     + y_deceleration;
     double y_r = yin_P[iprojectile] - y_deceleration;
     new_string->y_l_i = yin_T[itarget]    ;
     new_string->y_r_i = yin_P[iprojectile];
     new_string->y_l = y_l;
     new_string->y_r = y_r;
     
     double t_c = string_list[istring]->time;
     double z_c = string_list[istring]->z;
     double eta_c = log((t_c+z_c)/(t_c-z_c))*0.5;
     double tau_c = sqrt(t_c*t_c - z_c*z_c);
     new_string->eta_c = eta_c;
     new_string->tau_c = tau_c;

     double t_l_lab = t_c + delta_tau_f*(    (delta_tau_f*tension*sinh(yin_T[itarget]    )/(2.*nucleon_mass))
                                        + sqrt(pow(delta_tau_f*tension/(2.*nucleon_mass),2.) + 1.)*cosh(yin_T[itarget]))    ;
     double t_r_lab = t_c + delta_tau_f*(-1.*(delta_tau_f*tension*sinh(yin_P[iprojectile])/(2.*nucleon_mass))
                                        + sqrt(pow(delta_tau_f*tension/(2.*nucleon_mass),2.) + 1.)*cosh(yin_P[iprojectile]));
     double z_l_lab = z_c + delta_tau_f*(    (delta_tau_f*tension*cosh(yin_T[itarget]    )/(2.*nucleon_mass))
                                        + sqrt(pow(delta_tau_f*tension/(2.*nucleon_mass),2.) + 1.)*sinh(yin_T[itarget]))    ;
     double z_r_lab = z_c + delta_tau_f*(-1.*(delta_tau_f*tension*cosh(yin_P[iprojectile])/(2.*nucleon_mass))
                                        + sqrt(pow(delta_tau_f*tension/(2.*nucleon_mass),2.) + 1.)*sinh(yin_P[iprojectile]));

     new_string->eta_l = log((t_l_lab + z_l_lab)/(t_l_lab - z_l_lab));
     new_string->eta_r = log((t_r_lab + z_r_lab)/(t_r_lab - z_r_lab));

     string_final_list.push_back(new_string);
     yin_T[itarget]     = y_l;
     yin_P[iprojectile] = y_r;
   }
 }

 double Sample_y_loss(double yinT, double yinP, gsl_rng* random)
 {
   double boost = 0.5*(yinT + yinP);
   double yinT_lrf = yinT - boost;
   double yinP_lrf = yinP - boost;

   double yin_lrf = abs(yinT_lrf) + abs(yinP_lrf);

   double rannum = gsl_rng_uniform(random);

   double y_loss = rannum*(2.*cosh(yin_lrf) - 1.)*sinh(yin_lrf);
   y_loss -= sinh(2.*yin_lrf);
   y_loss = asinh(-1.*y_loss);
   y_loss = 2.*yin_lrf - y_loss;

   return y_loss;
 }

 double Calculate_total_final_energy(std::vector<std::shared_ptr<string_final>>& string_final_list)
 {
   double total_energy = 0.;
   for (int istring = 0; istring < string_final_list.size(); istring++)
   {
     double E_string = (nucleon_mass*cosh(string_final_list[istring]->y_l_i) + nucleon_mass*cosh(string_final_list[istring]->y_r_i)
                           - nucleon_mass*cosh(string_final_list[istring]->y_l) - nucleon_mass*cosh(string_final_list[istring]->y_r));
     total_energy += E_string;
   }
   return total_energy;
 }

 void Propagate_gaussians(std::vector<std::shared_ptr<binary_coll>>& binary_list, std::vector<std::shared_ptr<gaussians>>& gaussian_list)
 {

   for (int icoll = 0; icoll < binary_list.size(); icoll++)
   {
    std::shared_ptr<gaussians> new_gaussian(new gaussians);

    double sqrtsNN = binary_list[icoll]->sqrtsNN;
//    double gaussian_width = (3./10.)*log(sqrtsNN/(2.*nucleon_mass));
    double gaussian_width = sqrt(log(sqrtsNN/(2.*nucleon_mass)));
    new_gaussian->width_eta = gaussian_width;

//    double eta_c = atanh(binary_list[icoll]->z/binary_list[icoll]->time);// + binary_list[icoll]->y_com;
//    double tau = sqrt(pow(binary_list[icoll]->time,2.) - pow(binary_list[icoll]->z,2.)) + tau_thermalize;
    double eta_c = atanh(binary_list[icoll]->z/(binary_list[icoll]->time+tau_thermalize));
    double tau = sqrt(pow((binary_list[icoll]->time+tau_thermalize),2.) - pow(binary_list[icoll]->z,2.));// + tau_thermalize;
    new_gaussian->eta_c = eta_c;
    new_gaussian->tau = tau;
    new_gaussian->x_c = binary_list[icoll]->x;
    new_gaussian->y_c = binary_list[icoll]->y;

    cout << binary_list[icoll]->time << "  " << tau << "  " << binary_list[icoll]->z << "  " << eta_c << endl;

    double energy_milne =  cosh(eta_c)*binary_list[icoll]->energy - sinh(eta_c)*binary_list[icoll]->p_z;
    double pz_milne     = -sinh(eta_c)*binary_list[icoll]->energy + cosh(eta_c)*binary_list[icoll]->p_z;

//    new_gaussian->energy = binary_list[icoll]->energy;
//    new_gaussian->p_z    = binary_list[icoll]->p_z   ;
    new_gaussian->energy = energy_milne;
    new_gaussian->p_z    = pz_milne    ;

    gaussian_list.push_back(new_gaussian);
   }
   return;
 }

 void Propagate_remnants(nucleon* Target, nucleon* Projectile, std::vector<std::shared_ptr<remnant>>& remnant_list, double sqrtsNN,
                        int TargetA, int ProjectileA)
 {
   double y_beam = acosh(sqrtsNN/(2.*nucleon_mass));

   for (int inucleon = 0; inucleon < TargetA; inucleon++)
   {
     if(abs(Target[inucleon].rapidity) < y_beam)
     {
       std::shared_ptr<remnant> new_remnant(new remnant);

       new_remnant->x = Target[inucleon].x;
       new_remnant->y = Target[inucleon].y;
       
       double rapidity = Target[inucleon].rapidity;
  
       double tau_i = sqrt(pow(Target[inucleon].time,2.) - pow(Target[inucleon].z,2.));
       new_remnant->tau = tau_i + tau_thermalize;

       double velocity = tanh(Target[inucleon].rapidity);
       double eta_i = atanh(Target[inucleon].z/Target[inucleon].time);

       double b = (tau_i*sinh(eta_i) - velocity*tau_i*cosh(eta_i))/(tau_thermalize + tau_i);

       double eta_f = log((-1.*b - sqrt(1-velocity*velocity + b*b))/(velocity-1.));
       new_remnant->eta = eta_f;

       double energy = nucleon_mass*cosh(rapidity);
       double p_z    = nucleon_mass*sinh(rapidity);
       double energy_milne =  cosh(eta_f)*energy - sinh(eta_f)*p_z;
       double pz_milne     = -sinh(eta_f)*energy + cosh(eta_f)*p_z;

       new_remnant->energy = energy_milne;
       new_remnant->p_z    = pz_milne;

       remnant_list.push_back(new_remnant);
     }
   }
   for (int inucleon = 0; inucleon < ProjectileA; inucleon++)
   {
     if(abs(Projectile[inucleon].rapidity) < y_beam)
     {
       std::shared_ptr<remnant> new_remnant(new remnant);

       new_remnant->x = Projectile[inucleon].x;
       new_remnant->y = Projectile[inucleon].y;
       
       double rapidity = Projectile[inucleon].rapidity;
  
       double tau_i = sqrt(pow(Projectile[inucleon].time,2.) - pow(Projectile[inucleon].z,2.));
       new_remnant->tau = tau_i + tau_thermalize;

       double velocity = tanh(Projectile[inucleon].rapidity);
       double eta_i = atanh(Projectile[inucleon].z/Projectile[inucleon].time);

       double b = (tau_i*sinh(eta_i) - velocity*tau_i*cosh(eta_i))/(tau_thermalize + tau_i);

       double eta_f = log((-1.*b - sqrt(1-velocity*velocity + b*b))/(velocity-1.));
       new_remnant->eta = eta_f;

       double energy = nucleon_mass*cosh(rapidity);
       double p_z    = nucleon_mass*sinh(rapidity);
       double energy_milne =  cosh(eta_f)*energy - sinh(eta_f)*p_z;
       double pz_milne     = -sinh(eta_f)*energy + cosh(eta_f)*p_z;

       new_remnant->energy = energy_milne;
       new_remnant->p_z    = pz_milne;

       remnant_list.push_back(new_remnant);
     }
   }

   return;
 }
} 

