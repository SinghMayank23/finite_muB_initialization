#ifndef _SRC_STRUCT_H_
#define _SRC_STRUCT_H_

 struct nucleon
 {
   double x;
   double y;
   double z;
   double time;
   double rapidity;
 };

 struct remnant
 {
   double x;
   double y;
   double eta;
   double tau;
   double energy;
   double p_z;
 };

 struct binary_coll
 {
   double time;
   double x;
   double y;
   double z;
   int itarget;
   int iprojectile;
   double rapidity;
   double y_com;
   double energy;
   double p_z;
   double sqrtsNN;
 };

 struct gaussians
 {
   double tau;
   double eta_c;
   double width_eta;
   double x_c;
   double y_c;
   double energy;
   double p_z;
 };

 struct string_final
 {
   double tau_c;
   double eta_c;
   double x_c;
   double y_c;
   double delta_tau_f;
   double eta_l;
   double eta_r;
   double y_l;
   double y_r;
   double y_l_i;
   double y_r_i;
   int one_over_rem_l;
   int one_over_rem_r;
 };

#endif
