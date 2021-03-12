#ifndef _SRC_STRUCT_H_
#define _SRC_STRUCT_H_

 struct nucleon
 {
   double x;
   double y;
   double z;
 };

 struct string_initial
 {
   double time;
   double x;
   double y;
   double z;
   int itarget;
   int iprojectile;
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
