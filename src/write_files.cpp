#include "./write_files.h"

using namespace std;

namespace Write_files {
 void Write_nucleon_positions_and_binary_collisions(nucleon* Target, nucleon* Projectile, std::vector<std::shared_ptr<binary_coll>>& string_list, FILE *fileout1,
                                                    int TargetA, int ProjectileA)
 {
   fprintf(fileout1, "%s \n", "Target x, y, z");
   for (int inucleon = 0; inucleon < TargetA; inucleon++)
   {
     fprintf(fileout1, "%e %e %e\n", Target[inucleon].x, Target[inucleon].y, Target[inucleon].z);
   }
   fprintf(fileout1, "\n");
   fprintf(fileout1, "%s \n", "Projectile x, y, z");
   for (int inucleon = 0; inucleon < ProjectileA; inucleon++)
   {
     fprintf(fileout1, "%e %e %e \n", Projectile[inucleon].x, Projectile[inucleon].y, Projectile[inucleon].z);
   }
   fprintf(fileout1, "\n");

   int sizeofstringlist = string_list.size(); 
   fprintf(fileout1, "%s \n", "Binary Collisions");
   fprintf(fileout1, "%s \n", " time, x, y, z, itarget, iprojectile");
   for (int istring = 0; istring < sizeofstringlist; istring++)
   {
     fprintf(fileout1, "%e %e %e %e %d %d \n", string_list[istring]->time, string_list[istring]->x,
                                                  string_list[istring]->y, string_list[istring]->z,
                                   string_list[istring]->itarget, string_list[istring]->iprojectile); 
   }

   return;
 }

 void Write_final_terms(std::vector<std::shared_ptr<string_final>>& string_final_list, FILE *fileout2, double total_energy)
 {
   fprintf(fileout2, "%s %e \n", "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", total_energy);
   fprintf(fileout2, "%s \n", " "); 
 
   for (int istring = 0; istring < string_final_list.size(); istring++)
   {
     fprintf(fileout2, "%e %e %e %e %e %e \n", nucleon_mass, nucleon_mass/tension, 0.5,
                                         string_final_list[istring]->tau_c, string_final_list[istring]->eta_c,
                                         string_final_list[istring]->x_c  , string_final_list[istring]->y_c  , 
                                         string_final_list[istring]->x_c  , string_final_list[istring]->y_c  , 
                                         string_final_list[istring]->x_c  , string_final_list[istring]->y_c  ,
                                         string_final_list[istring]->eta_l, string_final_list[istring]->eta_r,
                                         string_final_list[istring]->y_l  , string_final_list[istring]->y_r  );
   }

   return;
 }

 void Write_gaussians(std::vector<std::shared_ptr<gaussians>>& gaussian_list, FILE *fileout1)
 {
   for(int igauss = 0; igauss < gaussian_list.size(); igauss++)
   {
     fprintf(fileout1, "%e %e %e %e %e %e %e \n", gaussian_list[igauss]->tau      , gaussian_list[igauss]->x_c,
                                                  gaussian_list[igauss]->y_c      , gaussian_list[igauss]->eta_c,
                                                  gaussian_list[igauss]->width_eta, gaussian_list[igauss]->energy,
                                                  gaussian_list[igauss]->p_z     );
   }
   return;
 }

 void Write_remnants(std::vector<std::shared_ptr<remnant>>& remnant_list, FILE *fileout2)
 {
   for(int irem = 0; irem < remnant_list.size(); irem++)
   {
     fprintf(fileout2, "%e %e %e %e %e %e \n",   remnant_list[irem]->tau, remnant_list[irem]->x,
                                              remnant_list[irem]->y     , remnant_list[irem]->eta,
                                              remnant_list[irem]->energy, remnant_list[irem]->p_z);
   }
   return;
 }
} 
