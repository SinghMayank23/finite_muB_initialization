#include "./write_file.h"

using namespace std;

namespace Write_files {
 void Write_nucleon_positions_and_binary_collisions(nucleon* Target, nucleon* Projectile, std::vector<std::shared_ptr<string_initial>>& string_list, FILE *fileout1)
 {
   fprintf(fileout1, "%s \n", "Target x, y, z; Projectile x, y, z");
   for (int inucleon = 0; inucleon < 197; inucleon++)
   {
     fprint(fileout1, "%e %e %e %e %e %e \n", Target[inucleon].x, Target[inucleon].y, Target[inucleon].z
                                  Projectile[inucleon].x, Projectile[inucleon].y, Projectile[inucleon].z);
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
   fprintf(fileout2, "%s \n", " ")  
 
   for (int istring = 0; istring < string_final_list.size(); istring++)
   {
     fprint(fileout2, "%e %e %e %e %e %e \n", nucleon_mass, nucleon_mass/tension, 0.5,
                                         string_final_list[istring]->tau_c, string_final_list[istring]->eta_c,
                                         string_final_list[istring]->x_c  , string_final_list[istring]->y_c  , 
                                         string_final_list[istring]->x_c  , string_final_list[istring]->y_c  , 
                                         string_final_list[istring]->x_c  , string_final_list[istring]->y_c  ,
                                         string_final_list[istring]->eta_l, string_final_list[istring]->eta_r,
                                         string_final_list[istring]->y_l  , string_final_list[istring]->y_r  ,
                                        
   }

   return;
 }

} 
       text_stream >> new_string->mass >> new_string->m_over_sigma
                    >> new_string->tau_form
                    >> new_string->tau_0 >> new_string->eta_s_0
                    >> new_string->x_perp >> new_string->y_perp
                    >> new_string->x_pl >> new_string->y_pl
                    >> new_string->x_pr >> new_string->y_pr
                    >> new_string->eta_s_left >> new_string->eta_s_right
                    >> new_string->y_l >> new_string->y_r
                    >> new_string->remnant_l >> new_string->remnant_r
                    >> new_string->y_l_i >> new_string->y_r_i
                    >> new_string->eta_s_baryon_left
                    >> new_string->eta_s_baryon_right
                    >> new_string->y_l_baryon
                    >> new_string->y_r_baryon
                    >> new_string->baryon_frac_l;
        if (!text_stream.eof()) {
            // read in the last element
            text_stream >> new_string->baryon_frac_r;

