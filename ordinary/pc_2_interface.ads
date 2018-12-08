
with pc_2_coeff_18; -- Order 18, highest stability
with pc_2_coeff_19;
with pc_2_coeff_20;
with pc_2_coeff_21;
with pc_2_coeff_22; -- Order 22, lowest stability

package pc_2_interface is

 --generic package Predictor_Corrector_Rules renames pc_2_coeff_18;
 --generic package Predictor_Corrector_Rules renames pc_2_coeff_19;

   generic package Predictor_Corrector_Rules renames pc_2_coeff_20;
   -- best general purpose. 
   -- Not much reason to use others except possibly pc_2_coeff_21.

 --generic package Predictor_Corrector_Rules renames pc_2_coeff_21;
 --generic package Predictor_Corrector_Rules renames pc_2_coeff_22;

end pc_2_interface;

