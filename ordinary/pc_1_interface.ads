

-- PACKAGE pc_1_interface
--
-- The 17th order coefficients are almost always the best choice.
--
-- The Predictor-Corrector method uses least-squares to fit a 17th 
-- order polynomial to the 33 previous values of F in dY/dt = F(t,Y).
-- The least-squares-fit polynomial is used to predict the next
-- values of F and Y.

with pc_1_coeff_16;
with pc_1_coeff_17;
with pc_1_coeff_18;

package pc_1_interface is
   
   --  predictor_Coeff_16 (16th order) gives bit better stability.
   --  predictor_Coeff_17 (17th order) is usually best for numerical 
   --                      accuracy, and best choice for general use.
   --  predictor_Coeff_18 (18th order) gives a bit better numerical
   --                      accuracy when larger stepsizes are used (and when
   --                      desired accuracy is well below machine's ultimate).

 --generic package Predictor_Corrector_Rules renames pc_1_coeff_16;
   generic package Predictor_Corrector_Rules renames pc_1_coeff_17; --use this
 --generic package Predictor_Corrector_Rules renames pc_1_coeff_18;

end pc_1_interface;

