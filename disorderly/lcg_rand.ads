

with mwc_constants;
use  mwc_constants;

generic

   type Parent_Random_Int is mod <>;

package LCG_Rand
   with Spark_Mode => On
is 

   pragma Pure (LCG_Rand);
 
   pragma Assert (Parent_Random_Int'Modulus > 2**63-1);
   pragma Assert (m0 < Parent_Random_Int'Last); 
   pragma Assert (m1 < Parent_Random_Int'Last); 
   
   subtype Valid_LCG_0_Range is Parent_Random_Int range 1 .. m0 - 1;
   subtype Valid_LCG_1_Range is Parent_Random_Int range 1 .. m1 - 1;
 
   procedure Get_Random_LCG_64_0 (X0 : in out Parent_Random_Int);
   pragma Inline (Get_Random_LCG_64_0);
   -- x0 = 0 gives period of 1; needs to be rejected as a seed.
 
   procedure Get_Random_LCG_64_1 (X1 : in out Parent_Random_Int);
   pragma Inline (Get_Random_LCG_64_1);
   -- x1 = 0 gives period of 1; needs to be rejected as a seed.
 
   procedure Get_Random_LCG_64_Combined
     (S0       : in out Parent_Random_Int;
      S1       : in out Parent_Random_Int;
      Random_x :    out Parent_Random_Int); -- result
   -- S1 = 0 and S2 = 0 give reduced periods; need to be rejected as a seeds.
 
end LCG_Rand;
