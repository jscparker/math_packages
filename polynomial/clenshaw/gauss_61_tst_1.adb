
with Text_IO; use Text_io;
with Gauss_Quadrature_61;
with Ada.Numerics.Generic_Elementary_Functions;

procedure gauss_61_tst_1 is

   type Real is digits 15;

   package Math is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use Math;

   package Quad is new Gauss_Quadrature_61 (Real);
   use Quad;

   package rio is new Float_io(Real);
   use rio;

   X_First, X_Last, True_Area, True_Error : Real;
   Area, Rough_Error : Real := 0.0;
   X_gauss   : Gauss_Values;
   F_val : Function_Values;

begin

   Gauss_61_Coeff_Test; -- look for mutated constants among the Gauss nodes.

   X_First := 0.1;
   X_Last  := 1.7777;
   new_line(1); put("Sin Test");
   new_line(1); 

   Sin_Test: for Pow in 1..30 loop

      X_Last := X_Last + 1.7;

      Find_Gauss_Nodes (X_First, X_Last, X_gauss);

      for I in Gauss_Index_61 loop
         F_val(I) := Sin(X_gauss(I));
      end loop;

      Get_Integral(F_val, X_First, X_Last, Area, Rough_Error);

      True_Area  := -Cos(X_Last) + Cos(X_First);

      True_Error := Area - True_Area;

      new_line(1); put ("True_Error, Rough Estimate::  ");
      put(True_Error);  put(" "); put(Rough_Error);

   end loop Sin_Test;

   X_First := -1.3;
   X_Last  := -0.7777;
   new_line(2); put("Exp Test");
   new_line(1); 

   Exp_Test: for Pow in 1..30 loop

      --X_Last  := X_Last  + 3.7;
      X_First := X_First - 1.7;

      Find_Gauss_Nodes (X_First, X_Last, X_gauss);

      for I in Gauss_Index_61 loop
         F_val(I) := Exp (X_gauss(I));
      end loop;

      Get_Integral(F_val, X_First, X_Last, Area, Rough_Error);

      True_Area  := Exp (X_Last) - Exp (X_First);

      True_Error := (Area - True_Area) / Area;

      new_line(1); put ("True_Error, Rough Estimate::  ");
      put(True_Error);  put(" "); put(Rough_Error);

   end loop Exp_Test;


   X_First := -1.3;
   X_Last  := -0.7777;
   new_line(2); put("Exp Test");
   new_line(1); 

   Exp2_Test: for Pow in 1..30 loop

      X_Last  := X_Last  + 1.7;
      --X_First := X_First - 4.7;

      Find_Gauss_Nodes (X_First, X_Last, X_gauss);

      for I in Gauss_Index_61 loop
         F_val(I) := Exp (X_gauss(I));
      end loop;

      Get_Integral(F_val, X_First, X_Last, Area, Rough_Error);

      True_Area  := Exp (X_Last) - Exp (X_First);
      True_Error := (Area - True_Area) / Area;

      new_line(1); put ("True_Error, Rough Estimate::  ");
      put(True_Error);  put(" "); put(Rough_Error);

   end loop Exp2_Test;

   X_First := 0.3;
   X_Last  := 0.7777;
   new_line(2); put("Log Test");
   new_line(1); 

   Log_Test: for Pow in 1..30 loop

      X_First := X_First / 1.2;

      Find_Gauss_Nodes (X_First, X_Last, X_gauss);

      for I in Gauss_Index_61 loop
         F_val(I) := 1.0 / (X_gauss(I));
      end loop;

      Get_Integral(F_val, X_First, X_Last, Area, Rough_Error);

      True_Area  := Log (X_Last) - Log (X_First);
      True_Error := (Area - True_Area) / Area;

      new_line(1); put ("True_Error, Rough Estimate::  ");
      put(True_Error);  put(" "); put(Rough_Error);

   end loop Log_Test;

   X_First := 1.0;
   X_Last  := 2.7;
   new_line(2); put("Reciprocal_Test");
   new_line(1); 

   Reciprocal_Test: for Pow in 1..30 loop

      X_Last  := X_Last + 2.0;

      Find_Gauss_Nodes (X_First, X_Last, X_gauss);

      for I in Gauss_Index_61 loop
         F_val(I) := 1.0 / (X_gauss(I));
      end loop;

      Get_Integral(F_val, X_First, X_Last, Area, Rough_Error);

      True_Area  := Log (X_Last) - Log (X_First);
      True_Error := (Area - True_Area) / Area;

      new_line(1); put ("True_Error, Rough Estimate::  ");
      put(True_Error);  put(" "); put(Rough_Error);

   end loop Reciprocal_Test;

end gauss_61_tst_1;
