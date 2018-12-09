
with Text_IO; use Text_io;
with Chebychev_Quadrature;
with Ada.Numerics.Generic_Elementary_Functions;

procedure cheby_quad_tst_1 is

   type Real is digits 15;

   package Maths is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use Maths;

   No_Gauss_Pts : constant Positive := 2**16;

   package Quad is new Chebychev_Quadrature (Real, No_Gauss_Pts, Cos, Sin);
   use Quad;

   package rio is new Float_io(Real);
   use rio;

   X_First, X_Last, True_Area, True_Error : Real;
   Area : Real := 0.0;
   X_gauss : Gauss_Values;
   F_val   : Function_Values;

   procedure Pause (s1,s2,s3,s4,s5,s6,s7,s8 : string := "") is
      Continue : Character := ' ';
   begin
      new_line;
      if S1 /= "" then put_line (S1); end if;
      if S2 /= "" then put_line (S2); end if;
      if S3 /= "" then put_line (S3); end if;
      if S4 /= "" then put_line (S4); end if;
      if S5 /= "" then put_line (S5); end if;
      if S6 /= "" then put_line (S6); end if;
      if S7 /= "" then put_line (S7); end if;
      if S8 /= "" then put_line (S8); end if;

      dialog: loop
         begin
            New_Line;
            Put ("Type a character to continue: ");
            Get_Immediate (Continue);
            exit dialog;
         exception
            when others => null;
         end;
      end loop dialog;
      new_line;
   end pause;

begin

   X_First := 0.1;
   X_Last  := 0.1777;

   New_Line;
   Pause ("Test Chebychev quadrature on sinusoids.");

   Sin_Test: for Pow in 1..30 loop

      Find_Gauss_Nodes (X_First, X_Last, X_gauss);

      for I in Gauss_Index loop
         F_val(I) := Sin (X_gauss(I));
      end loop;

      Get_Integral(F_val, X_First, X_Last, Area);

      True_Area  := -Cos(X_Last) + Cos(X_First);
      True_Error := Area - True_Area;

      new_line(1); put ("True_Error:  ");
      put(True_Error);

      X_Last := X_Last + 1.7;

   end loop Sin_Test;

   X_First := -1.3;
   X_Last  := -0.7777;

   New_Line;
   Pause ("Test Chebychev quadrature on the Exp function.");

   Exp_Test: for Pow in 1..30 loop

    --X_Last  := X_Last  + 3.7;
      X_First := X_First - 1.7;

      Find_Gauss_Nodes (X_First, X_Last, X_gauss);

      for I in Gauss_Index loop
         F_val(I) := Exp (X_gauss(I));
      end loop;

      Get_Integral(F_val, X_First, X_Last, Area);

      True_Area  := Exp (X_Last) - Exp (X_First);
      True_Error := (Area - True_Area) / Area;

      new_line(1); put ("True_Error:  ");
      put(True_Error);

   end loop Exp_Test;


   X_First := -1.3;
   X_Last  := -0.7777;

   New_Line;
   Pause ("Second test of Chebychev quadrature on the Exp function.");

   Exp2_Test: for Pow in 1..30 loop

      X_Last  := X_Last  + 1.7;
      --X_First := X_First - 4.7;

      Find_Gauss_Nodes (X_First, X_Last, X_gauss);

      for I in Gauss_Index loop
         F_val(I) := Exp (X_gauss(I));
      end loop;

      Get_Integral(F_val, X_First, X_Last, Area);

      True_Area  := Exp (X_Last) - Exp (X_First);
      True_Error := (Area - True_Area) / Area;

      new_line(1); put ("True_Error:  ");
      put(True_Error);

   end loop Exp2_Test;

   X_First := 0.3;
   X_Last  := 0.7777;
   new_line(2);

   New_Line;
   Pause ("Test Chebychev quadrature on the Log function.");

   Log_Test: for Pow in 1..30 loop

      X_First := X_First / 1.2;

      Find_Gauss_Nodes (X_First, X_Last, X_gauss);

      for I in Gauss_Index loop
         F_val(I) := 1.0 / (X_gauss(I));
      end loop;

      Get_Integral(F_val, X_First, X_Last, Area);

      True_Area  := Log (X_Last) - Log (X_First);
      True_Error := (Area - True_Area) / Area;

      new_line(1); put ("True_Error:  ");
      put(True_Error);

   end loop Log_Test;

   X_First := 1.0;
   X_Last  := 2.7;
   new_line(2);

   New_Line;
   Pause ("Second test of Chebychev quadrature on the Log function.");

   Log2_Test: for Pow in 1..30 loop

      X_Last := X_Last  + 2.0;

      Find_Gauss_Nodes (X_First, X_Last, X_gauss);

      for I in Gauss_Index loop
         F_val(I) := 1.0 / (X_gauss(I));
      end loop;

      Get_Integral(F_val, X_First, X_Last, Area);

      True_Area  := Log (X_Last) - Log (X_First);
      True_Error := (Area - True_Area) / Area;

      new_line(1); put ("True_Error:  ");
      put(True_Error); 

   end loop Log2_Test;
end cheby_quad_tst_1;
