with Orthogonal_Polys;
with text_io; use text_io;
with Ada.Numerics.Generic_Elementary_Functions;

procedure opolys_tst_2 is  -- fit data with a nth order polynomial

   type Real is digits 15;

   package Math is new Ada.Numerics.Generic_Elementary_Functions (Real);

   function Sin (X : Real) return Real renames Math.Sin;
   function Cos (X : Real) return Real renames Math.Cos;

   X_Start : constant Real := 0.0;
   X_End   : constant Real := 8.0;
   type array_Index is range 0..63;
   subtype Points_Index is array_Index range 0..32;
   --  This determines the maximum storage size of the Data array.
   --  the Data array will start at Points_Index'First and end
   --  at Points_Index'Last.
   --  Least squares fitting can be done on any sub interval of
   --  this array by setting the obvious parameters in the procedure
   --  below.

    type Data is array (array_Index) of Real;

    package LS is new Orthogonal_Polys (Real, array_Index, Data, 0.0, 1.0);
    --                          "*", "+", "-", "/", "<", "**", "-");
    use LS;

    ls_order : constant Coeff_Index := 24;

    X_axis  : Data;
    Y_Data  : Data;
    Weights : constant Data := (others => 1.0);

    Delta_X : Real;

    -- Variables for IO:

    X0, Y0: Real;
    X1, Y1 : Real;

    E              : Real;
    Best_Fit       : Poly_Data;
    Poly_Set       : Polynomials;
    C              : Poly_Sum_Coeffs;

    First   : constant Points_Index := Points_Index'First;
    Last    : constant Points_Index := Points_Index'Last;

    package rio is new Text_IO.FLoat_IO(Real);
    use rio;

    ---------
    -- Sum --    
    ---------
    
    function Sum (X : Real; 
		 Coeffs : Powers_Of_X_Coeffs) return Real is
      Result : Real := 0.0;
    begin
      Result := Coeffs(Coeff_Index'Last);
      for I in reverse Coeff_Index'First..Coeff_Index'Last-1 loop
	 Result := Result * X;
	 Result := Result + Coeffs(I);
      end loop;

      return Result;
   end Sum;

   -----------
   -- Pause --    
   -----------

   procedure Pause (s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11 : string := "") is
      Continue : Character := ' ';
   begin
      new_line;
      if S0 /= "" then put_line (S0); end if;
      if S1 /= "" then put_line (S1); end if;
      if S2 /= "" then put_line (S2); end if;
      if S3 /= "" then put_line (S3); end if;
      if S4 /= "" then put_line (S4); end if;
      if S5 /= "" then put_line (S5); end if;
      if S6 /= "" then put_line (S6); end if;
      if S7 /= "" then put_line (S7); end if;
      if S8 /= "" then put_line (S8); end if;
      if S9 /= "" then put_line (S9); end if;
      if S10 /= "" then put_line (S10); end if;
      if S11 /= "" then put_line (S11); end if;
      new_line;
      begin
	 Put ("Enter a character to continue: ");
	 Get_Immediate (Continue); New_Line;
      exception
	 when others => null;
      end;
   end pause;


 begin

    Delta_X := (X_end - X_Start) / Real (X_axis'Length - 1);
    for I in Points_Index loop
      X_axis(I) := X_Start + Real (I - Points_Index'First) * Delta_X;
    end loop;

    --  Want the wgts to go to zero at the points that are just before and
    --  just beyond the end points, so Diameter of the circle-wgt is greater
    --  than X_end - X_start by 2 * Delta_X.

    --     Diameter := X_end - X_start + 2.0 * Delta_X;
    --     Radius   := Diameter / 2.0;
    --     middle   := (X_end + X_start) / 2.0;
    --     for I in Points_Index loop
    --        Weights(I) := SQRT (Radius*Radius - (X_axis(I) - middle)**2);
    --     end loop;

    for I in Points_Index loop
       Y_Data(I) := Sin (X_axis(I));
    end loop;

    Poly_Fit (Y_Data, X_axis, Weights, First, Last, ls_order,
	     Best_Fit, C, Poly_Set, E);

    pause("Test 1: Make sinusoidal data, fit a high order polynomial to it,",
	  "and subtract the Integral of the fit from the Integral", 
	  "of the original data.");

     for I in Points_Index'first+1 .. Points_Index'last loop

	X0 := X_axis(I-1);
	Y0 := Poly_Integral (X0, C, Poly_Set);
	X1 := X_axis(I);
	Y1 := Poly_Integral (X1, C, Poly_Set);
	new_line; put (Y1 - Y0 - (-Cos (X1) + Cos (X0)) ); put ("  ");

     end loop;

end opolys_tst_2;
