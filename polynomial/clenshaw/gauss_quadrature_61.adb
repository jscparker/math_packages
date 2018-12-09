
-----------------------------------------------------------------------
-- package body Gauss_Quadrature_61. Gaussian quadrature.
-- Copyright (C) 2018 Jonathan S. Parker
--
-- Permission to use, copy, modify, and/or distribute this software for any
-- purpose with or without fee is hereby granted, provided that the above
-- copyright notice and this permission notice appear in all copies.
-- THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
-- WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
-- MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
-- ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
-- WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
-- ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
-- OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
---------------------------------------------------------------------------

with Gauss_Nodes_61, Gauss_Nodes_61_bak;

package body Gauss_Quadrature_61 is

  package Nodes is new Gauss_Nodes_61 (Gauss_Index_61, Real, Gauss_Values);
  use Nodes;

  package Test_Nodes is new Gauss_Nodes_61_bak (Gauss_Index_61, Real, Gauss_Values);
  use Test_Nodes;

  -------------------------
  -- Gauss_61_Coeff_Test --
  -------------------------

  --  Use redundant gauss coeffs to search for mutated coefficients.

  procedure Gauss_61_Coeff_Test is
     Error1, Error2, Err : Real := 0.0;
  begin 
     for I in Gauss_Index_61 loop
        Error1 := Abs (Gauss_Weights_61 (I) - Gauss_Weights_61_bak (I));
        Error2 := Abs (Gauss_Roots_61 (I) - Gauss_Roots_61_bak (I));
        if (abs(Error1) > Err) then
           Err := Abs (Error1);
        end if;
        if (abs(Error2) > Err) then
           Err := Abs (Error2);
        end if;
     end loop;
   
     for I in Gauss_Index_30 loop
        Error1 := Abs (Gauss_Weights_30 (I) - Gauss_Weights_30_bak (I));
        if (abs(Error1) > Err) then
           Err := Abs (Error1);
        end if;
     end loop;
     if (Err > Real'Epsilon * 2.0) then
        raise Program_Error with "Gaussian coefficients are wrong.";
     end if;
  end Gauss_61_Coeff_Test;

  -----------------------
  -- Find_Gauss_Nodes --
  -----------------------

  Gauss_Roots : Gauss_Values renames Gauss_Roots_61;

  procedure Find_Gauss_Nodes
    (X_Starting : in     Real;
     X_Final    : in     Real;
     X_gauss    :    out Gauss_Values)
  is
     Abs_Half_Delta_X : constant Real := Abs (X_Final - X_Starting) * 0.5;
     Central_X        : constant Real := (X_Final + X_Starting) * 0.5;
  begin
     for I in Gauss_Index_61 loop
	X_gauss(I) := Central_X + Gauss_Roots(I) * Abs_Half_Delta_X;
     end loop;
  end Find_Gauss_Nodes;

  ------------------
  -- Get_Integral --
  ------------------

  procedure Get_Integral
    (F_val       : in     Function_Values;
     X_Starting  : in     Real;
     X_Final     : in     Real;
     Area        :    out Real;
     Rough_Error :    out Real)
  is
     Area_30, Area_61, Area_abs : Real;
     Sum, Error_Est, Area_scaler, Scaled_Error, F_mean: Real;
     Half_Delta_X     : constant Real := (X_Final - X_Starting) * 0.5;
     Abs_Half_Delta_X : constant Real := Abs (Half_Delta_X);

     function Min(X,Y : Real) return Real is
     begin
       if X < Y then return X; else return Y; end if;
     end Min;

     function Max(X,Y : Real) return Real is
     begin
       if X > Y then return X; else return Y; end if;
     end Max;

  begin

     --  Get the 30-point Gaussian Quadrature:

     Sum := 0.0;
     for I in Gauss_Index_30 loop
        Sum := Sum + Gauss_Weights_30 (I) * F_val (2*I-1);
     end loop;
     Area_30 := Sum * Half_Delta_X;

     --  Calculate the 61-pt formula

     Sum := 0.0;
     for I in Gauss_Index_61 loop
        Sum := Sum + Gauss_Weights_61 (I) * F_val (I);
     end loop;
     F_mean  := Sum * 0.5;
     Area_61 := Sum * Half_Delta_X;

     --  Calculate Rough_Error:

     Sum := 0.0;
     for I in Gauss_Index_61 loop
        Sum := Sum + Gauss_Weights_61 (I) * Abs (F_val (I) - F_mean);
     end loop;
     Area_scaler := Sum * Abs_Half_Delta_X;

     Sum := 0.0;
     for I in Gauss_Index_61 loop
        Sum := Sum + Gauss_Weights_61 (I) * Abs (F_val (I));
     end loop;
     Area_abs := Sum * Abs_Half_Delta_X;

     Error_Est := Abs (Area_61 - Area_30);

     if Error_Est > 0.0 and Area_scaler > 0.0 then
        Scaled_Error := Error_Est / Area_scaler;
        Error_Est := Scaled_Error * Min (1.0, (Sqrt (200.0 * Scaled_Error))**3);
     end if;
     if Area_abs > Real'Small / (Real'Epsilon * 8.0) then
        Error_Est := Max (Error_Est, Area_abs * Real'Epsilon * 2.0);
     end if;
     --  Very rough fit to true error

     Rough_Error := Error_Est;
     Area        := Area_61;

  end Get_Integral;

begin

  Gauss_61_Coeff_Test; -- Search for typos in the tables of coefficients.

end Gauss_Quadrature_61;
