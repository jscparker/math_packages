
-----------------------------------------------------------------------
-- package body Chebychev_Quadrature. Coefficients for Chebychev-Gaussian quadrature.
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

package body Chebychev_Quadrature is

   Pii : constant := 3.14159_26535_89793_23846_26433_83279_50288_41971_694;
        
   Gauss_Root   : Gauss_Values := (others => 0.0);
   Gauss_Weight : Gauss_Values := (others => 0.0);
   --  Arrays initialized after the "begin" below.

   ---------------------------
   -- Construct_Gauss_Roots --
   ---------------------------

   procedure Construct_Gauss_Roots is
      Factor : constant Real := Pii / Real (Gauss_Index'Last);
   begin
      for i in Gauss_Index loop
         Gauss_Root (i) := -Cos (Factor * (Real (i) - 0.5));
      end loop;
   end Construct_Gauss_Roots;

   -----------------------------
   -- Construct_Gauss_Weights --
   -----------------------------

   procedure Construct_Gauss_Weights is
      Factor : constant Real := Pii / Real (Gauss_Index'Last);
   begin
      for i in Gauss_Index loop
         Gauss_Weight (i) := 0.5 * Factor * Abs (Sin (Factor * (Real (i) - 0.5)));
      end loop;
   end Construct_Gauss_Weights;

   ----------------------
   -- Find_Gauss_Nodes --
   ----------------------

   procedure Find_Gauss_Nodes
     (X_Starting : in     Real;
      X_Final    : in     Real;
      X_gauss    :    out Gauss_Values)
   is
      Half_Delta_X : constant Real := (X_Final - X_Starting) * 0.5;
      X_Center     : constant Real := X_Starting + Half_Delta_X;
   begin
      for i in Gauss_Index loop
         X_gauss(i) := X_Center + Gauss_Root(i) * Half_Delta_X;
      end loop;
   end Find_Gauss_Nodes;

   ------------------
   -- Get_Integral --
   ------------------
 
   procedure Get_Integral
     (F_val       : in     Function_Values;
      X_Starting  : in     Real;
      X_Final     : in     Real;
      Area        :    out Real)
   is
      Sum : Real;
      Delta_X : constant Real := (X_Final - X_Starting);
   begin
      Area := 0.0;
      Sum  := 0.0;
      for i in Gauss_Index loop
         Sum := Sum + Gauss_Weight (i) * F_val (i);
      end loop;
      Area := Sum * Delta_X;
   end Get_Integral;
 
begin

   Construct_Gauss_Roots;
   Construct_Gauss_Weights;

end Chebychev_Quadrature;
