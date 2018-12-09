
-----------------------------------------------------------------------
-- package body Spherical_Harmonics. Generates Spherical Harmonics.
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
--
with Clenshaw;
with A_Legendre;

package body Spherical_Harmonics is

   subtype Poly_Index is Base_Poly_Index range 0..Poly_Limit;
   
   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;

   package A_Legen is new A_Legendre (Real, Sqrt, Exp, Log, Base_Poly_Index);
 
   package Clen is new Clenshaw 
     (Real, Base_Poly_Index, Poly_Limit,
      A_Legen.Alpha, A_Legen.Beta, A_Legen.Q_0);
   use Clen;

   --
   -- Some global types: a table of Normalization
   -- factors, and present settings.
   --
   Norm_Factor_Table  : array (Poly_Index) of Real;
   m_Setting_Of_Table : Base_Poly_Index := Poly_Index'First - 1;
   Max_k_Of_Table     : Base_Poly_Index := Poly_Index'First - 1;
   
   -------------------------
   -- Make_Table_Of_Norms --
   -------------------------

   procedure Make_Table_Of_Norms 
     (m           : in Real;
      From_This_k : in Poly_Index;
      To_That_k   : in Poly_Index)
   is
   begin
     for k in From_This_k .. To_That_k loop
        Norm_Factor_Table (k) := A_Legen.Normalization_Factor (k, m);
     end loop;
   end Make_Table_Of_Norms;

   -------------------
   -- X_Lower_Bound --
   -------------------

   function X_Lower_Bound return Real is
   begin 
      return A_Legen.X_Lower_Bound;
   end;
   
   -------------------
   -- X_Upper_Bound --
   -------------------

   function X_Upper_Bound return Real is
   begin 
      return A_Legen.X_Upper_Bound;
   end;

   ----------------------
   -- Spherical_Harm_2 --
   ----------------------

   -- This is designed to be best when k changes a lot, and m
   -- changes rarely.  A table of Normalization factors for fixed
   -- m and varying k is made.
   --

   C0 : constant Coefficients := (others => Zero);  -- init. important

   function Spherical_Harm_2 
     (l          : Base_Poly_Index;
      m          : Real;
      X          : Real;
      Iterations : Positive := 1)      
      return Real 
   is
      Int_m : constant Poly_Index := Poly_Index (m);
      k     : constant Poly_Index := l - Int_m;
      Norm_Factor, Result : Real;
      Abs_m  : constant Real := Abs (m);
      Real_l : constant Real := +Real (l);
      C     : Coefficients; 
   begin
   
      if Abs_m > Real_l then  raise Constraint_Error;  end if;
      if Abs X > One    then  raise Constraint_Error;  end if;
    
      --  Designed to be efficient when X and k change a lot, and m changes rarely:
      
      if Int_m /= m_Setting_Of_Table then
         Make_Table_Of_Norms (m, From_This_k => Poly_Index'First, To_That_k => k);
         m_Setting_Of_Table := Int_m;
         Max_k_Of_Table     := k;
      end if;
 
      if k > Max_k_Of_Table then
         Make_Table_Of_Norms (m, From_This_k => Max_k_Of_Table, To_That_k => k);
         Max_k_Of_Table := k;
      end if;
 
      Norm_Factor := Norm_Factor_Table (k);
 
      C      := C0;
      C(k)   := One;
      Result := Norm_Factor * Sum (X, C, k, Abs_m, Iterations);
      --  Do have negative m in Spherical Harm., but
      --  don't have negative m in A_Legendre polys, so use Abs_m.
 
      return Result;

   end Spherical_Harm_2;
   
   --------------------
   -- Spherical_Harm --
   --------------------

   --  Uses calls to  Evaluate_Q, rather than to Sum.

   function Spherical_Harm
     (l          : Base_Poly_Index;
      m          : Real;
      X          : Real;
      Iterations : Positive := 1)      
      return Real 
   is
      Int_m : constant Poly_Index := Poly_Index (m);
      k     : constant Poly_Index := l - Int_m;
      Abs_m  : constant Real := Abs (m);
      Real_l : constant Real := +Real (l);
      Norm_Factor, Result : Real;
      Q      : Poly_Values; 
   begin
   
      if Abs_m > Real_l then  raise Constraint_Error;  end if;
      if Abs (X) > One  then  raise Constraint_Error;  end if;
    
      --  Designed to be efficient when X and k change a lot, and m changes rarely:
      
      if Int_m /= m_Setting_Of_Table then
         Make_Table_Of_Norms (m, From_This_k => Poly_Index'First, To_That_k => k);
         m_Setting_Of_Table := Int_m;
         Max_k_Of_Table     := k;
      end if;
 
      if k > Max_k_Of_Table then
         Make_Table_Of_Norms (m, From_This_k => Max_k_Of_Table, To_That_k => k);
         Max_k_Of_Table := k;
      end if;
 
      Norm_Factor := Norm_Factor_Table (k);
 
      Evaluate_Qs (X, Q, k, Abs_m, Iterations);
      Result := Norm_Factor * Q(k);
      --  Do have negative m in Spherical Harm., but
      --  don't have negative m in A_Legendre polys, so use Abs_m.

      return Result;

   end Spherical_Harm;

   -----------------
   -- Poly_Weight --
   -----------------
   
   function Poly_Weight (X : Real) return Real is
   begin
      return A_Legen.Poly_Weight(X);
   end Poly_Weight;

end Spherical_Harmonics;
