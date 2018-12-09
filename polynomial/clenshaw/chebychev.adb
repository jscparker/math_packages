
-----------------------------------------------------------------------
-- package body Chebychev. Data structure for generating Chebychev functions.
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
package body Chebychev is

   Half_Pi : constant Real := 0.5 * 3.14159_26535_89793_23846_26433_83279_502884;
   --  for norm.

   Two  : constant := 2.0;
   One  : constant := 1.0;
   Zero : constant := 0.0;
   Half : constant := 0.5;

   -------------------
   -- Log_of_1_plus --
   -------------------

   --  More accurate than using Log (1 + x) for Abs x << 1.

   function Log_of_1_plus
      (x : Real)
      return Real
   is
      u   : constant Real := One + x;
      Sqrt_Eps : constant Real := Two ** (-Real'Machine_Mantissa / 2 - 6);
   begin
      if u <= Zero then
         raise Constraint_Error;
      end if;

      -- Use Log(1+x) = x - x^2/2 + x^3/3 - ... = x*(1 - x/2 + x^2/3 ...)
      -- So if x is somewhat smaller than Sqrt (Real'Epsilon) then 1 + x^3/3 = 1:

      if Abs (u - One) < Sqrt_Eps then
         return x - Half * x * x;
      end if;

    --return Log(u) * x / (u - One); -- more accurate?  (u /= One; see above).
      return Log(u);  

   end Log_of_1_plus;

   function Alpha (k : Poly_ID_Integer; m : Real := 0.0; X : Real) return Real
   is
   begin 
      return Two * X;
   end;

   function Beta (k : Poly_ID_Integer; m : Real := 0.0; X : Real) return Real
   is
   begin 
      return -One;
   end;

   function Q_0 (X : Real; m : Real := 0.0) return Real
   is
   begin 
      return One;
   end;

   function Q_1 (X : Real; m : Real := 0.0) return Real 
   is
   begin 
      return Two * X;
   end;

   --  Sqrt (1.0 - X*X)

   function Poly_Weight (X : Real) return Real 
   is
      Result, Arg : Real;
   begin 

      if Abs (X) = One then return Zero;            end if;
      if Abs (X) > One then raise Constraint_Error; end if;

      Arg := 0.5 * Log_of_1_plus(-X*X); -- Arg always < 0

      if Abs (Arg) > 600.0 then  -- Exp(-600.0) ~  10^(-262)
         Result := Zero;
      else
         Result := Exp (Arg);
      end if;

      return Result;

   end Poly_Weight;

   function Norm (k : Poly_ID_Integer; m : Real := 0.0) return Real 
   is
   begin
       return Half_Pi;
   end;
end Chebychev;
