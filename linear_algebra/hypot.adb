
---------------------------------------------------------------------------
-- package body Hypot
-- Copyright (C) 2018 Jonathan S. Parker.
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

with Ada.Numerics;
with Ada.Numerics.Generic_Elementary_Functions;

package body Hypot is

   package Math is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use Math;

   Zero : constant Real := +0.0;
   Half : constant Real := +0.5;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;

   -- uses (-1 + sqrt(1+t*t) * (1 + sqrt(1+t*t) = t*t
   --
   -- sqrt(1+t*t) = 1  +  t*t / (1 + sqrt(1+t*t)
   --
   -- if t = a/b, and b > a, then
   --
   -- sqrt(a^2 + b^2) = |b| + a*a / (b + sqrt(a^2 + b^2))
   -- sqrt(a^2 + b^2) = |b| + correction
   --
   -- Hypot              =       sqrt(a^2 + b^2) = |b|   + correction
   -- Hypot_Plus_Max_Arg = |b| + sqrt(a^2 + b^2) = 2*|b| + correction
   --
   -- Min_Arg_Over_Hypot := a/b - a*(a**2/(b*hypot*Hypot_Plus_Max_Arg));
   -- (not gd near a~b, where a/Hypot is better. But a/Hypot is better always.)
   --
   -- Idea is to get Hypot_Plus_Max_Arg to a little better precision if
   -- possible using the 2*|b|.
   --
   -- Let a = |max|, b = |min|.
   -- Need better method when a ~ b.
   -- Expand about
   --
   --    Sqrt(0.5*(Max-Min)**2 + 0.5*(Min + Max)**2) = Sqrt(0.5*Diff**2 + 0.5*Sum**2)
   --
   -- and formulate it so you add correction terms to Min:
   --
   --    Hypot := Min*Sqrt_2 +
   --               Dif*Reciprocal_Sqrt_2 + Dif**2 / (Sqrt_2*Sum + 2.0*Hypot)))
   --
   -- Below we add small terms to each other first in order to reduce roundoff,
   -- (by using Sqrt_2 = 1 + Sqrt_2_minus_1).
   --
   -- To get Min / Hypot when min ~ max:
   --
   --    (1/sqrt(1+e**2) - 1/Sqrt_2) * (1/sqrt(1+e**2) - 1/Sqrt_2) =
   --          (1 - e**2)/(2*sqrt(1+e**2) + Sqrt_2*((1+e**2)))
   --
   -- where e = max / min
   --
   --  1/sqrt(1+e**2) = 0.5 + ((1/Sqrt_2-0.5) + (min-max)*(min+max)/(2max*hypot + Sqrt_2*hypot**2))
   --

   procedure Get_Hypotenuse
     (a, b                       : in     Real;
      Hypot                      :    out Real;
      Min_Arg_Over_Hypot         :    out Real;
      Max_Arg_Over_Hypot_Minus_1 :    out Real)
   is
      Emin : constant Integer := Real'Machine_Emin;
      Emax : constant Integer := Real'Machine_Emax;

      Max, Min, Min_Squared, Max_Squared : Real;
      M_Exp                              : Integer;
      Scale, Un_Scale                    : Real := One;
      Arg, Correction                    : Real;
      Three_Quarters_Max, Reciprocal_Min : Real;
   begin
      if (not a'Valid) or else (not b'Valid) then
         raise Constraint_Error with "Invalid input data in Hypot";
      end if;

      if abs b >= abs a then       -- need >= !!!
         Max := abs b;
         Min := abs a;
      else
         Max := abs a;
         Min := abs b;
      end if;

      -- defaults for out params

      Hypot                      := Max;
      Min_Arg_Over_Hypot         := Zero;
      Max_Arg_Over_Hypot_Minus_1 := Zero;

      -- force Max to underflow to zero at Two**Emin (~1.0E-307 for 15 digits)
      -- Won't get the best answer here, but it prevents some NaNs.

      if Max < Two**Emin then
         return; -- use defaults: hypot = Max
      end if;

      -- force Min to underflow to zero at Two**Emin (~1.0E-307 for 15 digits)
      -- Won't get the best answer here, but it prevents some NaNs.

      if Min < Two**Emin then    -- set Min = 0
         return; -- use defaults: hypot = Max
      end if;

      -- Case: Sqrt (Max**2 + Min**2) = Max.  (i.e. (1 + (Min/Max)**2)**0.5 = 1)
      -- The following also deals with Min = 0.
      -- (Dealt with Max ~ 0 above, so can divide by Max.)

      if Max > Two**Emin * Two**(+Real'Machine_Mantissa / 2 + 10) then
      if Min <       Max * Two**(-Real'Machine_Mantissa / 2 - 10) then
         Hypot                      :=  Max;
         Min_Arg_Over_Hypot         :=  Min / Hypot;
         Max_Arg_Over_Hypot_Minus_1 := -Half * Min_Arg_Over_Hypot**2;
         -- - Min**2 / (Max*(Max+Hypot)) = - Min**2 / (2*max**2)
         return;
      end if;
      end if;

      -- Scale the arguments so we can square them w/o under/overflow.
      -- Above we made sure that Min and Max have similar exponents.

      if Max > Two**(Emax / 4)  or else  Min < Two**(Emin / 4) then
         M_Exp    := Real'Exponent (Max) - 2;
         Scale    := Two**(-M_Exp);
         Un_Scale := Two**(M_Exp);

         Max := Scale * Max;
         Min := Scale * Min;
      else
         Scale    := One;
         Un_Scale := One;
      end if;

      -- Get hypotenuse: formula Sqrt(Min**2 + Max**2) not bst if fpmath=387

      Min_Squared        := Min**2;
      Max_Squared        := Max**2;
      Arg                := Min_Squared + Max_Squared;
      Reciprocal_Min     := One / Min;
      Three_Quarters_Max := 0.75 * Max;

      Get_Hypot : declare
         Sq             : constant Real := Sqrt (Arg);
         Two_Max        : constant Real := Two * Max;
         Sum, Lost_Bits : Real          := Zero;
      begin

         Correction := Min_Squared / (Max + Sq);  -- Sq is first estimate of Hypot.
         Hypot      := Max + Correction;
         -- Hypot should be an improved Sq, so replace Sq with Hypot and iterate:
         --     Correction := Min_Squared / (Max + Max + Correction);
         --
         -- Max_x_Corr := Min_Squared / (2.0 + Min_Squared / (Max_Squared + Hypot*Max));

         if 1.5 * Min > Max then -- iterate
            Sum        := Two_Max + Correction;
            Lost_Bits  := Correction - (Sum - Two_Max);
            Correction := Min_Squared / Sum;
            Correction := Correction - Lost_Bits * (Correction * Reciprocal_Min)**2;
            Hypot      := Max + Correction;
         end if;

         -- On 387 the following is better:
         -- Sqrt(Min**2 + Max**2) - Sqrt((Max + Min**2/(2Max))**2) with
         -- Sqrt(x) - Sqrt(y) = (x - y) /  (Sqrt(x) + Sqrt(y)) to get,
         --
         --Term       := Half * Min**2 / Max;
         --Correction := Term - Term**2 / ((Max + Term) + Sq);
         --   expand out hypot = max + correction:
         -- Max_x_Corr := 0.5*Min**2 - 0.25*(Min**2)**2/(((Max**2 + 0.5*Min**2) + Sq*Max));
         -- Correction := Max_x_Corr/Max;
         -- Hypot      := Max + Correction;

      end Get_Hypot;

      --
      -- Cos_Minus_1 = Max_Arg_Over_Hypot_Minus_1 = |Max| / Sqrt(Max**2 + Min**2) - 1
      --
      -- Cos_Minus_1 = -Min2_Over_Hypot_Max = - Min**2 / (Hypot * (Hypot + Max));
      -- Cos_Minus_1 = -Min2_Over_Hypot_Max = - Min**2 / (Min**2 + 2*Max**2 + Max*Correction);
      --   (uses Hypot = Max + Correction)

      Get_Cos_Minus_1 : declare
         Two_Max_Squared          : constant Real := Two * Max_Squared;
         Max_x_Corr, Sum, Denom   : Real;
         Lost_Bits_2, Lost_Bits_1 : Real;
         Min2_Over_Hypot_Max      : Real;
      begin

         --  Max*Correction  <  Min**2  <  2*Max**2

         Max_x_Corr := Max * Correction;

         Sum         := Min_Squared + Two_Max_Squared;
         Lost_Bits_1 := Min_Squared - (Sum - Two_Max_Squared);
         -- Sum = Small + Large.  Lost_Bits = Small - (Sum - Large)

         Denom       := Sum + (Max_x_Corr + Lost_Bits_1);
         Lost_Bits_2 := Max_x_Corr - (Denom - Sum);

         Min2_Over_Hypot_Max := Min_Squared / Denom;

         Min2_Over_Hypot_Max := Min2_Over_Hypot_Max -
           (Min2_Over_Hypot_Max * Reciprocal_Min)**2 * (Lost_Bits_2 + Lost_Bits_1);

         Max_Arg_Over_Hypot_Minus_1 := -Min2_Over_Hypot_Max;

      end Get_Cos_Minus_1;

      -- Use Hypot to get  Min / Hypot

      Get_Min_over_Hypot : declare
         Sqrt_2 : constant Real := 1.41421_35623_73095_04880_16887_24210;
         Sqrt_Half_Minus_Half : constant Real := 0.20710_67811_86547_52440_08443_62105;
         Sum, Dif  : Real;
         Lost_Bits : Real;
      begin

         if Min > Three_Quarters_Max then    -- diff becomes small.

            Sum                := Max + Min;
            Dif                := Max - Min;
            Min_Arg_Over_Hypot := Half +
              (Sqrt_Half_Minus_Half - (Dif * Sum) / (Two * Min * Hypot + Sqrt_2 * Arg));

         else

            Min_Arg_Over_Hypot := Min / Hypot;   -- now unscaled

            -- Following helps on sse if Max/Min ~ ... 1-1.25, 2-2.5, 4-5, 8-10, ...
            -- Uses Hypot = Max + Correction
            -- Sum = Small + Large.  Lost_Bits = Small - (Sum - Large)

            Lost_Bits := Correction - (Hypot - Max);

            Min_Arg_Over_Hypot :=
              Min_Arg_Over_Hypot -
              Min_Arg_Over_Hypot**2 * Reciprocal_Min * Lost_Bits;

         end if;

      end Get_Min_over_Hypot;

      Hypot := Un_Scale * Hypot;

   end Get_Hypotenuse;

   function Hypotenuse (a, b : in Real) return Real is
      Emin            : constant Integer := Real'Machine_Emin;
      Emax            : constant Integer := Real'Machine_Emax;
      Max, Min        : Real;
      M_Exp           : Integer;
      Scale, Un_Scale : Real := One;
      Result          : Real := Zero;
      Min_Squared, Max_Squared, Arg : Real;
      Correction, Reciprocal_Min : Real := 0.0;
   begin
      if (not a'Valid) or else (not b'Valid) then
         raise Constraint_Error with "Invalid input data in Hypot";
      end if;

      if abs b >= abs a then       -- need >= !!!
         Max := abs b;
         Min := abs a;
      else
         Max := abs a;
         Min := abs b;
      end if;

      -- force Max to underflow to zero near the threshold (~ 1.0E-307 for 15 digits)
      -- Won't get best answer here, but it prevents some NaNs.

      if Max < Two**Emin then
         return Max; -- use defaults: hypot = Max
      end if;

      -- force Min to underflow to zero near the threshold (~ 1.0E-307 for 15 digits)

      if Min < Two**Emin then    -- set Min = 0
         return Max; -- use defaults: hypot = Max
      end if;

      -- Case: Sqrt (Max**2 + Min**2) = Max.  (i.e. (1 + (Min/Max)**2)**0.5 = 1)

      if Max > Two**Emin * Two**(+Real'Machine_Mantissa / 2 + 10) then
      if Min <       Max * Two**(-Real'Machine_Mantissa / 2 - 10) then
         return Max;
      end if;
      end if;

      -- Scale the arguments so we can square them w/o under/overflow.
      -- Above we made sure that Min and Max have similar exponents.

      if Max > Two**(Emax / 4)  or else  Min < Two**(Emin / 4) then
         M_Exp    := Real'Exponent (Max) - 2;
         Scale    := Two**(-M_Exp);
         Un_Scale := Two**(M_Exp);

         Max := Scale * Max;
         Min := Scale * Min;
      else
         Scale    := One;
         Un_Scale := One;
      end if;

      Min_Squared    := Min**2;
      Max_Squared    := Max**2;
      Arg            := Min_Squared + Max_Squared;
      Reciprocal_Min := One / Min;

      Get_Hypot : declare
         Sq             : constant Real := Sqrt (Arg);
         Two_Max        : constant Real := Two * Max;
         Sum, Lost_Bits : Real          := Zero;
      begin

         Correction := Min_Squared / (Max + Sq);  -- Sq is first estimate of Hypot.
         Result     := Max + Correction;
         -- sq ~ hypot = max+correction, so iteration improves the estimate.

         if 1.5 * Min > Max then -- iterate
            Sum        := Two_Max + Correction;
            Lost_Bits  := Correction - (Sum - Two_Max);
            Correction := Min_Squared / Sum;
            Correction := Correction - Lost_Bits * (Correction * Reciprocal_Min)**2;
            Result     := Max + Correction;
         end if;

         Result := Un_Scale * Result;

      end Get_Hypot;

      return Result;
   end Hypotenuse;

end Hypot;
