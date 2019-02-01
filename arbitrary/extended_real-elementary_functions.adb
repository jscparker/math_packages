
-----------------------------------------------------------------------
-- package body Extended_Real.Elementary_Functions, extended precision functions
-- Copyright (C) 2008-2018 Jonathan S. Parker
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

package body Extended_Real.Elementary_Functions is

   Max_Available_Bits : constant e_Integer
                       := Desired_No_Of_Bits_in_Radix * e_Real_Machine_Mantissa;
   -- This equals:  Bits_per_Mantissa = Bits_per_Digit * No_of_Digits_per_Mantissa
   -- (e_Real_Machine_Mantissa = Mantissa'Length = No_of_Digits_per_Mantissa).

   --  The following are frequently used global constants:

   Radix : constant Real := e_Real_Machine_Radix;

   Half_Digit   : constant E_Digit := Make_E_Digit (0.5);
   Two_Digit    : constant E_Digit := Make_E_Digit (2.0);
   Three_Digit  : constant E_Digit := Make_E_Digit (3.0);
   Four_Digit   : constant E_Digit := Make_E_Digit (4.0);
   Nine_Digit   : constant E_Digit := Make_E_Digit (9.0);
   Twelve_Digit : constant E_Digit := Make_E_Digit (12.0);

 --One   : constant e_Real := +1.0; -- in Extended_Real spec
 --Zero  : constant e_Real := +0.0;
   Two   : constant e_Real := +2.0;
   Three : constant e_Real := +3.0;
   Half  : constant e_Real := +0.5;
   Three_Quarters    : constant e_Real := +(0.75);
   Less_Than_Half_Pi : constant e_Real := +(3.14159265 / 2.0);


   --  Global memory for important constants.  They're initialized
   --  on the first call to the functions that return them, and
   --  there after are simply returned from memory:
   --  Pi, Sqrt_2, etc

   type Pi_mem is record
     Val : e_Real; -- automatically initialized to Zero.
     Initialized : Boolean := False;
   end record;
   Pi_Memory : Pi_Mem;

   type Inverse_Pi_mem is record
     Val : e_Real; -- automatically initialized to Zero.
     Initialized : Boolean := False;
   end record;
   Inverse_Pi_Memory : Inverse_Pi_Mem;

   type Quarter_Pi_mem is record
     Val : e_Real; -- automatically initialized to Zero.
     Initialized : Boolean := False;
   end record;
   Quarter_Pi_Memory : Quarter_Pi_Mem;

   type Inverse_Sqrt_2_mem is record
     Val : e_Real;
     Initialized : Boolean := False;
   end record;
   Inverse_Sqrt_2_memory : Inverse_Sqrt_2_mem;

   type Half_Inverse_Log_2_mem is record
     Val : e_Real;
     Initialized : Boolean := False;
   end record;
   Half_Inverse_Log_2_memory : Half_Inverse_Log_2_mem;

   type Log_2_mem is record
     Val : e_Real;
     Initialized : Boolean := False;
   end record;
   Log_2_memory : Log_2_mem;

   ------------
   -- Arcsin --
   ------------

   --  The result of the Arcsin function is in the quadrant containing
   --  the point (1.0, x), where x is the value of the parameter X. This
   --  quadrant is I or IV; thus, the range of the Arcsin function is
   --  approximately -Pi/2.0 to Pi/2.0 (-Cycle/4.0 to Cycle/4.0, if the
   --  parameter Cycle is specified).
   --  Argument_Error is raised by Arcsin when the absolute
   --  value of the parameter X exceeds one.
   --
   --  Uses Newton's method: Y_k+1 = Y_k + (A - Sin(Y_k)) / Cos(Y_k)
   --  to get Arcsin(A).
   --  Requires call to Arcsin (Real) and assumes that this call gets
   --  the first two radix digits correct: 48 bits usually.  (It almost
   --  always gets 53 bits correct.)
   --
   --  Arcsin(x) = x + x^3/6 + 3*x^5/40 ...
   --     (so Arctan(x) = x if x <  e_Real_Model_Epsilon)

   function Arcsin (X : e_Real)
   return e_Real
   is
      Y_0     : e_Real;
      X_Abs : e_Real := Abs (X);
      No_Correct_Bits   : E_Integer;
      Sign_Is_Negative  : Boolean := False;
      Scaling_Performed : Boolean := False;
   begin

      if Are_Equal (X_Abs, Positive_Infinity) then
         raise E_Argument_Error;
      end if;

      if One < X_Abs then
         raise E_Argument_Error;
      end if;

      if X_Abs < e_Real_Model_Epsilon then
         return X; -- series solution: arcsin = x + x^3/6 + ...
      end if;

      Sign_Is_Negative := False;
      if X < Zero then
         Sign_Is_Negative := True;
      end if;

      if Are_Equal (X_Abs, One) then
         Y_0 := Two_Digit * e_Quarter_Pi;
         if Sign_Is_Negative then Y_0 := -Y_0; end if;
         return Y_0;
      end if;

      -- STEP 2.  We may have to scale the argument if it's near 1.  Newton_Raphson
      -- doesn't do well there.  So we use identity:
      --
      --   arcsin(x) - arcsin(y) = arcsin(x*Sqrt(1-y*y) - y*Sqrt(1-x*x)),
      --
      -- so setting y = 1 (get the arccos(x) = Pi/2 - Arcsin(x)):
      --
      --   arcsin(x)  = -arcsin(Sqrt(1-x*x)) + Pi/2.
      --
      -- Well, we can't use Sqrt(1-x*x) because there's too much error near
      -- X = small.  We can use Sqrt(1-X) * Sqrt(1+X) but don't want the 2 sqrts.
      -- So!  we want Arcsin (Sqrt(1-X) * Sqrt(1+X)) = Y, or
      -- Sin(Y) = 2 * Sqrt((1-X)/2) * Sqrt((1+X)/2). Then Sin(Y) = 2*Sin(Z)Cos(Z)
      -- where Z = Arcsin(Sqrt((1-X)/2)).  So, Y = Arcsin (Sin(Z + Z)) = 2*Z.
      -- The final answer is Arcsin(X) = Pi/2 - 2*Arcsin(Sqrt((1-X)/2)).
      --
      -- IMPORTANT: We can't scale if X_Abs <= 0.5 because we use this
      -- routine with an argument of 0.5 to calculate Pi, and scaling
      -- requires a knowledge of Pi.

      Scaling_Performed := False;
      if Three_Quarters < X_Abs then
         X_Abs := Sqrt (Half - Half_Digit * X_Abs);
         Scaling_Performed := True;
      end if;

      -- Need starting value for Newton's iteration of Arcsin (X_scaled).
      -- Assume positive arg to improve likelihood that
      -- external Arcsin fctn. does what we want.

      Y_0 := Make_Extended (Arcsin (Make_Real (X_Abs)));

      No_Correct_Bits := Real'Machine_Mantissa - 2;

      loop

       --Y_0 := Y_0 + (X_Abs - Sin(Y_0)) / Cos(Y_0);
         Y_0 := Y_0 + Divide ((X_Abs - Sin(Y_0)) , Cos(Y_0));
	 -- Cos is never near 0.

         No_Correct_Bits := No_Correct_Bits * 2;
         exit when No_Correct_Bits >= Max_Available_Bits;

      end loop;

      if Scaling_Performed then
         Y_0 := Two_Digit * (e_Quarter_Pi - Y_0);
      end if;

      if Sign_Is_Negative then
         Y_0 := -Y_0;
      end if;

      return Y_0;

   end Arcsin;

   ------------
   -- Arccos --
   ------------

   --  The result of the Arccos function is in the quadrant containing
   --  the point (x, 1.0), where x is the value of the parameter X. This
   --  quadrant is I or II; thus, the Arccos function ranges from 0.0 to
   --  approximately Pi (Cycle/2.0, if the parameter Cycle is
   --  specified).
   --
   --  Argument_Error is raised by Arccos when the absolute
   --  value of the parameter X exceeds one.
   --
   --  In all cases use Arccos(X) = Pi/2 - Arcsin(X).
   --  When Abs(X) < 0.5 we use  Pi/2 - Arcsin(X).  When
   --  Abs(X) > 0.5 it's better to use the following formula for Arcsin:
   --      Arcsin(X) =  Pi/2 - 2*Arcsin(Sqrt((1-|X|)/2)).   (X > 0.0)
   --      Arcsin(X) = -Pi/2 + 2*Arcsin(Sqrt((1-|X|)/2)).   (X < 0.0)
   --
   function Arccos (X : e_Real)
   return e_Real
   is
      Result : e_Real;
      X_Abs  : constant e_Real := Abs (X);
   begin

      if X_Abs > One then
         raise Constraint_Error;
      end if;

      if Are_Equal (X, Zero) then
         return Two_Digit * e_Quarter_Pi;
      end if;

      if Are_Equal (X, One) then
         return Zero;
      end if;

      if Are_Equal (X, -One) then
         return e_Pi;
      end if;

      if Are_Equal (X, Positive_Infinity) then
         raise E_Argument_Error;
      end if;

      if Are_Equal (X, Negative_Infinity) then
         raise E_Argument_Error;
      end if;

      if X > Half then
         Result :=        Two_Digit * Arcsin (Sqrt (Half - Half_Digit*X_Abs));
      elsif X < -Half then
         Result := Four_Digit*
           (e_Quarter_Pi - Half_Digit*Arcsin (Sqrt (Half - Half_Digit*X_Abs)));
      else
         Result := (e_Quarter_Pi - Arcsin(X)) + e_Quarter_Pi;
      end if;
      --  e_Quarter_Pi is stored with more correct binary digits than
      --  E_Half_Pi because it's slightly less than 1.0.

     return Result;

   end Arccos;

   ------------
   -- Arctan --
   ------------

   --  Newton's method avoids divisions:
   --
   --     Result := Result + Cos(Result) * (Cos(Result) * Arg - Sin(Result))
   --
   --  Arctan(x) =  Pi/2 -  Arctan(1/x)
   --  Arctan(x) = -Pi/2 +  Arctan(1/x)  (x<0)
   --
   --  Arctan(x) = x - x^3/3 + x^5/5 - x^7/7 ...
   --     (so Arctan(x) = x if x <  e_Real_Model_Epsilon)
   --
   --  Arctan (X) =  Arcsin(X / Sqrt(1 + X**2)).
   --  
   --  Not really Ada95-ish for Arctan.
   --  
   function Arctan
     (X : e_Real)
      return e_Real
   is
      Y_0, Arg, Cos_Y_0, Sin_Y_0 : e_Real;
      X_Abs : constant e_Real := Abs (X);
      X_real : Real;
      Sign_Is_Negative : Boolean := False;
      Argument_Reduced : Boolean := False;
      No_Correct_Bits : E_Integer;
   begin
      if X_Abs < e_Real_Model_Epsilon then
         return X; -- series solution: arctan = x - x^3/3 + ...
      end if;

      Sign_Is_Negative := False;
      if X < Zero then
         Sign_Is_Negative := True;
      end if;

      --  returns Pi/2 at +/- inf.  (Note Reciprocal returns 1/inf as Zero.)
      --  inf is regarded as large finite number.  Raising exceptions
      --  in these cases may also be good idea.  Which is right?

      if Are_Equal (X_Abs, Positive_Infinity) then
         Y_0 := Two_Digit * e_Quarter_Pi;
         if Sign_Is_Negative then Y_0 := -Y_0; end if;
         return Y_0;
      end if;

      -- function Make_Real underflows to 0.0 as desired for small X_Abs.

      if X_abs < Two then
         Argument_Reduced := False;
	 Arg := X_abs;
      else
         --  use: Arctan(x) = Pi/2 - Arctan(1/x)
         Argument_Reduced := True;
         Arg := Reciprocal (X_abs);
      end if;

      X_real := Make_Real (Arg);
      Y_0    := Make_Extended (Arctan (X_real));

      No_Correct_Bits := Real'Machine_Mantissa - 2;

      loop

         Sin_Y_0 := Sin (Y_0);
         Cos_Y_0 := Cos (Y_0);                    -- bst accuracy.
       --Cos_Y_0 := Sqrt (One - Sin_Y_0*Sin_Y_0); -- fstr, not too bad accuracy-wise.

         Y_0 := Y_0 + (Arg*Cos_Y_0 - Sin_Y_0) * Cos_Y_0;

         No_Correct_Bits := No_Correct_Bits * 2;
         exit when No_Correct_Bits >= Max_Available_Bits;

      end loop;

      if Argument_Reduced then
         Y_0 := Two_Digit * (e_Quarter_Pi - Half_Digit*Y_0);
      end if;
      --  Lose a whole guard digit of precision when we double e_Quarter_Pi.
      --  (It becomes slightly > 1, pushes a digit off the end of the mantissa.)
      --  So do the subtraction 1st.

      if Sign_Is_Negative then
         Y_0 := -Y_0;
      end if;

      return Y_0;

   end Arctan;

   ---------
   -- Log --
   ---------

   --  Uses Newton's method: Y_k+1 = Y_k + (A - Exp(Y_k)) / Exp(Y_k)
   --  to get Log(A) = Natural Log, the inverse of Exp.
   --  Requires call to Log (Real) and assumes that this call gets
   --  the first two radix digits correct: 48 bits usually.  (It almost
   --  always gets 53 bits correct.)  Argument reduction due to Brent:
   --  1st step is to get a rough approximate value of Log(X), called
   --  Log_X_approx. Then get Y = Log(X/exp(Log_X_approx)).
   --  The final answer is Log(X) = Y + Log_X_approx.  (Disabled at present.)
   --  This gets Y very near 0, which greatly speeds up subsequent
   --  calls to Exp in the Newton iteration above.  Actually,
   --  the scaling business is commented out below.  If you uncomment it,
   --  the routine runs 60% faster, in some tests, but is less accurate.
   --  We'll err on the side of accuracy and use an unscaled X.  But
   --  if you use extended floating point with 2 guard digits, it would
   --  be sensible to use the scaled X, because the loss in accuracy is
   --  negligable compared to the extra precision of another guard digit.
   --  Test for X = One: must return Zero.
   --
   --  The exception Constraint_Error is raised, signaling a pole of the
   --  mathematical function (analogous to dividing by zero), in the following
   --  cases, provided that Float_Type'Machine_Overflows is True:
   --  by the Log, Cot, and Coth functions, when the value of the
   --  parameter X is zero;

   function Log (X : e_Real)
      return e_Real
   is
      Result, X_scaled, Y_0 : e_Real;
      X_Exp : E_Integer;
      No_Correct_Bits   : E_Integer;
      Log_X_approx_real : Real;
   begin

      if X < Zero then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, Zero) then
          raise Constraint_Error;
      end if;

      if Are_Equal (X, Positive_Infinity) then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, Negative_Infinity) then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, One) then
          return Zero;
      end if;

      -- STEP 1.  First argument reduction step.  Get Log_X_approx_real using
      -- a call to log(real).  If X=0 or inf, then X_scaled=0..get exception.

      X_Exp    := Exponent (X);
      X_scaled := Fraction (X);  -- not zero or infinity.

      Log_X_approx_real := Log (Make_Real(X_scaled)) + Log (Radix) * Real (X_Exp);
      X_scaled          := X;

      -- Log_X_approx := Make_Extended (Log_X_approx_real);
      -- X_scaled     := X / Exp (Log_X_approx);
      --  The above Scaling is the fastest, since then Exp does no scaling.
      --  It is clearly less accurate, tho' tolerably so, especially if you
      --  use two guard digits instead of 1.  We use the unscaled X here,
      --  because it (for example) returns a value of Log(2.0)
      --  with an error that is about 50 times smaller than the above.


      -- STEP 2.  Need starting value for Newton's iteration of Log (X_scaled).

    --Y_0 := Make_Extended (Log (Make_Real (X_scaled)));  -- slightly > Zero

      Y_0 := Make_Extended (Log_X_approx_real);


      -- STEP 3.  Start the iteration.  Calculate the number of iterations
      -- required as follows.  num correct digits doubles each iteration.
      -- 1st iteration gives 4 digits, etc.  Each step set desired precision
      -- to one digit more than that we expect from the Iteration.

      No_Correct_Bits := Real'Machine_Mantissa - 2;

      loop

         Y_0 := Y_0 + (X_scaled * Exp(-Y_0) - One);

         No_Correct_Bits := No_Correct_Bits * 2;
         exit when No_Correct_Bits >= Max_Available_Bits;

      end loop;

    --Result := Y_0 + Log_X_approx; -- for use with scaled version (Step 2)
      Result := Y_0;

      return Result;

   end Log;

   --------------------------
   -- Log (arbitrary base) --
   --------------------------

   --  The exception E_Argument_Error is raised, signaling a parameter
   --  value outside the domain of the corresponding mathematical function,
   --  in the following cases:
   --  by the Log function with specified base, when the value of the
   --  parameter Base is zero, one, or negative;
   --
   --  The exception Constraint_Error is raised, signaling a pole of the
   --  mathematical function (analogous to dividing by zero), in the following
   --  cases, provided that Float_Type'Machine_Overflows is True:
   --  by the Log, Cot, and Coth functions, when the value of the
   --  parameter X is zero.
   --
   --  Struggling to remember: Z = Log(X, Base) implies X = Base**Z
   --  or X = Exp (Log(Base)*Z) which implies Log(X) = Log(Base) * Z, so
   --  Log(X, Base) = Z = Log(X) / Log(Base).
   --
   function Log (X : e_Real; Base : e_Real) return e_Real is
      Result : e_Real;
   begin

     if Are_Equal (Base,Zero) then
        raise E_Argument_Error;
     end if;

     if Base < Zero then
        raise E_Argument_Error;
     end if;

     if Are_Equal (Base,One) then
        raise E_Argument_Error;
     end if;

     if Are_Equal (Base,Two) then

        Result := Two_Digit * (Log(X) * e_Half_Inverse_Log_2);
        --  Divide by e_log_2. Multiply by 0.5/Log(2) not 1/log(2), because
        --  0.5/Log(2) is slightly less than 1, hence contains more correct
        --  digits.  (And multiplication is preferred for efficiency.)

     else

        Result := Log(X) / Log(Base);

     end if;

     return Result;

   end Log;

   ----------
   -- "**' --
   ----------

   --  Say X**N = Exp (Log (X) * N).
   --
   --  Exponentiation by a zero exponent yields the value one.
   --  Exponentiation by a unit exponent yields the value of the left
   --  operand.  Exponentiation of the value one yields the value one.
   --  Exponentiation of the value zero yields the value zero.
   --  The results of the Sqrt and Arccosh functions and that of the
   --  exponentiation operator are nonnegative.
   --
   --  Argument_Error is raised by "**" operator, when the value of the left
   --  operand is negative or when both operands have the value zero;
   --
   function   "**" (Left : e_Real; Right : e_Real) return e_Real is
      Result : e_Real;
   begin

     -- Errors:

     if Left < Zero then
        raise E_Argument_Error;
     end if;

     if Are_Equal (Left, Zero) and then Are_Equal (Right, Zero) then
        raise E_Argument_Error;
     end if;

     -- Special Cases.  We now know that they aren't both Zero:

     if Are_Equal (Right, Zero) then  -- Left is not Zero
        return One;
     end if;

     if Are_Equal (Left, Zero) then   -- Right is not Zero
        return Zero;
     end if;

     if Are_Equal (Right, One) then
        return Left;
     end if;

     if Are_Equal (Left, One) then    -- Still OK if Right = Zero
        return One;
     end if;

     -- Should we optimize for integer N?

     Result := Exp (Log (Left) * Right);

     return Result;

   end "**";

   ---------
   -- Exp --
   ---------

   --  Sum Taylor series for Exp(X).
   --  Actually, we sum series for Exp(X) - 1 - X, because scaling makes
   --  X small, and Exp(X) - 1 - X has more correct digits for small X.
   --  [get max arg size and test for it.]
   --
   function Exp
     (X : e_Real)
      return e_Real
   is
      Order          : Real      := 0.0;
      Delta_Exponent : E_Integer := 0;
      Next_Term, Sum : e_Real;

      X_Scaled_2, X_scaled_1 : e_Real;
      Total_Digits_To_Use : E_Integer;

      N        : Integer;
      N_e_Real : e_Real;

      J              : constant Integer := 11;
      Scale_Factor_2 : constant Integer := 2**J;
      --  Must be power of 2 for function call to Make_E_Digit.
      --  The optimal value increases with desired precision.  But higher
      --  order terms in the series are cheap, so it's not *too* important.

      Inverse_Two_To_The_J : constant E_Digit
                               := Make_E_Digit (1.0 / Real (Scale_Factor_2));

      We_Flipped_The_Sign_Of_X_scaled : Boolean := False;
      First_Stage_Scaling_Performed   : Boolean := False;
      Second_Stage_Scaling_Performed  : Boolean := False;
   begin

      if Are_Equal (X, Zero) then
          return One;
      end if;

      if Are_Equal (X, Positive_Infinity) then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, Negative_Infinity) then
          raise E_Argument_Error;
      end if;

      -- STEP 1.  Reduce argument in 2 stages: 1st Remainder(X, Log(2)),
      -- then divide by 2**J.
      -- So X_Scaled = Remainder (X, Log_2), or approximately:
      -- X_Scaled = X - Unbiased_Rounding (X / Log(2)) * Log(2) = X - N * Log(2)
      -- Then Exp(X) = Exp(X_Scaled) * Exp(N*Log(2)) = Exp(X_Scaled) * 2**N.
      --
      -- Second stage of argument reduction: divide X_Scaled by 2**J:
      -- Exp(X) = Exp(X_Scaled/2**J)**2**J * 2**N.
      --
      -- E_log_2 is calculated by recursively calling this routine, but
      -- with an argument very near 0.0, so First stage scaling is not
      -- performed.  (It also calls with arg of approx. log(2) = 0.69, so
      -- must not allow 1st stage scaling for args that small.)

      N                             := 0;
      X_Scaled_1                    := X;
      First_Stage_Scaling_Performed := False;

      if Three_Quarters < Abs (X) then

        N_e_Real := Unbiased_Rounding (Two_Digit * (X_Scaled_1 * E_Half_Inverse_Log_2));
        X_Scaled_1 := Remainder (X_Scaled_1, E_Log_2);

        -- Use X_Scaled_1 := X_Scaled_1 - N_e_Real * E_Log_2; ?
        -- Not much faster.  Somewhat less accurate.  Seems OK for small args.

        if Make_Extended (Real (Integer'Last)) < Abs (N_e_Real) then
           raise E_Argument_Error with "Argument too large in Exp.";
        end if;

        N := Integer (Make_Real (N_e_Real));

        First_Stage_Scaling_Performed := True;

      end if;

      -- STEP 1b. We want to get Exp() slightly less than one, to maximize
      -- the precision of the calculation. So make sure arg is negative.

      We_Flipped_The_Sign_Of_X_scaled := False;
      if not (X_scaled_1 < Zero) then
         We_Flipped_The_Sign_Of_X_scaled := True;
         X_scaled_1 := -X_scaled_1;
      end if;


      -- STEP 2.  2nd stage of argument reduction.  Divide X_scaled by 2**J.
      -- Don't scale if arg is already small to avoid complications due
      -- to underflow of arg to zero.  Arg may already be 0.  It's OK.

      if Exponent (X_Scaled_1) >= -2 then             -- it'll do, Zero OK here
         X_Scaled_2 := Inverse_Two_To_The_J * X_Scaled_1;
         Second_Stage_Scaling_Performed := True;
      else
         X_scaled_2 := X_scaled_1;
         Second_Stage_Scaling_Performed := False;
      end if;


      -- STEP 3.  Start the sum.  Calculate Exp(X_Scaled) - (1 + X_Scaled),
      -- because for small X_Scaled, this contains more correct digits than Exp.
      -- Start summing the series at order 2 instead of order 0.
      -- We have verified above that X_scaled /= Zero.

      Order     := 2.0;
      Next_Term := Half_Digit * X_Scaled_2 * X_Scaled_2;
      Sum       := Next_Term;

      loop

         Order := Order + 1.0;
         if Order = Radix then
            raise E_Argument_Error with "Too many terms needed in Exp taylor sum.";
         end if;

         --  Use relative Exponents of Sum and Next_Term to check convergence
         --  of the sum.  Exponent doesn't work for args of 0, so check.
         --  Abs (Next_Term) <= Abs (Sum), so we need only check Next_Term.

         if Are_Equal (Next_Term, Zero) then exit; end if;

         Delta_Exponent := Exponent (Sum) - Exponent (Next_Term);
         Total_Digits_To_Use := e_Real_Machine_Mantissa - Delta_Exponent + 1;
         exit when Total_Digits_To_Use <= 0;

         Next_Term := (X_Scaled_2 * Next_Term) / Make_E_Digit (Order);

         Sum := Sum + Next_Term;
         --  Sum can overflow to infinity?  Not with our scaled arguments.

      end loop;

      -- STEP 4.  Undo effect of 2nd stage of argument scaling.  Recall we
      -- divided the arg by 2**J, and found Exp(X_Scaled/2**J).  Now to get
      -- Exp(X_Scaled), must take Exp(X_Scaled/2**J)**2**J, which means
      -- repeated squaring of Exp(X_Scaled/2**J) (J times).  It's more complicated
      -- than that because we calculated G(X) = Exp(X) - 1 - X (since G contains
      -- more correct digits than Exp, expecially for small X.)   So we
      -- use G(2X) = Exp(2X) - 1 - 2X = (G + (1 + X))*(G + (1 + X)) - 1 - 2X
      --           = G*G + 2*G*(1+X) + X*X
      --     G(2X) = (G(X)+X)**2 + 2G(X).
      --     G(4X) = (G(2X)+2X)**2 + 2G(2X).
      -- Repeat J times to unscale G.  The following also returns X_scaled*2**J.

      if Second_Stage_Scaling_Performed then
         for I in 1..J loop
            Sum        := (Sum + X_scaled_2)**2 + Two_Digit * Sum;
            X_Scaled_2 := Two_Digit * X_Scaled_2;
         end loop;
      end if;

      -- DO the following whether or not Second_Stage or First_Stage
      -- scaling was performed (because the series sum neglected the
      -- the 1 and the X.  If there were no  extra guard digit for
      -- subtraction ((X_scaled_1 + Sum) is negative) then it would be best
      -- to use (0.5 + (Sum + Xscaled)) + 0.5.  Following is OK though to
      -- get a number slightly less than one with full precision.
      -- Recover Exp = G(X) + 1 + X = Sum + 1 + X = (Sum + X_scaled) + 1:

      Sum := (Sum + X_scaled_1) + One;

      -- Second stage unscaling.  We now have Y = Exp(-|X_scaled_1|), which is
      -- slightly less than 1.0.  Keep
      -- in Y < 1.0 form as we unscale: might preserve more precision that
      -- way, cause we lose much precision if invert a number that's slightly
      -- less than one.

      if First_Stage_Scaling_Performed then
         if We_Flipped_The_Sign_Of_X_scaled then
            Sum := Sum * Two**(-N);
         else
            Sum := Sum * Two**(N);
         end if;
      end if;

      if We_Flipped_The_Sign_Of_X_scaled then
         Sum := Reciprocal (Sum);
         --  X_scaled was positive. We flipped its sign so must invert the result.
      end if;

      return Sum;

   end Exp;

   ---------------------------
   -- Sin (arbitrary cycle) --
   ---------------------------

   --  The exception E_Argument_Error is raised, signaling a parameter
   --  value outside the domain of the corresponding mathematical function,
   --  in the following cases:
   --  by any forward or inverse trigonometric function with specified
   --  cycle, when the value of the parameter Cycle is zero or negative;
   --  The results of the Sin, Cos, Tan, and Cot functions with
   --  specified cycle are exact when the mathematical result is zero;
   --  those of the first two are also exact when the mathematical
   --  result is +/-1.0.

   function Sin 
     (X     : e_Real;
      Cycle : e_Real) 
      return e_Real 
   is
      Fraction_Of_Cycle, Result : e_Real;
   begin
     -- The input parameter X is units of Cycle.  For example X = Cycle
     -- is same as X = 2Pi, so could call Sin (2 * Pi * X / Cycle), but we
     -- want to apply the remainder function here to directly meet certain
     -- requirements on returning exact results.
     -- Recall:   Remainder = X - Round (X/Cycle) * Cycle = X - N * Cycle
     -- which is in the range  -Cycle/2 .. Cycle/2.  The formula will be
     --
     --     Sin (X, Cycle) = Sin (2 * Pi * X / Cycle)
     --                    = Sin (2 * Pi * (X - N * Cycle) / Cycle)
     --                    = Sin (2 * Pi * Remainder(X,Cycle) / Cycle)

     if Are_Equal (Cycle, Zero) then
        raise E_Argument_Error;
     end if;

     if Cycle < Zero then
        raise E_Argument_Error;
     end if;

     Fraction_Of_Cycle := Remainder (X, Cycle) / Cycle;

     if Are_Equal (Fraction_Of_Cycle, Zero) then
        return Zero;
     end if;

     if Are_Equal (Abs (Fraction_Of_Cycle), Half) then
        return Zero;
     end if;

     if Are_Equal (Fraction_Of_Cycle, Make_Extended(0.25)) then
        return One;
     end if;

     if Are_Equal (Fraction_Of_Cycle, Make_Extended(-0.25)) then
        return -One;
     end if;

     Result := Sin (Make_E_Digit(8.0) * (e_Quarter_Pi * Fraction_Of_Cycle));
     --  Pi/4 is used instead of Pi/2, because it contains more correct
     --  binary digits.

     return Result;

   end Sin;

   ---------
   -- Sin --
   ---------

   --  Sum Taylor series for Sin (X).
   --
   --  Max argument is at present set by requirement:
   --      Exponent(X) < Present_Precision-1

   function Sin
     (X : e_Real)
      return e_Real
   is
      Half_Sqrt_Of_Radix : constant Real := 2.0**(Desired_No_Of_Bits_in_Radix/2-1);
      Order          : Real      := 0.0;
      Delta_Exponent : E_Integer := 0;
      Next_Term, Sum : e_Real;

      X_Scaled_1, X_Scaled_2, X_Scaled_2_Squared : e_Real;

      Total_Digits_To_Use : E_Integer;

      N_e_Real, Half_N : e_Real;

      J : constant Integer := 8;
      Three_To_The_J : constant E_Digit := Make_E_Digit (3.0**J);

      Factorial_Part : E_Digit;

      Sign_Of_Term_Is_Pos : Boolean := True;
      Arg_Is_Negative     : Boolean := False;
      N_Is_Odd            : Boolean := False;

      First_Stage_Scaling_Performed  : Boolean := False;
      Second_Stage_Scaling_Performed : Boolean := False;
   begin
      if Are_Equal (X, Zero) then
         return Zero;
      end if;

      if Are_Equal (X, Positive_Infinity) then
        raise E_Argument_Error;
      end if;

      if Are_Equal (X, Negative_Infinity) then
         raise E_Argument_Error;
      end if;

      Arg_Is_Negative := False;
      if X < Zero then
         Arg_Is_Negative := True;
      end if;

      if Exponent(X) >= e_Real_Machine_Mantissa-1 then
         raise E_Argument_Error;
      end if;
      --  Can't figure out N below. N_e_Real has to be
      --  integer valued: 0, 1, ... 2**(Radix*e_Real_Machine_Mantissa) - 1
      --  This determines Max allowed Argument.
      --  If Exponent(N_e_Real) is too large, then can't tell if N is even or odd.


      -- STEP 1.  First argument reduction: Get modulo Pi by Remainder(X, Pi).
      -- Get X in range -Pi/2..Pi/2.
      --     X_scaled_1 := X - N_e_Real * e_Pi;

      if Less_Than_Half_Pi < Abs(X) then
         X_scaled_1   := Remainder (Abs(X), e_Pi);
         N_e_Real     := Unbiased_Rounding ((Abs(X)-X_scaled_1) * e_Inverse_Pi);
         First_Stage_Scaling_Performed := True;
      else
         X_Scaled_1            := Abs(X);
         N_e_Real              := Zero;
         First_Stage_Scaling_Performed := False;
      end if;

      -- Need to know if N is even or odd.  N is Positive.
      -- If Exponent(N_e_Real) is too large, then we can't tell if N is even or
      -- odd.  So raise Arg error.  This determines Max Arg.

      N_Is_Odd := False;
      if not Are_Equal (N_e_Real, Zero) then
         Half_N := Half_Digit * N_e_Real;
         if Truncation (Half_N) < Half_N then
            N_Is_Odd := True;
         end if;
      end if;

      -- STEP 2.  Second stage of argument reduction.  Divide by 3**5 = 243
      -- to get Arg less than about 0.01?.  Call this 3**J in what follows.
      -- Later we recursively use J repetitions of the formula
      -- Sin(3*Theta) =  Sin(Theta)*(3 - 4*Sin(Theta)**2), to get Sin(Theta*3**J)
      -- Cos(3*Theta) = -Cos(Theta)*(3 - 4*Cos(Theta)**2), to get Cos(Theta*3**J)
      -- to get Sin (Original_Arg).
      --
      --  MUST avoid underflow to Zero in this step.  So only if X_Scaled is big.
      --  Actually, X_scaled = 0 may pass through, but it's OK to start out 0.

      if Exponent (X_Scaled_1) >= -2 then             -- it'll do
         X_Scaled_2 := X_scaled_1 / Three_To_The_J;
         Second_Stage_Scaling_Performed := True;
      else
         X_Scaled_2 := X_scaled_1;
         Second_Stage_Scaling_Performed := False;
      end if;

      -- STEP 3.  Start the sum.  Terms are labeled Order = 1, 2, 3
      -- but the series is X - X**3/3! + X**5/5! + ...+- X**(2*Order-1)/(2*Order-1)!.
      -- Summed G(X) = Sin(X) - X, which contains more correct digits at
      -- the end.  We need these extra digits when we unscale the result.

      X_Scaled_2_Squared  := X_scaled_2 * X_Scaled_2;

      Order               := 2.0;
      Sign_Of_Term_Is_Pos := False;
      Next_Term           := X_Scaled_2 * X_Scaled_2_Squared / Make_E_Digit (6.0);
      Sum                 := -Next_Term;
      --  Above we make the 1st term in G, and begin the sum.

      loop

         Sign_Of_Term_Is_Pos := not Sign_Of_Term_Is_Pos;
         Order               := Order + 1.0;

         --  Can't make Factorial part if, roughly, 2*Order-1.0 > Radix-1,
         --  Because max argument of Make_E_Digit is Radix-1.

         if Order >= Radix / 2.0 then
            raise E_Argument_Error with "Too many terms needed in Exp taylor sum.";
         end if;

         --  Use relative Eponents of Sum and Next_Term to check convergence
         --  of the sum.  Exponent doesn't work for args of 0, so check.
         --  Abs (Next_Term) <= Abs (Sum), so we need only check Next_Term.

         if Are_Equal (Next_Term, Zero) then exit; end if;

         Delta_Exponent := Exponent (Sum) - Exponent (Next_Term);

         Total_Digits_To_Use := e_Real_Machine_Mantissa - Delta_Exponent + 1;

         exit when Total_Digits_To_Use <= 0;

         if Order < Half_Sqrt_Of_Radix then

            Factorial_Part := Make_E_Digit ((2.0*Order-1.0)*(2.0*Order-2.0));
            Next_Term := (X_Scaled_2_Squared * Next_Term) / Factorial_Part;

         else

            Factorial_Part := Make_E_Digit ((2.0*Order-1.0));
            Next_Term      := (X_Scaled_2_Squared * Next_Term) / Factorial_Part;
            Factorial_Part := Make_E_Digit ((2.0*Order-2.0));
            Next_Term      := Next_Term / Factorial_Part;

         end if;

         if Sign_Of_Term_Is_Pos then
            Sum := Sum + Next_Term;
         else
            Sum := Sum - Next_Term;
         end if;

      end loop;

      -- STEP 4.  Scale the result iteratively. Recall we divided the arg by 3**J,
      -- so we recursively use J repetitions of the formula
      -- Sin(3*X) = Sin(X)*(3 - 4*Sin(X)**2), to get Sin(X*3**J).
      -- Actually, we summed G(X) = Sin(X) - X.  So the formula for G(X) is
      -- G(3X) = S(3X) - 3X = S(X)*(3 - 4S(X)**2) - 3X,
      --                    = (G+X)*(3 - 4(G+X)**2) - 3X,
      --                    = 3G     - 4(G+X)**3, (Cancel out the 3X),
      --              G(3X) = 3G(X)  - 4(G(X)+X)**3.
      --              G(9X) = 3G(3X) - 4(G(3X)+3X)**3, etc.
      -- Still requires only 2 (full) mults per loop, just like the original formula.
      -- Notice below that we output X_scaled * 3**J, which is required next step.

      if Second_Stage_Scaling_Performed then
         for I in 1..J loop
            Sum        := Three_Digit * Sum - Four_Digit * (Sum + X_scaled_2)**3;
            X_scaled_2 := Three_Digit * X_scaled_2;
         end loop;
      end if;


      -- STEP 5.  We have Sin(X - N * Pi).  Want Sin(X).  If N is odd, then
      -- flip sign of Sum.  Next, flip sign again if the
      -- original argument is neg: Arg_Is_Negative = True.
      -- Remember, we summed for Sum = G = Sin - X, whether or not scaling
      -- was performed.  (X is called X_scaled, no matter what.)
      -- So we recover Sin = G + X

      Sum := Sum + X_scaled_1;

      if First_Stage_Scaling_Performed then
         if N_Is_Odd then
            Sum := -Sum;
         end if;
      end if;

      if Arg_Is_Negative then
         Sum := -Sum;
      end if;

      return Sum;

   end Sin;

   ---------------------------
   -- Cos (arbitrary cycle) --
   ---------------------------

   --  The exception E_Argument_Error is raised, signaling a parameter
   --  value outside the domain of the corresponding mathematical function,
   --  in the following cases:
   --  by any forward or inverse trigonometric function with specified
   --  cycle, when the value of the parameter Cycle is zero or negative;
   --  The results of the Sin, Cos, Tan, and Cot functions with
   --  specified cycle are exact when the mathematical result is zero;
   --  those of the first two are also exact when the mathematical
   --  result is +/-1.0.
   --
   function Cos (X : e_Real; Cycle : e_Real) return e_Real is
      Fraction_Of_Cycle, Result : e_Real;
   begin

     -- The input parameter X is units of Cycle.  For example X = Cycle
     -- is same as X = 2Pi, so could use Cos (2 * Pi * X / Cycle), but we
     -- want to apply the remainder function here to directly meet certain
     -- requirements on returning exact results.
     -- Recall:   Remainder = X - Round (X/Cycle) * Cycle = X - N * Cycle
     -- which is in the range  -Cycle/2 .. Cycle/2.  The formula will be
     --
     --     Cos (X, Cycle) = Cos (2 * Pi * X / Cycle)
     --                    = Cos (2 * Pi * (X - N * Cycle) / Cycle)
     --                    = Cos (2 * Pi * Remainder(X,Cycle) / Cycle)

     if Are_Equal (Cycle, Zero) then
        raise E_Argument_Error;
     end if;

     if Cycle < Zero then
        raise E_Argument_Error;
     end if;

     --  Now get twice the fraction of the cycle, and handle special cases:

     Fraction_Of_Cycle := Remainder (X, Cycle) / Cycle;

     if Are_Equal (Fraction_Of_Cycle, Zero) then
        return One;
     end if;

     if Are_Equal (Abs (Fraction_Of_Cycle), Half) then
        return -One;
     end if;

     if Are_Equal (Abs (Fraction_Of_Cycle), Make_Extended(0.25)) then
        return Zero;
     end if;

     Result := Cos (Make_E_Digit(8.0) * (e_Quarter_Pi * Fraction_Of_Cycle));
     --  Use Pi/4 becase it contains more correct binary digits than Pi.

     return Result;

   end Cos;

   ---------
   -- Cos --
   ---------

   --  Sum Taylor series for Cos (X).  Actually sum series for G = Cos(X) - 1.
   --  Reason is, G contains more correct digits than Cos, which is
   --  required when we undo effects of argument reduction.
   --  Max argument is at present
   --  set by requirement: Exponent(X) < Present_Precision-1
   --
   function Cos
     (X : e_Real)
      return e_Real
   is
      Half_Sqrt_Of_Radix : constant Real := 2.0**(Desired_No_Of_Bits_in_Radix/2-1);
      Order              : Real          := 0.0;
      Delta_Exponent : E_Integer := 0;
      Next_Term, Sum, Sum_2, Sum_3 : e_Real;

      X_Scaled_1, X_Scaled_Squared : e_Real;

      Total_Digits_To_Use : E_Integer;

      N_e_Real, Half_N : e_Real;

      J : constant Integer := 8;
      Three_To_The_J : constant E_Digit := Make_E_Digit (3.0**J);

      Factorial_Part : E_Digit;

      Sign_Of_Term_Is_Pos : Boolean := True;
      N_Is_Odd            : Boolean := False;

      First_Stage_Scaling_Performed  : Boolean := False;
      Second_Stage_Scaling_Performed : Boolean := False;
   begin
      if Are_Equal (X, Zero) then
          return One;
      end if;

      if Are_Equal (X, Positive_Infinity) then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, Negative_Infinity) then
          raise E_Argument_Error;
      end if;

      if Exponent(X) >= e_Real_Machine_Mantissa-1 then
         raise E_Argument_Error;
      end if;
      --  Can't figure out N below. N_e_Real has to be
      --  integer valued: 0, 1, ... 2**(Radix*e_Real_Machine_Mantissa) - 1
      --  This determines Max allowed Argument.


      -- STEP 1.  First stage argument reduction.
      -- Take X modulo Pi by Remainder(X, Pi).  Get X in range -Pi/2..Pi/2.
      --     X_scaled_1 := X - N_e_Real * e_Pi;

      if Less_Than_Half_Pi < Abs(X) then
         X_scaled_1   := Remainder (Abs(X), e_Pi);
         N_e_Real     := Unbiased_Rounding ((Abs(X)-X_scaled_1) * e_Inverse_Pi);
         First_Stage_Scaling_Performed := True;
      else
         X_Scaled_1                    := Abs(X);
         N_e_Real                      := Zero;
         First_Stage_Scaling_Performed := False;
      end if;

      -- Need to know if N is even or odd.  N is Positive.

      N_Is_Odd := False;
      if not Are_Equal (N_e_Real, Zero) then
         Half_N := Half_Digit * N_e_Real;
         if Truncation (Half_N) < Half_N then
            N_Is_Odd := True;
         end if;
      end if;

      -- STEP 1b. If X_Scaled is nearing Pi/2 then it's too big for Taylor's.
      -- Must call Sin (Pi/2 - X_Scaled).  IMPORTANT: only do this for
      -- X_scaled > Pi/6, because Arcsin(X => 0.5) = Pi/6 is used to calculate
      -- Pi, and would get infinite recursive calls when Arcsin calls Cos
      -- which calls Sin, while Sin and Cos call e_Quarter_Pi, which
      -- calls Arcsin.  e_Quarter_Pi is used instead of the less accurate Pi/2.

      if One < X_Scaled_1 then

         Sum := Sin ((e_Quarter_Pi - X_Scaled_1) + e_Quarter_Pi);
         if N_Is_Odd then
            Sum := -Sum;
         end if;
         return Sum;

      end if;

      -- STEP 2.  Second stage of argument reduction.  Divide by 3**8 = 81*81
      -- to get argument less than about .02.  Call this 3**J in what follows.
      -- Later we recursively use J repetitions of the formula
      -- Sin(3*Theta) =  Sin(Theta)*(3 - 4*Sin(Theta)**2), to get Sin(Theta*3**J)
      -- Cos(3*Theta) = -Cos(Theta)*(3 - 4*Cos(Theta)**2), to get Cos(Theta*3**J).
      --
      -- It's OK if X_scaled is 0 at this point, but not if the following
      -- forces it to underflow to 0.  Therefore only scale large args:

      if Exponent (X_Scaled_1) >= -2 then             -- it'll do
          X_Scaled_1 := X_scaled_1 / Three_To_The_J;
          Second_Stage_Scaling_Performed := True;
      else
          Second_Stage_Scaling_Performed := False;
      end if;


      -- STEP 3.  Start the sum.  Terms are labeled Order = 0, 1, 2, 3
      -- but the series is 1 - X**2/2! + X**4/4! + ...+- X**(2*Order)/(2*Order)!.
      -- Below we actually calculate Cos(X) - 1.
      -- Start summing the series at order 1 instead of order 0.

      Order               := 1.0;
      X_Scaled_Squared    := X_scaled_1 * X_Scaled_1;
      Next_Term           := Half_Digit * X_Scaled_Squared;
      Sum                 := -Next_Term;
      Sign_Of_Term_Is_Pos := False;

      loop

         Sign_Of_Term_Is_Pos := not Sign_Of_Term_Is_Pos;
         Order               := Order + 1.0;

         -- Can't make Factorial part if, roughly, 2*Order > Radix-1.

         if Order >= (Radix-1.0) / 2.0 then
            raise E_Argument_Error with "Too many terms needed in Exp taylor sum.";
         end if;

         --  Use relative Eponents of Sum and Next_Term to check convergence
         --  of the sum.  Exponent doesn't work for args of 0, so check.
         --  Abs (Next_Term) <= Abs (Sum), so we need only check Next_Term.
         --  If Next_Term is 0, we are finished anyway.

         if Are_Equal (Next_Term, Zero) then exit; end if;

         Delta_Exponent := Exponent (Sum) - Exponent (Next_Term);
         Total_Digits_To_Use := e_Real_Machine_Mantissa - Delta_Exponent + 1;
         exit when Total_Digits_To_Use <= 0;

         if Order < Half_Sqrt_Of_Radix then

           Factorial_Part := Make_E_Digit ((2.0*Order)*(2.0*Order-1.0));
           Next_Term      := (X_Scaled_Squared * Next_Term) / Factorial_Part;

         else    -- Do it the slow way. (Should rarely happen.)

           Factorial_Part := Make_E_Digit (2.0*Order);
           Next_Term      := (X_Scaled_Squared * Next_Term) / Factorial_Part;
           Factorial_Part := Make_E_Digit (2.0*Order-1.0);
           Next_Term      := Next_Term / Factorial_Part;

         end if;

         if Sign_Of_Term_Is_Pos then
            Sum := Sum + Next_Term;
         else
            Sum := Sum - Next_Term;
         end if;

      end loop;

      -- STEP 4.  Scale the result iteratively. Recall we got Cos(Arg/3**J).  Now
      -- we want Cos(Arg).  So we use J repetitions of the formula
      -- Cos(3*Theta) = -Cos(Theta)*(3 - 4*Cos(Theta)**2), to get Cos(Theta*3**J).
      -- Recall we summed for Cos(X) - 1, because we retain more correct digits
      -- this way for small X.  (The 1 would have shifted correct digits off the
      -- array.) So we actually have is Sum = G(X) = Cos(X) - 1.  So the formula
      -- for Cos(3X) is (1+G)*(4(G+1)**2 - 3) = (1+G)*(1 + 8G + 4G**2).  Then
      -- G(3X) = Cos(3X) - 1 = 9G + 12G*2 + 4G**3.  Next, unscale G:

      if Second_Stage_Scaling_Performed then
         for I in 1..J loop
          --Sum := Sum * (Four_Digit * Sum * Sum - Three);
            Sum_2 := Sum*Sum;
            Sum_3 := Sum*Sum_2;
            Sum := Nine_Digit * Sum + Twelve_Digit * Sum_2 + Four_Digit * Sum_3;
         end loop;
      end if;

      -- STEP 5.  We have Cos(X - N * Pi).  Want Cos(X).  If N is odd, then
      -- flip sign of Sum.  First remember we summed for G = Cos - 1, whether or
      -- not scaling was performed.  Must recover Cos next:

      Sum := Sum + One; -- Get Cos(X) = G(X) + 1 = Sum + 1.

      if First_Stage_Scaling_Performed then
         if N_Is_Odd then
            Sum := -Sum;
         end if;
      end if;

      return Sum;

   end Cos;

   -------------------------
   -- Reciprocal_Nth_Root --
   -------------------------

   --  Uses Newton's method to get A**(-1/N):
   --      Y_k+1 = Y_k + (1 - A * Y_k**N) * Y_k / N.
   --  Requires call to log(A) and assumes that this call gets
   --  the first two radix digits correct: 48 bits usually.  (It almost
   --  always gets more correct.)
   --  REMEMBER, this calculates to the max possible precision;
   --  Does not reflect dynamic precision floating point.
   --  N must be less than Radix - 1, which is  usually 2**24 - 1.

   function Reciprocal_Nth_Root 
     (X : e_Real; 
      N : Positive) 
   return e_Real is
      Exponent_Of_X, Scaled_Exp_Of_X, Exp_Mod_N : E_Integer;
      Result, X_Fraction, Scaled_X_Fraction, Y_0 : e_Real;

      Shift_Sign, Ratio : Real := 0.0;
      Log_Of_Scaled_X_Fraction, Real_X_Fraction : Real := 0.0;

      No_Correct_Bits : E_Integer;
      E_Digit_N, Inverse_E_Digit_N : E_Digit; -- already initialized

      Optimization_Is_Possible : Boolean := False;
   begin
      if Are_Equal (X, Zero) then
         raise E_Argument_Error;
      end if;

      if Are_Equal (X, Positive_Infinity) then
         raise E_Argument_Error;
      end if;

      if Are_Equal (X, Negative_Infinity) then
         raise E_Argument_Error;
      end if;

      if Real(N) > (Radix-1.0) then
         raise E_Argument_Error;
      end if;

      -- STEP 0b.  An optimization.  If N is a power of two, we can speed
      -- up the calculation by multiplying by 1/N rather than dividing by
      -- N in the newton iteration.

      Optimization_Is_Possible := False;
      for I in 0..Desired_No_Of_Bits_In_Radix-2 loop
        if 2**I = N then
          Optimization_Is_Possible := True;
        end if;
      end loop;

      if Optimization_Is_Possible then
         Inverse_E_Digit_N := Make_E_Digit (1.0 / Real(N));
      else
         E_Digit_N := Make_E_Digit (Real(N));
      end if;

      -- STEP 1.  Argument reduction.  Break X into Fraction and Exponent.
      -- Choose to decrement Abs(Exponent) by (Exponent MOD N),
      -- and multiply fraction by Radix ** (Exponent MOD N).
      -- The reason is of course, we want to divide decremented Exponent by N,
      -- so we want Scaled_Exp_Of_X to be an integral multiple of -N.

      Exponent_Of_X := Exponent (X);  -- What if X is 0, inf, etc???
      X_Fraction    := Fraction (X);

      Exp_Mod_N :=   Abs (Exponent_Of_X) MOD E_Integer(N); -- N never < 1.

      -- Make sure that Scaled_Exp_Of_X is in range of e_Real, and also
      -- make sure that it is scaled to an integral multiple of N.
      if Exponent_Of_X < 0 then
         Scaled_Exp_Of_X := Exponent_Of_X + Exp_Mod_N;
         Shift_Sign := +1.0;
      else
         Scaled_Exp_Of_X := Exponent_Of_X - Exp_Mod_N;
         Shift_Sign := -1.0;
      end if;

      --  Scale the fraction to compensate for the above shift in the Exponent:
      if Exponent_Of_X < 0 then
         Scaled_X_Fraction := Scaling (X_Fraction, - Exp_Mod_N);
      else
         Scaled_X_Fraction := Scaling (X_Fraction, + Exp_Mod_N);
      end if;

      -- STEP 2.   Get starting value for Newton's iteration.
      -- Want Real number estimate of Scaled_X_Fraction**(-1/N).
      -- Get the first 2 radix digits correct (48 bits usually), and call it Y_0.
      -- Must worry about exponents too large for type Real in the value
      -- Scaled_X_Fraction prior to the **(-1/N) operation, so we do it indirectly.
      -- Arg ** (-1/N): use exp (log (Arg**(-1/N))) = exp ( (-1/N)*log (Arg) ).
      -- Need estimate: [ X_Fraction * Radix**(-Shift_Sign * Exp_Mod_N) ]**(-1/N).
      -- First want natural Log(X_Fraction * Radix**(-Shift_Sign * Exp_Mod_N))
      -- which equals Log(X_Fraction) - Shift_Sign * Log(Radix) * Exp_Mod_N.
      -- Next divide this quantity by -N, and take exp():
      -- Exp (-Log (X_Fraction) / N + Shift_Sign * Log(Radix) * Exp_Mod_N / N).
      -- This is the estimate we want, and because Exp_Mod_N / N is always
      -- < 1.0, the arguments should be well within range of Log and Exp,
      -- because Exp (Shift_Sign * Log(Radix) * Exp_Mod_N / N) is less than Radix.

      Real_X_Fraction := Make_Real (X_Fraction);
      Ratio           := Real (Exp_Mod_N) / Real(N);
      Log_Of_Scaled_X_Fraction
          := -Log (Real_X_Fraction) / Real(N) + Shift_Sign * Ratio * Log (Radix);
      Y_0 := Make_Extended (Exp (Log_Of_Scaled_X_Fraction));
      --  Starting val in Newton's iteration.


      -- STEP 3.  Start the iteration.  Calculate the number of iterations
      -- required as follows.  num correct digits doubles each iteration.
      -- 1st iteration gives 4 digits, etc.  Each step set desired precision
      -- to one digit more than that we expect from the Iteration.
      --
      -- It is important to remember that e_Real_Machine_Mantissa includes
      -- the 1 or 2 guard digits.  The last of these may have a lot of error in
      -- the end, the first of these may have some error.  That's why they're
      -- there.   Also remember that when Set_No_Of_Digits_To_Use is
      -- called, the precision it sets includes the 2 guard digits, both
      -- of which may be wrong, so we add 2 to the setting below, just
      -- in case.

      No_Correct_Bits := Real'Machine_Mantissa - 2;

      loop

         if Optimization_Is_Possible then  --  multiply by inverse of N:

            Y_0 := Y_0 + Inverse_E_Digit_N * ((One - Scaled_X_Fraction * Y_0**N) * Y_0);

         else  --  divide by N:

            Y_0 := Y_0 + (One - Scaled_X_Fraction * Y_0**N) * Y_0 / E_Digit_N;

         end if;

         No_Correct_Bits := No_Correct_Bits * 2;
         exit when No_Correct_Bits > Max_Available_Bits;

      end loop;

      Result := Scaling (Y_0, (-Scaled_Exp_Of_X) / E_Integer (N));
      --  Product of Y_0 = Scaled_X_Fraction**(-1/N) with
      --  Radix**(-Scaled_Exponent(X) / N) equals X**(-1/N).

      return Result;

   end Reciprocal_Nth_Root;

   ------------
   -- Divide --
   ------------

   --  Uses Newton's method: Y_k+1 = Y_k + (1 - A*Y_k) * Y_k to get Z / A.
   --  Requires call to 1 / Make_Real(A).

   function Divide (Z, X : e_Real)
      return e_Real
   is
      Exponent_Of_X           : E_Integer;
      Result, X_Fraction, Y_0 : e_Real;
      No_Correct_Bits       : E_Integer;
   begin

      if Are_Equal (X, Zero) then
         raise E_Argument_Error;
      end if;

      if Are_Equal (X, Positive_Infinity) then -- like underflow
         return Zero;
      end if;

      if Are_Equal (X, Negative_Infinity) then -- like underflow
         return Zero;
      end if;

      -- Argument reduction.   Break X into Fraction and Exponent.
      -- Iterate to get inverse of fraction.  Negate to get inverse of Exp.

      Exponent_Of_X := Exponent (X);
      X_Fraction    := Fraction (X);

      -- Get the first 2 radix digits correct (48 bits usually).  Remember that the
      -- Newton's iteration here produces 1/(X_Fraction).  The result will be
      -- the product of the newton's iteration and Radix to the power Exp_Scale_Val.

      Y_0 := Make_Extended (1.0 / Make_Real (X_Fraction));
      --  Starting val in Newton's iteration.

      No_Correct_Bits := Real'Machine_Mantissa - 2;

      -- Iterate:

      loop

       --Y_0:= Y_0 *(Two_Digit + (-X_Fraction) * Y_0);  -- faster, much less accurate
         Mult (Y_0, (Two - X_Fraction * Y_0)); -- faster, much less accurate

         No_Correct_Bits := No_Correct_Bits * 2;
         exit when No_Correct_Bits >= Max_Available_Bits / 2 + 1;
         --  final correction is outside the loop.

      end loop;

      Result := Z * Y_0; -- Z / X

      Result := Result + (Z - X_Fraction * Result) * Y_0; --bst so far
     --  Y_0 := Y_0 + (One - X_Fraction * Y_0) * Y_0;
     --  The iteration for Y_0 is the final step for 1/X. Multiplied by Z to get Result.

      Result := Scaling (Result, -Exponent_Of_X);
      --  Product of 1/Fraction(X) with Radix**(-Exponent(X)) equals 1/X.

      return Result;

   end Divide;

   ----------------
   -- Reciprocal --
   ----------------

   --  Uses Newton's method: Y_k+1 = Y_k + (1 - A*Y_k) * Y_k to get 1 / A.

   function Reciprocal (X : e_Real)
      return e_Real
   is
      Exponent_Of_X           : E_Integer;
      Result, X_Fraction, Y_0 : e_Real;
      No_Correct_Bits       : E_Integer;
   begin

      if Are_Equal (X, Zero) then
         raise E_Argument_Error;
      end if;

      if Are_Equal (X, Positive_Infinity) then -- like underflow
         return Zero;
      end if;

      if Are_Equal (X, Negative_Infinity) then -- like underflow
         return Zero;
      end if;

      -- Argument reduction.   Break X into Fraction and Exponent.
      -- Iterate to get inverse of fraction.  Negate to get inverse of Exp.

      Exponent_Of_X := Exponent (X);
      X_Fraction    := Fraction (X);

      -- Newton's iteration here produces 1/(X_Fraction).  Final result will be
      -- the product of the newton's iteration and Radix to the power Exp_Scale_Val.

      Y_0 := Make_Extended (1.0 / Make_Real (X_Fraction));
      --  Starting val in Newton's iteration.

      No_Correct_Bits := Real'Machine_Mantissa - 2;

      -- Iterate:

      loop

       --Y_0:= Y_0 *(Two_Digit + (-X_Fraction) * Y_0);  -- faster, much less accurate
         Mult (Y_0, (Two - X_Fraction * Y_0)); -- faster, much less accurate

         No_Correct_Bits := No_Correct_Bits * 2;
         exit when No_Correct_Bits > Max_Available_Bits / 2 + 1;
         --  final correction is below.

      end loop;

      Y_0 := Y_0 + (One - X_Fraction * Y_0) * Y_0; -- accurate final step.

      Result := Scaling (Y_0, -Exponent_Of_X);
      --  Product of 1/Fraction(X) with Radix**(-Exponent(X)) equals 1/X.

      return Result;

   end Reciprocal;

   ---------------------
   -- Reciprocal_Sqrt --
   ---------------------

   --  Uses Newton's method: Y_k+1 = Y_k + (1 - A*Y_k**2) * Y_k / 2
   --  to get 1 / sqrt(A).  Multiply by A to get desired result; then refine.

   function Reciprocal_Sqrt (X : e_Real)
      return e_Real
   is
      Result          : e_Real;
      X_scaled, Y_0   : e_Real;
      Exp_Scale_Val   : E_Integer;
      No_Correct_Bits : E_Integer;
   begin

      if X < Zero then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, Positive_Infinity) then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, Negative_Infinity) then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, Zero) then
          raise E_Argument_Error;
      end if;


      -- Break X into Fraction and Exponent.  If Exponent is
      -- odd then add or subtract 1.  (Increase X if X < 1, decrease otherwise.)
      -- Only important thing is scale X
      -- down to somewhere near 1, and to scale X by an even power of Radix.
      -- We break X up because the Newton's method works better on X_scaled than
      -- than on X in general.  Also we use  Sqrt(Make_Real(X_scaled)) to start
      -- things off for Newton's method, so we want X_scaled in range of Sqrt(Real).

      Exp_Scale_Val := Exponent (X); -- what if X is 0, inf, etc..???

      --  Exp is odd.  Make it even, but keep things in range of e_Real:

      if Abs (Exp_Scale_Val) mod 2 /= 0 then
         if Exp_Scale_Val < 0 then
            Exp_Scale_Val :=  Exp_Scale_Val + 1;
         else
            Exp_Scale_Val :=  Exp_Scale_Val - 1;
         end if;
      end if;

      X_scaled := Scaling (X, -Exp_Scale_Val);


      -- Take Sqrt by dividing the even Exp_Scale_Val by 2, and by taking
      -- the SQRT of X_scaled.  Start the iteration off with a call to SQRT
      -- in the standard library for type Real.

      Exp_Scale_Val  := Exp_Scale_Val / 2;
      Y_0            := Make_Extended (Sqrt (1.0 / Make_Real (X_scaled)));
      --  Starting val in Newton's iteration.

      No_Correct_Bits := Real'Machine_Mantissa - 2;

      loop

       --Y_0:= Y_0 + Half_Digit * ((One - X_scaled * (Y_0 * Y_0)) * Y_0);
       --Y_0:= Y_0*(Half_Digit * (Three - X_scaled * (Y_0 * Y_0))); --inaccurate
         Mult (Y_0, Half_Digit * (Three - X_scaled * (Y_0 * Y_0)));

         No_Correct_Bits := No_Correct_Bits * 2;
         exit when No_Correct_Bits > Max_Available_Bits / 2 + 1;
         --  final correction is below.

      end loop;

    -- both work:
    --Y_0 := Y_0 + Half_Digit * ((One - X_scaled * (Y_0 * Y_0)) * Y_0);
      Y_0 := Y_0 - Half_Digit * (X_scaled * (Y_0 * Y_0) - One) * Y_0;
      Result := Scaling (Y_0, -Exp_Scale_Val);

      return Result;

   end Reciprocal_Sqrt;

   ----------
   -- Sqrt --
   ----------

   --  Uses Newton's method: Y_k+1 = Y_k + (1 - A*Y_k**2) * Y_k / 2
   --  to get 1 / sqrt(A).  Multiply by A to get desired result; then refine.
   --  Requires call to Sqrt(Real) and assumes that this call gets
   --  the first radix digit correct. (It almost
   --  always gets 53 bits correct.)

   function Sqrt (X : e_Real)
      return e_Real
   is
      Result, X_scaled, Y_0 : e_Real;
      Exp_Scale_Val         : E_Integer;
      No_Correct_Bits       : E_Integer;
   begin

      if X < Zero then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, Positive_Infinity) then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, Negative_Infinity) then
          raise E_Argument_Error;
      end if;

      if Are_Equal (X, Zero) then
          return Zero;
      end if;

      --if Are_Equal (X, One) then  -- Ada9X.
          --return One;
      --end if;

      -- Break X into Fraction and Exponent.  If Exponent is
      -- odd then add or subtract 1.  (Increase X if X < 1, decrease otherwise.)
      -- Only important thing is scale X
      -- down to somewhere near 1, and to scale X by an even power of Radix.
      -- We break X up because the Newton's method works better on X_scaled than
      -- than on X in general.  Also we use  Sqrt(Make_Real(X_scaled)) to start
      -- things off for Newton's method, so we want X_scaled in range of Sqrt(Real).

      Exp_Scale_Val := Exponent (X); -- what if X is 0, inf, etc..???
      -- This Exp is powers of Radix = 2**30 or 2**29.

      --  Exp is odd.  Make it even, but keep things in range of e_Real:

      if Abs (Exp_Scale_Val) mod 2 /= 0 then
         if Exp_Scale_Val < 0 then
            Exp_Scale_Val :=  Exp_Scale_Val + 1;
         else
            Exp_Scale_Val :=  Exp_Scale_Val - 1;
         end if;
      end if;

      X_scaled := Scaling (X, -Exp_Scale_Val);

      Exp_Scale_Val  := Exp_Scale_Val / 2;
      Y_0            := Make_Extended (Sqrt (1.0 / Make_Real (X_scaled)));
      --  Starting val in Newton's iteration.

      No_Correct_Bits := Real'Machine_Mantissa - 2;

      loop

       --Y_0:= Y_0 + Half_Digit * ((One - X_scaled * Y_0 * Y_0) * Y_0);
       --Y_0:= Y_0*(Half_Digit * (Three - X_scaled * Y_0 * Y_0)); --inaccurate
         Mult (Y_0, Half_Digit * (Three - X_scaled * Y_0 * Y_0));

         No_Correct_Bits := No_Correct_Bits * 2;
         exit when No_Correct_Bits > Max_Available_Bits / 2 + 1;
         --  Have the final correction outside the loop.

      end loop;

      Result := Y_0 * X_scaled;   -- now it's  SQRT(X_scaled); Y_0 = 1/SQRT(X_scaled)

      Result := Result + Half_Digit * (X_scaled - Result*Result) * Y_0; --important

      Result := Scaling (Result, Exp_Scale_Val);   -- equals   SQRT(X).

      return Result;

   end Sqrt;

  -------------
  -- Make_Pi --
  -------------

  --  This is an important independent test of Arcsin (hence Sin and Arcos).
  --  Once we verify that Arcsin (hence Sin) is correct, can test other
  --  arguments with (eg) Sin (A + B) = etc.
  --  Has much greater error than 4*Arcsin(1/Sqrt(2)).
  --  Need Pi to full precision for the trigonometic functions.
  --  Here is the Salamin-Brent algorithm.
  --  A_0 = 1.0,  B_0 = 1.0/Sqrt(2.0),  and D_0 = Sqrt(2.0) - 0.5.
  --
  --  A_k = (A_{k-1} + B_{k-1}) / 2
  --  B_k = Sqrt (A_{k-1} * B_{k-1})
  --  D_k = D_{k-1} - 2**k * (A_k - B_k)**2
  --
  --  Then P_k =  (A_k + B_k)**2 / D_k converges quadratically to
  --  Pi.  All steps must be done at full precision.
  --
  --    function Make_Pi return e_Real is
  --       A_0, B_0, D_0 : e_Real;
  --       A_1, B_1, D_1 : e_Real;
  --       C, Result : e_Real;
  --       Two_To_The_k  : E_Digit;
  --       Two_To_The_20 : constant E_Digit := Make_E_Digit (2.0**20);
  --       We_Are_Finished : Boolean := False;
  --       No_Of_Correct_Digits : E_Integer;
  --       Old_Precision : constant E_Integer := Present_Precision;
  --    begin
  --
  --      A_0 := One;
  --      B_0 := E_Inverse_Sqrt_2;
  --      D_0 := Two * E_Inverse_Sqrt_2 - Half;
  --  --  this give Pi_0 = 3.1877, or error of about 1.47 %.  This is smaller
  --  --  1 part in 64, so 6 bits correct.  There follows (k=1) 12, (k=2) 24.
  --  --  So requires two more iterations to get one digit, 3 to get 2
  --  --  digits.  Seems to work better than this estimate.
  --
  --      No_Of_Correct_Digits := 1;
  --
  --  --  The following loop should get us up to half a million digits.  In
  --  --  the unlikely case you need more, then another loop follows.
  --  --  k in 1..7 gives you 33 Radix 2**24 digits.
  --
  --      for k in 1..20 loop
  --
  --        Two_To_The_k := Make_E_Digit (2.0**k);
  --
  --        A_1 := Half_Digit * (A_0 + B_0);
  --        B_1 := Sqrt (A_0 * B_0);
  --        C   := (A_1 - B_1);
  --        D_1 := D_0 - Two_To_The_k * (C * C);
  --
  --        if k >= 3 then
  --       --  We did 3rd iteration to get 2 correct digits.
  --       --  No_Correct.. was initialized to 1.
  --           No_Of_Correct_Digits := No_Of_Correct_Digits * 2;
  --        end if;
  --    --  Should be OK overflow-wise here.  Range of E_Integer is 4 times
  --    --  the limit set by Max_Available_Precision.
  --
  --        if No_Of_Correct_Digits > e_Real_Machine_Mantissa then
  --          We_Are_Finished := True;
  --          exit;
  --        end if;
  --
  --        A_0 := A_1; B_0 := B_1; D_0 := D_1;
  --
  --
  --      end loop;
  --
  --  --  We want to optimize the calculation of D_1 above by multiplying
  --  --  by an E_Digit on the left (Two_To_The_k) instead of an e_Real.
  --  --  Stop doing this at Two_To_The_k = 2**20 to stay in range of E_Digit.
  --  --  Below we finish up if necessary by multiplying twice..still much
  --  --  more efficient than e_Real*e_Real.
  --
  --    if  not We_Are_Finished then    -- keep trying
  --      for k in 21..40 loop
  --
  --        Two_To_The_k := Make_E_Digit (2.0**(k-20));
  --
  --        A_1 := Half_Digit * (A_0 + B_0);
  --        B_1 := Sqrt (A_0 * B_0);
  --        C   := (A_1 - B_1);
  --        D_1 := D_0 - Two_To_The_k * (Two_To_The_20 * (C * C));
  --
  --        No_Of_Correct_Digits := No_Of_Correct_Digits * 2;
  --        exit when No_Of_Correct_Digits > e_Real_Machine_Mantissa;
  --
  --        A_0 := A_1; B_0 := B_1; D_0 := D_1;
  --
  --      end loop;
  --    end if;
  --
  --      C := (A_1 + B_1);
  --      Result :=  C * C / D_1;
  --
  --      Set_No_Of_Digits_To_Use (Old_Precision); -- Restore precision.
  --
  --      return Result;
  --
  --    end Make_Pi;

  ----------------
  -- Make_Log_2 --
  ----------------

  --  Important independent test of Log(X).  Verify that Log(X) is correct
  --  at X = 2, and use (eg) Log(XY) = Log(X) + Log(Y) to test other vals.
  --  This is for testing other routines: Has greater error than Log(X).
  --  Log_2, hopefully, for testing purposes.
  --  A_0 = 1.0,  B_0 = Two**(2-M);
  --
  --  A_k = (A_{k-1} + B_{k-1}) / 2
  --  B_k = Sqrt (A_{k-1} * B_{k-1})
  --
  --  Then Log(2) =  Pi / (2 * B_k * m)  ???
  --
  -- Here M = N/2+1 where N = 29*e_Real_Machine_Mantissa = number of bits desired.
  --
  --    function Make_Log_2 return e_Real is
  --       A_0, B_0 : e_Real;
  --       A_1, B_1 : e_Real;
  --       Result   : e_Real;
  --       We_Are_Finished : Boolean := False;
  --       No_Of_Correct_Digits : E_Integer;
  --       N : Integer := 24 * Integer(e_Real_Machine_Mantissa);    -- Upper estimate.
  --       M : Integer := N/2 + 24;     -- Need only N/2 + 1
  --    begin
  --
  --      A_0 := One;
  --      B_0 := Two**(2-M);    -- clean this up with the scaling ftcn.
  --
  --      No_Of_Correct_Digits := 1;
  --
  --    --  The following loop should get us up to half a million digits.  In
  --    --  the unlikely case you need more, then another loop follows.
  --
  --      for k in 1..16 loop    -- I suspect far fewer than 20 iterations required.
  --
  --        A_1 := Half_Digit * (A_0 + B_0);
  --        B_1 := Sqrt (A_0 * B_0);
  --
  --        A_0 := A_1; B_0 := B_1;
  --
  --      end loop;
  --
  --      Result := Half_Digit * e_Pi / (B_1 * Make_Extended(Real(M)));
  --
  --      Set_No_Of_Digits_To_Use (Old_Precision);   --  Restore precision.
  --
  --      return Result;
  --
  --    end Make_Log_2;

  ----------
  -- e_Pi --
  ----------

  --  Returns Pi to Max Available Precision.  Use Arcsin, cause it has
  --  much lower error than Make_Pi.
  --  Used for scaling trig functions, etc.

  function e_Pi return e_Real is
  begin

    if not Pi_memory.Initialized then
     --Pi_memory.Val := Make_Pi;
       Pi_memory.Val := Four_Digit * e_Quarter_Pi;
       --  Only works because arg is so small no scaling by E_pi is done.

       Pi_memory.Initialized := True;

    end if;

    return Pi_memory.Val;

  end e_Pi;

  ------------------
  -- e_Inverse_Pi --
  ------------------

  --  Returns Pi to Max Available Precision.  Use Arcsin, cause it has
  --  lower error than Make_Pi.
  --  Used for scaling trig functions, etc.

  function e_Inverse_Pi return e_Real is
  begin

    if not Inverse_Pi_memory.Initialized then
       Inverse_Pi_memory.Val         := (+0.25) / e_Quarter_Pi;
       Inverse_Pi_memory.Initialized := True;
    end if;

    return Inverse_Pi_memory.Val;

  end e_Inverse_Pi;

  ------------------
  -- e_Quarter_Pi --
  ------------------

  --  Returns Pi/4 to Max Available Precision.
  --  Used for scaling trig functions, etc.

  function e_Quarter_Pi return e_Real is
  begin

    if not Quarter_Pi_memory.Initialized then
       Quarter_Pi_memory.Val := (+1.5) * Arcsin (Half);
       Quarter_Pi_memory.Initialized := True;
    end if;

    return Quarter_Pi_memory.Val;

  end e_Quarter_Pi;

  ----------------------
  -- e_Inverse_Sqrt_2 --
  ----------------------

  --  Returns 1/Sqrt(2.0) to Max Available Precision.
  --  Used for making Pi.

  function e_Inverse_Sqrt_2 return e_Real is
  begin

    if not Inverse_Sqrt_2_memory.Initialized then
       Inverse_Sqrt_2_memory.Val := Reciprocal_Nth_Root (Two, 2);
       Inverse_Sqrt_2_memory.Initialized := True;
    end if;

    return Inverse_Sqrt_2_memory.Val;

  end e_Inverse_Sqrt_2;

  --------------------------
  -- e_Half_Inverse_Log_2 --
  --------------------------

  --  Returns Exp(1.0) to Max Available Precision.
  --  Used for scaling arguments of Exp.

  function e_Half_Inverse_Log_2 return e_Real is
  begin

    if not Half_Inverse_Log_2_memory.Initialized then
       Half_Inverse_Log_2_memory.Val := Half / E_Log_2;
       Half_Inverse_Log_2_memory.Initialized := True;
    end if;

    return Half_Inverse_Log_2_memory.Val;

  end e_Half_Inverse_Log_2;

  --------------
  -- e_Log_2  --
  --------------

  --  Returns Log(2.0) to Max Available Precision.
  --  Used for scaling Exp(X).

  function e_Log_2 return e_Real is
  begin

    if not Log_2_memory.Initialized then
       Log_2_memory.Val := Log (Two);
       Log_2_memory.Initialized := True;
    end if;

    return Log_2_memory.Val;

  end e_Log_2;

end Extended_Real.Elementary_Functions;

