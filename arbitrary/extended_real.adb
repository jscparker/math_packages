
-----------------------------------------------------------------------
-- package body Extended_Real, extended precision floating point arithmetic
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

-- Internally the extended numbers are stored in such a way that the
-- value of e_Real number X  is
--
--                    Max
--   X = Radix**Exp * SUM {Radix**(-I) * Digit(I)}.
--                    I=0
--
-- Externally, the user sees e_Real (via the Exponent, and Fraction
-- attribute functions) as tho' it were normalized.  In other words, the
-- value of X is
--
--                     Max
--  X = Radix**Exp_n * SUM {Radix**(-I-1) * Digit(I)}
--                     I=0
--
-- Exp_n is called the "normalized" exponent.  If Exp_n is the normalized exponent
-- then, say, a binary number would be written:
--
--   0.111011010001 * 2**(Exp_n).
--
-- In other words the first binary digit in the mantissa is of power 2**(-1).
-- It is important to know this because the function Real'Exponent(x) returns
-- the *normalized* exponent, and the function Real'Fraction(x) returns
-- x * 2**(-Exp_n) where Exp_n is the normalized exponent.  So in the above case,
-- 'Fraction would return 0.111011010001.
-- Also, in normalized form, the first binary digit of the mantissa is always
-- non-zero.


package body Extended_Real is

   Disable_Program_Error_Tests : constant Boolean := False;
   --  A few assertions and the like.

   ----------------------------------------
   -- Shift_Right_2x_No_of_Bits_in_Radix --
   ----------------------------------------

   function Shift_Right_2x_No_of_Bits_in_Radix
     (Digit : in Digit_Type)
      return Digit_Type
   is
   begin
    --  Works in general, but use only for flt types:
    --return Digit_Type (Real'Floor (Real (Digit) * Inverse_Radix_Squared));
      return Digit / 2**(2*No_of_Bits_in_Radix);
   end Shift_Right_2x_No_of_Bits_in_Radix;

   pragma Inline (Shift_Right_2x_No_of_Bits_in_Radix);

   -------------------------------------
   -- Shift_Right_No_of_Bits_in_Radix --
   -------------------------------------

   function Shift_Right_No_of_Bits_in_Radix
     (Digit : in Digit_Type)
      return Digit_Type
   is
   begin
    --  Works in general, but use only for flt types:
    --return Digit_Type (Real'Floor (Real (Digit) * Inverse_Radix));
      return Digit / 2**No_of_Bits_in_Radix;
   end Shift_Right_No_of_Bits_in_Radix;

   pragma Inline (Shift_Right_No_of_Bits_in_Radix);


   -----------------
   -- Digit_Floor --
   -----------------

   --  Used by "/" and Make_Extended:

   function Digit_Floor (X : Real) return Digit_Type is
   begin
      return  Digit_Type (Real'Floor (X));
   end Digit_Floor;

   pragma Assert (Real'Machine_Mantissa > No_Of_Bits_In_Radix); -- ie 53 > 30

   pragma Inline (Digit_Floor);

   ----------------------------------
   -- Minimum_No_Of_Digits_Allowed --
   ----------------------------------

   --  Return setting of the constant Min_No_Of_Digits.
   --
   function Minimum_No_Of_Digits_Allowed return e_Integer is
   begin return Min_No_Of_Digits;
   end Minimum_No_Of_Digits_Allowed;

   ----------------------------
   -- Number_Of_Guard_Digits --
   ----------------------------

   --  Return setting of the constant No_Of_Guard_Digits.
   --
   function Number_Of_Guard_Digits return e_Integer is
   begin return No_Of_Guard_Digits;
   end Number_Of_Guard_Digits;


-- SECTION IV.
--
-- Ada94 attributes.  More information on the machine model is given
-- above in the introduction.


   ---------------------------
   -- e_Real_Machine_Rounds --
   ---------------------------

   function e_Real_Machine_Rounds return Boolean is
   begin return False;
   end e_Real_Machine_Rounds;

   ------------------------------
   -- e_Real_Machine_Overflows --
   ------------------------------

   function e_Real_Machine_Overflows return Boolean is
   begin return False;
   end e_Real_Machine_Overflows;

   -------------------------
   -- e_Real_Signed_Zeros --
   -------------------------

   function e_Real_Signed_Zeros return Boolean is
   begin return False;
   end e_Real_Signed_Zeros;

   -------------------
   -- e_Real_Denorm --
   -------------------

   function e_Real_Denorm return Boolean is
   begin return False;
   end e_Real_Denorm;

   -------------------------
   -- e_Real_Machine_Emax --
   -------------------------

   function e_Real_Machine_Emax return e_Integer is
   begin return Max_Exponent;
   end e_Real_Machine_Emax;

   -------------------------
   -- e_Real_Machine_Emin --
   -------------------------

   function e_Real_Machine_Emin return e_Integer is
   begin return Min_Exponent;
   end e_Real_Machine_Emin;

   -----------------------------
   -- e_Real_Machine_Mantissa --
   -----------------------------

   -- Number of digits in machine mantissa.

   function e_Real_Machine_Mantissa 
      return e_Integer is
   begin 
      return  Mantissa'Length;
   end e_Real_Machine_Mantissa;

   --------------------------
   -- e_Real_Machine_Radix --
   --------------------------

   function e_Real_Machine_Radix return Real is
   begin return Real_Radix;
   end e_Real_Machine_Radix;

   ----------------------------
   -- e_Real_Machine_Epsilon --
   ----------------------------

   function e_Real_Machine_Epsilon return e_Real is
   begin
      return e_Real_Model_Epsilon_1;
   end e_Real_Machine_Epsilon;

   --------------------------
   -- e_Real_Model_Epsilon --
   --------------------------

   function e_Real_Model_Epsilon return e_Real is
   begin
      return e_Real_Model_Epsilon_2;
   end e_Real_Model_Epsilon;

   ----------------------------
   -- e_Real_Model_Epsilon_1 --
   ----------------------------

   --  1 unit in the larger of the 2 guard digits.

   function e_Real_Model_Epsilon_1 return e_Real is
      Result : e_Real; -- equals Zero
   begin
      Result.Is_Zero   := False;
      Result.Digit (0) := Digit_One; -- 1 here with Exp=0 => Eps=radix**(-1)
      Result.Exp       := -(e_Real_Machine_Mantissa-1);
      return Result;
   end e_Real_Model_Epsilon_1;

   ----------------------------
   -- e_Real_Model_Epsilon_2 --
   ----------------------------

   --  Guard_Digits = 2 always; assume neither of them is correct.
   --  if there's 3 digits of Radix 2^30 then eps is 2^(-30).
   --  if there's 4 digits of Radix 2^30 then eps is 2^(-60).
   --  if there's 5 digits of Radix 2^30 then eps is 2^(-90) or ~ 10**(-30).

   function e_Real_Model_Epsilon_2 return e_Real is
      Result : e_Real; -- equals Zero
   begin
      Result.Is_Zero   := False;
      Result.Digit (0) := Digit_One; -- 1 here with Exp=0 => Eps=radix**(-1)
      Result.Exp       := -(e_Real_Machine_Mantissa-2);
      return Result;
   end e_Real_Model_Epsilon_2;

   --------------
   -- Exponent --
   --------------

   -- For every value x of a floating point type T, the normalized exponent
   -- of x is defined as follows:
   --     - the normalized exponent of zero is (by convention) zero;
   --     - for nonzero x, the normalized exponent of x is the unique
   --
   --                                         k-1                        k
   --      integer k such that T'Machine_Radix   <= |x| < T'Machine_Radix .
   --
   -- For example, if x = 0.1101011 * 2**1 then k = 1 since
   -- 2**0 <= x < 2**1.  If X = 0.100000 * 2**1 then 2**0 = x < 2**1.
   --
   -- BY CONVENTION, the normalized exponent of 0 is 0.  (Internally it ain't.)

   function Exponent (X : e_Real) return e_Integer is
      Normalized_Exponent : constant e_Integer := X.Exp + 1;
   begin
      if X.Is_Zero then
         return 0;
      else
         return Normalized_Exponent;
      end if;
      -- proper choice because internal form of e_Real is not
      -- normalized; internally, Exp is smaller than the normalized exp
      -- by 1, because the internal mantissa is larger by a factor
      -- of Radix than the normalized mantissa.
   end Exponent;

   ------------------
   -- Leading_Part --
   ------------------

   --  Let k be the normalized exponent of x.  The function Leading_Part(x)
   --  yields the value
   --
   --                             k-r                 k-r
   --         - |x/T'Machine_Radix   |*T'Machine_Radix   ,
   --         when x is nonnegative and r is positive;
   --                             k-r                 k-r
   --         - |x/T'Machine_Radix   |*T'Machine_Radix   ,
   --         when x is negative and r is positive.
   --
   --  Constraint_Error is raised when r is zero or negative.
   --  A zero result, which can only occur when x is zero, has
   --  the sign of X.
   --  Returns the leading digits of X in the range 0..Radix_Digits-1,
   --  regardless of sign of X.  Other digits are set to 0.0.  The exponent
   --  is unchanged. Only X = Zero returns Zero.  Not sure about infinity.
   --  Makes sense that Leading_Part of inf is inf.
   --  Notice assume Digit_Index'Last = Ultimate digit.  The results are
   --  consistant with a dynamic Digit_Index'Last, just a little slower.
   --
   function Leading_Part
     (X            : e_Real;
      Radix_Digits : e_Integer)
      return e_Real
   is
      No_Of_Digits : constant e_Integer := Ultimate_Digit + 1;
      Result : e_Real := X;
   begin
      -- By convention:
      if Radix_Digits <= 0 then
         raise Constraint_Error with "Must have Radix_Digits > 0.";
      end if;

      if X.Is_Zero then
         return Zero;
      end if;

      if X.Is_Infinite then
         return X;
      end if;

      if Radix_Digits >= No_Of_Digits then
         return X;
      end if;

      -- So now Radix_Digits < No_Of_Digits  which implies that
      --        Radix_Digits <= Ultimate_Digit

      for I in Digit_Index range Radix_Digits .. Ultimate_Digit loop
         Result.Digit (I) := Digit_Zero;
      end loop;
      --  The above uses the fact that Digit_Index starts at 0.

      return Result;

   end Leading_Part;

   --------------
   -- Fraction --
   --------------

   --                                                 -k
   --  The function yields the value x*T'Machine_Radix  , where
   --  k is the normalized exponent of x.  A zero result, which
   --  can only occur when x is zero, has the sign of X.
   --  Not sure about inf, so raise contraint error, because the
   --  user may make assumptions about the size of the exponent returned
   --  by Exponent.
   --
   function Fraction (X : e_Real) return e_Real is
      X2 : e_Real := X;
   begin
      if X.Is_Zero then
         return Zero;
      end if;

      if X.Is_Infinite then
         raise Constraint_Error with "Cannot take Fraction of inf.";
      end if;

      X2.Exp := -1;
      --  Proper choice because format of e_Real
      --  is not normalized.  The effect of the -1 is
      --  to shift the internal format down to the normalized format.

      return X2;

   end Fraction;

   -------------
   -- Compose --
   -------------

   --   S'Compose (X, Exp) = Fraction (X) * Machine_Radix ** Exp
   --
   --                                       e-k
   --   Let v be the value X*T'Machine_Radix   , where k is the
   --   normalized exponent of X.  If v is a machine number of
   --   the type T, or if |v|GT'Model_Small, the function yields
   --   v; otherwise, it yields either one of the machine
   --   numbers of the type T adjacent to v.    Constraint_Error
   --   is optionally raised if v is outside the base range of
   --   S. A zero result has the sign of Fraction when S'Signed_
   --   Zeros is True.
   --
   function Compose
     (Fraction : e_Real;
      Exponent : e_Integer)
      return e_Real
   is
      X2 : e_Real := Fraction;
   begin
      if Fraction.Is_Zero then
         return Zero;
      end if;

      if Fraction.Is_Infinite then
         raise Constraint_Error with "Cannot compose inf.";
      end if;

      X2.Exp := Exponent - 1;
      --  The minus 1 comes from the Fraction(X) operation.

      if X2.Exp < Min_Exponent or else X2.Exp > Max_Exponent then
         raise Constraint_Error with "Exponent out of range in Compose operation.";
      end if;

      return X2;
   end Compose;

   -------------
   -- Scaling --
   -------------

   --  S'Scaling (X, Exp)
   --                                      Exp
   --  Let v be the value X*T'Machine_Radix  .  If v is a
   --  machine number of the type T, or if |v|GT'Model_Small,
   --  the function yields v; otherwise, it yields either one
   --  of the machine numbers of the type T adjacent to v.
   --  Constraint_Error is optionally raised if v is outside
   --  the base range of S. A zero result has the sign of X
   --  when S'Signed_Zeros is True.
   --
   function Scaling
     (X : e_Real;
      Adjustment : e_Integer)
      return e_Real
   is
      X2 : e_Real := X;
   begin
      if X.Is_Zero then
         return Zero;
      end if;

      if X.Is_Infinite then
         raise Constraint_Error with "Cannot scale inf.";
      end if;

      X2.Exp := X.Exp + Adjustment;

      if X2.Exp < Min_Exponent or else X2.Exp > Max_Exponent then
         raise Constraint_Error with "Exp out of range in Scaling operation.";
      end if;

      return X2;
   end Scaling;

   ----------------
   -- Truncation --
   ----------------

   --  Strip off fraction for both + and - numbers.
   --  This sets all digits beyond the decimal point to 0.0.
   --  The function yields the value [x] when X > 0.0.  ([x] by defn
   --  is floor: largest int less than or equal to x.)  When x is zero,
   --  the result has the sign of X; a zero result otherwise has a
   --  positive sign.  When X < 0.0, then it's [x]+1 unless [x] = x.
   --  Notice assume Digit_Index'Last = Ultimate_Digit.

   function Truncation
     (X : e_Real)
      return e_Real
   is
      Result : e_Real;
      First_Zeroed_Out_Digit : Digit_Index;
      No_Of_Digits         : constant e_Integer := Ultimate_Digit + 1;
      No_Of_Digits_To_Keep : constant e_Integer := X.Exp + 1;
      --  Normalized Exponent of X.  Remember, if X.Exp = 0,then X has one
      --  non-zero digit that's >= 1.0.  So (Exp + 1)=1 means keep 1 digit 
   begin

      if X.Is_Zero then
         return Zero;
      end if;

      if No_Of_Digits_To_Keep < 1 then
         return Zero;
      end if;
      --  According to the internal format of X, X.Exp = 0 means X >= 1.0.
      --  If X.Exp < 0 then X < 1.0.

      if X.Is_Infinite then
         return X;
         --Print_Text ("Cannot truncate inf.");
         --raise Constraint_Error;
      end if;

      if No_Of_Digits_To_Keep >= No_Of_Digits then
         return X;
      end if;

      -- So now No_Of_Digits_To_keep < No_Of_Digits  which implies that
      --        No_Of_Digits_To_keep <= Ultimate_Digit
      -- Remember, Digit_Index starts at 0, and is a subtype of e_Integer.

      Result := X;
      First_Zeroed_Out_Digit := Digit_Index'First + No_Of_Digits_To_keep;

      for I in First_Zeroed_Out_Digit .. Ultimate_Digit loop
         Result.Digit (I) := Digit_Zero;
      end loop;

      return Result;

   end Truncation;

   -------------
   -- Ceiling --
   -------------

   --  Let [x] = floor(x).
   --  Function yields the value [x]+1, unless [x] = x, in which case,
   --  function returns x.  When x is zero, the
   --  result has the sign of x; a zero result otherwise has a
   --  negative sign when S'Signed_Zeros is True.

   function Ceiling (X : e_Real) 
      return e_Real 
   is
      Result             : e_Real;
      Lowest_Order_Digit : e_Integer;
   begin
      if X.Is_Zero then
         return Zero;
      end if;

      if X.Is_Infinite then
         return X;
      end if;

      --  Step 1.
      --  Special case, X.Exp < 0 so X < 1.0. Have to round 0.0 up to 1.0.
      --  We know X /= 0.0.  if X.Exp < 0, then know |X| < 1.0.  If
      --  X > 0.0, then result must be 1.0;  if X < 0.0 then result
      --  is 0.0.  So can save some time:
      if X.Exp < 0 then
         if X.Is_Positive then
            return One;
         else
            return Zero;
         end if;
      end if;

      --  Step 2.
      --  Now know that X.Exp >= 0.  Have to do some work:

      Result := Truncation (X);

      if not X.Is_Positive then  -- we're done by defn: Ceiling = Trunc.
         return Result;
      end if;

      --  Now know that X is positive.  Must add one if Trunc(x) /= x.
      --  We also know that Result.Exp >= 0, by defn of Trunc.
      Lowest_Order_Digit := Result.Exp;

      if (Lowest_Order_Digit <= Ultimate_Digit) and then Result < X then
         if Result.Digit (Lowest_Order_Digit) < Digit_Radix_Minus_1 then
            Result.Digit (Lowest_Order_Digit)
                          := Result.Digit (Lowest_Order_Digit) + Digit_One;
         else
            Result := Result - One;
         end if;
      end if;

      return Result;

   end Ceiling;

   -----------
   -- Floor --
   -----------

   --  Function yields the value [x].  ([x] by defn is floor: largest int
   --  less than or equal to x.)  This equals trunc(x) when x > 0.0,
   --  and = trunc(x)-1, when X < 0.0, unless Trunc(x) = x.  When x = 0.0,
   --  result has the sign of x; a zero result otherwise has a
   --  negative sign when S'Signed_Zeros is True.

   function Floor (X : e_Real) return e_Real is
      Result             : e_Real;
      Lowest_Order_Digit : e_Integer;
   begin

      if X.Is_Zero then
         return Zero;
      end if;

      if X.Is_Infinite then
        return X;
      end if;

      --  Step 1.
      --  Special case, X.Exp < 0.
      --  Know X /= 0.0.  if X.Exp < 0, then know |X| < 1.0.  If then
      --  X > 0.0, then result must be 0.0;  if X < 0.0 then result
      --  is -1.0.  So can save some time:
      if X.Exp < 0 then
         if X.Is_Positive then
            return Zero;
         else
            return -One;
         end if;
      end if;

      --  Step 2.
      --  Now know that X.Exp >= 0.  Have to do some work:

      Result := Truncation (X);

      if X.Is_Positive then  -- we're done by defn: Floor = Trunc.
         return Result;
      end if;

      --  Now know that X is negative.  Must subtract one if Trunc(x) > x.
      --  We also know that Result.Exp >= 0, by defn of Trunc, unless trunc
      --  returned Zero.  (But it didn't cause X.Exp started out >= 0).

      Lowest_Order_Digit := Result.Exp;

      if (Lowest_Order_Digit <= Ultimate_Digit) and then Result > X then
         if Result.Digit (Lowest_Order_Digit) > Digit_Zero then
            Result.Digit (Lowest_Order_Digit)
                             := Result.Digit (Lowest_Order_Digit) - Digit_One;
         else
            Result := Result - One;
         end if;
      end if;

      return Result;

   end Floor;

   -------------------------------------
   -- Round_Away_Smallest_Guard_Digit --
   -------------------------------------

   --  Round away  Digit_Index'Last

   function Round_Away_Smallest_Guard_Digit (X : e_Real) return e_Real is
      Result : e_Real := X;
      Penultimate : constant Digit_Index := Digit_Index'Last-1;
   begin

      if X.Digit(Digit_Index'Last) <= Half_Radix then
         Result.Digit (Digit_Index'Last) := Digit_Zero;
      else  
         -- X is not Zero.
         Result.Digit (Digit_Index'Last) := Digit_Zero;
         if X.Digit(Penultimate) < Digit_Radix_Minus_1 then
            Result.Digit(Penultimate) := Result.Digit(Penultimate) + Digit_One;
         else
            if Result.Is_Positive then
               Result := Result + e_Real_Model_Epsilon_1; --1 unit in penultimate digit.
            else
               Result := Result - e_Real_Model_Epsilon_1;
            end if;
         end if;
      end if;

      return Result;

   end Round_Away_Smallest_Guard_Digit;

   -------------
   -- Machine --
   -------------

   --  Rounds away smallest of the 2 guard digits of X.

   function Machine (X : e_Real) return e_Real is
      Y : e_Real;
   begin
      Y := Round_Away_Smallest_Guard_Digit (X);
      return Y;
   end Machine;

   ---------------
   -- Copy_Sign --
   ---------------

   --  S'Copy_Sign denotes a function with the following specification:
   --            function S'Copy_Sign (Value, Sign : T) return T
   --        If the value of Value is nonzero, the function yields a
   --        result whose magnitude is that of Value and whose sign
   --        is that of Sign; otherwise, it yields the value zero.
   --        Constraint_Error is optionally raised if the result is
   --        outside the base range of S. A zero result has the sign
   --        of Sign when S'Signed_Zeros is True.
   --
   function Copy_Sign (Value, Sign : e_Real) return e_Real is
      Result : e_Real := Value;
   begin
      if Value.Is_Zero then
         return Zero;
      end if;

      -- following holds even if Value is inf:
      Result.Is_Positive := Sign.Is_Positive;
      return Result;

   end Copy_Sign;

   --------------
   -- Adjacent --
   --------------

   --  If t=x, the function yields x; otherwise, it yields the
   --  machine number of the type T adjacent to x in the
   --  direction of t, if that machine number exists.    If the
   --  result would be outside the base range of S, Constraint_
   --  Error is raised.  When T'Signed_Zeros is True, a zero
   --  result has the sign of X.  When t is zero, its sign has
   --  no bearing on the result.
   --
   --function Adjacent (X, Towards : e_Real)  return e_Real is
   --begin
   --end Adjacent;

   --------------
   -- Rounding --
   --------------

   --  The function yields the integral value nearest to x,
   --  rounding away from zero if x lies exactly halfway
   --  between two integers.  A zero result has the sign of X
   --  when S'Signed_Zeros is True.
   --
   function Rounding (X : e_Real)  return e_Real is
      Result : e_Real;
      Half : constant e_Real := +0.5;
      Del : e_Real;
   begin

      if X.Is_Zero then
        return Zero;
      end if;

      if X.Is_Infinite then
        return X;
      end if;

      Result := Truncation (X);
      Del := Abs (X - Result);

      if not (Del < Half) then
         if X.Is_Positive then
            Result := Result + One;  -- because Trunc chopped it toward zero
         else
            Result := Result - One;
         end if;
      end if;

      return Result;

   end Rounding;

   -----------------------
   -- Unbiased_Rounding --
   -----------------------

   --  The function yields the integral value nearest to x,
   --  rounding toward the even integer if x lies exactly
   --  halfway between two integers.  A zero result has the
   --  sign of X when S'Signed_Zeros is True.
   --
   function Unbiased_Rounding (X : e_Real)  return e_Real is
      Result : e_Real;
      Half   : constant e_Real := +0.5;
      Del    : e_Real;
      Least_Significant_Digit : e_Integer;
   begin

      if X.Is_Zero then
        return Zero;
      end if;

      if X.Is_Infinite then
        return X;
      end if;

      Result := Truncation (X);
      Del := Abs (X - Result);

    --if Del < Half then -- result is unmodified

      if Del > Half then
         if X.Is_Positive then
            Result := Result + One;  -- because Trunc chopped it toward zero
         else
            Result := Result - One;
         end if;
      end if;

      if Are_Equal (Del, Half) then

        --  Must find out if Result (= Truncation (X) = int) is even or not.
        --  If it's not even, then add (or subtract) One as above.
        --  To find out, must examine lowest order bit of least significant
        --  digit.

        Least_Significant_Digit := Exponent (Result) - 1;
        --  If Least_Significant_Digit = 0 then Exponent (Result) = +1 cause
        --  it returns the normalized Exp.

        if (Least_Significant_Digit REM 2 = 1) then  -- it's odd
            if X.Is_Positive then
              Result := Result + One;  -- because Trunc chopped it toward zero
            else
              Result := Result - One;
            end if;
        end if;

      end if;

      return Result;

   end Unbiased_Rounding;

   ---------------
   -- Remainder --
   ---------------

   --  Code, algorithm and comments from Ken Dritz. (Ada 9X Language study note.)
   --  Modified to use e_Real's, so errors my own.
   --  It does seem to work well .. much more accurate than the simple alternatives
   --  I tried for Sin, Cos arguments, (but that is all the testing I've done.)
   --  For nonzero y, let v be the value x-n*y, where n is the
   --  integer nearest to the exact value of x/y; if
   --  |n-x/y|=1/2, then n is chosen to be even.  If v is a
   --  machine number of the type T, the function yields v;
   --  otherwise, it yields zero.    Constraint_Error is raised
   --  if y is zero.  A zero result always has a positive sign.

   function Remainder
     (X, Y : e_Real) 
      return e_Real 
   is
      Residue, Temp, Reducer, Reducer_Head, Reducer_Tail, N : e_Real;
      Abs_Y : constant e_Real := Abs (Y);
      -- See comments above about the possibility of overflow here on
      -- radix-complement machines.
      Scaled_Up, Negated : Boolean;
      CONST1 : constant e_Real
           := Scaling (One, e_Real_Machine_Emin + e_Real_Machine_Mantissa - 2);
      CONST2 : constant e_Real
           := Scaling (One, e_Real_Machine_Emin + e_Real_Machine_Mantissa - 1);
   begin
      if Y.Is_Zero then
        raise Constraint_Error with "Must have Y /= 0 in function Remainder.";
      end if;
 
      Residue := X;
      Negated := False;
 
      loop
 
       -- This loop is executed at most once more than the difference between
       -- the exponents of X and Y.
 
       if Copy_Sign (One, Residue) < Zero then
          -- The following two statements are to be executed when the sign of
          -- Residue is negative, that is, when Residue is less than zero or is
          -- a negative zero.  Simply comparing Residue to zero is not good
          -- enough when T'Signed_Zeros is True.  An actual implementation might
          -- instead examine the sign bit.  In an implementation in which
          -- T'Signed_Zeros is False, the condition above can be simplified to
          -- Residue < 0.0.
          Residue := -Residue;
          Negated := not Negated;
       end if;
       -- At this point, Residue has a positive sign, and Negated records the
       -- parity of sign flippings that Residue has undergone.
 
       exit when Residue < Abs_Y;
 
       -- At this point, Residue is greater than or equal to Abs_Y.  Its
       -- exponent is the same as, or greater than, that of Abs_Y.

       Reducer :=  Compose (Abs_Y,  Exponent(Residue));

       -- Reducer now has the fraction of Abs_Y and the exponent of Residue.
       -- Thus, it is a (possibly large) exact multiple of Abs_Y.
 
       if Reducer > Residue then
 
          -- Reducer is greater than Residue only when
          -- T'Fraction(Abs_Y) is greater than T'Fraction(Residue).
          -- Reduce its exponent by one.
          Reducer := Scaling(Reducer, -1);
          -- It can be proved that underflow cannot occur in the above scaling.
 
          -- At this point, 1.0 < Residue/Reducer < e_Real(T'Machine_Radix).
          N := Unbiased_Rounding (Residue / Reducer);
          -- Thus, 1.0 <= N <= e_Real (Machine_Radix).
 
          -- Now basically want to subtract N*Reducer from Residue exactly,
          -- but the product may have one too many digits to be represented
          -- exactly.  That occurs when the exponent of N*Reducer exceeds that
          -- of Reducer; in the present case, that can happen for N as small as
          -- two.
 
          -- The following almost works:
          --    Reducer_Head := T'Leading_Part(Reducer, 1);
          --    Reducer_Tail := Reducer - Reducer_Head;
          --    Residue := (Residue - N*Reducer_Head) - N*Reducer_Tail;
          -- It fails only when Reducer is so small that underflow occurs when
          -- subtracting Reducer_Head from it.  Note that this is only a problem
          -- when T'Denorm is False; when T'Denorm is True, the above suffices.
 
          if Reducer < CONST1 then
             -- Reducer is near the underflow threshold, and of course Residue
             -- is near Reducer.  Scale both of them up by a sufficient amount
             -- to prevent underflow when subtracting Reducer_Head from Reducer;
             -- scale back down later.
             Residue := Scaling(Residue, e_Real_Machine_Mantissa);
             Reducer := Scaling(Reducer, e_Real_Machine_Mantissa);
             Scaled_Up := True;
          else
             Scaled_Up := False;
          end if;
 
          Reducer_Head := Leading_Part (Reducer, 1);
          Reducer_Tail := Reducer - Reducer_Head;
          -- Because cancellation occurs in the above subtraction, the result is
          -- exact.
 
          -- Now the subtraction can be performed in two stages.
          Residue := (Residue - N*Reducer_Head) - N*Reducer_Tail;
          -- In the present case, although N*Reducer can have too many digits to
          -- be representable, it cannot overflow.
 
          if Scaled_Up then
             -- Scale back down.  Note that underflow can occur in rare
             -- circumstances here (i.e., when  T'Denorm is False and the
             -- remainder is less than the underflow threshold, which requires
             -- that Y be near the underflow threshold and X be near a multiple
             -- of Y).  The specification calls for zero to be returned, but
             -- T'Scaling might not return a zero when it underflows.  If it
             -- does, and the zero is properly signed, the if-then-else below
             -- can be replaced by the else part (or by the equivalent
             -- multiplication or division, if it yields a properly signed
             -- zero on underflow).
             if Abs (Residue) < CONST2 then
                Residue := Copy_Sign (Zero, Residue);
             else
                Residue := Scaling (Residue, -E_Real_Machine_Mantissa);
             end if;
          end if;
 
       else
 
          -- This case is for Reducer <= Residue.
 
          -- At this point, 1.0 <= Residue/Reducer < e_Real(T'Machine_Radix).

          N := Unbiased_Rounding (Residue / Reducer);

          -- Thus, 1.0 <= N <= e_Real(T'Machine_Radix).
 
          -- Here, the technique for subtracting N*Reducer exactly from Residue
          -- is different.  In the present case, N*Reducer may have one too many
          -- digits to be represented exactly only when the rounding was upward,
          -- hence (N-1.0)*Reducer must necessarily be representable.  Also,
          -- N*Reducer could even overflow, but (N-1.0)*Reducer cannot.
 
          if N > One then
             -- The optimization represented by the above test is probably
             -- worthwhile.
             Residue := Residue - (N - One) * Reducer;
          end if;

          Residue := Residue - Reducer;
          -- The above subtraction can underflow when T'Denorm is False, in
          -- which case the desired result is zero.  It is assumed that when
          -- subtraction underflows, it underflows to zero.
 
       end if;
 
       -- Residue may now be negative, but its absolute value is less than or
       -- equal to half of Reducer.
 
    end loop;
 
    -- At this point, Residue has a positive sign and a magnitude less than that
    -- of Abs_Y.  If Residue is greater than half of Abs_Y, correct it by
    -- subtracting Abs_Y one more time.  We do this without computing half of
    -- Abs_Y, which could underflow or be inexact.
 
    Temp := Residue - Abs_Y;
    -- The above can underflow.  It is assumed here that underflow produces a
    -- zero result.  Note that Temp now has a negative sign (a zero produced on
    -- underflow is presumably a negative zero when T'Signed_Zeros is True).

    if Temp + Residue > Zero then
       -- True only if Residue is greater than half of Abs_Y, or if the
       -- computation of Temp underflowed to zero.  Note that the condition
       -- might, on some machines, be more efficiently evaluated as
       -- -Temp < Residue, or even as abs Temp < Residue.
       Residue := Temp;
    end if;
 
    -- The above step might even be slightly more efficiently evaluated as
    -- follows (here Temp is the negation of the value computed above and 
    -- hence always has a positive sign):
    --    Temp := Abs_Y - Residue;
    --    if Temp < Residue then
    --       Residue := -Temp;
    --    end if;
    -- This version, which is clearly equivalent but harder to motivate, is
    -- used in the binary case at the end of this LSN.
 
    -- The desired remainder is now Residue, with a possible sign flip 
    -- (i.e., if Negated is True at this point).
 
    if Negated then
       return -Residue;
    else
       return Residue;
    end if;
 
   end Remainder;


-- SECTION III.
--
-- Routines for conversion from Real to e_Real and back again.


   No_Of_Usable_Bits_In_Real : constant := 52;

   pragma Assert (Real'Machine_Mantissa >= No_Of_Usable_Bits_In_Real);

   No_Of_e_Digits_Per_Real : constant
                   := (No_Of_Usable_Bits_In_Real-1) / No_Of_Bits_In_Radix + 1;
   --  Used only by:  Make_Extended to convert a Real object into an e_Real.
   --  Equals:  Ceiling (No_Of_Usable_Bits_In_Digit / No_Of_Bits_In_Radix)


   -------------------
   -- Make_Extended --
   -------------------

   function Make_Extended
     (X : Real)
      return e_Real
   is
      Result         : e_Real; -- Initialized to zero (important).
      Abs_X          : Real          := Abs (X);
      Exponent       : e_Integer     := 0;
   begin

      if not X'Valid then
         raise Constraint_Error with "Failure in routine Make_Real: Input is inf or NaN.";
      end if;

      if X = 0.0 then
         return Zero;
      elsif X = -0.0 then
         return Zero;
      else
         Result.Is_Zero := False;
      end if;

      -- Abs_X = Abs(X):

      if Abs_X > X then
         Result.Is_Positive := False;
      else
         Result.Is_Positive := True;
      end if;

      -- Get power of 2 exponent and 1st digit.  This is not usually in an
      -- inner loop, so do it the brute force way.  If Abs(X) < 1.0
      -- then keep multiplying by Radix until 1.0 <= X <= Radix-1.
      -- Strip off the fraction to get the first digit.

      if Abs_X < Real_One then  -- mult. by Radix until >= 1.0:

         for I in e_Integer loop
            Abs_X    := Abs_X * Real_Radix;
            Exponent := Exponent - 1;
            if Abs_X >= Real_One then
               exit;
            end if;
         end loop;
         -- Abs_X is now 1.0 <= Abs_X < Radix.  When strip off the fractional
         -- part with Real_Floor, then Abs_X will be in 1.0..Radix-1.

      elsif Abs_X >= Real_Radix then -- divide by Radix until just right:

         for I in e_Integer loop
            Abs_X    := Abs_X * Inverse_Radix;
            Exponent := Exponent + 1;
            if Abs_X < Real_Radix then
               exit;
            end if;
         end loop;
         -- Abs_X is now 1.0 <= Abs_X < Radix.  When strip off the fractional
         -- part with Real_Floor, then Abs_X will be in 1.0..Radix-1.

      else -- Abs_X is in desired range:

         Exponent := 0;

      end if;

      Result.Exp := Exponent;

      -- Now we've got Result.Exp, Result.Is_Positive, Result.Is_Zero all set.
      -- Is_Infinite is initialized to False.  Next get the first digit:

      Result.Digit(0) := Digit_Floor (Abs_X);
      Abs_X           := Abs_X - Real (Result.Digit(0));

      -- Now just the Abs_X < 1.0 fraction remains in Abs_X.
      -- Optimization: if Abs_X = 0.0 then return early.  Result is
      -- already initialized to zero...no need to get next digits.

      if Abs_X = 0.0 then
         return Result;
      end if;

      -- Get subsequent digits.  These digits are in the range
      -- 0.0 <= X <= Radix-1. (Run the loop longer by 1 for safety.)

      for I in Digit_Index range Digit_Index'First+1..No_Of_e_Digits_Per_Real loop
         Abs_X           := Abs_X * Real_Radix;
         Result.Digit(I) := Digit_Floor (Abs_X);
         Abs_X           := Abs_X - Real (Result.Digit(I));
      end loop;

      if Abs_X > Real_One then
         raise Constraint_Error with "Error in Make_Extended.  Probably bad input.";
      end if;

      return Result;

   end Make_Extended;

   ---------
   -- "+" --
   ---------

   --  Only works in range of Real (15 digits usually).
   --  So X = 2**62 raises Constraint_Error if Real'Digits = 15.

   function "+" (X : Integer) return e_Real
   is
      X_Real : constant Real := Real (X);
   begin
      if Abs X_Real > 2.0**(Real'Machine_Mantissa-1) then
         raise Constraint_Error with "Can't make extended. Argument too large.";
      end if;
      return Make_Extended (X_Real);
   end "+";

   ------------------
   -- Make_e_Digit --
   ------------------

   function Make_e_Digit
     (X : Real)
      return e_Digit
   is
      Ext_Result : e_Real;  -- Initialized to zero.
      Result     : e_Digit; -- Initialized to zero.
   begin

      Ext_Result := Make_Extended (X);

      if Ext_Result.Digit(1) /= Digit_Zero or Ext_Result.Digit(2) /= Digit_Zero then
         raise Constraint_Error with "Error in Make_e_Digit: arg not in range.";
      end if;

      if Ext_Result.Is_Zero then
         return Result;
      end if;

      Result.Exp         := Ext_Result.Exp;
      Result.Is_Positive := Ext_Result.Is_Positive;
      Result.Is_Zero     := Ext_Result.Is_Zero;
      Result.Digit       := Ext_Result.Digit(0);

      return Result;

   end Make_e_Digit;

   ---------
   -- "-" --
   ---------

   -- Same as Make_Extended, but changes the sign.

   function "-" (X : Real) return e_Real is
      Z : e_Real := Make_Extended (X);
   begin
      Z.Is_Positive := not Z.Is_Positive;

      if Z.Is_Zero then -- get sign right again
         Z.Is_Positive := True;
      end if;
      return Z;
   end;

   ---------------
   -- Make_Real --
   ---------------

   --  Most of the arithmetic here is between powers of the machine_radix.
   --  Results should be exact out to the last place of Real. But
   --  can't be guaranteed.

   function Make_Real (X : e_Real)
      return Real
   is
      Result           : Real    := Real_Zero;
      Mantissa         : Real    := Real_Zero;
   begin

      if X.Is_Zero then
         return Real_Zero;
      end if;

      if X.Is_Infinite then
         raise Constraint_Error with "Failure in routine Make_Real: Number is infinite.";
      end if;

      -- Here is the general case.  It produces a Mantissa that is a Factor
      -- of Radix larger than the Normalized Fraction that appears in
      -- Value_Of_Real = Normalized Fraction * Radix**Normalized_Exponent.
      --
      --Power_Of_Radix := 1.0;
      --Result         := Real (X.Digit(0));
      --for I in Digit_Index range 1..No_Of_e_Digits_Per_Real-1 loop
         --Power_Of_Radix := Power_Of_Radix * Inverse_Radix;
         --Result         := Result + Power_Of_Radix * Real (X.Digit(I));
      --end loop;
      --
      --  The following is the usual case.  This is the inner product form.
      --  This sometimes gives the best results because it is more often
      --  done in the machine's extended arithmetic, if that's available.
      --  The following produces a Mantissa that is a Factor
      --  of Radix larger than the Normalized_Fraction that appears in
      --  Value_Of_Real = Normalized_Fraction * Radix**Normalized_Exponent.
      --  Recall that X.Exp is less than the Normalized exponents by 1.

      Mantissa := Real (X.Digit(0))
                + Real (X.Digit(1)) * Inverse_Radix
                + Real (X.Digit(2)) * Inverse_Radix_Squared
                + Real (X.Digit(3)) * Inverse_Radix_Squared * Inverse_Radix
                + Real (X.Digit(4)) * Inverse_Radix_Squared * Inverse_Radix_Squared;


      --  Possible overflows are left to the external float package to raise.
      --  underflows to Zero are done explicitly.

      if Integer(X.Exp) * No_Of_Bits_In_Radix < Real'Machine_Emin then
         Result := 0.0;
      else
         Result := Mantissa * Real_Radix**Integer(X.Exp);
      end if;

      -- Here is the Ada94 way:
      --
      --Real_Exponent_Shift := No_Of_Bits_In_Radix * Integer (X.Exp - 1.0);
      --Result := Real_Scaling (Mantissa, Real_Exponent_Shift);

      -- The scaling function multiplies
      -- Mantissa by Real'Machine_Radix ** Real_Exponent_Shift.

      -- At present leave it up to Real floating point to raise the
      -- constraint errors if they exist.

      if X.Is_Positive then
         null;
      else
         Result := -Result;
      end if;

      return Result;

   end Make_Real;


-- SECTION II.
--
-- Standard arithmetic operators.


   ---------
   -- Abs --
   ---------

   function "Abs" (X : e_Real) return e_Real is
     X2 : e_Real := X;
   begin
     X2.Is_Positive := True;
     return X2;
   end "Abs";

   ---------
   -- "-" --
   ---------

   function "-" (X : e_Real) return e_Real is
     X2 : e_Real := X;
   begin

     X2.Is_Positive := not X2.Is_Positive;

     if X2.Is_Zero then
         X2 := Zero;
     end if;
     return X2;
   end "-";

   -------------------
   -- Abs_Is_Lesser --
   -------------------

   -- Is Abs(X) less than Abs(Y)?
   -- This performs the comparison for all digits: 0..Digit_Index'Last.
   -- The user is epected to call "Round" first if he wants the comparison
   -- in the range 0..Last_Correct_Digit  ==  0..Digit_Index'Last-No_Of_Guard_Digits.

   function Abs_Is_Lesser (X, Y : e_Real)
      return Boolean
   is
   begin

      -- Step 0. Handle the infinities. inf < inf raises c.e. but not here.

      if X.Is_Infinite and not Y.Is_Infinite then
         return False; -- |X| > |Y|
      elsif not X.Is_Infinite and Y.Is_Infinite then
         return True;  -- |X| < |Y|
      end if;

      -- Step 0b. Handle the Zeros.  Another case where the Exp does not tell
      -- us the magnitude of the number.

      if X.Is_Zero and not Y.Is_Zero then
         return True;   -- |X| < |Y|
      elsif not X.Is_Zero and Y.Is_Zero then
         return False;  -- |X| > |Y|
      elsif X.Is_Zero and Y.Is_Zero then
         return False;  -- |X| = |Y|
      end if;

      -- Step 1. Find the lesser number, Exponent-wise.  Must have the
      -- number normalized or the following is false.  Also must have filtered
      -- out the special cases in which Exp is unrelated to the size of
      -- numbers: Zero and Infinity.

      if X.Exp < Y.Exp then
         return True;
      elsif X.Exp > Y.Exp then
         return False;
      end if;

      -- Step 2. If got this far, then the Exponents are equal.  Find the
      -- the first unequal digit.  The following makes use of the fact that
      -- the digits are essentially INTEGER values (all zeros beyond the
      -- the decimal point.)

      for I in Digit_Index'First .. Digit_Index'Last loop
         if X.Digit(I) < Y.Digit(I) then
            return True;
         elsif X.Digit(I) > Y.Digit(I) then
            return False;
         end if; -- if got this far, then the digits are equal
                 -- so continue on to the next digit and try again.
      end loop;

      -- If got this far, then the numbers are equal.
      return False;

   end Abs_Is_Lesser;

   --------------------
   -- Abs_Is_Greater --
   --------------------

   -- Is Abs(X) greater than Abs(Y)?

   function Abs_Is_Greater (X, Y : e_Real) return Boolean is
   begin
      -- Step 0. Handle the infinities.

      if X.Is_Infinite and not Y.Is_Infinite then
         return True;  -- |X| > |Y|
      elsif not X.Is_Infinite and Y.Is_Infinite then
         return False; -- |Y| > |X|
      end if;

      -- Step 0b. Handle the Zeros.  Another case where the Exp does not tell
      -- us the magnitude of the number.

      if X.Is_Zero and not Y.Is_Zero then
         return False;   -- |X| < |Y|
      elsif not X.Is_Zero and Y.Is_Zero then
         return True;  -- |X| > |Y|
      elsif X.Is_Zero and Y.Is_Zero then
         return False;  -- |X| = |Y|
      end if;

      -- Step 1b. Find the larger number, Exponent-wise.  Must have the
      -- number normalized or the following is false.

      if X.Exp > Y.Exp then
         return True;
      elsif X.Exp < Y.Exp then
         return False;
      end if;

      -- Step 2. If got this far, then the Exponents are equal.  Find the
      -- the first unequal digit.  The following makes use of the fact that
      -- the digits are essentially INTEGER valued (all zeros beyond the
      -- the decimal point.)

      for I in Digit_Index'First .. Digit_Index'Last loop
         if X.Digit(I) > Y.Digit(I) then
            return True;
         elsif X.Digit(I) < Y.Digit(I) then
            return False;
         end if; -- if got this far, then the digits are equal
                 -- so continue on to the next digit and try again.
      end loop;

      -- If got this far, then the numbers are equal up to Digit_Index'Last.
      return False;

   end Abs_Is_Greater;

   ---------------
   -- Are_Equal --
   ---------------

   --  X equals Y?  Checks all digits except the digits beyond Digit_Index'Last.
   --  No rounding is performed.
   --  Can use this routine recognize Positive_Infinity and Zero.

   function Are_Equal (X, Y : e_Real)
      return Boolean
   is
   begin
      --  both zero, go home early:

      if X.Is_Zero AND Y.Is_Zero then
         return True;
      end if;

      --  one is zero and the other isn't go home early:

      if X.Is_Zero XOR Y.Is_Zero then
         return False;
      end if;

      --  Check signs. Return False if they have different signs.
      -- We've already checked for Zero's.

      if X.Is_Positive XOR Y.Is_Positive then
         return False;
      end if;

      -- both infinite but have different signs, then
      -- the above step already returned false.
      -- Make inf = inf so one can use this functions to recognize inf.
      -- Another reasonable option would be to make it false.

      if X.Is_Infinite AND Y.Is_Infinite then
         return True;
      end if;

      --  One is infinite, the other not:

      if X.Is_Infinite XOR Y.Is_Infinite then
         return False;
      end if;

      -- ANDing and XORing Is_Zero and Is_Infinite now know
      -- that the neither of the numbers is Zero or Infinite.
      -- Check equality, Exponent-wise.  Must have both
      -- numbers normalized or the following doesn't work.  Remember that the
      -- the only unnormalized nums are Zero, and the 2 infinities.  If got
      -- this far then neither X nor Y is one of those three.

      if X.Exp /= Y.Exp then
         return False;
      end if;

      -- got this far, then the Exponents are equal.  Find the
      -- the first unequal digit.  Makes use of the fact that the digits are
      -- essentially integer valued.

      for I in Digit_Index loop
         if X.Digit(I) /= Y.Digit(I) then
            return False;
         end if; -- if got this far, then the digits are equal
                 -- so continue onto the next digit and try again.
      end loop;

      --If got this far, then digits, exponent, and sign are equal.
      return True;

   end Are_Equal;

   -------------------
   -- Are_Not_Equal --
   -------------------

   function Are_Not_Equal (X, Y : e_Real) return Boolean is
   begin
      return NOT Are_Equal (X, Y);
   end Are_Not_Equal;

   ---------
   -- ">" --
   ---------

   -- Is X > Y?
   function ">" (X, Y : e_Real) return Boolean is
   begin


      -- Step 0. Check Zeros.

      if X.Is_Zero AND Y.Is_Zero then
         return False;
      end if;


      -- Step 0b. Need some optimizations for the common case in which one
      -- attempts to determine Positivity by X > Zero or signs by Zero > Y.
      -- The following lets infinities through for the non-zero part.

      if (not X.Is_Zero) AND Y.Is_Zero then
         if X.Is_Positive then
            return True;   -- X is pos. but not zero, then (X > Zero).
         else
            return False;  -- X is neg. but not zero, then not (X > Zero).
         end if;
      end if;

      if X.Is_Zero AND (not Y.Is_Zero) then
         if Y.Is_Positive then
            return False;  -- Y is pos. but not zero, then not (Zero > Y).
         else
            return True;   -- Y is neg. but not zero, then (Zero > Y).
         end if;
      end if;

      -- Step 1.  Now do things more systematically.
      -- Check signs.  Notice that these give us efficient way to
      -- check sign of a number.  If X is negative, this is fast because
      -- Zero is classified as positive.  So go home early if:

      if X.Is_Positive and not Y.Is_Positive then
         return True;
      elsif not X.Is_Positive and Y.Is_Positive then
         return False;
      end if;

      -- Step 1b. Now they are either both positive or both negative.
      -- If they are both inf, then raise ce, since don't know:

      if X.Is_Infinite AND Y.Is_Infinite then
         raise Constraint_Error with "Constraint_Error in routine >. Arguments are inf.";
      end if;

      -- Step 2. Now they are either both positive or both negative.
      -- If they are both neg. return true if Abs X < Abs Y.

      if X.Is_Positive and Y.Is_Positive then
         return Abs_Is_Greater (X, Y);  -- Abs X > Abs Y
      else
         return Abs_Is_Lesser (X, Y);   -- Abs X < Abs Y
      end if;

   end ">";


   ---------
   -- "<" --
   ---------

   function "<" (X, Y : e_Real) return Boolean is
   begin


      -- Step 0. Check Zeros.

      if X.Is_Zero AND Y.Is_Zero then
         return False;
      end if;


      -- Step 0b. Need some optimizations for the common case in which one
      -- attempts to determine signs by X < Zero or positivity by Zero < Y.
      -- The following lets infinities through for the non-zero part.

      if (not X.Is_Zero) AND Y.Is_Zero then
         if X.Is_Positive then
            return False; -- X is pos. but not zero, then not (X < Zero).
         else
            return True;  -- X is neg. but not zero, then (X < Zero).
         end if;
      end if;

      if X.Is_Zero AND (not Y.Is_Zero) then
         if Y.Is_Positive then
            return True;  -- Y is pos. but not zero, then (Zero < Y).
         else
            return False; -- Y is neg. but not zero, then not (Zero < Y).
         end if;
      end if;


      -- Step 1. Now do things more sytematically.
      -- Check signs.  Notice that these give us efficient way to
      -- check sign of a number.  If X is negative, this is fast because
      -- Zero is classified as positive.  (If want to find if it's Pos., use
      -- not (X < Zero). Since they aren't both 0 go home early if:

      if X.Is_Positive and not Y.Is_Positive then
         return False;
      elsif not X.Is_Positive and Y.Is_Positive then
         return True;
      end if;


      -- Step 1b. Now they are either both positive or both negative.
      -- If they are both inf, then raise ce:

      if X.Is_Infinite AND Y.Is_Infinite then
         raise Constraint_Error with "Error in routine <. Arguments are inf.";
      end if;


      -- Step 2. Now they are either both positive or both negative.
      -- If they are both neg. return true if Abs X > Abs Y.

      if X.Is_Positive and Y.Is_Positive then
         return Abs_Is_Lesser (X, Y);  -- Abs X < Abs Y
      else
         return Abs_Is_Greater (X, Y); -- Abs X > Abs Y
      end if;

   end "<";

   ----------
   -- ">=" --
   ----------

   function ">="
     (X, Y : e_Real)
      return Boolean
   is
   begin
      return (not (X < Y));
   end ">=";

   ----------
   -- "<=" --
   ----------

   function "<="
     (X, Y : e_Real)
      return Boolean
   is
   begin
      return (not (X > Y));
   end "<=";

   ----------
   --  Add --
   ----------

   -- Add the numbers.  The individual digits may overflow the 0..Radix-1 range
   -- but not the range of the base floating point number used to represent the
   -- digit.  Carrying is done later.

   function Add
     (Larger      : e_Real;
      Smaller     : e_Real;
      Digit_Shift : Digit_Index)
      return e_Real
   is
      Z : e_Real;
   begin
      if Digit_Shift > Digit_Index'First then
         for I in Digit_Index'First .. Digit_Shift-1 loop
            Z.Digit(I) := Larger.Digit(I);
         end loop;
      end if;
      for I in Digit_Shift .. Digit_Index'Last loop
         Z.Digit(I) := Larger.Digit(I) + Smaller.Digit(I - Digit_Shift);
      end loop;

      return Z;

   end Add;

   --pragma Inline (Add);

   --------------
   -- Subtract --
   --------------

   -- Subtract the smaller from the larger.  If they have the same
   -- exponent, (ie Digit_Shift = 0),
   -- then use the quantity "First_Unequal_Digit" to optimize
   -- the subtraction, by inserting zeros for all of the equal digits.
   -- We already verified that Larger and Smaller are not equal.

   procedure Subtract
     (Larger              : in  e_Real;
      Smaller             : in  e_Real;
      Digit_Shift         : in  Digit_Index;
      First_Unequal_Digit : in  Digit_Index;
      Result              : out e_Real;
      Extra_Guard_Digit   : out Digit_Type)
   is
      I : Digits_Base;
   begin

      if Digit_Shift = 0 then  -- We can make use of First_Unequal_Digit.

         if First_Unequal_Digit > 0 then
            for I in 0 .. First_Unequal_Digit-1 loop
               Result.Digit(I) := Digit_Zero;
            end loop;
         end if;
         for I in First_Unequal_Digit .. Digit_Index'Last loop
            Result.Digit(I) := Larger.Digit(I) - Smaller.Digit(I);
         end loop;

         Extra_Guard_Digit := Digit_Zero;  -- important initialization.

      else

         for I in 0 .. Digit_Shift-1 loop
            Result.Digit(I) := Larger.Digit(I);
         end loop;
         for I in Digit_Shift .. Digit_Index'Last loop
            Result.Digit(I) := Larger.Digit(I) - Smaller.Digit(I - Digit_Shift);
         end loop;

         I := Digit_Index'Last + 1;
         Extra_Guard_Digit := -Smaller.Digit(I - Digit_Shift);
         --  Here the Larger.Digit(I) = 0.0 because it ran out of digits.

      end if;

   end Subtract;

   --pragma Inline (Subtract);

   ----------------------------
   -- Borrow_For_Subtraction --
   ----------------------------

   -- Do the Borrowing.  This can have much overhead in "-" or "+", so some
   -- optimizations are performed.  If the digits have become less than zero
   -- then must borrow from the next higher order digit: subtract 1 from that
   -- digit, and add Radix to digit in question.  Start
   -- at the least significant digit, Digit_Index'Last, work up to point at which
   -- the subtraction began: Digit_Shift.  Digit_Shift is the first digit of
   -- Z on which a subtraction was performed.  After that, unless a borrow is
   -- performed, the process may end.  Borrowing should be extremely
   -- rare, so don't do unless necessary.
   -- The Extra_Guard_Digit is a virtual digit at index Digit_Index'Last+1.  Very
   -- important in improving precision in some cases: particularly subtracting
   -- a small number from 1.0.

   procedure Borrow_For_Subtraction
     (Z                 : in out e_Real;
      Extra_Guard_Digit : in     Digit_Type;
      Digit_Shift       : in     Digit_Index)
   is
      First_Nonzero_Digit     : Digits_Base := Digit_Index'Last+1;
      All_The_Digits_Are_Zero : Boolean     := True;

      Borrow      : constant Digit_Type := -Digit_One;
      Guard_Digit : Digit_Type := Extra_Guard_Digit;
      I           : Digits_Base;
   begin


      -- This is the general case in which many numbers were subtracted:
      -- Here it is possible that Borrow > 1.0.
      -- if Digit_Shift < Digit_Index'Last then
      -- for I in reverse Digit_Shift+1..Digit_Index'Last loop
      --   if Z.Digit(I) < 0.0 then
      --      Borrow       := Real_Floor (Z.Digit(I) * Inverse_Radix);
      --      Z.Digit(I)   := Z.Digit(I)   - Borrow * Radix; -- Borrow is < 0.0.
      --      Z.Digit(I-1) := Z.Digit(I-1) + Borrow;
      --   end if;
      -- end loop;
      -- end if;


      -- We are subtracting only 2 numbers, so borrow at most 1 digit.

      --  Special case: the extra guard digit at Digit_Index'Last+1:
    --Borrow := -Digit_One;
      I := Digit_Index'Last+1;
      if Guard_Digit < Digit_Zero then
         Guard_Digit  := Guard_Digit  + Digit_Radix;
         Z.Digit(I-1) := Z.Digit(I-1) + Borrow;
      end if;

      if Digit_Shift < Digit_Index'Last then
    --Borrow := -Digit_One;
      for I in reverse Digit_Shift+1 .. Digit_Index'Last loop
         if Z.Digit(I) < Digit_Zero then
            Z.Digit(I)   := Z.Digit(I)   + Digit_Radix;
            Z.Digit(I-1) := Z.Digit(I-1) + Borrow;
         end if;
      end loop;
      end if;


      -- Step 1. Do everything between the 2nd digit and Digit_Shift.
      -- If no borrowing is performed, then are done, since these are the
      -- digits on which no subtractions were performed initially. (With the
      -- exception of digit Digit_Shift: still must check that it
      -- is not < 0.0.)

      -- The general case:
      -- Borrow_Loop:
      -- for I in reverse Digit_Index'First+1 .. Digit_Shift loop
      --    if Z.Digit(I) < 0.0 then
      --      Borrow       := Real_Floor (Z.Digit(I) * Inverse_Radix);
      --      Z.Digit(I)   := Z.Digit(I)   - Borrow * Radix;
      --      Z.Digit(I-1) := Z.Digit(I-1) + Borrow;
      --    else
      --       exit Borrow_Loop;
      --    end if;
      -- end loop Borrow_Loop;


      -- We are subtracting only 2 numbers, so borrow at most 1 digit.

    --Borrow := -Digit_One;
      Borrow_Loop:
      for I in reverse Digit_Index'First+1..Digit_Shift loop
         if Z.Digit(I) < Digit_Zero then
            Z.Digit(I)   := Z.Digit(I)   + Digit_Radix;
            Z.Digit(I-1) := Z.Digit(I-1) + Borrow;
         else
            exit Borrow_Loop;
         end if;
      end loop Borrow_Loop;


      -- Step 2. If Z.Digit(0) < 0.0 then the result is < 0.0, which means
      -- a failure in the "+" or "-" routines below.

      if not Disable_Program_Error_Tests then
      if Z.Digit(0) < Digit_Zero then
          raise Program_Error with "Some error in Borrow_For_Subtraction.";
      end if;
      end if;


      -- Step 3. Normalize the result if the highest order Digit is zero.
      -- Shift the exponent accordingly.  Recall that Z should not be Zero;
      -- checked for that possibility before subtracting.
      -- So shift the entire mantissa left by the number of leading zeros,
      -- and decrement the exponent by the same amount.  If do any left-shifts,
      -- then put at the end of the mantissa the extra guard digit dragged
      -- along just for this event.

      First_Nonzero_Digit     := Digit_Index'Last + 1;
      All_The_Digits_Are_Zero := True;
      for I in Digit_Index loop
         if Z.Digit(I) /= Digit_Zero then
             First_Nonzero_Digit     := I;
             All_The_Digits_Are_Zero := False;
             exit;
         end if;
      end loop;

      if All_The_Digits_Are_Zero then
         if Guard_Digit = Digit_Zero then
            --  But we checked equality of X and Y. 
	    --  Only time this happened it was a compiler bug.
            -- but maybe not an err?
            Z := Zero;
            if not Disable_Program_Error_Tests then
               raise Program_Error with "Might be a bug in Borrow_For_Subtraction.";
            end if;
         else
         --  This is certainly possible:
            Z.Digit(0) := Guard_Digit;
            Z.Exp      := Z.Exp - e_Integer (First_Nonzero_Digit);
         end if;
      end if;

      if not All_The_Digits_Are_Zero then -- First_Nonzero_Digit < Max_Index+1

      if First_Nonzero_Digit > 0 then  -- shift the mantissa left by this amount

         for I in 0 .. Digit_Index'Last-First_Nonzero_Digit loop
            Z.Digit(I) := Z.Digit(I + First_Nonzero_Digit);
         end loop;
         --  Shift the mantissa left by this amount.

         for I in Digit_Index'Last-First_Nonzero_Digit+1 .. Digit_Index'Last loop
            Z.Digit(I) := Digit_Zero;
         end loop;
         --  Set the rest of the mantissa to 0.0.

         Z.Digit(Digit_Index'Last-First_Nonzero_Digit+1) := Guard_Digit;
         --  Even tho' set Digit to 0.0 above, set it right now.

         Z.Exp := Z.Exp - e_Integer (First_Nonzero_Digit);

      end if;

      end if;

      -- If First_Nonzero_Digit = 0, the usual case, then are done.

   end Borrow_For_Subtraction;

   --pragma Inline (Borrow_For_Subtraction);

   ------------------------
   -- Carry_For_Addition --
   ------------------------

   -- Do the carrying.  This can have much overhead, so some
   -- optimizations are performed.  If the digits have become larger
   -- than Radix-1 then must break the digit into 2 parts and add the larger
   -- to a higher order digit. This carry distance is at most one digit,
   -- rather than the 2 possible in the multiplication routine.  Start
   -- at the least significant digit, Digit_Index'Last, work up to point at which
   -- the addition began: Digit_Shift.  Digit_Shift is the first digit of
   -- Z on which an addition was performed.  After that, unless a carry is
   -- performed, the process may be ended.  Carrying should be extremely
   -- rare, so don't do unless necessary.

   procedure Carry_For_Addition
     (Z           : in out e_Real;
      Digit_Shift : in     Digit_Index)
   is
      Digit_Minus_1   : Digit_Type := Digit_Zero;
      We_Are_Finished : Boolean    := False;
      Must_Normalize  : Boolean    := False;
      Carry           : constant Digit_Type := Digit_One;
   begin
      -- Step 1.  Do the carrying among the digits that have been added to each
      -- other (Digit_ID in Digit_Shift+1..Digit_Index'Last).  Actually,
      -- Digit_ID = Digit_Shift had an addition performed to it also..that's
      -- dealt with in step 2.
      --
      -- This is the general case.  Useful for an optimized Sum(many numbers).
      -- if Digit_Shift < Digit_Index'Last then
      -- for I in reverse Digit_Shift+1..Digit_Index'Last loop
      --   if Z.Digit(I) > Radix_Minus_1 then
      --      Carry        := Real_Floor (Z.Digit(I) * Inverse_Radix);
      --      Z.Digit(I)   := Z.Digit(I)   - Carry * Radix;
      --      Z.Digit(I-1) := Z.Digit(I-1) + Carry;
      --   end if;
      -- end loop;
      -- end if;
      --
      -- We're summing at most 2 numbers, so Carry is at most 1.0.

      if Digit_Shift < Digit_Index'Last then
      for I in reverse Digit_Shift+1 .. Digit_Index'Last loop
         if Z.Digit(I) > Digit_Radix_Minus_1 then
          --Carry        := Digit_One;
            Z.Digit(I)   := Z.Digit(I)   - Digit_Radix;
            Z.Digit(I-1) := Z.Digit(I-1) + Carry;
         end if;
      end loop;
      end if;
      We_Are_Finished := False; -- We have at least Digit(Digit_Shift) to check.


      -- Step 2. Do everything between the 2nd digit and Digit_Shift.
      -- If no carry is performed, then we're done, since these are the
      -- digits on which no additions were performed initially. (With the
      -- exception of digit Digit_Shift: still must check that it
      -- is not larger than Radix-1.)
      --
      -- Carry_Loop:
      -- for I in reverse Digit_Index'First+1 .. Digit_Shift loop
      --    if Z.Digit(I) > Radix_Minus_1 then
      --       Carry        := Real_Floor (Z.Digit(I) * Inverse_Radix);
      --       Z.Digit(I)   := Z.Digit(I) - Carry * Radix;
      --       Z.Digit(I-1) := Z.Digit(I-1) + Carry;
      --    else
      --       We_Are_Finished := True;
      --       exit Carry_Loop;
      --    end if;
      -- end loop Carry_Loop;

      -- When summing at most 2 numbers:

      Carry_Loop:
      for I in reverse Digit_Index'First+1..Digit_Shift loop
         if Z.Digit(I) > Digit_Radix_Minus_1 then
          --Carry        := Digit_One;
            Z.Digit(I)   := Z.Digit(I)   - Digit_Radix;
            Z.Digit(I-1) := Z.Digit(I-1) + Carry;
         else
            We_Are_Finished := True;
            exit Carry_Loop;
         end if;
      end loop Carry_Loop;

      -- Step 3. If left the carry_loop early, then go home now.
      -- No need to normalize the result (i.e. make sure that the first
      -- digit is not 0). No need to increment Z.Exp.  This should be the usual
      -- case.  First however, a debugging test:

      if We_Are_Finished then
        if not Disable_Program_Error_Tests then
        if Z.Digit(0) <= Digit_Zero then
           raise Program_Error with "Some error in Carrying for + operator.";
        end if;
        end if;
        return;
     end if;

     -- Step 4. Perform the final carry if Z.Digit(0) > Radix-1 (to a digit
     -- that doesn't exist yet, called Digit_Minus_1.)

     -- Must_Normalize := False;
     -- if Z.Digit(0) > Radix_Minus_1 then
     --    Carry          := Real_Floor (Z.Digit(0) * Inverse_Radix);
     --    Z.Digit(0)     := Z.Digit(0) - Carry * Radix;
     --    Digit_Minus_1  := Carry;
     --    Must_Normalize := True;
     -- end if;

     Must_Normalize := False;
     if Z.Digit(0) > Digit_Radix_Minus_1 then
       --Carry          := Digit_One;
         Z.Digit(0)     := Z.Digit(0) - Digit_Radix;
         Digit_Minus_1  := Carry;  -- Digit_Minus_1 is initially 0.0.
         Must_Normalize := True;
     end if;

     -- Step 5. Normalize the result if Digit(0) was > Radix-1,
     -- hence a carry occurred to a larger digit.
     -- Is it possible that Digit(0) is 0 and Digit_minus_1 is also 0?
     -- No.  To get Digit(0) to zero, it would have to = Radix..then a
     -- carry to Digit_Minus_1 would make it zero.  But then Digit_Minus_1
     -- would be non-zero.

     if Must_Normalize then

         Z.Exp := Z.Exp + 1;

         for I in reverse Digit_Index'First+1 .. Digit_Index'Last loop
            Z.Digit(I) := Z.Digit(I-1);
         end loop;
         Z.Digit(0) := Digit_Minus_1;

      end if;

      -- Test for failure in algorithm:

      if not Disable_Program_Error_Tests then
      if Z.Digit(0) <= Digit_Zero then
         raise Program_Error with "Some error in Carrying for + operator.";
      end if;
      end if;

   end Carry_For_Addition;

   --pragma Inline (Carry_For_Addition);

   ---------
   -- "+" --
   ---------

   -- Just the identity operator, so you can type  A := +1.23E+02
   -- and the statement will be accepted whether or not A is e_Real.

   function "+" (X : e_Real) return e_Real is
   begin
      return X;
   end "+";

   ---------
   -- "+" --
   ---------

   -- Surprizingly complicated.

   function "+"(X, Y : e_Real) return e_Real is
      Z           : e_Real;
      Delta_Exp   : e_Integer;
      Digit_Shift : Digit_Index := 0;
      First_Unequal_Digit : Digit_Index := 0;

      Extra_Guard_Digit : Digit_Type := Digit_Zero; -- Important init. (for subtraction)

      Final_Sign_Is_Positive : Boolean;
      Mantissas_Are_Equal : Boolean;

      type Max_Info is (X_Is_Max, Y_Is_Max);
      Max_Num_ID  : Max_Info := X_Is_Max;

      type Add_Choice is (Add_Them, Subtract_Y_From_X, Subtract_X_From_Y);
      Add_Code : Add_Choice;
   begin

      -- Step 0. If Either of the numbers is 0.0, then return the other.

      if Y.Is_Zero then
         return X;
      elsif X.Is_Zero then
         return Y;
      end if;

      -- Step 0b. If one is infinite, but not the other, return the infinity,
      -- sign of the inf unchanged. If both are inf, say inf + inf = +inf.
      -- And say inf - inf raises c.e., because don't know what it is.

      if Y.Is_Infinite and not X.Is_Infinite then
         return Y;
      elsif not Y.Is_Infinite and X.Is_Infinite then
         return X;
      end if;

      if Y.Is_Infinite and X.Is_Infinite then
         if not X.Is_Positive and not Y.Is_Positive then
            return Negative_Infinity;
         elsif X.Is_Positive and Y.Is_Positive then
            return Positive_Infinity;
         else
            raise Constraint_Error with "Subtraction of inf by inf is undefined.";
         end if;
      end if;


      -- Step 1. Find the larger number, Exponent-wise, and return it if it is
      -- so much larger than the other that there is no addition to be done.
      -- If they are equal, exponent-wise, then say X is the larger.

      if Y.Exp > X.Exp then
         Max_Num_ID  := Y_Is_Max;
      else
         Max_Num_ID  := X_Is_Max;
      end if;

      Delta_Exp := Abs (X.Exp - Y.Exp);
      if Delta_Exp > e_Integer(Digit_Index'Last) then -- ie, Delta_Exp >= No_Of_Digits
         case Max_Num_ID is
         when X_Is_Max =>
            return X;
         when Y_Is_Max =>
            return Y;
         end case;
      end if;


      -- Step 2. When the exponents are equal, and subtraction is going to be
      -- done (ie, one of the numbers is negative, the other pos., X>0 XOR Y>0),
      -- need more information about which number is smaller.
      -- Now correctly find the larger (Abs-wise) even if the Exp's are equal.
      -- Only do this in the subtraction case; if then the numbers turn
      -- out to be equal, then return Zero.  It is important to handle that
      -- special case here.

      First_Unequal_Digit := Digit_Index'First;

      if  (X.Is_Positive XOR Y.Is_Positive) and then X.Exp = Y.Exp  then

        -- Find the first non-equal word in the mantissas:
        Mantissas_Are_Equal := True;
        for I in Digit_Index'First .. Digit_Index'Last loop
           if X.Digit(I) /= Y.Digit(I) then
               Mantissas_Are_Equal := False;
               First_Unequal_Digit := I;
               exit;
           end if;
        end loop;

        -- We're finished if the Exp's are equal, the Mantissas are equal
        -- and we're subtracting one from the other:

        if Mantissas_Are_Equal then
           return Zero;
        end if;

        -- Find the larger of the Two (Absolute values of course):

        if X.Digit(First_Unequal_Digit) > Y.Digit(First_Unequal_Digit) then
           Max_Num_ID  := X_Is_Max;
        else
           Max_Num_ID  := Y_Is_Max;
        end if;

      end if;


      -- Step 3. Do add or subtract?  Depends on their signs.

      if X.Is_Positive and Y.Is_Positive then
         Add_Code := Add_Them;
         Final_Sign_Is_Positive := True;
      end if;

      if (not X.Is_Positive) and (not Y.Is_Positive) then
         Add_Code := Add_Them;    -- add 2 neg nums as tho they were pos.
         Final_Sign_Is_Positive := False;
      end if;

      if (not X.Is_Positive) and Y.Is_Positive then
         case Max_Num_ID is
         when X_Is_Max =>
            Add_Code := Subtract_Y_From_X; -- I mean Abs(X) - Abs(Y)
            Final_Sign_Is_Positive := False;
         when Y_Is_Max =>
            Add_Code := Subtract_X_From_Y;
            Final_Sign_Is_Positive := True;
         end case;
      end if;

      if X.Is_Positive and (not Y.Is_Positive) then
         case Max_Num_ID is
         when X_Is_Max =>
            Add_Code := Subtract_Y_From_X;  -- I mean Abs(X) - Abs(Y)
            Final_Sign_Is_Positive := True;
         when Y_Is_Max =>
            Add_Code := Subtract_X_From_Y;  -- I mean Abs(Y) - Abs(X)
            Final_Sign_Is_Positive := False;
         end case;
      end if;


      -- Step 4. We're now ready to do the adding (or subtracting).
      -- The adding and subtracting are separated from the
      -- carrying/borrowing/normalizing process, because i) the
      -- adding can then be vectorized or otherwise optimized and
      -- ii) in other versions many numbers will be summed (not just X and Y),
      -- in which case it really pays off to do the normalizing
      -- and carrying just once, after many additions, because the
      -- overhead of carrying can be higher than summing.

      Digit_Shift := Digit_Index (Delta_Exp);

      case Add_Code is
      when Add_Them =>

         case Max_Num_ID is
         when X_Is_Max =>
            Z := Add (Larger => X, Smaller => Y, Digit_Shift => Digit_Shift);
         when Y_Is_Max =>
            Z := Add (Larger => Y, Smaller => X, Digit_Shift => Digit_Shift);
         end case;

      when Subtract_Y_From_X =>

         Subtract (Larger => X, Smaller =>  Y, Digit_Shift =>  Digit_Shift,
                   First_Unequal_Digit => First_Unequal_Digit,
                   Result => Z, Extra_Guard_Digit => Extra_Guard_Digit);

      when Subtract_X_From_Y =>

         Subtract (Larger => Y, Smaller =>  X, Digit_Shift =>  Digit_Shift,
                   First_Unequal_Digit => First_Unequal_Digit,
                   Result => Z, Extra_Guard_Digit => Extra_Guard_Digit);

      end case;


      -- Step 5.  Do everything except the carrying or borrowing.
      -- We were careful about subtracting the smaller from the larger, so
      -- the following steps can be done before or after the carrying, (provided
      -- only 2 numbers are being summed).

      if Final_Sign_Is_Positive then
        Z.Is_Positive := True;
      else
        Z.Is_Positive := False;
      end if;

      -- Set Z.Exp but remember it still may be raised by 1 (if carrying occurs),
      -- or lowered by 1 if borrowing occurs.
      case Max_Num_ID is
         when X_Is_Max =>
            Z.Exp := X.Exp;
         when Y_Is_Max =>
            Z.Exp := Y.Exp;
      end case;

      -- The Z = 0 case has already been considered.
      Z.Is_Zero := False;


      -- Do the carrying or borrowing, as the case may be.  These also
      -- normalize the number, and produce an additional correction to the Exp.

      case Add_Code is
      when Add_Them =>

          Carry_For_Addition (Z, Digit_Shift);

      when Subtract_Y_From_X | Subtract_X_From_Y =>

          Borrow_For_Subtraction (Z, Extra_Guard_Digit, Digit_Shift);

      end case;


      -- Step 6. Handle over and under flow.  Here underflow goes to Zero,
      -- overflow to Infinity.  This analysis is all isolated to the end
      -- of the arithmetic routines so that it is more easily modified to
      -- raise exceptions if that is what is desired.   In order to do it
      -- here, must assume that the parmeters Min_Exponent and Max_Exponent
      -- limit the dynamic range of the Exp to about 1/4 of that allowed
      -- by the base type used to represent the exponent.  This is
      -- checked in the spec with an assertion.  (The reason is, the above
      -- code will go outside the accepted range of Exp with out being
      -- checked till down here.)  This limit is OK because the base type
      -- allows excessively large exponents anyway.

      if Z.Exp < Min_Exponent then
         Z := Zero;
      end if;

      if Z.Exp > Max_Exponent then
         if Z.Is_Positive then
            Z := Positive_Infinity;
         else
            Z := Negative_Infinity;
         end if;
      end if;

      return Z;

   end "+";

   ---------
   -- "-" --
   ---------

   -- Subtract Y from X:

   function "-"(X, Y : e_Real) return e_Real is
      Y2 : e_Real := Y;
   begin
      Y2.Is_Positive := not Y.Is_Positive;
      return X + Y2;
   end "-";

   ------------------------------------
   -- Do_Carrying_For_Multiplication --
   ------------------------------------

   -- The digits have probably overflowed their allotted range of 0..Radix-1.
   -- Break the digits into three parts:
   --   digit = d2 * 2**(Radix*2) + d1 * 2**(Radix*1) + d0 * 2**(Radix*0)
   -- where d_n is in the range 0..Radix-1.
   -- Carry d2 to two digits to the left of Digit; carry d1 one digit to the
   -- left of Digit.  So Carry_Minus_1 = d1, Carry_Minus_2 = d2.

   procedure Do_Carrying_For_Multiplication
     (Z             : in out e_Real;
      Digit_Minus_1 : in out Digit_Type;
      Digit_Minus_2 : in out Digit_Type)
   is
      Carry_Minus_1, Carry_Minus_2 : Digit_Type := Digit_Zero; -- Essential init.
      The_Digit : Digit_Type;
      I : Digit_Index;
   begin

      --*******************************************************************
      -- The real valued digits are effectively integers in the range
      -- 0..2**63-1. (Std. vers.) Break them into three words: 1st 30 bits (d0),
      -- 2nd 30 bits (d1), and last 3 bits (d2).  (If Radix is 2**29, then these
      -- will be 29 bit, 29 bit, 5 bit word repectively.)  The 3rd (5 bit) word
      -- will be Carried 2 digits to the left (Carry_Minus_2).  The second
      -- 29 bit word will be carried 1 digit to the left.  The remaining
      -- lowest order 29 bit word will be the desired digit.
      --
      -- Overhead is large from carrying so don't calculate and
      -- carry the 3rd (5 bit) word unless necessary.  This carry's usually
      -- not necessary for the lowest order digits, because usually
      -- the bits in the words are random, so half the time the word is
      -- > Radix/2, half smaller.  Or, <X> = Radix / 2.
      -- Assume <X*Y> = <X><Y> (assume independence).  Then
      -- < SUM(X*Y) > = SUM <(X*Y)> = SUM <X><Y> = SUM <X>**2.  So after
      -- 4 sums (I = 3), you break the Radix**2 barrier, on the average.
      -- So the optimization matters when the total number of digits is
      -- small.  When there are many digits, then it doesn't help much
      -- but the overhead from trying is small enough that might as well
      -- optimize for the common case: relatively small number of total
      -- digits.
      --*******************************************************************

      for I in reverse Digit_Index'First+2 .. Digit_Index'Last loop
         The_Digit := Z.Digit(I);
         if The_Digit >= Digit_Radix_Squared then
            Carry_Minus_2 := Shift_Right_2x_No_of_Bits_in_Radix (The_Digit);
            The_Digit     := The_Digit - Carry_Minus_2 * Digit_Radix_Squared;
            Z.Digit(I-2)  := Z.Digit(I-2) + Carry_Minus_2;
         end if;
         Carry_Minus_1 := Shift_Right_No_of_Bits_in_Radix (The_Digit);
         Z.Digit(I)    := The_Digit - Carry_Minus_1 * Digit_Radix;
         Z.Digit(I-1)  := Z.Digit(I-1) + Carry_Minus_1;
      end loop;

      -- Special case I = Digit_Index'First + 1 = 1.

      I := Digit_Index'First + 1;
      The_Digit := Z.Digit(I);
      if The_Digit >= Digit_Radix_Squared then
         Carry_Minus_2 := Shift_Right_2x_No_of_Bits_in_Radix (The_Digit);
         The_Digit     := The_Digit - Carry_Minus_2 * Digit_Radix_Squared;
         Digit_Minus_1 := Digit_Minus_1 + Carry_Minus_2;
      end if;
      Carry_Minus_1 := Shift_Right_No_of_Bits_in_Radix (The_Digit);
      Z.Digit(I)    := The_Digit - Carry_Minus_1 * Digit_Radix;
      Z.Digit(I-1)  := Z.Digit(I-1) + Carry_Minus_1;

      -- Special case I = Digit_Index'First = 0

      I := Digit_Index'First;
      The_Digit := Z.Digit(I);
      if The_Digit >= Digit_Radix_Squared then
         Carry_Minus_2 := Shift_Right_2x_No_of_Bits_in_Radix (The_Digit);
         The_Digit     := The_Digit - Carry_Minus_2 * Digit_Radix_Squared;
         Digit_Minus_2 := Digit_Minus_2 + Carry_Minus_2;
      end if;
      Carry_Minus_1 := Shift_Right_No_of_Bits_in_Radix (The_Digit);
      Z.Digit(I)    := The_Digit - Carry_Minus_1 * Digit_Radix;
      Digit_Minus_1 := Digit_Minus_1 + Carry_Minus_1;

      -- Special case I = Digit_Index'First - 1

      if Digit_Minus_1 > Digit_Radix_minus_1 then
         Carry_Minus_1 := Shift_Right_No_of_Bits_in_Radix (Digit_Minus_1);
         Digit_Minus_1 := Digit_Minus_1 - Carry_Minus_1 * Digit_Radix;
         Digit_Minus_2 := Digit_Minus_2 + Carry_Minus_1;
      end if;

   end Do_Carrying_For_Multiplication;

   --pragma Inline (Do_Carrying_For_Multiplication);

   ---------------
   -- Normalize --
   ---------------

   -- Normalize the result if the highest order (negative Index) digits are
   -- non-zero.  (The usual case.)
   -- Shift Mantissa and the exponent accordingly.  (The canonical form
   -- requires that the first digit is non-zero.)

   procedure Normalize
     (Z             : in out e_Real;
      Digit_Minus_1 : in     Digit_Type;
      Digit_Minus_2 : in     Digit_Type)
   is
      First_Nonzero_Digit : Digit_Index := Digit_Index'First; -- Init. essential
      First_Nonzero_Digit_Is_Minus_1 : Boolean := False;   -- Init. essential
      First_Nonzero_Digit_Is_Minus_2 : Boolean := False;   -- Init. essential
      All_Digits_Are_Zero : Boolean := False;              -- Init. essential
   begin
     -- Step 0. Infinities and Zero.

     if Z.Is_Infinite then
        return;
     end if;

     Z.Is_Zero := False; --  Will be toggled if all digits 0.

     -- Step 1. Find the first non-zero digit:

     if Digit_Minus_2 > Digit_Zero then
        First_Nonzero_Digit_Is_Minus_2 := True;
     elsif Digit_Minus_1 > Digit_Zero then
        First_Nonzero_Digit_Is_Minus_1 := True;

     else

       All_Digits_Are_Zero := True;
       for I in Digit_Index loop
          if Z.Digit(I) /= Digit_Zero then
             First_Nonzero_Digit := I; -- So First_Nonzero_Digit <= Digit_Index'Last.
             All_Digits_Are_Zero := False;
             exit;
          end if;
       end loop;

     end if;

     if All_Digits_Are_Zero then
        Z := Zero;
        return;
     end if;

     -- Step 2. Shift the array to the right if the Minus_N digits are
     -- non Zero.  Shift the array to the left if necessary (shouldn't be).

     if First_Nonzero_Digit_Is_Minus_2 then    -- Shift right by 2:

        for I in reverse Digit_Index'First+2 .. Digit_Index'Last loop
           Z.Digit(I) := Z.Digit(I - 2);
        end loop;
        Z.Digit(1) := Digit_Minus_1;
        Z.Digit(0) := Digit_Minus_2;
        Z.Exp      := Z.Exp + 2;

     elsif First_Nonzero_Digit_Is_Minus_1 then -- Shift right by 1:

        for I in reverse Digit_Index'First+1 .. Digit_Index'Last loop
           Z.Digit(I) := Z.Digit(I - 1);
        end loop;
        Z.Digit(0) := Digit_Minus_1;
        Z.Exp      := Z.Exp + 1;

     elsif First_Nonzero_Digit > Digit_Index'First then
     --  Shift left by val of First_Non...:

        for I in 0 .. Digit_Index'Last-First_Nonzero_Digit loop
           Z.Digit(I) := Z.Digit(I + First_Nonzero_Digit);
        end loop;
        for I in Digit_Index'Last-First_Nonzero_Digit+1 .. Digit_Index'Last loop
           Z.Digit(I) := Digit_Zero;
        end loop;

        Z.Exp := Z.Exp - e_Integer (First_Nonzero_Digit); --assumes Digit_Index'First=0

     end if;

     -- Test for failure in algorithm:
     if not Disable_Program_Error_Tests then
     if Z.Digit(0) > Digit_Radix_Minus_1 or Z.Digit(0) <= Digit_Zero then
        raise Program_Error with "Some error in Normalization for * operator.";
     end if;
     end if;

   end Normalize;

   --pragma Inline (Normalize);

   ----------------------------
   -- General_Multiplication --
   ----------------------------

   function General_Multiplication
     (X, Y : e_Real)
      return e_Real
   is
      Z : e_Real;
      Digit_Minus_1, Digit_Minus_2 : Digit_Type := Digit_Zero; -- Essential init

      Sum : Digit_Type := Digit_Zero;
      No_Of_Digits, No_Of_Segments, Remaining_Sums : Digits_Base;
      Starting_Digit, Ending_k : Digit_Index;
      Start_Sum, End_Sum       : Digit_Index;
      Allowed_Digits_Per_Carry : constant := Sums_per_Carry + 1;
      --  If you sum 9 numbers you only do 8 sums.
   begin

      -- Step 0. Handle the Zeros. We'll say 0 * infinity is 0, since inf is
      -- really just a large finite number.

      if X.Is_Zero or Y.Is_Zero then
         return Zero;
      end if;

      Z.Is_Zero := False; -- toggled below if underflow.

      -- Step 1. If one or more is infinite..Notice inf * inf = inf here.

      if Y.Is_Infinite or X.Is_Infinite then
         if X.Is_Positive xor Y.Is_Positive then -- opposite signs.
            return Negative_Infinity;
          else
            return Positive_Infinity;
          end if;
      end if;

      -- Step 2. Handle the signs and exponents:

      Z.Is_Positive := not (X.Is_Positive XOR Y.Is_Positive);

      Z.Exp := X.Exp + Y.Exp; -- Will be further adjusted by Carry/Normalize.


      No_Of_Digits   := Digit_Index'Last + 1;
      No_Of_Segments := No_Of_Digits  /  Allowed_Digits_Per_Carry;
      Remaining_Sums := No_Of_Digits REM Allowed_Digits_Per_Carry;


      --  Z has been initialized to Zero: essential if Remaining_Sums = 0

      -- First do the stragglers, digits of index (k in 0..Remaining_Sums-1):

      if  Remaining_Sums > 0 then
      for Digit_ID in Digit_Index loop

         Ending_k := Digit_Index'Min (Digit_ID, Remaining_Sums-1);

         Sum := Digit_Zero;
         for k in Digit_Index'First .. Ending_k loop
            Sum := Sum + X.Digit(k) * Y.Digit(Digit_ID - k);
         end loop;

         Z.Digit(Digit_ID) := Sum; -- init Z.

      end loop;

      Do_Carrying_For_Multiplication (Z, Digit_Minus_1, Digit_Minus_2);
      end if;

      -- Now do the segments of length (up to) (usually) 32:

      if No_Of_Segments > 0 then
      for Segment in 0 .. No_Of_Segments-1 loop

        Start_Sum := (Segment+0) * Allowed_Digits_Per_Carry + Remaining_Sums;
        End_Sum   := (Segment+1) * Allowed_Digits_Per_Carry + Remaining_Sums - 1;
      --End_Sum   := Start_Sum + (Allowed_Digits_Per_Carry - 1);
        Starting_Digit := Start_Sum;

        for Digit_ID in Starting_Digit .. Digit_Index'Last loop

           Ending_k := Digit_Index'Min (Digit_ID, End_Sum);

           Sum := Digit_Zero;
           for k in Start_Sum .. Ending_k loop
              Sum := Sum + X.Digit(k) * Y.Digit(Digit_ID - k);
           end loop;

           Z.Digit(Digit_ID) := Z.Digit(Digit_ID) + Sum;
           -- Z.Digit(Digit_ID) is close enough to 0 (ie, < Radix-1) that this
           -- does not count as a sum in the Sums_Per_Carry rule.
           -- That's why   Allowed_Sums_Per_Carry = Sums_Per_Carry+1   here.

         end loop;

         Do_Carrying_For_Multiplication (Z, Digit_Minus_1, Digit_Minus_2);

      end loop;
      end if;

      -- Must Normalize: shift digit array to make sure that the first digit
      -- is non-zero: if Infinity or Zero gets this far, then problem occurs.
      -- Should catch in normalize.

      Normalize (Z, Digit_Minus_1, Digit_Minus_2);

      -- Step 4. Handle over and under flow.  Here underflow goes to Zero,
      -- overflow to Infinity. In order to do it
      -- here, must assume that the parmeters Min_Exponent and Max_Exponent
      -- limit the dynamic range of the Exp to about 1/4 of that allowed
      -- by the base type used to represent the exponent.  This is
      -- checked in the spec with an assertion.  (The reason is, the above
      -- code will go well outside the accepted range of Exp with out being
      -- checked till down here.)  This limit is OK because the base type
      -- allows excessively large exponents anyway, up to 2**31-1.

      if Z.Exp < Min_Exponent then
         Z := Zero;
      end if;

      if Z.Exp > Max_Exponent then
         if Z.Is_Positive then
            Z := Positive_Infinity;
         else
            Z := Negative_Infinity;
         end if;
      end if;

      return Z;

   end General_Multiplication;

   ------------
   -- Square --
   ------------

   -- Suppose Z, and X are composed of n digits 0..n-1.   Then X*X is
   --
   -- Xn-1 := X0*Xn-1 + X1*Xn-2 + .. + Xn-1*X0
   -- ...
   -- X2   := X0*X2   + X1*X1   + X2*X0
   -- X1   := X0*X1   + X1*X0
   -- X0   := X0*X0
   --
   -- Now follows the upper half of the table, which produces words beyond
   -- the precision of the two number X and X that are being multiplied. These
   -- don't calculate.  (Just make N larger if more precision is needed).
   --
   procedure Square
     (X : in out e_Real)
   is
      Digit_Minus_1, Digit_Minus_2 : Digit_Type := Digit_Zero; -- Essential init

      Ultimate_No_of_Digits : constant := Ultimate_Digit + 1;
      --  equals: Digit_Index'Last - Digit_Index'First + 1  =  Mantissa'Length

      -------------------------------------
      -- Product_if_digits_fewer_than_17 --
      -------------------------------------

      procedure Product_if_digits_fewer_than_17
      is
         pragma Assert (Ultimate_No_of_Digits >= 5); --because no if-then for 1st 5.
         pragma Assert (Ultimate_No_of_Digits < 17);
         pragma Assert (No_Of_Bits_In_Radix <= 30);
	 --  (not a) or b  is same as a implies b.
	 --  no need to carry until the end if following evaluate to true:
         pragma Assert (not (No_Of_Bits_In_Radix = 30) or Ultimate_No_of_Digits <= 8);
         pragma Assert (not (No_Of_Bits_In_Radix = 29) or Ultimate_No_of_Digits <= 32);
         pragma Suppress (Index_Check);
         A      : Mantissa renames X.Digit;
         S      : constant Digit_Index := Digit_Index'First;
      begin
         --  Ultimate_No_of_Digits is named number, so optimizer should
	 --  eliminate unused blocks of code below...not that it matters.

         --  need to do intermediate carrying if No_Of_Bits_In_Radix >= 30
	 --  and more than 8 digits.

         --  ! Must be done in the following order !

         if Ultimate_No_of_Digits >= 16 then
            A (S+15) :=(A(S+0)*A(S+15) + A(S+1)*A(S+14) + A(S+2)*A(S+13) +
                        A(S+3)*A(S+12) + A(S+4)*A(S+11) + A(S+5)*A(S+10)
                      + A(S+6)*A(S+9)  + A(S+7)*A(S+8))*Digit_Two;
         end if;

         if Ultimate_No_of_Digits >= 15 then
            A (S+14) :=(A(S+0)*A(S+14) + A(S+1)*A(S+13) + A(S+2)*A(S+12) +
                        A(S+3)*A(S+11) + A(S+4)*A(S+10) + A(S+5)*A(S+9)
                      + A(S+6)*A(S+8))*Digit_Two
		      + A(S+7)*A(S+7);
         end if;

         if Ultimate_No_of_Digits >= 14 then
            A (S+13) :=(A(S+0)*A(S+13) + A(S+1)*A(S+12) + A(S+2)*A(S+11) +
                        A(S+3)*A(S+10) + A(S+4)*A(S+9)  + A(S+5)*A(S+8)
                      + A(S+6)*A(S+7))*Digit_Two;
         end if;

         if Ultimate_No_of_Digits >= 13 then
            A (S+12) :=(A(S+0)*A(S+12) + A(S+1)*A(S+11) + A(S+2)*A(S+10) +
                        A(S+3)*A(S+9)  + A(S+4)*A(S+8)  + A(S+5)*A(S+7))*Digit_Two
                      + A(S+6)*A(S+6);
         end if;

         if Ultimate_No_of_Digits >= 12 then
            A (S+11) :=(A(S+0)*A(S+11) + A(S+1)*A(S+10) + A(S+2)*A(S+9) +
                        A(S+3)*A(S+8)  + A(S+4)*A(S+7)  + A(S+5)*A(S+6))*Digit_Two;
         end if;

         if Ultimate_No_of_Digits >= 11 then
            A (S+10) := (A(S+0)*A(S+10) + A(S+1)*A(S+9) + A(S+2)*A(S+8)
                       + A(S+3)*A(S+7)  + A(S+4)*A(S+6)) * Digit_Two
                       + A(S+5)*A(S+5);
         end if;

         if Ultimate_No_of_Digits >= 10 then
            A (S+9) := (A(S+0)*A(S+9) + A(S+1)*A(S+8) + A(S+2)*A(S+7)
                      + A(S+3)*A(S+6) + A(S+4)*A(S+5)) * Digit_Two;
         end if;

         if Ultimate_No_of_Digits >= 9 then
            A (S+8) := (A(S+0)*A(S+8) + A(S+1)*A(S+7) + A(S+2)*A(S+6)
                      + A(S+3)*A(S+5)) * Digit_Two
	              + A(S+4)*A(S+4);
         end if;

         if Ultimate_No_of_Digits >= 8 then
            A (S+7) := (A(S+0)*A(S+7) + A(S+1)*A(S+6) + A(S+2)*A(S+5)
                      + A(S+3)*A(S+4)) * Digit_Two;
         end if;

         if Ultimate_No_of_Digits >= 7 then
            A (S+6) := (A(S+0)*A(S+6) + A(S+1)*A(S+5) + A(S+2)*A(S+4)) * Digit_Two
                      + A(S+3)*A(S+3);
         end if;

         if Ultimate_No_of_Digits >= 6 then
            A (S+5) := (A(S+0)*A(S+5) + A(S+1)*A(S+4) + A(S+2)*A(S+3)) * Digit_Two;
         end if;

         A (S+4) := A(S+2)*A(S+2) +
                   (A(S+0)*A(S+4) + A(S+1)*A(S+3)) * Digit_Two;

         A (S+3) :=(A(S+0)*A(S+3) + A(S+1)*A(S+2)) * Digit_Two;

         A (S+2) := A(S+1)*A(S+1) +
                    A(S+0)*A(S+2) * Digit_Two;

         A (S+1) := A(S+0)*A(S+1) * Digit_Two;

         A (S+0) := A(S+0)*A(S+0);

      end Product_if_digits_fewer_than_17;

   begin

      -- Step 0. Handle the Zeros. We'll say 0 * infinity is 0, since inf is
      -- really just a large finite number.

      if X.Is_Zero then
         return;
      end if;

      -- Step 2. Do the multiplication.  We can only sum 32 elements of the sum
      -- before the carrys must be done (in one stnd setting).

      -- We use the inner product version (inner loop is an
      -- dot-product.)  To get the outer product version (BLAS routine DAXPY),
      -- interchange the order of the loops.
      --
      -- Here's the idea:
      -- for Digit_ID in Digit_Index loop
      --    Sum := 0.0;
      --    for k in 0..Digit_ID loop
      --       Sum := Sum + X(k) * Y(Digit_ID-k);
      --    end loop;
      --    Z(Digit_ID) := Sum;
      -- end loop;
      --
      -- Break sum into segments of 32 each in index k.  Perform a carry
      -- after each sum of 32 elements.

      if (Ultimate_No_of_Digits < 9  and then No_Of_Bits_In_Radix <= 30) or
         (Ultimate_No_of_Digits < 17 and then No_Of_Bits_In_Radix <= 29)   then

         Product_if_digits_fewer_than_17;

         X.Exp := X.Exp + X.Exp; -- Will be further adjusted by "Normalize.." routine.

         Do_Carrying_For_Multiplication (X, Digit_Minus_1, Digit_Minus_2);

         -- Must Normalize: shift digit array to make sure that the first digit
         -- is non-zero: if Infinity or Zero gets this far, then problem occurs.
         -- Should catch in normalize.

         Normalize (X, Digit_Minus_1, Digit_Minus_2);

         X.Is_Positive := True;

      else

         X := General_Multiplication (X, X);
         return;

      end if;


      -- Step 4. Handle over and under flow.  Here underflow goes to Zero,
      -- overflow to Infinity.  This is all isolated to the end
      -- of the arithmetic routines so that it is easily modified to
      -- raise exceptions if that's what is desired.   In order to do it
      -- here, must assume that the parmeters Min_Exponent and Max_Exponent
      -- limit the dynamic range of the Exp to about 1/4 of that allowed
      -- by the base type used to represent the exponent.  This is
      -- checked in the spec with an assertion.  (The reason is, the above
      -- code will go well outside the accepted range of Exp with out being
      -- checked till down here.)  This limit is OK because the base type
      -- allows excessively large exponents anyway, up to 2**31-1.

      if X.Exp < Min_Exponent then
         X := Zero;
      end if;

      if X.Exp > Max_Exponent then
         X := Positive_Infinity;
      end if;

   end Square;

   -------------------
   -- Multiply_Stnd --
   -------------------

   -- Suppose Z, X, and Y are composed of n digits 0..n-1.   Then X*Y is
   --
   -- Z0   := X0*Y0
   -- Z1   := X0*Y1   + X1*Y0
   -- Z2   := X0*Y2   + X1*Y1   + X2*Y0
   -- ...
   -- Zn-1 := X0*Yn-1 + X1*Yn-2 + .. + Xn-1*Y0
   --
   -- Now follows the upper half of the table, which produces words beyond
   -- the precision of the two number X and Y that are being multiplied. These
   -- don't calculate.  (Just make N larger if more precision is needed).
   -- If n is minimal for efficiency reasons, then it would be more efficient
   -- in some calculations to make use the following instead of increasing the
   -- number of digits to get the required precision.
   --
   -- Zn    := Xn-1*Y1    + ... +  X1*Yn-1
   -- ...
   -- Z2n-2 := Xn-1*Yn-1
   --
   function Multiply_Stnd
     (X, Y : e_Real)
      return e_Real
   is
      Z   : e_Real := X;
      Digit_Minus_1, Digit_Minus_2 : Digit_Type := Digit_Zero; -- Essential init

      Ultimate_No_of_Digits : constant := Ultimate_Digit + 1;
      --  equals: Digit_Index'Last - Digit_Index'First + 1  =  Mantissa'Length

      -------------------------------------
      -- Product_if_digits_fewer_than_17 --
      -------------------------------------

      procedure Product_if_digits_fewer_than_17
      is
         pragma Assert (Ultimate_No_of_Digits >= 5); --because no if-then's for 1st 5.
         pragma Assert (Ultimate_No_of_Digits < 17);
         pragma Assert (No_Of_Bits_In_Radix <= 30);
	 --  (not a) or b  is same as a implies b.
	 --  no need to carry until the end if following evaluate to true:
         pragma Assert (not (No_Of_Bits_In_Radix = 30) or Ultimate_No_of_Digits <= 8);
         pragma Assert (not (No_Of_Bits_In_Radix = 29) or Ultimate_No_of_Digits <= 32);
         pragma Suppress (Index_Check);
         A      : Mantissa renames Z.Digit;
         B      : Mantissa renames Y.Digit;
         S      : constant Digit_Index := Digit_Index'First;
      begin
         --  would need to do intermediate carrying if
	 --  No_Of_Bits_In_Radix >= 30 and Ultimate_No_of_Digits > 5.

         --  Ultimate_No_of_Digits is named number, so optimizer should
	 --  eliminate unused blocks of code below...

         if Ultimate_No_of_Digits >= 16 then
            A(S+15) := A(S+0)*B(S+15) + A(S+1)*B(S+14) + A(S+2)*B(S+13)
                     + A(S+3)*B(S+12) + A(S+4)*B(S+11) + A(S+5)*B(S+10)
                     + A(S+6)*B(S+9)  + A(S+7)*B(S+8)  + A(S+8)*B(S+7)
                     + A(S+9)*B(S+6)  + A(S+10)*B(S+5) + A(S+11)*B(S+4)
                     + A(S+12)*B(S+3) + A(S+13)*B(S+2) + A(S+14)*B(S+1)
                     + A(S+15)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 15 then
            A(S+14) := A(S+0)*B(S+14) + A(S+1)*B(S+13) + A(S+2)*B(S+12)
                     + A(S+3)*B(S+11) + A(S+4)*B(S+10) + A(S+5)*B(S+9)
                     + A(S+6)*B(S+8)  + A(S+7)*B(S+7)  + A(S+8)*B(S+6)
                     + A(S+9)*B(S+5)  + A(S+10)*B(S+4) + A(S+11)*B(S+3)
                     + A(S+12)*B(S+2) + A(S+13)*B(S+1) + A(S+14)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 14 then
            A(S+13) := A(S+0)*B(S+13) + A(S+1)*B(S+12) + A(S+2)*B(S+11)
                     + A(S+3)*B(S+10) + A(S+4)*B(S+9)  + A(S+5)*B(S+8)
                     + A(S+6)*B(S+7)  + A(S+7)*B(S+6)  + A(S+8)*B(S+5)
                     + A(S+9)*B(S+4)  + A(S+10)*B(S+3) + A(S+11)*B(S+2)
                     + A(S+12)*B(S+1) + A(S+13)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 13 then
            A(S+12) := A(S+0)*B(S+12) + A(S+1)*B(S+11) + A(S+2)*B(S+10)
                     + A(S+3)*B(S+9)  + A(S+4)*B(S+8)  + A(S+5)*B(S+7)
                     + A(S+6)*B(S+6)  + A(S+7)*B(S+5)  + A(S+8)*B(S+4)
                     + A(S+9)*B(S+3)  + A(S+10)*B(S+2) + A(S+11)*B(S+1)
                     + A(S+12)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 12 then
            A(S+11) := A(S+0)*B(S+11) + A(S+1)*B(S+10) + A(S+2)*B(S+9)
                     + A(S+3)*B(S+8)  + A(S+4)*B(S+7)  + A(S+5)*B(S+6)
                     + A(S+6)*B(S+5)  + A(S+7)*B(S+4)  + A(S+8)*B(S+3)
                     + A(S+9)*B(S+2)  + A(S+10)*B(S+1) + A(S+11)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 11 then
            A(S+10) := A(S+0)*B(S+10) + A(S+1)*B(S+9) + A(S+2)*B(S+8)
                     + A(S+3)*B(S+7)  + A(S+4)*B(S+6) + A(S+5)*B(S+5)
                     + A(S+6)*B(S+4)  + A(S+7)*B(S+3) + A(S+8)*B(S+2)
                     + A(S+9)*B(S+1)  + A(S+10)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 10 then
            A(S+9) := A(S+0)*B(S+9) + A(S+1)*B(S+8) + A(S+2)*B(S+7)
                    + A(S+3)*B(S+6) + A(S+4)*B(S+5) + A(S+5)*B(S+4)
                    + A(S+6)*B(S+3) + A(S+7)*B(S+2) + A(S+8)*B(S+1)
                    + A(S+9)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 9 then
            A(S+8) := A(S+0)*B(S+8) + A(S+1)*B(S+7) + A(S+2)*B(S+6)
                    + A(S+3)*B(S+5) + A(S+4)*B(S+4) + A(S+5)*B(S+3)
                    + A(S+6)*B(S+2) + A(S+7)*B(S+1) + A(S+8)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 8 then
            A(S+7) := A(S+0)*B(S+7) + A(S+1)*B(S+6) + A(S+2)*B(S+5)
                    + A(S+3)*B(S+4) + A(S+4)*B(S+3) + A(S+5)*B(S+2)
                    + A(S+6)*B(S+1) + A(S+7)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 7 then
            A(S+6) := A(S+0)*B(S+6) + A(S+1)*B(S+5) + A(S+2)*B(S+4)
                    + A(S+3)*B(S+3) + A(S+4)*B(S+2) + A(S+5)*B(S+1)
                    + A(S+6)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 6 then
            A(S+5) := A(S+0)*B(S+5) + A(S+1)*B(S+4) + A(S+2)*B(S+3)
                    + A(S+3)*B(S+2) + A(S+4)*B(S+1) + A(S+5)*B(S+0);
         end if;

         A(S+4) := A(S+0)*B(S+4) + A(S+1)*B(S+3) + A(S+2)*B(S+2)
                 + A(S+3)*B(S+1) + A(S+4)*B(S+0);

         A(S+3) := A(S+0)*B(S+3) + A(S+1)*B(S+2) + A(S+2)*B(S+1)
                 + A(S+3)*B(S+0);

         A(S+2) := A(S+0)*B(S+2) + A(S+1)*B(S+1) + A(S+2)*B(S+0);

         A(S+1) := A(S+0)*B(S+1) + A(S+1)*B(S+0);

         A(S+0) := A(S+0)*B(S+0);

      end Product_if_digits_fewer_than_17;

   begin
      -- Step 0. Handle the Zeros. We'll say 0 * infinity is 0, since inf is
      -- really just a large finite number.

      if X.Is_Zero or Y.Is_Zero then
         return Zero;
      end if;

      Z.Is_Zero := False; -- toggled below if underflow.

      -- Step 1. If one or more is infinite..Notice inf * inf = inf here.

      if Y.Is_Infinite or X.Is_Infinite then
         if X.Is_Positive xor Y.Is_Positive then -- opposite signs.
            return Negative_Infinity;
          else
            return Positive_Infinity;
          end if;
      end if;

      -- Step 2. Handle the signs and exponents:

      Z.Is_Positive := not (X.Is_Positive XOR Y.Is_Positive);

      Z.Exp := X.Exp + Y.Exp; -- Will be further adjusted by "Carry.." routine.


      -- Step 3. Do the multiplication.  Can only sum (say) 32 elements of
      -- the sum before the carrys must be done.

      -- We use the inner product version (inner loop is an
      -- dot-product.)  To get the outer product version (BLAS routine DAXPY),
      -- interchange the order of the loops.
      --
      -- Here's the idea:
      -- for Digit_ID in Digit_Index loop
      --    Sum := 0.0;
      --    for k in 0..Digit_ID loop
      --       Sum := Sum + X(k) * Y(Digit_ID-k);
      --    end loop;
      --    Z(Digit_ID) := Sum;
      -- end loop;
      --
      -- Break sum into segments of 32 each in index k.  Perform a carry
      -- after each sum of 32 elements.

      if (Ultimate_No_of_Digits < 9  and then No_Of_Bits_In_Radix <= 30) or
         (Ultimate_No_of_Digits < 17 and then No_Of_Bits_In_Radix <= 29) then

         Product_if_digits_fewer_than_17; -- Z = Z * Y  (Z has been initialized to X).

         Do_Carrying_For_Multiplication (Z, Digit_Minus_1, Digit_Minus_2);
         --  the procedure requires initialized (0) Digit_Minus_1, ...

      else

         Z := General_Multiplication (X, Y);
         return Z;

      end if;

      -- Must Normalize: shift digit array to make sure that the first digit
      -- is non-zero: if Infinity or Zero gets this far, then problem occurs.
      -- Should catch in normalize.

      Normalize (Z, Digit_Minus_1, Digit_Minus_2);

      -- Step 4. Handle over and under flow.  Here underflow goes to Zero,
      -- overflow to Infinity. In order to do it
      -- here, must assume that the parmeters Min_Exponent and Max_Exponent
      -- limit the dynamic range of the Exp to about 1/4 of that allowed
      -- by the base type used to represent the exponent.  This is
      -- checked in the spec with an assertion.  (The reason is, the above
      -- code will go well outside the accepted range of Exp with out being
      -- checked till down here.)  This limit is OK because the base type
      -- allows excessively large exponents anyway, up to 2**31-1.

      if Z.Exp < Min_Exponent then
         Z := Zero;
      end if;

      if Z.Exp > Max_Exponent then
         if Z.Is_Positive then
            Z := Positive_Infinity;
         else
            Z := Negative_Infinity;
         end if;
      end if;

      return Z;

   end Multiply_Stnd;

   pragma Inline (Multiply_Stnd);

   ----------
   -- Mult --
   ----------

   -- Suppose Z, X, and Y are composed of n digits 0..n-1.   Then X*Y is
   --
   -- Z0   := X0*Y0
   -- Z1   := X0*Y1   + X1*Y0
   -- Z2   := X0*Y2   + X1*Y1   + X2*Y0
   -- ...
   -- Zn-1 := X0*Yn-1 + X1*Yn-2 + .. + Xn-1*Y0
   --
   -- Now follows the upper half of the table, which produces words beyond
   -- the precision of the two number X and Y that are being multiplied. These
   -- don't calculate.  (Just make N larger if more precision is needed).
   -- If n is minimal for efficiency reasons, then it would be more efficient
   -- in some calculations to make use the following instead of increasing the
   -- number of digits to get the required precision.
   --
   -- Zn    := Xn-1*Y1    + ... +  X1*Yn-1
   -- ...
   -- Z2n-2 := Xn-1*Yn-1
   --
   procedure Mult
     (X : in out e_Real;
      Y : in     e_Real)
   is
      Digit_Minus_1, Digit_Minus_2 : Digit_Type := Digit_Zero; -- Essential init

      Ultimate_No_of_Digits : constant := Ultimate_Digit + 1;
      --  equals: Digit_Index'Last - Digit_Index'First + 1  =  Mantissa'Length

      -------------------------------------
      -- Product_if_digits_fewer_than_17 --
      -------------------------------------

      procedure Product_if_digits_fewer_than_17
      is 
         pragma Assert (Ultimate_No_of_Digits >= 5); --because no if-then's for 1st 5.
         pragma Assert (Ultimate_No_of_Digits < 17);
         pragma Assert (No_Of_Bits_In_Radix <= 30);
	 --  (not a) or b  is same as a implies b.
	 --  no need to carry until the end if following evaluate to true:
         pragma Assert (not (No_Of_Bits_In_Radix = 30) or Ultimate_No_of_Digits <= 8);
         pragma Assert (not (No_Of_Bits_In_Radix = 29) or Ultimate_No_of_Digits <= 32);
         pragma Suppress (Index_Check);
         A : Mantissa renames X.Digit;
         B : Mantissa renames Y.Digit;
         S : constant Digit_Index := Digit_Index'First;
	 --Result : Mantissa;
      begin
         --  would need to do intermediate carrying if
	 --  No_Of_Bits_In_Radix >= 30 and Ultimate_No_of_Digits > 5.

         --  Ultimate_No_of_Digits is named number, so optimizer should
	 --  eliminate unused blocks of code below...

         --  ! Must be done in the following order !

         if Ultimate_No_of_Digits >= 16 then
            A (S+15) := A(S+0)*B(S+15) + A(S+1)*B(S+14) + A(S+2)*B(S+13)
                      + A(S+3)*B(S+12) + A(S+4)*B(S+11) + A(S+5)*B(S+10)
                      + A(S+6)*B(S+9)  + A(S+7)*B(S+8)  + A(S+8)*B(S+7)
                      + A(S+9)*B(S+6)  + A(S+10)*B(S+5) + A(S+11)*B(S+4)
	              + A(S+12)*B(S+3) + A(S+13)*B(S+2) + A(S+14)*B(S+1)
                      + A(S+15)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 15 then
            A (S+14) := A(S+0)*B(S+14) + A(S+1)*B(S+13) + A(S+2)*B(S+12)
                      + A(S+3)*B(S+11) + A(S+4)*B(S+10) + A(S+5)*B(S+9)
                      + A(S+6)*B(S+8)  + A(S+7)*B(S+7)  + A(S+8)*B(S+6)
                      + A(S+9)*B(S+5)  + A(S+10)*B(S+4) + A(S+11)*B(S+3)
                      + A(S+12)*B(S+2) + A(S+13)*B(S+1) + A(S+14)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 14 then
            A (S+13) := A(S+0)*B(S+13) + A(S+1)*B(S+12) + A(S+2)*B(S+11)
                      + A(S+3)*B(S+10) + A(S+4)*B(S+9)  + A(S+5)*B(S+8)
                      + A(S+6)*B(S+7)  + A(S+7)*B(S+6)  + A(S+8)*B(S+5)
                      + A(S+9)*B(S+4)  + A(S+10)*B(S+3) + A(S+11)*B(S+2)
                      + A(S+12)*B(S+1) + A(S+13)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 13 then
            A (S+12) := A(S+0)*B(S+12) + A(S+1)*B(S+11) + A(S+2)*B(S+10)
                      + A(S+3)*B(S+9)  + A(S+4)*B(S+8)  + A(S+5)*B(S+7)
                      + A(S+6)*B(S+6)  + A(S+7)*B(S+5)  + A(S+8)*B(S+4)
                      + A(S+9)*B(S+3)  + A(S+10)*B(S+2) + A(S+11)*B(S+1)
                      + A(S+12)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 12 then
            A (S+11) := A(S+0)*B(S+11) + A(S+1)*B(S+10) + A(S+2)*B(S+9)
                      + A(S+3)*B(S+8)  + A(S+4)*B(S+7)  + A(S+5)*B(S+6)
                      + A(S+6)*B(S+5)  + A(S+7)*B(S+4)  + A(S+8)*B(S+3)
                      + A(S+9)*B(S+2)  + A(S+10)*B(S+1) + A(S+11)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 11 then
            A (S+10) := A(S+0)*B(S+10) + A(S+1)*B(S+9) + A(S+2)*B(S+8)
                      + A(S+3)*B(S+7)  + A(S+4)*B(S+6) + A(S+5)*B(S+5)
                      + A(S+6)*B(S+4)  + A(S+7)*B(S+3) + A(S+8)*B(S+2)
                      + A(S+9)*B(S+1)  + A(S+10)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 10 then
            A (S+9) := A(S+0)*B(S+9) + A(S+1)*B(S+8) + A(S+2)*B(S+7)
                     + A(S+3)*B(S+6) + A(S+4)*B(S+5) + A(S+5)*B(S+4)
                     + A(S+6)*B(S+3) + A(S+7)*B(S+2) + A(S+8)*B(S+1)
                     + A(S+9)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 9 then
            A (S+8) := A(S+0)*B(S+8) + A(S+1)*B(S+7) + A(S+2)*B(S+6)
                     + A(S+3)*B(S+5) + A(S+4)*B(S+4) + A(S+5)*B(S+3)
                     + A(S+6)*B(S+2) + A(S+7)*B(S+1) + A(S+8)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 8 then
            A (S+7) := A(S+0)*B(S+7) + A(S+1)*B(S+6) + A(S+2)*B(S+5)
                     + A(S+3)*B(S+4) + A(S+4)*B(S+3) + A(S+5)*B(S+2)
                     + A(S+6)*B(S+1) + A(S+7)*B(S+0);
         end if;

         if Ultimate_No_of_Digits >= 7 then
            A (S+6) := A(S+0)*B(S+6) + A(S+1)*B(S+5) + A(S+2)*B(S+4)
                     + A(S+3)*B(S+3) + A(S+4)*B(S+2) + A(S+5)*B(S+1)
                     + A(S+6)*B(S+0);
         end if;


         if Ultimate_No_of_Digits >= 6 then
            A (S+5) := A(S+0)*B(S+5) + A(S+1)*B(S+4) + A(S+2)*B(S+3)
                     + A(S+3)*B(S+2) + A(S+4)*B(S+1) + A(S+5)*B(S+0);
         end if;


         A (S+4) := A(S+0)*B(S+4) + A(S+1)*B(S+3) + A(S+2)*B(S+2)
                  + A(S+3)*B(S+1) + A(S+4)*B(S+0);

         A (S+3) := A(S+0)*B(S+3) + A(S+1)*B(S+2) + A(S+2)*B(S+1)
                  + A(S+3)*B(S+0);

         A (S+2) := A(S+0)*B(S+2) + A(S+1)*B(S+1) + A(S+2)*B(S+0);

         A (S+1) := A(S+0)*B(S+1) + A(S+1)*B(S+0);

         A (S+0) := A(S+0)*B(S+0);

      end Product_if_digits_fewer_than_17;

   begin

      -- Step 0. Handle the Zeros. We'll say 0 * infinity is 0, since inf is
      -- really just a large finite number.

      if X.Is_Zero or Y.Is_Zero then
         X := Zero;
         return;
      end if;

      -- Step 1. If one or more is infinite..Notice inf * inf = inf here.

      if Y.Is_Infinite or X.Is_Infinite then
         if X.Is_Positive xor Y.Is_Positive then -- opposite signs.
            X := Negative_Infinity;
         else
            X := Positive_Infinity;
         end if;
         return;
      end if;

      -- Step 3. Do the multiplication.  We can only sum 32 elements of the sum
      -- before the carrys must be done.

      -- We use the inner product version (inner loop is an
      -- dot-product.)  To get the outer product version (BLAS routine DAXPY),
      -- interchange the order of the loops.
      --
      -- Here's the idea:
      -- for Digit_ID in Digit_Index loop
      --    Sum := 0.0;
      --    for k in 0..Digit_ID loop
      --       Sum := Sum + X(k) * Y(Digit_ID-k);
      --    end loop;
      --    Z(Digit_ID) := Sum;
      -- end loop;
      --
      -- Break sum into segments of 32 each in index k.  Perform a carry
      -- after each sum of 32 elements.

      if (Ultimate_No_of_Digits < 9  and then No_Of_Bits_In_Radix <= 30) or
         (Ultimate_No_of_Digits < 17 and then No_Of_Bits_In_Radix <= 29) then

         --  Notice we only change X here, not above, in case the else is used.

         X.Is_Positive := not (X.Is_Positive XOR Y.Is_Positive);

         X.Exp := X.Exp + Y.Exp; -- Will be further adjusted by "Normalize.." routine.

         Product_if_digits_fewer_than_17;

         Do_Carrying_For_Multiplication (X, Digit_Minus_1, Digit_Minus_2);

      else

         X := General_Multiplication (X, Y);
         return;

      end if;

      Normalize (X, Digit_Minus_1, Digit_Minus_2);

      -- Step 4. Handle over and under flow.  Here underflow goes to Zero,
      -- overflow to Infinity.  This analysis is all isolated to the end
      -- of the arithmetic routines so that it is more easily modified to
      -- raise exceptions if that's what is desired.   In order to do it
      -- here, must assume that the parmeters Min_Exponent and Max_Exponent
      -- limit the dynamic range of the Exp to about 1/4 of that allowed
      -- by the base type used to represent the exponent.  This is
      -- checked in the spec with an assertion.  (The reason is, the above
      -- code will go well outside the accepted range of Exp with out being
      -- checked till down here.)  This limit is OK because the base type
      -- allows excessively large exponents anyway, up to 2**31-1.

      if X.Exp < Min_Exponent then
         X := Zero;
      end if;

      if X.Exp > Max_Exponent then
         if X.Is_Positive then
            X := Positive_Infinity;
         else
            X := Negative_Infinity;
         end if;
      end if;

   end Mult;

   -----------------------
   -- Multiply_In_Place --
   -----------------------

   function Multiply_In_Place 
     (X, Y : e_Real)
      return e_Real
   is
      Z   : e_Real := X;
   begin
      Mult (Z, Y);
      return Z;
   end Multiply_In_Place;

   function "*" (X : e_Real; Y : e_Real) return e_Real renames Multiply_In_Place;
 --function "*" (X : e_Real; Y : e_Real) return e_Real renames Multiply_stnd;
 --  These 2 about the same. Multiply_In_Place calls procedure Mult.

   ---------
   -- "*" --
   ---------

   -- Multiply the first and only digit of X by Y.

   function "*"(X : e_Digit;
                Y : e_Real) return e_Real is

      Z   : e_Real;
      Digit_Minus_1 : Digit_Type := Digit_Zero; --Essential init
      Digit_Minus_2 : constant Digit_Type := Digit_Zero; --Essential init
      Carry_Minus_1 : Digit_Type := Digit_Zero; --Essential init
      I : Digit_Index;

   begin


     -- Step 0.  Sign etc.
     -- Notice this assumes 0 * inf = 0.

     if X.Is_Zero or Y.Is_Zero then
        return Zero;
     end if;


     -- Step 1.  Infinities.  Digit can't be inf.

     Z.Is_Positive := not (X.Is_Positive XOR Y.Is_Positive);

     if Y.Is_Infinite then
        if Z.Is_Positive then
           return Positive_Infinity;
        else
           return Negative_Infinity;
        end if;
     end if;

     Z.Is_Zero     := False;
     Z.Exp         := X.Exp + Y.Exp;

     for I in Digit_Index'First .. Digit_Index'Last loop
        Z.Digit(I) := X.Digit * Y.Digit(I);
     end loop;


     -- Step 2.  Do the carries.  This is simpler than the more general version
     -- because you carry at most one digit.  (The Max is of Digit*Digit is
     -- Radix**2 - 2*Radix + 1. This results in a carry of approx. Radix.  Add
     -- this the the next higher order digit to get a max of Radix**2 - Radix.
     -- So still don't overflow the Range over which a single carry is
     -- all that's needed: 0..Radix**2 - 1.)
     -- Carry_Minus_2 is always zero, and Digit_Minus_2 is always 0.

     for I in reverse Digit_Index'First+1 .. Digit_Index'Last loop
        Carry_Minus_1      := Shift_Right_No_of_Bits_in_Radix (Z.Digit(I));
        Z.Digit(I)         := Z.Digit(I) - Carry_Minus_1 * Digit_Radix;
        Z.Digit(I-1)       := Z.Digit(I-1) + Carry_Minus_1;
     end loop;

     -- Special case I = Digit_Index'First = 0

        I := 0;
        Carry_Minus_1      := Shift_Right_No_of_Bits_in_Radix (Z.Digit(I));
        Z.Digit(I)         := Z.Digit(I) - Carry_Minus_1 * Digit_Radix;
        Digit_Minus_1      := Carry_Minus_1;


     -- Step 3.  Must Normalize: shift digit array to make sure that the first
     -- digit is non-zero.

     Normalize (Z, Digit_Minus_1, Digit_Minus_2);


     -- Step 4. Handle over and under flow.  Here underflow goes to Zero,
     -- overflow to Infinity.

     if Z.Exp < Min_Exponent then
        Z := Zero;
     end if;

     if Z.Exp > Max_Exponent then
        if Z.Is_Positive then
           Z := Positive_Infinity;
        else
           Z := Negative_Infinity;
        end if;
     end if;

     return Z;

   end "*";

   -------------
   -- Scaling --
   -------------

   --  S'Scaling (X, Exp)
   --                                      Exp
   --  Let v be the value X*T'Machine_Radix  .  If v is a
   --  machine number of the type T, or if |v|GT'Model_Small,
   --  the function yields v; otherwise, it yields either one
   --  of the machine numbers of the type T adjacent to v.
   --  Constraint_Error is optionally raised if v is outside
   --  the base range of S. A zero result has the sign of X
   --  when S'Signed_Zeros is True.
   --
   function Scaling (X : e_Digit;
                     Adjustment : e_Integer)       return e_Digit is
     X2 : e_Digit := X;
   begin
     if X.Is_Zero then
        return X2;
     end if;

     X2.Exp := X.Exp + Adjustment;

     if X2.Exp < Min_Exponent or else X2.Exp > Max_Exponent then
        raise Constraint_Error with "Exp out of range in Scaling e_Digit operation.";
     end if;

     return X2;
   end Scaling;

   -------------------
   -- Make_Extended --
   -------------------

   --  Turn Digit into Extended:

   function Make_Extended (X : e_Digit) return e_Real is
      Z : e_Real; -- initialized to Zero. Import.
   begin

      Z.Digit(0)    := X.digit;
      Z.Is_Positive := X.Is_Positive;
      Z.Exp         := X.Exp;
      Z.Is_Zero     := X.Is_Zero;
      Z.Is_Infinite := False;
      return Z;

   end Make_Extended;

   ------------
   -- Sum_Of --
   ------------

   --  Optimized SUM routine for e_Digit + e_Real -> e_Real.

   function Sum_Of (X : e_Digit; Y : e_Real) return e_Real is
      Z           : e_Real := Y;
      Delta_Exp   : e_Integer;

      New_Digit_1   : Digit_Type    := Digit_Zero;
      New_Digit_2   : Digit_Type    := Digit_Zero;
      Need_To_Carry : Boolean := False;

      type Max_Info is (X_Is_Max, Y_Is_Max);
      Max_Num_ID  : Max_Info := X_Is_Max;

   begin

      -- Step 0. If Either of the numbers is 0.0, then return the other.

      if X.Is_Zero then
         return Y;
      elsif Y.Is_Zero then
         return Make_Extended (X);
      end if;

      -- Step 0b. If X is infinite, it doesn't matter what do.

      if Y.Is_Infinite then
         return Y;
      end if;

      -- Step 0c. If one is positive and the other neg., then are lazy.
      -- Do it the slow way.

      if X.Is_Positive XOR Y.Is_Positive then
         return (Y + Make_Extended (X));
      end if;

      -- Step 1. Now they either both Positive or both Neg.  Sum them if it's
      -- easy and return with the right sign.  Start by finding the larger
      -- number, Exponent-wise. *Notice Y.Exp = X.Exp is classified as Y_Is_Max.*

      if X.Exp > Y.Exp then
         Max_Num_ID  := X_Is_Max;
      else
         Max_Num_ID  := Y_Is_Max;
      end if;

      Delta_Exp := Abs (X.Exp - Y.Exp);

      if Delta_Exp > e_Integer(Digit_Index'Last) then -- ie, Delta_Exp >= No_Of_Digits
         case Max_Num_ID is
         when Y_Is_Max =>
            return Y;
         when X_Is_Max =>
            return Make_Extended (X);
         end case;
      end if;

      -- Step 2. If the digit X has the smaller exponent, then try
      -- an optimization.  Otherwise just use the full scale "+".
      -- We verified above that Delta_Exp is in range of index of X.Digit.
      -- The optimization covers the most common case by far in the "/" routine
      -- and most other uses of this function: Y > X, and same sign.
      -- We allow at most one carry..if more are required, give up trying.

      if Max_Num_ID = Y_Is_Max then

         New_Digit_1 := X.digit + Y.Digit(Delta_Exp);
         if New_Digit_1 > Digit_Radix_Minus_1 then
            Need_To_Carry := True;
         end if;

         if Need_To_Carry and then (Delta_Exp = 0) then -- would have to normalize.
            goto Abort_Optimization;
         end if;

         if Need_To_Carry then
            New_Digit_1 := New_Digit_1 - Digit_Radix;
            New_Digit_2 := Y.Digit(Delta_Exp-1) + Digit_One; -- Carry.
            if New_Digit_2 > Digit_Radix_Minus_1 then        -- give up trying.
               goto Abort_Optimization;
            end if;
         end if;

         --  If got this far, then are going through with it.
         --  just change the 1 or 2 digits of X, call it Z, and return it:

         -- Z := Y;
         -- Z is initialized to Y:

         Z.Digit(Delta_Exp) := New_Digit_1;
         if Need_To_Carry then
            Z.Digit(Delta_Exp-1) := New_Digit_2;
         end if;

         --  Recall that X and Y have same sign. This must be the
         --  sign of Z.  OK, the Z := Y establishes this.

         return Z;

      end if;

      <<Abort_Optimization>>

      -- Step 3. If got this far, the optimization failed.  Do it the
      -- slow way.

      return (Y + Make_Extended (X));

   end Sum_Of;

   ---------
   -- "/" --
   ---------

   -- Calculate   Quotient = X / Y.
   -- Schoolboy algorithm.  Not fast.  Get Quotient a digit at a time.
   -- Estimate trial Next_Digit by dividing first few digits of X by first few
   -- digits of Y.  Make this ratio into the proper range for a digit.
   -- Multiply this Next_Digit by Y to get Product.  (This step is why do it
   -- a digit at a time.  The above multiplication by Y can be highly optimized
   -- then.)  Subtract Product from X to get Remainder.  Increment Quotient
   -- by Next_Digit.  Repeat the above steps with X replaced by Remainder.
   -- Continue until Remainder has an exponent so small that subsequent
   -- Next_Digit's are too small to contribute to Quotient.
   --
   -- Few remarks: Next_Digit may not really be the next digit in Quotient.
   -- It should usually be, but sometimes it's too small, and very rarely, it's
   -- too large.  Nevertheless, you increment Quotient by Next_Digit and
   -- converge on the answer.  Also, when Next_Digit is too large, Remainder
   -- becomes negative, and one *subtracts* the next Next_Digit from quotient,
   -- and *adds* the next Product to Remainder, and continues doing this as
   -- long as Remainder is negative.  In other words can converge on the
   -- quotient from either above or below.  So are calculating each
   -- iteration another pseudo-digit, DigitJ, such that, for example.
   --
   --  Y * (Digit0  + Digit1 + Digit3 - Digit4 + Digit5 ...) = X,
   --
   -- Add up the digits to get the approximation of X / Y.
   -- Some are negative (rare), some positive, but they add up to Quotient.
   -- More precisely,
   --
   --  X - Y * (Digit0  + Digit1 + Digit3 - Digit4 + Digit5 ...) = Remainder,
   --
   -- where X.Exp - Remainder.Exp  > Digit_Index'Last. Further
   -- Digits contribute negligably to Quotient when this inequality holds.
   -- For example, if there are 2 digits, then Digit_Index'Lasts = 1, so that if
   -- X.Exp - Remainder.Exp > Digit_Index'Last = 1 then the difference in
   -- exponents guarantees that subsequent iterations contribute to digits
   -- beyond Digit_Index'Last in the Quotient.

   function "/" (X, Y : e_Real)
      return e_Real
   is
      Delta_Remainder     : e_Real;
      Remainder_Of_X      : e_Real := X;    -- init important.
      Quotient            : e_Real;  -- Initialized to zero (important).
      Next_Digit          : e_Digit; -- Initialized to zero (important).
      Next_Digit_Val      : Digit_Type;
      Real_Next_Digit_Val : Real;
      Real_Y, Inverse_Real_Y, Real_Remainder_Of_X : Real;
      Decrement_Digit_Exp  : Boolean := False;

      Count : e_Integer := 0;
      Max_Allowed_Iterations : constant e_Integer := Digit_Index'Last * 2;
   begin
      -- Step 0.  Handle Zeros and Infinities.  Do signs at very end.  Until
      -- the very end pretend that Quotient is positive.

      if Y.Is_Zero then
         raise Constraint_Error with "Division by zero.";
      end if;

      if X.Is_Zero then
         return Zero;
      end if;

      if X.Is_Infinite and Y.Is_Infinite then
         raise Constraint_Error with "Division of inf by inf is undefined.";
      end if;

      if X.Is_Infinite then
         if  (X.Is_Positive xor Y.Is_Positive) then
            return Negative_Infinity;
         else
            return Positive_Infinity;
         end if;
      end if;

      if Y.Is_Infinite and (not X.Is_Infinite) then
         return Zero;
      end if;

      -- Step 1.  We initialized Remainder_Of_X to X.  Below make it positive,
      -- and assume all quantities are positive.  Signs are set at end.
      -- Inverse_Real_Y      is in range [1.0/(Radix-Eps) .. 1.0].
      -- Real_Remainder_Of_X is in range [1.0 .. (Radix-Eps)].
      -- Next_Digit is Real_Floor (Real_Remainder_Of_X*Inverse_Real_Y) where
      -- the product has been multiplied by Radix if it's less than 1.0.
      -- Possible range of Next_Digit: Max is Floor (Radix-eps / 1.0), which
      -- is Radix-1.  Min is Floor (Radix * (1.0 / (Radix-eps))) = 1.0.

      Real_Y  :=   Real (Y.Digit(0))
                 + Real (Y.Digit(1)) * Inverse_Radix
                 + Real (Y.Digit(2)) * Inverse_Radix_Squared;

      Inverse_Real_Y  := Real_One / Real_Y;

      --  Pretend all quantities are positive till end.
      Quotient.Is_Positive            := True;
      Remainder_Of_X.Is_Positive      := True; -- This may go negative.

    --  Iterate until the remainder is small enough:

    Iteration: for Iter in e_Integer range 1..Max_Allowed_Iterations loop

      -- Initialized to X, so it's not zero.  Must lead with a nonzero digit.
      -- Important here to sum ALL three leading digits of Remainder_Of_X.

      Real_Remainder_Of_X :=   Real (Remainder_Of_X.Digit(0))
                             + Real (Remainder_Of_X.Digit(1))*Inverse_Radix
                             + Real (Remainder_Of_X.Digit(2))*Inverse_Radix_Squared;

      Real_Next_Digit_Val := Real_Remainder_Of_X * Inverse_Real_Y;

      -- Need Next_Digit in proper range for Mult_For_Div fctn: 1..Radix-1:

      Decrement_Digit_Exp := False;
      if Real_Next_Digit_Val < Real_One then
         Real_Next_Digit_Val := Real_Next_Digit_Val * Real_Radix;
         Decrement_Digit_Exp := True;
      end if;

      Next_Digit_Val := Digit_Floor (Real_Next_Digit_Val);

      if Next_Digit_Val > Digit_Radix_Minus_1 then
         Next_Digit_Val := Digit_Radix_Minus_1;
      end if;

      --if Next_Digit_Val < Digit_One then   -- Should never happen, but test it.
      --   Next_Digit_Val := Digit_One;
      --end if;


      -- Step 2.  We now have Next_Digit, which is approx. Remainder_Of_X / Y.
      -- We are ready to make Next_Digit an e_Digit number.
      -- (so "*" can be optimized.)  It must be in
      -- the range 1..Radix-1, so truncate first.  It has the exponent of
      -- the Remainder_Of_X.Exp - Y.Exp, but decremented by 1.0 if multiplied
      -- Next_Digit by Radix above.
      -- Below remember that Next_Digit was initialized to Zero above,
      -- and all that is changed in this loop is the digit value, and the Exp.
      -- Also want the product of Next_Digit and Y to be >0, regardless of Y.
      -- The conversion from Real to e_Digit could be done by the
      -- function Make_e_Digit, but lets optimize:

      Next_Digit.Digit       := Next_Digit_val;
      Next_Digit.Exp         := Remainder_Of_X.Exp - Y.Exp;
      Next_Digit.Is_Zero     := False;
      Next_Digit.Is_Positive := True;

      if Decrement_Digit_Exp then
         Next_Digit.Exp := Next_Digit.Exp - 1;
      end if;


      -- Step 3. Make the trial product of the next digit with the divisor..
      -- this will be subtracted from the remainder.

      Delta_Remainder             := Next_Digit * Y;
      Delta_Remainder.Is_Positive := True;


      -- Step 3.  Calculate the new "Quotient" and the new "Remainder_Of_X".
      -- Add Extended_Next_Digit to Quotient, if Old Remainder was > 0.
      -- Subtract Extended_Next_Digit from Quotient if Old Remainder was < 0.
      -- If the Old Remainder = 0 then are done; that was checked on
      -- previous pass of loop by line below.  It was checked initially by
      -- by checking if X = 0.  (Remainder was initialized to X, and then made
      -- positive).

      if Remainder_Of_X.Is_Positive then
          -- Add Next_Digit to Quotient, and subtract Delta_Remainder from Remainder.
          Delta_Remainder.Is_Positive := False;
          Next_Digit.Is_Positive      := True;
      else
          -- Subtract Next_Digit from Quotient, and add Delta_Remainder to Remainder.
          Delta_Remainder.Is_Positive := True;
          Next_Digit.Is_Positive      := False;
      end if;

      Remainder_Of_X  := Remainder_Of_X + Delta_Remainder;
      Quotient        := Sum_Of (Next_Digit, Quotient);


      -- Step 4.  Are finished?
      -- Remainder_Of_X.Exp started at X.Exp.  Have made it small enough?
      -- Remember that the calls above may have returned Infinity.  Or if they
      -- returned Zero, then the Exponent is 0. (None of these should happen tho.)

      if Remainder_Of_X.Is_Zero then exit Iteration; end if;

      if Remainder_Of_X.Is_Infinite then exit Iteration; end if;

      if (X.Exp-Remainder_Of_X.Exp) > Digit_Index'Last then exit Iteration; end if;

      Count := Iter;

    end loop Iteration;

      --  Max_Allowed_Iterations is twice the number usually required.  I've
      --  never seen it happen, but just in case:

      if Count = Max_Allowed_Iterations then  -- Should raise error?
          raise Program_Error with "Convergence problem in division routine.";
      end if;
      --Print_Text (Integer'Image(Count));
      --Print_Text (" Iterations.");
      --text_io.new_line;


      -- Step 5. Set Sign of Quotient.  Handle over and under flow.
      -- Here underflow goes to Zero, overflow to Infinity.
      -- We can do this here because allowed room to maneuver in range of Exp.

      Quotient.Is_Positive := not (X.Is_Positive xor Y.Is_Positive);

      if Quotient.Exp < Min_Exponent then
         Quotient := Zero;
      end if;

      if Quotient.Exp > Max_Exponent then
         if Quotient.Is_Positive then
            Quotient := Positive_Infinity;
         else
            Quotient := Negative_Infinity;
         end if;
      end if;

      return Quotient;

   end "/";

   ---------
   -- "*" --
   ---------

   -- Multiply the first and only digit of X by the first and only digit
   -- of Y.  the reduction in Carrying is the big win here.

   function "*"(X : e_Digit;
                Y : e_Digit) return e_Real is

      Z   : e_Real;  -- Init important.
      Digit_Minus_1 : Digit_Type := Digit_Zero; -- Init important.
      Digit_Minus_2 : constant Digit_Type := Digit_Zero; -- Init important.
      Carry_Minus_1 : Digit_Type := Digit_Zero; -- Init important.
      I : Digit_Index;

   begin


     -- Step 0.  Zeros.

     if X.Is_Zero or Y.Is_Zero then
        return Zero;
     end if;

     Z.Is_Zero  := False;


     -- Step 1.  Sign and Infinities.  Digits can't be inf.

     Z.Is_Positive := not (X.Is_Positive XOR Y.Is_Positive);

     Z.Exp      := X.Exp + Y.Exp;

     Z.Digit(0) := X.Digit * Y.Digit;


     -- Step 2.  Do the carries.  This is simpler than the more general version
     -- because you carry at most one digit.  (The Max is of Digit*Digit is
     -- Radix**2 - 2*Radix + 1. This results in a carry of approx. Radix.  Add
     -- this the the next higher order digit to get a max of Radix**2 - Radix.
     -- So still don't overflow the Range over which a single carry is
     -- all that's needed: 0..Radix**2 - 1.)
     -- Carry_Minus_2 is always zero, and Digit_Minus_2 is always 0.

      I := 0;
      Carry_Minus_1      := Shift_Right_No_of_Bits_in_Radix (Z.Digit(I));
      Z.Digit(I)         := Z.Digit(I) - Carry_Minus_1 * Digit_Radix;
      Digit_Minus_1      := Carry_Minus_1;


     -- Step 3.  Must Normalize: shift digit array to make sure that the first
     -- digit is non-zero.

     Normalize (Z, Digit_Minus_1, Digit_Minus_2);


     -- Step 4. Handle over and under flow.  Here underflow goes to Zero,
     -- overflow to Infinity.

     if Z.Exp < Min_Exponent then
        Z := Zero;
     end if;

     if Z.Exp > Max_Exponent then
        if Z.Is_Positive then
           Z := Positive_Infinity;
        else
           Z := Negative_Infinity;
        end if;
     end if;

     return Z;

   end "*";


   ---------
   -- "/" --
   ---------

   -- The only difference between this "/" and the general "/" is that
   -- the (e_Digit * e_Real) operation has been replaced
   -- with an (e_Digit * e_Digit) operation.  This produces
   -- big savings in time.  This is used in elementary math function
   -- routines, where the efficiency is important.

   function "/" (X : e_Real; Y : e_Digit) return e_Real is
     Delta_Remainder     : e_Real;
     Remainder_Of_X      : e_Real := X;    -- init important.
     Quotient            : e_Real;  -- Initialized to zero (important).
     Next_Digit          : e_Digit; -- Initialized to zero (important).
     Real_Next_Digit_Val : Real;
     Next_Digit_Val      : Digit_Type;
     Real_Y, Inverse_Real_Y, Real_Remainder_Of_X : Real;
     Decrement_Digit_Exp  : Boolean := False;

     Count : e_Integer := 0;
     Max_Allowed_Iterations : constant e_Integer := Digit_Index'Last * 2;
   begin


     -- Step 0.  Handle Zeros and Infinities.  Do signs at very end.  Until
     -- the very end pretend that Quotient is positive.

     if Y.Is_Zero then
        raise Constraint_Error with "Division by zero.";
     end if;

     if X.Is_Zero then
        return Zero;
     end if;

     if X.Is_Infinite then
        if  (X.Is_Positive xor Y.Is_Positive) then
           return Negative_Infinity;
        else
           return Positive_Infinity;
        end if;
     end if;


     -- Step 1.  We initialized Remainder_Of_X to X.  Below make it positive,
     -- and assume all quantities are positive.  Signs are set at end.
     -- Possible range of Next_Digit: Max is Floor (Radix-eps / 1.0), which
     -- is Radix-1.  Min is Floor (Radix * (1.0 / (Radix-eps))) = 1.0.
     -- Remember, X and Y are normalized and we've tested for 0's, so
     -- both are in the range 1..Radix-1.

     Real_Y         := Real (Y.Digit);        -- >= 1.0, since it's quantized.
     Inverse_Real_Y := Real_One / Real_Y;

     --  Pretend all quantities are positive till end.
     Quotient.Is_Positive            := True;
     Remainder_Of_X.Is_Positive      := True; -- This may go negative.

   --  Iterate until the remainder is small enough.

   Iteration: for Iter in e_Integer range 1..Max_Allowed_Iterations loop

     Real_Remainder_Of_X :=
        Real (Remainder_Of_X.Digit(0))
      + Real (Remainder_Of_X.Digit(1))*Inverse_Radix
      + Real (Remainder_Of_X.Digit(1))*Inverse_Radix_Squared;

     Real_Next_Digit_Val := Real_Remainder_Of_X * Inverse_Real_Y;

     -- Need Next_Digit in proper range for Mult_For_Div fctn: 1..Radix-1:

     Decrement_Digit_Exp := False;
     if Real_Next_Digit_Val < Real_One then
        Real_Next_Digit_Val := Real_Next_Digit_Val * Real_Radix;
        Decrement_Digit_Exp := True;
     end if;

     Next_Digit_Val := Digit_Floor (Real_Next_Digit_Val);

     if Next_Digit_Val > Digit_Radix_Minus_1 then
        Next_Digit_Val := Digit_Radix_Minus_1;
     end if;

     --if Next_Digit_Val < Digit_One then   -- Should never happen. Perform tests w/o
     --   Next_Digit_Val := Digit_One;
     --end if;


     -- Step 2.  We now have Next_Digit, which is approx. Remainder_Of_X / Y.
     -- We are ready to make Next_Digit an e_Digit number.
     -- (so "*" can be optimized.)  It must be in
     -- the range 1..Radix-1, so truncate first.  It has the exponent of
     -- the Remainder_Of_X.Exp - Y.Exp, but decremented by 1.0 if multiplied
     -- Next_Digit by Radix above.
     -- Below remember that Next_Digit was initialized to Zero above,
     -- and all that is changed in this loop is the digit value, and the Exp.
     -- Also want the product of Next_Digit and Y to be >0, regardless of Y.
     -- The conversion from Real to e_Digit could be done by the
     -- function Make_e_Digit, but lets optimize:

     Next_Digit.Digit       := Next_Digit_Val;
     Next_Digit.Exp         := Remainder_Of_X.Exp - Y.Exp;
     Next_Digit.Is_Zero     := False;  -- Is this always true?
     Next_Digit.Is_Positive := True;

     if Decrement_Digit_Exp then
        Next_Digit.Exp := Next_Digit.Exp - 1;
     end if;


     -- Step 3. Make the trial product of the next digit with the divisor..
     -- this will be subtracted from the remainder.

     Delta_Remainder             := Next_Digit * Y;
     Delta_Remainder.Is_Positive := True;


     -- Step 3.  Calculate the new "Quotient" and the new "Remainder_Of_X".
     -- Add Extended_Next_Digit to Quotient, if Old Remainder was > 0.
     -- Subtract Extended_Next_Digit from Quotient if Old Remainder was < 0.
     -- If the Old Remainder = 0 then are done; that was checked on
     -- previous pass of loop by line below.  It was checked initially by
     -- by checking if X = 0.  (Remainder was initialized to X, and then made
     -- positive).

     if Remainder_Of_X.Is_Positive then
         --  Add Next_Digit to Quotient, and
         --  subtract Delta_Remainder from Remainder.
         Delta_Remainder.Is_Positive := False;
         Next_Digit.Is_Positive      := True;
     else
         --  Subtract Next_Digit from Quotient, and
         --  add Delta_Remainder to Remainder.
         Delta_Remainder.Is_Positive := True;
         Next_Digit.Is_Positive      := False;
     end if;

     Remainder_Of_X  := Remainder_Of_X + Delta_Remainder;
     Quotient        := Sum_Of (Next_Digit, Quotient);


     -- Step 4.  Are finished?
     -- Remainder_Of_X.Exp started at X.Exp.  Have made it small enough?
     -- Remember that the calls above may have returned Infinity.  Or if they
     -- returned Zero, then the Exponent is 0. (None of these should happen tho.)

     if Remainder_Of_X.Is_Zero then exit Iteration; end if;

     if Remainder_Of_X.Is_Infinite then exit Iteration; end if;

     if (X.Exp-Remainder_Of_X.Exp) > Digit_Index'Last then exit Iteration; end if;

     Count := Iter;

   end loop Iteration;


   --  Max_Allowed_Iterations is twice the number usually required.  I've
   --  never seen it happen, but just in case:

   if Count = Max_Allowed_Iterations then  -- Should raise error?
      raise Program_Error with "Convergence problem in division routine.";
   end if;
   --Print_Text (Integer'Image(Count));
   --Print_Text (" Iterations.");
   --text_io.new_line;


   -- Step 4. Set Sign of Quotient.  Handle over and under flow.
   -- Here underflow goes to Zero, overflow to Infinity.

   Quotient.Is_Positive := not (X.Is_Positive xor Y.Is_Positive);

   if Quotient.Exp < Min_Exponent then
      Quotient := Zero;
   end if;

   if Quotient.Exp > Max_Exponent then
      if Quotient.Is_Positive then
         Quotient := Positive_Infinity;
      else
         Quotient := Negative_Infinity;
      end if;
   end if;

   return Quotient;

   end "/";

   ----------
   -- "**" --
   ----------

   -- Standard algorithm.  Write N in binary form:
   --
   --   N = 2**m0 + 2**m1 + 2**m2 + ... + 2**mn,
   --
   -- where the m0, m1, m2...are the indices of the nonzero binary digits of N.
   -- Then
   --
   --   X**N = X**(2**m0) * X**(2**m1) * X**(2**m2) * ... * X**(2**mn).
   --
   -- We have written X as product of Powers of X.
   -- Powers of X are obtained by squaring X, squaring the square of X,
   -- etc.
   --
   -- 0.0**0 is defined to be 1.0.  Anything to the 0 is defined to be 1.
   --
   function "**"
     (X : e_Real;
      N : Integer)
      return e_Real
   is
      Power_Of_X, Product : e_Real;
      Last_Bit_ID         : constant Integer := Integer'Size - 2;
      subtype Bit_Index is Integer range 0..Last_Bit_ID;
      type Bit_Array is array(Bit_Index) of Integer;
      Bit         : Bit_Array := (others => 0);
      Exp         : Integer   := Abs (N);
      Final_Bit   : Bit_Index;
      Exp_Is_Even : Boolean;
   begin
      -- The following seems to be what the lrm implies, even if X=0.0:
      -- If the exponent N is zero then return 1.  (If X is inf this makes
      -- sense since inf is really just a large finite number.)

      if N = 0 then
         return One;
      end if;

      -- Next, if the argument X is zero, set the result to zero.
      -- If the Exponent is negative then raise constraint error.

      if X.Is_Zero then
         if N > 0 then
            return Zero;
         else
            raise Constraint_Error with "Error in ** operation, division by 0.0.";
         end if;
      end if;

      -- If the argument X is inf, set the result to inf, or 0 if N < 0.
      -- If X < 0 then sign depends on whether Exp is even or odd.
      -- (IS THIS TRUE if N=0?  Will assume NO.)

      if X.Is_Infinite then
         Exp_Is_Even := ((Exp rem 2) = 0);

         if X.Is_Positive then
            Product := Positive_Infinity;
         else
            if Exp_Is_Even then
               Product := Positive_Infinity;
            else
               Product := Negative_Infinity;
            end if;
         end if;

         if N < 0 then
            return Zero;
         else
            return Product;
         end if;
      end if;

      -- Should try to avoid possible run-time errors (if, for example,
      -- N * X.Exp > Max_Exponent)?  No: the following
      -- algorithm uses extended precision "*", which will overflow to
      -- inf and short-circuit the process efficiently.
      --
      -- Get binary form of the exponent Exp = Abs(N)

      for I in Bit_Index loop
         Bit(I)    := Exp REM 2;
         Exp       := Exp  /  2;
         Final_Bit := I;
         if Exp = 0 then
            exit;
         end if;
      end loop;

      -- Do the arithmetic.

      Product    := One; -- X**0
      Power_Of_X := X;   -- X**(2**0).  Preserves sign of X if N is odd.
      if Bit(Bit_Index'First) /= 0 then
       --Product := Product * Power_Of_X;
         Mult (Product, Power_Of_X);
      end if;
      for I in Bit_Index'First+1 .. Final_Bit loop
       --Power_Of_X := Power_Of_X * Power_Of_X;     -- X**(2**I)
         Square (Power_Of_X);                       -- X**(2**I), good speed-up
         if Bit(I) = 1 then
          --Product := Product * Power_Of_X;
            Mult (Product, Power_Of_X);
         end if;
      end loop;

      --  Under flow to zero.  THe "/" operator should correctly do this, but
      --  it's important to other routines, so make sure it's done right here:

      if Product.Is_Infinite and N < 0 then
         return Zero;
      end if;

      if N < 0 then
         Product := One / Product; -- notice we've already checked for X=0.
      end if;

      -- need to do the final over/under flow check?  No.  The "/" and
      -- "*" did that for us above.

      return Product;

   end "**";

end Extended_Real;

