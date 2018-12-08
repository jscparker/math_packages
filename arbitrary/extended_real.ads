
-- PACKAGE Extended_Real
--
-- package Extended_Real provides:
--    An arbitrary precision floating-point data type:  e_Real.
--
--    Lower limit on precision is 28 decimals. No upper limit is
--    enforced. All internal arithmetic is done on 64-bit Integers,
--    so its most efficient on 64-bit CPU's.  The package is Pure. 
--    Floating point attributes (Ada 95) are implemented as function
--    calls. The package exports standard floating point operators:
--    "*", "+", "/", "**", "Abs", "<", ">", "<=" , ">=", etc. 
--    The standard operators make it easy to modify existing code to
--    use extended precision arithmetic. Procedure calls Mult(X,Y) and
--    Square(X) are also provided. They do multiplication "in-place",
--    (overwrite X with the result) and are somewhat faster than the 
--    equivalent X := X*Y, and X := X*X.
--
-- To set the precision search below for:  
--
--    Desired_Decimal_Digit_Precision  
--
-- package Extended_Real.Elementary_Functions provides:
--    Sin, Cos, Sqrt, Arcsin, Arccos, Arctan, Log, Exp, Reciprocal (1/x),
--    Reciprocal_Nth_Root (x to the power of -1/N), Divide, and "**" for
--    e_Real arguments and e_Real exponents. Routines are Ada 95'ish.
--
-- package e_Derivs provides:
--    Extended precision routines for taking high order derivatives of
--    functions.  Functions made from  "*", "+", "/", "**", Sin, Cos, 
--    Sqrt, Arcsin, Arccos, Arctan, Log, Exp, Compose = f(g(x)),
--    and Reciprocal can be differentiated to order specified by user.
--
-- package Extended_Real.IO provides:
--    Text to extended-precision e_Real translation routines, and
--    e_Real to Text translation routines.
--
-- package Extended_Real.Rand provides:
--    a (very) basic Random Number Generator. Really just for making 
--    test vectors.
--
-- procedure e_real_demo_1.adb is:
--    an introductory routine that demonstrates use of Extended_Real.
--
-- procedure e_function_demo_1.adb is:
--    an introductory routine that demonstrates use of 
--    Extended_Real.Elementary_Functions.
--
-- procedure e_jacobi_eigen_demo_1.adb demonstrates:
--    extended-precision eigen-decomposition on Hilbert's matrix using
--    package e_Jacobi_Eigen.
--
-- package e_Jacobi_Eigen is:
--    a Jabobi iterative eigen-decomposition routine
--    that shows how easy it is to upgrade a floating point routine
--    to extended precision.  e_Jacobi_Eigen uses package Extended_Real.
--
-- good optimization on gcc/GNAT:
--    -gnatNp -O3 -march="your machine architecture" -funroll-loops -ffast-math
--    (sometimes:-ftree-vectorize  -funroll-all-loops  -falign-loops=4, 
--     -falign-loops=3, or -frename-registers are worth trying.)
--
-- latest GNAT (gcc 4.3) try
--    gnatmake -gnatNp -O3 -march=native -mtune=native -funroll-loops -ffast-math
--
-- Always do a preliminary run which exercizes Assertions, and other Checks:
--    -gnato -gnatV -gnata
--
--
-- Because precision is arbitrary, Extended_Real is not specially
-- optimized for any particular precision. The chosen design works best
-- in the limit of 100's of decimal digits.  If the package had been
-- designed for 32 decimal digits of precision, then almost every feature
-- of the design would have been different. On the other hand, performance
-- seems to be respectable on 64-bit CPU's even at the lower limit (eg
-- 28 or 38 decimal digits). (Comparison is with Intel's machine-optimized
-- 32 decimal-digit floating point: i.e. Intel Fortran Real*16 on an Intel
-- 64-bit CPU.) 32 digit floating point is probably the most often used 
-- (and most often needed) extended precision floating point.
-- Most Fortrans (including the gcc Fortran) don't offer anything higher 
-- than 18 digit floating point.
--
-- Common applications:
-- 0. Estimation of error in lower precision floating-point calculations.
-- 1. Evaluation of constants for math routines and Table-driven algorithms.
-- 2. Evaluation of series solutions of special function, especially when
--    the terms are constructed of large factorials and exponentials.
-- 3. Evaluation of recurrance relations for special functions.
--
-- Generics greatly reduce the work you have to do in modifying programs
-- to use extended floating point:
-- 
-- 1. place generic formal declarations
--    of the required extended arithmetic functions at the the top of the
--    package or subprogram to be modified.
--
-- 2. use the unary "-" and "+" routines that convert Real to Extended:
--
--    so that declarations
--      Number : Generic_Formal_Type := +1.234;
--    and statements like
--      Z  := (+4.567834E+012) * X;
--    will be acceptible to both Real and Extended types.
--
-- Underflows to Zero.  Overflows to Positive or Negative infinity.  I am
-- still trying to decide if the Infinities are worth the trouble, but the
-- trouble isn't great and there seem to be benefits.  Sometimes you can
-- put off worrying about overflow in intermediate calculation and test
-- for it at the end by testing for Positive_Infinity.  There are no NaNs.
--
-- At the expense of purity, error messages using text_io can be 
-- re-enabled in the body - see top of body of Extended_Real.
--
--***************************************************************************
--
-- SECTION I.
--
-- Constants and overflow/underflow/constraint_error conventions.
--
-- To test an arbitrary X : e_Real to see if X is Zero or infinity use the
-- function Are_Equal (X, Zero) etc.  Its written to make the test efficiently.
--
-- Underflows are to (unsigned) Zero; overflows to (signed) infinity:
--
-- Infinity here means a finite number that is too large to represent 
-- in the floating point, and whose inverse is too small to represent in 
-- floating point.  The following conventions seemed sensible. Treat inf's
-- as Constraint Errors if uneasy with them. Assuming X is a positive e_Real:
-- 0*inf = 0, inf * inf = inf, X / inf = 0, |X| * -inf = -inf,
-- inf + inf = inf, X + inf = inf, -inf * inf = -inf, X - inf = -inf,
-- inf > X = True, -inf < X = True, inf > -inf = True,
-- (inf = -inf) = False
--
-- Constraint_Error:
--
-- The following ops have no sensible meaning, so Constraint_Error (ce) is
-- raised.
-- inf - inf => ce, inf / inf => ce, X / 0 => ce, inf / 0 => ce,
-- inf < inf => ce.
--
--***************************************************************************
-- SECTION II.
--
-- Standard arithmetic operators.
--
-- The arithmetic performed by these routines is supposed to be correct out
-- to the number of decimals specified by Desired_Decimal_Digit_Precision,
-- which you type in at the beginning of the spec.  But the arithmetic is
-- actually performed on digits well beyond this limit in order to guarantee
-- this level of accuracy. The values
-- held by these extra digits (guard digits) are usually almost correct.
-- None of the following operators rounds away these guard digits.
-- Rounding is done explicitly by calling Round_Away_Smallest_Guard_Digit().
-- In particular, none of the following comparison operators ("<", ">=", etc.)
-- rounds away guard digits of operands before performing the comparison.
-- All of them perform their comparisons out to the final guard digit.
-- To reduce much confusion, I decided to leave rounding entirely up to the
-- user, with Round_Away_....  Whether its best to round or not
-- depends on the algorithm.  In a lot of algorithms its better to
-- round before you use the "Are_Equal" operator, and better not to round
-- when you use the "<" and ">" operators.
--
--  "=" or Are_Equal(X, Zero) is the most efficient way to find out if X
--  is Zero. Same for Infinity. X < Zero is the efficient way find Sign of X.
--  X > Zero is efficient way to test positivity of X.  Zero < X and
--  Zero > X are also handled efficiently.
--
--***************************************************************************
-- SECTION III.
--
-- Routines for conversion from Real to e_Real and back again.
--
-- Makes it easy to write generics that can be instantiated with either
-- conventional floating point of this extended floating point.
-- The unary + and - are here to make it easier to convert programs from
-- ordinary floating point to extended, by making it easy to replace
--
--   X : Generic_Float_Type := 1.2345;     --here instantiate w/ Float.
--   X                      := 4.5678 * Y; --here instantiate w/ Float.
--
-- with
--
--   X : Generic_Float_Type := +1.2345; --here instantiate w/ e_Real or Float.
--   X                      := (+4.5678) * Y; --same here.
--
-- Now you can instantiate the generic with either e_Real or Float,
-- (but you have to add the unary "+" to the list of generic formals.)
--
--***************************************************************************
-- SECTION V.
--
-- Real * Extended operations.
--
-- More efficient operations.  The "Real" is not your ordinary real, but
-- something in the range 0.0 .. Radix-1, and integer valued, though it
-- can have an negative or positive exponent.  So its not very appropriate
-- for general use; 
--
-- The Real * Extended operations can be particularly efficient if
-- the Real number is in the same range as a Digit, ie, 0..Radix-1.
-- So we define a type e_Digit, a single real number with an
-- exponent.  These mixed multiplication "*" and "/" ops are used by the
-- ascii to real_extended and real_extended to ascii translators,
-- and by Newton's method calculations of elementary functions.
-- This efficiency only comes if the real number can be represented by
-- a single digit: integer values in the range 0..Radix-1, (times
-- an exponent in a power-of-2 Radix.  e.g. 0.5 is OK, 1.0/3.0 is not.)
-- Make_e_Digit will raise a constraint error if the range of the
-- intended real number is wrong.
--
--**********************************************************************
-- INTERNAL FORMAT OF e_Real
--
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
--                   Max
--  X = Radix**Exp * SUM {Radix**(-I-1) * Digit(I)}
--                   I=0
--
-- Exp is called the "normalized" exponent.  If Exp is the normalized exponent
-- then, say, a binary number would be written:
--
--   0.111011010001 * 2**(Exp).
--
-- In other words the first binary digit in the mantissa is of power 2**(-1).
-- It is important to know this because the function Real'Exponent(x) returns
-- the *normalized* exponent, and the function Real'Fraction(x) returns
-- x * 2**(-Exp) where Exp is the normalized exponent.  So in the above case,
-- 'Fraction would return 0.111011010001.
-- Also, in normalized form, the first binary digit of the mantissa is always
-- non-zero.
--***************************************************************************

generic

   type Real is digits <>;
   --  Make it 15 digits or greater.  This is checked.
   --  This is usually the type you are going to replace with e_Real.

package Extended_Real is

   pragma Pure (Extended_Real);

   pragma Assert (Real'Digits >= 15);

   type e_Real is private; -- The extended precision floating pt type.

   --  Instructions:
   --  The only things that need to be adjusted by the user are
   --
   --        Desired_Decimal_Digit_Precision
   --  and
   --        Desired_No_Of_Bits_In_Radix
   --
   --  The 2 parameters follow next, along with instructions.

   Desired_Decimal_Digit_Precision : constant := 28;
   --  If you request       28 Decimal Digits, you usually get 29 or more.
   --  If you request 29 to 37 Decimal Digits, you usually get 38 or more.
   --  If you request 38 to 46 Decimal Digits, you usually get 47 or more.
   --  If you request 47 to 55 Decimal Digits, you usually get 56 or more.
   --  (And so on, in jumps of 9.  Assumes Desired_No_Of_Bits_In_Radix = 30.)
   --
   --  The simple operators "*", "+", "/" usually give the best precision.
   --  They should get the 1st guard digit right, and by themselves:
   --  If you request 28 Decimal Digits, they're good to about 37 Decimal Digits.
   --  If you request 37 Decimal Digits, they're good to about 46 Decimal Digits.
   --
   --  Large complicated floating point computations will usually get both
   --  guard digits wrong and additional error will accumulate, so:
   --  If you request 28 Decimal Digits, ultimately expect <28 Decimal Digits.
   --
   --  Lower limit on Desired_Decimal_Digit_Precision is 28.

   pragma Assert (Desired_Decimal_Digit_Precision >= 28);


   Desired_No_Of_Bits_In_Radix : constant := 30;
   --  Anything under 31 works, but should be adjusted for best performance:
   --  30 is good if Desired_Decimal_Digit_Precision is 28 to 55.
   --  29 is good standard setting (use it when Desired_Decimal_Digit_Precision > 55).
   --  28 is good if Desired_Decimal_Digit_Precision >> 200. (But experiment.)
   --
   --  30 is necessary if you want the minimum decimal digits setting: 28.
   --  (If you choose 29 bits in Radix, you will get more decimals than you expect.)

   pragma Assert (Desired_No_Of_Bits_In_Radix <= 30);


   type e_Int is range -2**31+1 .. 2**31-1;
   subtype e_Integer is e_Int'Base;
   --  Type of Exponent.  Also takes the place of Universal_Integer in the
   --  "attribute" functions defined below for e_Real.
   --  Keep it 32-bit. Smallest usually fastest.

   pragma Assert (e_Integer'Size <= 32); 
   --  try fit e_Reals into small space; not essential, but Larger is slower.


   Zero              : constant e_Real;
   One               : constant e_Real;
   Positive_Infinity : constant e_Real;
   Negative_Infinity : constant e_Real;
   --  To test an arbitrary X : e_Real to see if X is Zero or infinity use the
   --  function:  Are_Equal (X, Zero), or X = Zero etc. Testing for Zero is fast.
   --  Infinity here means a finite number that is too large to represent in the
   --  floating point, and whose inverse is too small to represent in floating
   --  point Zero is positive.



   -- SECTION II.  Standard operators.
   --
   -- To reduce much confusion, rounding is entirely up to the
   -- user, with Round_Away_Guard_Digits().  Whether its best to round or not
   -- depends on the algorithm.  For example, in some cases it is better to
   -- round before you use the "Are_Equal" operator, and better not to round
   -- when you use the "<" and ">" operators.  (see intro.)


   function "*" (X, Y : e_Real) return e_Real; -- inline can slow it down.

   function "+" (X, Y : e_Real) return e_Real; -- inline can slow it down.

   function "-" (X, Y : e_Real) return e_Real;

   function "+" (X    : e_Real) return e_Real;
   function "-" (X    : e_Real) return e_Real;

   function "/" (X, Y : e_Real) return e_Real;

   function "**"(X : e_Real;
                 N : Integer)   return e_Real;

   procedure Square (X : in out e_Real);
   --  Same as X := X * X; (but usually faster if < 120 decimal digits.)

   procedure Mult 
     (X : in out e_Real;
      Y : in     e_Real);
   --  Same as X := X * Y; (but usually faster if < 120 decimal digits.)

   function "Abs" (X : e_Real)  return e_Real;

   function Are_Equal (X, Y : e_Real) return Boolean;
   --  Return true only if
   --  equality is exact in the cases of Zero and the 2 infinities.
   --  Are_Equal(X, Zero) is the most efficient way to find out if X is Zero.
   --  Same for Infinity. X < Zero is the efficient way find Sign of X.
   --  X > Zero is efficient way to test positivity of X.  Zero < X etc. OK too.

   function "<"  (X, Y : e_Real) return Boolean;
   function "<=" (X, Y : e_Real) return Boolean;
   function ">"  (X, Y : e_Real) return Boolean;
   function ">=" (X, Y : e_Real) return Boolean;
   function  "=" (X, Y : e_Real) return Boolean renames Are_Equal;

   function Are_Not_Equal (X, Y : e_Real) return Boolean; -- not Are_Equal

   -- SECTION III. Conversions between Real to e_Real. (see intro.)


   function Make_Real (X : e_Real) return Real;

   function Make_Extended (X : Real) return e_Real;
   function "+" (X : Real) return e_Real renames Make_Extended;
   function "-" (X : Real) return e_Real;
   --  The above 3 functions are identical, except "-" changes sign of X.
   --  Makes it easy to write generics that can be instantiated with either
   --  conventional floating point of this extended floating point, via:
   --    X : Generic_Float_Type := +1.2345;

   function "+" (X : Integer) return e_Real;
   --  Only works in range of Real (15 digits usually).
   --
   --     raises Constraint_Error
   --
   --  if X is greater than about 10**Real'Digits.
   --  So X = 2**62 raises Constraint_Error if Real'Digits = 15.
   --  Its really just for making e_Reals out of small ints: +7.


   -- SECTION IV.  Ada9X oriented attributes.
   --
   -- Below: Machine attributes and the function calls (like Truncation).
   -- More information on the machine model is given in the introduction.
   -- The machine model is static, so none of the Machine oriented attributes,
   -- and none of the functions reflect varying precision.  (see intro.)
   --
   -- Written in the spirit of the Ada attributes, but the fit is never
   -- perfect.

   function Remainder (X, Y : e_Real) return e_Real;

   function Copy_Sign (Value, Sign : e_Real) return e_Real;


   function e_Real_Machine_Rounds return Boolean;
   function e_Real_Machine_Overflows return Boolean;
   function e_Real_Signed_Zeros return Boolean;
   function e_Real_Denorm return Boolean;
   --  These functions always return False.

   function e_Real_Machine_Emax return e_Integer;

   function e_Real_Machine_Emin return e_Integer;

   function e_Real_Machine_Mantissa return e_Integer;
   --  Always returns  Mantissa'Length: all the digits including guards.
   --
   --  NOT binary digits, NOT decimal digits.


   function e_Real_Machine_Radix return Real;
   --  Usually 2.0**29 or 2.0**30 for Integer digits; 2.0**24 for Flt pt digits.
   --  Returns: No_of_Bits_in_Radix  (as a Real type).

   function Leading_Part (X            : e_Real;
                          Radix_Digits : e_Integer) return e_Real;
   --  Example: to set to zero all but the first digit of X  use
   --  First_Digit := Leading_Part (X, 1);

   function Exponent (X : e_Real)  return e_Integer;
   --  By convention return 0 for Zero. Else return nomalized Expon.
   --  Returns Max_Exponent+2 for the 2 infinities.
   --  NOT decimal, and NOT binary Exponent.

   function Fraction (X : e_Real)  return e_Real;

   function Compose (Fraction : e_Real;
                     Exponent : e_Integer) return e_Real;

   function Scaling (         X : e_Real;
                     Adjustment : e_Integer) return e_Real;


   --  Chop off fractional parts.
   --
   --  Rounding, Unbiased_Rounding, Ceiling, Floor return e_Real
   --  with Zero fractions.

   function Rounding (X : e_Real)  return e_Real;

   function Unbiased_Rounding (X : e_Real) return e_Real;

   function Truncation (X : e_Real) return e_Real;

   function Ceiling (X : e_Real) return e_Real;

   function Floor (X : e_Real) return e_Real;


   --  Round away guard digits.  
   --
   --  function Machine rounds away the smallest Guard digit.
   --  There's no one right way to round away Guard Digits or choose
   --  Model_Epsilon's. Doesn't follow the Ada95 model for rounding
   --  64-bit floats. That model doesn't seem to fit very well.

   function e_Real_Model_Epsilon return e_Real;
   --  At present this calls:  e_Real_Model_Epsilon_2 which is
   --  1 unit in the 3rd smallest digit. (The 3rd smallest digit
   --  is the 1st digit that is larger than the 2 guard digits.)

   function e_Real_Machine_Epsilon return e_Real;
   --  At present this calls:  e_Real_Model_Epsilon_1 which is
   --  1 unit in the 2nd smallest digit. (The 2nd smallest digit
   --  is the larger of the 2 guard digits.)

   function Machine (X : e_Real) return e_Real;
   --  This calls:  
   --               Round_Away_Smallest_Guard_Digit
   --

   function Round_Away_Smallest_Guard_Digit (X : e_Real) return e_Real;

   function e_Real_Model_Epsilon_1 return e_Real;
   --  One unit in the 2nd smallest digit.

   function e_Real_Model_Epsilon_2 return e_Real;
   --  One unit in the 3rd smallest digit. (The smallest digit that
   --  is *not* a Guard_Digit.
   -- 
   --  Guard_Digits = 2 always; assume neither of them is correct:
   --  if there's 3 digits of Radix 2^30 then eps_2 is 2^(-30).
   --  if there's 4 digits of Radix 2^30 then eps_2 is 2^(-60).
   --  if there's 5 digits of Radix 2^30 then eps_2 is 2^(-90) or about 10**-27.
   --
   --  So Eps_2 is the smallest number s/t eps_2+.999999999999 /= .999999999999  
   --  when you remove both guard digits. 


   -- SECTION V.  Digit * Extended operations.
   --
   -- More efficient operations.  "Digit" is not your ordinary real, but
   -- something in the range 0.0 .. Radix-1.0, and integral valued, though it
   -- can have a negative exponent.  So the following is not very appropriate
   -- for general use; in the '83 version we export it so that it can be used by
   -- elementary function packages.


   type e_Digit is private;

   function "*" (X : e_Digit; Y : e_Real) return e_Real;

   function "/" (X : e_Real; Y : e_Digit) return e_Real;

   function Sum_Of (X : e_Digit; Y : e_Real) return e_Real;

   function "+" (X : e_Digit; Y : e_Real) return e_Real renames Sum_Of;

   function Scaling (X : e_Digit; Adjustment : e_Integer) return e_Digit;
   --  Multiply X by Radix**N where N = Adjustment.

   function Make_Extended (X : e_Digit) return e_Real;

   function Make_e_Digit (X : Real) return e_Digit;
   --  X must be a whole number: 0.0, 1.0, 2.0 etc. in the range 0..Radix-1,
   --  times some integer power of the Radix.  So 0.5 is OK, but not 1/3.

   function Number_Of_Guard_Digits return e_Integer;
   --  Constant.  To get number of digits that are being correctly calculated
   --  (by conservative estimates) use
   --  No_Correct_Digits = Present_Precision - Number_Of_Guard_Digits.

   function Minimum_No_Of_Digits_Allowed return e_Integer;
   --  Constant.  Includes guard digits.


private

   --
   -- SECTION VII.  Make the Data structure for e_Real.
   --

   --  Using 32-bit ints for the Digits: (don't do it)

   --No_Of_Usable_Bits_In_Digit : constant := 31; -- bad idea; lots of trouble.
   --No_Of_Bits_In_Radix        : constant := 13; -- can't use 14 or >

   --  Using 64-bit floats for the Digits: (don't bother)

   --No_Of_Usable_Bits_In_Digit : constant := 53; -- if using flt pt Mantissa (slow)
   --No_Of_Bits_In_Radix        : constant := 24;

   --  Using 64-bit ints for the Digits:

   No_Of_Usable_Bits_In_Digit : constant := 63; -- Integer; must allow neg. vals
   No_Of_Bits_In_Radix        : constant := Desired_No_Of_Bits_In_Radix;

   --  30 is good if Desired_Decimal_Digit_Precision is 28 to 55.
   --  29 is good standard setting (especially: Desired_Decimal_Digit_Precision > 55).
   --  28 is good if Desired_Decimal_Digit_Precision is in the 100's. (But experiment.)

   Sums_Per_Carry : constant := 2**(No_Of_Usable_Bits_In_Digit-2*No_Of_Bits_In_Radix)-1;
   --  Sums_Per_Carry : This is number of sums you can accumulate during
   --  multiplication before the carrys need to be performed.
   --
   --  You can do a large number of X*Y < Radix*Radix products, and then sum
   --  (Sums_Per_Carry+1) of them before a Carry is necessary in multiplication.

   -- a implies b    is same as   not (a) or b:
   pragma Assert (not (No_Of_Bits_In_Radix = 30) or Sums_Per_Carry <= 8-1);
   pragma Assert (not (No_Of_Bits_In_Radix = 29) or Sums_Per_Carry <= 32-1);
   pragma Assert (not (No_Of_Bits_In_Radix = 28) or Sums_Per_Carry <= 128-1);
   pragma Assert (not (No_Of_Bits_In_Radix = 27) or Sums_Per_Carry <= 512-1);
   pragma Assert (not (No_Of_Bits_In_Radix = 26) or Sums_Per_Carry <= 2048-1);

   --
   -- Now that we know:   No_Of_Bits_In_Radix,
   --
   -- get number of binary digits and extended digits needed to make e_Real.
   -- Use the following formula for the number of Binary digits needed
   -- to meet Desired Decimal Digit precision D:
   --
   --            Binary_Digits  >=  Ceiling (D * Log_Base_2_Of_10) + 1
   --
   -- where  D = Desired_Decimal_Digit_Precision, and
   -- where  Log_Base_2_Of_10 = 3.321928094887362.
   -- Ceiling of Real numbers with static declarations? 
   -- 
   --       Ceiling (3.321928094887 * D) <= Ceiling (3.322  * D)
   --                                     = Ceiling((3322.0 * D) / 1000.0)
   --                                     = (3322 * D - 1) / 1000 + 1
   -- D is integer valued, so use integer Ceiling (A / B) = (A - 1) / B + 1.
   -- (for positive A). The above
   -- steps give us the number of binary digits required: No_Of_B_Digits.
   -- Next: min number of Radix 2.0**No_Of_Bits_In_Radix digits: No_Of_e_Digits.
   -- To get No_Of_e_Digits divide by No_Of_Bits_In_Radix and take the Ceiling.
   --

   -- B is for binary, E for extended:

   ILog_Base_2_Of_10_x_1000 : constant := 3322; -- Round UP.
   D              : constant := Desired_Decimal_Digit_Precision - 2;
   No_Of_B_Digits : constant := (ILog_Base_2_Of_10_x_1000 * D - 1) / 1000 + 2;
   No_Of_e_Digits : constant := (No_Of_B_Digits - 1) / No_Of_Bits_In_Radix + 1;

   --
   -- The following parameter settings give us 2 more words in the Mantissa
   -- than required.  These two are essential in getting the full desired
   -- precision in extensive floating pt calculation, and also in IO, and in
   -- functions that are evaluated by Newton's method.
   -- (At least one such guard digit is essential anyway, to compensate for
   -- leftward shift of the mantissa during normalization.)
   -- The index of the digits is a subtype of the Exponent type because
   -- there are frequent conversions between the two.
   --
   -- An assertion verifies that there are 2 guard digits.
   --
   -- e_Integer is used as the type of the index of the extended digits
   -- because e_Integer is the type of the exponent (defined below).
   -- There's a close relationship between the exponent of the number and the
   -- index of the digits of the number.  (They are often scaled
   -- simultaneously, there is a relationship between their ultimate ranges,
   -- so they are given the same type here.)
   --    Log_Base_2_Of_10 : constant := 3.321928094887362;
   --

   No_Of_Guard_Digits : constant := 2;
   --  Guard_Digits are extra digits of precision at the end of the mantissa.
   --
   --  The 2nd Guard_Digit makes
   --  the Elementary Math Functions full precision (or almost full).
   --  Also the IO routines need 2 Guard_Digits.

   pragma Assert (No_Of_Guard_Digits = 2);


   --  The following are not decimal digits.

   Ultimate_Correct_Digit : constant e_Integer := No_Of_e_Digits - 1;
   Ultimate_Digit : constant e_Integer := No_Of_Guard_Digits + Ultimate_Correct_Digit;

   subtype Digits_Base is e_Integer range 0..Ultimate_Digit+1;
   subtype Digit_Index is Digits_Base range 0..Ultimate_Digit;

   pragma Assert (Digit_Index'First = 0); 
   -- some of the arithmetic in "+" assumes this.

   --  The following are not decimal digits.

   Min_No_Of_Correct_Digits : constant := 3;
   --  Things stop working if this is less than 3.

   Min_No_Of_Digits : constant := Min_No_Of_Correct_Digits + 2;
   --  The 2 is the min number of guard digits.

   pragma Assert (Ultimate_Digit >= Min_No_Of_Digits - 1);
   pragma Assert (Ultimate_Digit >= Min_No_Of_Correct_Digits+No_Of_Guard_Digits-1);
   --  Digits go from  0..Ultimate_Digit


   Max_Exponent : constant e_Integer :=  2 ** (e_Integer'Size - 5);
   Min_Exponent : constant e_Integer := -Max_Exponent;
   --  The exponent is usually 16 or 32 bit int.  Limits on its range are set
   --  below what the base type allows: no more than 1/4 the dynamic range
   --  of the base type. If we use 1/8 of that limit, it allows us to delay
   --  overflow check to end of most routines (except "**"). If we use 1/32
   --  of that limit, it allows us to do IO more simply. So to make
   --  IO work, at present the requirement is 2 ** (e_Integer'Size - 5).

   pragma Assert (Max_Exponent <= 2 ** (e_Integer'Size - 5));

   --  function  e_Real_Machine_Emin    returns    Min_Exponent
   --  function  e_Real_Machine_Emax    returns    Max_Exponent


   --subtype Digit_Type is Real;
   --  Can use Floats with 53 bit mantissas as Digit_Type.  Make 2 changes
   --  above (search for No_Of_Usable_Bits_In_Digit and follow instructions)
   --  and 2 changes in body (compiler will tell you where). Also comment
   --  out next 3 statements.  Amazingly, it worked nicely last time I did it.
   --  Its slow, and it only makes sense when 64 bit ints are bad or missing.


   type D_Type is range -2**63+1 .. 2**63-1;
   subtype Digit_Type is D_Type'Base;
   --  Must allow negative digits. Use 64 bit Integer.

   pragma Assert (Digit_Type'Last = 2**(Digit_Type'Size-1)-1);
   pragma Assert (Digit_Type'Size-1 >= No_Of_Usable_Bits_In_Digit);

   Digit_Zero  : constant Digit_Type := Digit_Type (0);
   Digit_One   : constant Digit_Type := Digit_Type (1);
   Digit_Two   : constant Digit_Type := Digit_Type (2);

   Digit_Radix         : constant Digit_Type := Digit_Two**No_Of_Bits_In_Radix;
   Half_Radix          : constant Digit_Type := Digit_Two**(No_Of_Bits_In_Radix-1);
   Digit_Radix_Squared : constant Digit_Type := Digit_Radix * Digit_Radix;
   Digit_Radix_Minus_1 : constant Digit_Type := Digit_Radix - Digit_One;

   type Mantissa is array (Digit_Index) of Digit_Type;

   type e_Real is record
      Digit       : Mantissa  := (others => Digit_Zero);
      Exp         : e_Integer := 0;
      Is_Zero     : Boolean   := True;
      Is_Positive : Boolean   := True;
      Is_Infinite : Boolean   := False;
   end record;

   --for e_Real'Size use (Digit_Type'Size*Mantissa'Length + e_Integer'Size*2);
   --  Make e_Real'Size Integer number of 64-bit words. Usually doesn't matter.
   --  Only for integer Digit_Type. Comment out for Float. pt. Digit types.

   Zero : constant e_Real
         := e_Real' ((others => Digit_Zero), 0, True, True, False);

   One : constant e_Real
         := e_Real' ((0 => Digit_One, others => Digit_Zero), 0, False, True, False);

   Positive_Infinity : constant e_Real
         := e_Real' ((others => Digit_Zero), Max_Exponent+4, False, True, True);

   Negative_Infinity : constant e_Real
         := e_Real' ((others => Digit_Zero), Max_Exponent+4, False, False, True);


   -- For efficiency, we need an optimized (Real * Extended)
   -- operation.  So define type e_Digit, a single real number with
   -- an exponent.  Its a real number that's restricted  to integral values
   -- in the range to 0..Radix-1.

   type e_Digit is record
      Digit       : Digit_Type := Digit_Zero;
      Exp         : e_Integer  := 0;
      Is_Zero     : Boolean    := True;
      Is_Positive : Boolean    := True;
   end record;

   -- Constants used in body. Real is used for easy communication with e_Real.

   Real_Zero       : constant Real := Real (0.0);
   Real_One        : constant Real := Real (1.0);

   Real_Radix      : constant Real := 2.0**No_Of_Bits_In_Radix;
   Radix_Minus_1   : constant Real := 2.0**No_Of_Bits_In_Radix - 1.0;
   Radix_Squared   : constant Real := 2.0**(2*No_Of_Bits_In_Radix);
   Inverse_Radix   : constant Real := 2.0**(-No_Of_Bits_In_Radix);
   Inverse_Radix_Squared : constant Real := Inverse_Radix * Inverse_Radix;

end Extended_Real;
