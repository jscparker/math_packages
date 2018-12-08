--
-- Simple tests of Extended_Real.  Search for catastrophic
-- errors or portability problems.
--
with Extended_Real;
with text_io; use text_io;
procedure e_real_tst_1 is

  type Real is digits 15; 
  
  package Extend is new Extended_Real (Real);
  use Extend;
  package rio is new Text_IO.float_io(Real);
  use rio;
  package iio is new Text_IO.Integer_io(E_Integer);
  use iio;
  
  Radix                : constant Real := E_Real_Machine_Radix;  
  Radix_Minus_1        : constant Real := Radix - 1.0;
  Radix_Minus_1_Squ    : constant Real := Radix_Minus_1 * Radix_Minus_1; 
  
  Z1, Z2, Z3, Z4 : E_Real;
  
  --  Attempt to make Junk constants of order one with non-zero's all over.
  Junk1 : constant E_Real := One / (+3.0);
  Junk2 : constant E_Real := (+17.0) / (+19.0);
  
  Test_Vector_Seed : Real;
  
   --**********************************************************************
  procedure Test_Rem (Test_Vector_Seed : Real) is 
    Dig, Digit : Real;
    Min, DeltaExp   : E_Integer;
    Difference1 : E_Real;
    Z1, Z2 : E_Real := Zero;
    Z4 : E_Real;
  begin
    Min   := E_Integer'Last;
    Dig   := 0.0;
    
   new_line;
    for I in 1..100 loop
      Z1 := Z1 + Junk1;
      Z2 := Z2 + (+Test_Vector_Seed) * Junk2;
      
      -- try to reproduce Z1:
      
      Z4 := Remainder (Z1, Z2) + Z2 * Unbiased_Rounding (Z1 / Z2);
      -- Z4 should equal Z1.
      
      Difference1 := (Z4 - Z1) / Z4;

      --put (Make_Real (Remainder (Z1, Z2) / Z2)); new_line;
      --   should always be in range [-0.5..0.5]
      --put (Make_Real ((Z2*Unbiased_Rounding (Z1 / Z2)) / Z1)); New_Line;

      IF Are_Not_Equal (Difference1, Zero) THEN
       
        DeltaExp   := Exponent(Difference1);
        Digit      := Make_Real (Fraction (Leading_Part (Difference1,1)));
        IF Abs(DeltaExp) < Abs(Min) THEN
           Min := DeltaExp;
           Dig := Digit;
        END IF;
      
      END IF;
      
   end loop;
   
   IF Min = E_Integer'Last THEN
      Dig := 0.0;
      Min := -E_Integer'Last + 1;
   END IF;
      
   -- The above uses the normalized Exponent, which assumes a fractional 1st
   -- Digit.  Below we (effectively) multiply the first non-zero digit by Radix,
   -- and subtract 1 from the "power of Radix" exponent.
   put ("Approximate error ="); 
   put (Dig*Radix); put (" * Radix**("); put (E_Integer'Image(Min-1)); put(")");
   new_line;
   
  end Test_Rem;

  
  procedure Test_Digit_Mult_and_Div (Test_Vector_Seed : Real) is 
    Dig, Digit : Real;
    Min, DeltaExp   : E_Integer;
    Difference1 : E_Real;
    Real_Digit : Real := Radix + 1999.0;
    Digit1 : E_Digit;
    Z1, Z2 : E_Real := Zero;
  begin
    Min  := E_Integer'Last;
    Z1   := (+Test_Vector_Seed) * Junk1;
    
    for I in 1..1000 loop

      Real_Digit := Real_Digit - 2000.0;
      Digit1 := Make_E_Digit (Real_Digit);

      Z2          := Digit1 * (Z1 / Digit1);
      Difference1 := (Z1 - Z2) / Z1;

      IF Are_Not_Equal (Difference1, Zero) THEN
       
        DeltaExp   := Exponent(Difference1);
        Digit      := Make_Real (Fraction (Leading_Part (Difference1,1)));
        IF Abs(DeltaExp) < Abs(Min) THEN
           Min := DeltaExp;
           Dig := Digit;
        END IF;
      
      END IF;
      
   end loop;
   
   IF Min = E_Integer'Last THEN
      Dig := 0.0;
      Min := -E_Integer'Last + 1;
   END IF;
      
   -- The above uses the normalized Exponent, which assumes a fractional 1st
   -- Digit.  Below we (effectively) multiply the first non-zero digit by Radix,
   -- and subtract 1 from the "power of Radix" exponent.
   put ("Approximate error ="); 
   put (Dig*Radix); put (" * Radix**("); put (E_Integer'Image(Min-1)); put(")");
   new_line;
   
  end Test_Digit_Mult_and_Div;
  
  -----------------------
  -- Test_Mult_and_Add --
  -----------------------

   -- uses Newton's method to calculate 1/(X).  Square and invert to
   -- compare with X; square and multiply with with X to compare w. 1.0.

  procedure Test_Mult_and_Add (Test_Vector_Seed : Real) is 
    Dig, Digit : Real;
    Min, DeltaExp   : E_Integer;
    Difference1 : E_Real;
    
    --
    function Reciprocal (X : E_Real) return E_Real is
      X_isqr     : E_Real;
      X_start    : constant Real := 1.0 / Make_Real(X);
      Iterations : constant Integer := 9;
    begin
      X_isqr := Make_Extended (X_Start);
      for I in 1..Iterations loop
        X_isqr := X_isqr + (One - X * X_isqr) * X_isqr;
      end loop;
      return X_isqr;
    end Reciprocal;
    
  begin
    Min  := E_Integer'Last;
    Dig  := 0.0;
    Z1 := Zero;
    
    for I in 1..20 loop
      Z1 := Z1 + (+Test_Vector_Seed) * Junk1;
      
      Z2          := Reciprocal (Z1);
      Difference1 := One - Z1 * Z2;

      IF Are_Not_Equal (Difference1, Zero) THEN
       
        DeltaExp   := Exponent(Difference1);
        Digit      := Make_Real (Fraction (Leading_Part (Difference1,1)));
        IF Abs(DeltaExp) < Abs(Min) THEN
           Min := DeltaExp;
           Dig := Digit;
        END IF;
      
      END IF;
      
   end loop;
   
   IF Min = E_Integer'Last THEN
      Dig := 0.0;
      Min := -E_Integer'Last + 1;
   END IF;
      
   put ("Approximate error ="); 
   put (Dig*Radix); put (" * Radix**("); put (E_Integer'Image(Min-1)); put(")");
   new_line;
   
  end Test_Mult_and_Add;

   --**********************************************************************
  procedure Test_Mult_And_Div (Test_Vector_Seed : Real) is 
    Dig, Digit : Real;
    Min, DeltaExp   : E_Integer;
    Difference : E_Real;
  begin
    Min  := E_Integer'Last;
    Dig  := 0.0;
    for I in 1 .. 10_000 loop
      Z1  := Z1 + Junk2 * (+Test_Vector_Seed);
      Z2  := (Z1 * (One / Z1));
      Difference := One - Z2;
      DeltaExp   := Exponent(Difference);
      Digit      := Make_Real (Fraction (Leading_Part (Difference,1)));
      IF Abs(DeltaExp) < Abs(Min) THEN
         Min := DeltaExp;
         Dig := Digit;
      END IF;
   end loop;
   
   put ("Approximate error =");
   put (Dig*Radix); put (" * Radix**("); put (E_Integer'Image(Min-1)); put(")");
   new_line;
  end Test_Mult_And_Div;
  
  --**********************************************************************
  procedure Test_Make_Extended (Test_Vector_Seed : Real) is 
    Max : Real;
    Difference   : Real;
    Junk1, Junk2 : Real;
  begin
    Junk1 := 0.0;
    Max   := 0.0;
    for I in 1..1000 loop
      Junk1   := Junk1 + Test_Vector_Seed;
      Junk2  := Make_Real (Make_Extended (Junk1));
      Difference := Junk2 - Junk1;
      IF Abs(Difference) > Max THEN
         Max := Difference;
      END IF;
   end loop;
   put ("Max error in (X - Make_Real (Make_Extended(X))) = "); 
   put (Max);
   new_line;
  end Test_Make_Extended; 
   
  --**********************************************************************
begin

  -- Simple test of "*" and "-":
  Z1 := Make_Extended (Radix_Minus_1);
  Z2 := Z1 * Z1;
  Z3 := Z1 * Z2;
  Z4  := Z2 - Make_Extended (Radix_Minus_1_Squ);
  put("  Simple test of * and - : ");  new_line;
  put(Make_Real(Z1)); put("  This should be Radix - 1.");  new_line;
  put(Make_Real(Z2)); put("  This should be Radix - 1 squared."); new_line; 
  put(Make_Real(Z3)); put("  This should be Radix - 1 cubed."); new_line;
  put(Make_Real(Z4)); put("  This should be near zero, but not quite."); new_line; 
  new_line;
  
  
  new_line;
  put("  Simple test of Remainder:");
  new_line;
      
  Test_Vector_Seed := 1.2345678912345678E+77;
  Test_Rem (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+14;
  Test_Rem (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+8;
  Test_Rem (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-8;
  Test_Rem (Test_Vector_Seed);
    
  Test_Vector_Seed := 1.2345678912345678E-14;
  Test_Rem (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-35;
  Test_Rem (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-77;
  Test_Rem (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-123;
  Test_Rem (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-171;
  Test_Rem (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-201;
  Test_Rem (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-231;
  Test_Rem (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-271;
  Test_Rem (Test_Vector_Seed);
  
  
  new_line;
  put("  Simple test of mixed digit, extended operations:");
  new_line;
  
  Test_Vector_Seed := 1.2345678912345678E+8;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  new_line;
  put("  Simple test of + and *:");
  new_line;
    
  Test_Vector_Seed := 1.2345678912345678E+131;
  Test_Mult_and_Add (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+81;
  Test_Mult_and_Add (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+31;
  Test_Mult_and_Add (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+8;
  Test_Mult_and_Add (Test_Vector_Seed);
    
  Test_Vector_Seed := 1.2345678912345678E-8;
  Test_Mult_and_Add (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-31;
  Test_Mult_and_Add (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-81;
  Test_Mult_and_Add (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-131;
  Test_Mult_and_Add (Test_Vector_Seed);
  

  --**********************************************************************
  new_line; put ("Maximum available number of Digits is: ");
  put (e_Real_Machine_Mantissa); new_line;

  new_line(2); put ("Some tests of +,*: "); new_line;
  Test_Vector_Seed := 1.2345678912345678E-4;
  for I in 1..10 loop
     Test_Vector_Seed := Test_Vector_Seed + 1.2345678912345678E-4;
     Test_Mult_and_Add (Test_Vector_Seed);
  end loop;
  Test_Vector_Seed := 1.2345678912345678E+4;
  for I in 1..10 loop
     Test_Vector_Seed := Test_Vector_Seed * 1.2345678912345678E+4;
     Test_Mult_and_Add (Test_Vector_Seed);
  end loop;

  new_line(2); put ("Some tests of /,*: "); new_line;
  Test_Vector_Seed := 1.2345678912345678E-4;
  for I in 1..10 loop
     Test_Vector_Seed := Test_Vector_Seed + 1.2345678912345678E-4;
     Test_Mult_and_Div (Test_Vector_Seed);
  end loop;
  Test_Vector_Seed := 1.2345678912345678E+4;
  for I in 1..10 loop
     Test_Vector_Seed := Test_Vector_Seed * 1.2345678912345678E+4;
     Test_Mult_and_Div (Test_Vector_Seed);
  end loop;

  new_line(2); put ("Some tests of Make_Extended: "); new_line;
  Test_Vector_Seed := 1.2345678912345678E-4;
  for I in 1..10 loop
     Test_Vector_Seed := Test_Vector_Seed + 1.2345678912345678E-4;
     Test_Make_Extended (Test_Vector_Seed);
  end loop;
  Test_Vector_Seed := 1.2345678912345678E+4;
  for I in 1..10 loop
     Test_Vector_Seed := Test_Vector_Seed * 1.2345678912345678E+4;
     Test_Make_Extended (Test_Vector_Seed);
  end loop;
   
  --**********************************************************************   
  new_line;
  put("  Simple test of /:");
  new_line;
    
  Test_Vector_Seed := 1.2345678912345678E+31;
  Test_Mult_And_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+8;
  Test_Mult_And_Div (Test_Vector_Seed);
    
  Test_Vector_Seed := 1.2345678912345678E-8;
  Test_Mult_And_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-31;
  Test_Mult_And_Div (Test_Vector_Seed); 
  
  new_line;
  put("  Simple test of Make_Extended:");
  new_line;
  Test_Vector_Seed := 1.234567891234567891E+31;
  Test_Make_Extended (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.234567891234567891E+8;
  Test_Make_Extended (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.234567891234567891E-8;
  Test_Make_Extended (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.234567891234567891E-31;
  Test_Make_Extended (Test_Vector_Seed);

  --**********************************************************************   
  new_line;
  put("  Simple test of exponentiation:");  new_line;
  Z1 := (Make_Extended (7.2345)) ** (-277); 
  put (Make_Real(Z1)); put(" should be: "); put (7.2345**(-277)); new_line; 
  Z1 := (Make_Extended (7.2345)) ** 277; 
  put (Make_Real(Z1)); put(" should be: "); put (7.2345**277); new_line; 
  Z1 := (Make_Extended (1.2345)) ** 177; 
  put (Make_Real(Z1)); put(" should be: "); put (1.2345**177); new_line; 
  Z1 := (One + One + One) ** 97; 
  put (Make_Real(Z1)); put(" should be: "); put (3.0**97); new_line; 
  Z1 := (One + One) ** 67; 
  put (Make_Real(Z1)); put(" should be: "); put (2.0**67); new_line; 
  Z1 := (One + One) ** (-67); 
  put (Make_Real(Z1)); put(" should be: "); put (2.0**(-67)); new_line; 
  Z1 := (One + One) ** (0); 
  put (Make_Real(Z1)); put(" should be: "); put (2.0**(0)); new_line(1); 

  new_line;
  put("  Quick test of Machine (rounds away smallest guard digit) and Model_Epsilon:");
  new_line;
  Z1 := (One + e_Real_Model_Epsilon) - One;
  Z1 := Z1 / e_Real_Model_Epsilon;
  put(" should be 1.0:"); put (Make_Real(Z1)); new_line;

  Z1 := Machine (One + e_Real_Machine_Epsilon) - One;
  put(" should be 0.0:"); put (Make_Real(Z1)); new_line;

  Z2 := Make_Extended(0.99999999999999);
  Z1 := Machine (Z2 + e_Real_Model_Epsilon) - Z2;
  Z1 := Z1 / e_Real_Model_Epsilon;
  put(" should be 1.0:"); put (Make_Real(Z1)); new_line;

end;
