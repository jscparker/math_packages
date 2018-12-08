
-- Simple test of "*", "/", "+", and "-" using random arguments.
-- Doesn't test endpoints very well.

with Extended_Real;
with Extended_Real.e_Rand; 
with Extended_Real.IO; 
with Text_IO; use Text_IO;

procedure  e_real_tst_2 is

  type Real is digits 15; 
  
  package ext is new Extended_Real (Real);
  use ext;
  package rnd is new ext.e_Rand;
  use rnd;
  package eio is new ext.IO;
  use eio;
  package rio is new Text_IO.Float_IO(Real);
  use rio;
  package iio is new Text_IO.Integer_IO(e_Integer);
  use iio;
  
  X, Y, Z1            : e_Real;
  Delta_X, Big_Digits : e_Real;
  
  Frac_Part  : Real;  
  Exp_Part   : e_Integer;  
  Frac_Error : Real;  
  Exp_Error  : e_Integer; 
  Exp_First, Exp_Last : e_Integer; 
  Y_Exp_First, Y_Exp_Last : e_Integer; 
  Half_Digit : constant e_Digit := Make_e_Digit (0.5); 
  X_min : constant e_Real 
                := Scaling (Make_Extended(1.0), e_Real_Machine_Emin);
     
  type Integer32 is range -2**31+1..2**31-1;
  
  Mult_Limit : constant Integer32 := 100_000_000; -- usually > N * 10^6
  --  Number of iterations of multiplication div test.  The +/- tests
  --  are 8 times this numbers.
  
  Some_Seed : constant Integer := 9795178;
  --  Start at different rand stream.
  
  ---------
  -- Min --
  ---------
  
  function Min (X, Y : e_Integer) return e_Integer is
  begin
     return e_Integer'Min (X, Y);
  end;
  
  ---------
  -- Max --
  ---------
  
  function Max (X, Y : e_Integer) return e_Integer is
  begin
     return e_Integer'Max (X, Y);
  end;
  
  --------------------
  -- Get_Random_Exp --
  --------------------
  
  --  Make Random e_Integer in range Exp_First..Exp_Last
  --
  function Random_Exp 
    (Exp_First : in e_Integer;
     Exp_Last  : in e_Integer)
     return e_Integer
  is              
     Exp : e_Integer;
     X : Real;
  begin
    
     --  returns random Real in range [0, 1):

     X := Real (rnd.Next_Random_Int mod 2**7) * 2.0**(-7);

     Exp := Exp_First + e_Integer (X * (Real (Exp_Last) - Real (Exp_First)));

     return Exp;
     
  end;
  
  -------------
  -- Get_999 --
  -------------
  
  --  Make an e_real full of digits that have the max value that any
  --  any digit can have.
  --
  function Get_999 return e_Real is
 
     Max_Digit  : constant e_Digit   := Make_e_Digit (e_Real_Machine_Radix-1.0);
     Digit_Last : constant e_Integer := e_Real_Machine_Mantissa;
     
     Next_Digit : e_Digit;
     Delta_Exp  : e_Integer;
     Result     : e_Real;     -- Init to 0 important
     
  begin
  
     for I in e_Integer range 0..Digit_Last loop
     
        Delta_Exp  := -I - 1;
        Next_Digit := Scaling (Max_Digit, Delta_Exp);
        Result     := Next_Digit + Result;
        
     end loop;
      
     return Result;
     
  end;
  
  --------------------------
  -- Get_Normalized_Delta --
  --------------------------
  
  --  Want something like Delta_X / (X + X_Min), but want to avoid
  --  extended division and avoid possible exceptions.  Here X_min is 
  --  say about Radix**(e_Real_Machine_Emin-1). 
  --  RETURNS Abs (Frac_Part).
  --
  procedure Get_Normalized_Delta 
    (X : in e_Real;
     Delta_X : in e_Real;
     Frac_Part : out Real;
     Exp_Part  : out e_Integer)
  is
     X_tmp : e_Real := X; 
  begin
  
     IF (Are_Equal (Delta_X, Zero)) THEN
        Frac_Part := 0.0;
        Exp_Part  := e_Integer'First;
        return;
     END IF;

     --  OK, so Delta_X is not Zero.
     --  Find out if X is:
     
     IF (Are_Equal (X, Zero)) THEN
        X_tmp := X_min;
     END IF;
     
     --  OK, so Delta_X and X are not Zero.
     --  Now attempt an approximation of Delta_X / (X + X_min):
  
     Exp_Part  := Exponent (Delta_X) - Exponent (X_tmp);
     Frac_Part := Abs (Make_Real (Fraction (Delta_X)) 
                     / Make_Real (Fraction (X_tmp)));
     
  end Get_Normalized_Delta;

  -----------------------
  -- Print_e_Real_data --
  -----------------------

  --  Want a fast and very approximate translation to decimal, without
  --  linking to IO package.

  procedure Print_e_Real_data
    (Frac_Part : in Real;
     Exp_Part  : in e_Integer)
  is
     Decimal_Digits_in_Mantissa : constant Real := 
        Real (e_Real_Machine_Mantissa*Desired_No_Of_Bits_In_Radix)/3.322;
     Exp_02, Exp_10  : e_Integer;
     Frac_Correction : Real;
  begin

     -- go from radix**Exp_Part to 2**(30*Exp_10 + Exp_02):

     if Abs Exp_Part >= (e_Integer'Last-1) / 30 then -- Bits_In_Radix<=30 always.
       Exp_10 := 0;
       Exp_02 := 0;
     else
       Exp_10 := (Exp_Part * Desired_No_Of_Bits_In_Radix)  /  30;
       Exp_02 := (Exp_Part * Desired_No_Of_Bits_In_Radix) rem 30;
       if Exp_02 < 0 then
         Frac_Correction := 0.5 ** (Integer (Abs Exp_02));
       else
         Frac_Correction := 2.0 ** (Integer (Abs Exp_02));
       end if;
     end if;

     put(Frac_Part * Frac_Correction, aft => 3); 
     put(" * "); put("10**("); put(e_Integer'Image (Exp_10*9)); put(")");


     put("   or:  ");
     put(Frac_Part, aft => 3); put(" * ");  
     put("Radix**("); put(e_Integer'Image (Exp_Part)); put(")");

  exception 
     when others =>
     put("  Some printing error here."); new_line;

  end Print_e_Real_data;

  ----------------------------------
  -- Print_Extended_Real_Settings --
  ----------------------------------
 
  procedure Print_Extended_Real_Settings
  is
     Decimal_Digits_in_Mantissa : constant Real := 
        Real (e_Real_Machine_Mantissa*Desired_No_Of_Bits_In_Radix)/3.322;
  begin
     new_line(1); 
     put ("    Desired_Decimal_Digit_Precision =");
     put (Integer'Image(Desired_Decimal_Digit_Precision)); 
     new_line(2); 
     put ("Currently using "); 
     put (Integer'Image(Desired_Decimal_Digit_Precision));
     put ("  (or more) decimal digits of precision.");
     new_line(2); 
     put ("Number of digits in use (including 2 guard digits):"); 
     put (e_Integer'Image(e_Real_Machine_Mantissa));
     new_line(1); 
     put ("The digits are not decimal; they have Radix:        2**(");
     put (e_Integer'Image(Desired_No_Of_Bits_In_Radix)); put(")");
     new_line(1); 
     put ("In other words, each digits is in the range:        0 .. 2**(");
     put (e_Integer'Image(Desired_No_Of_Bits_In_Radix)); put(")"); put (" - 1.");
     new_line(2);
     put("If a number is correct to all digits except the last digit, expect:");
     new_line(2);
     put("Fractional part of error: 1.0000000000E+00");
     put("  Error's exponent: "); put(e_Integer'Image(-e_Real_Machine_Mantissa+1));
     new_line(2); 
     put("An error of this size is of the order:  10**(");
     put(Integer'Image (Integer  (-Decimal_Digits_in_Mantissa+9.0))); put(")");
     new_line(2); 
     put("If you lose two digits of accuracy (i.e. both guard digits) instead of 1");
     new_line(1); 
     put("(as in the above case) then you lose another 9 decimal digits of");
     new_line(1); 
     put("accuracy.  In this case expect errors of the order:");
     new_line(2);
     put("Fractional part of error: 1.0000000000E+00");
     put("  Error's exponent: "); put(e_Integer'Image(-e_Real_Machine_Mantissa+2));
     new_line(2); 
     put("An error of this size is of the order:  10**(");
     put(Integer'Image (Integer (-Decimal_Digits_in_Mantissa+18.0))); put(")"); 
     new_line(2); 
  end Print_Extended_Real_Settings;
  
begin

  rnd.Reset (Some_Seed);

  Print_Extended_Real_Settings;

  put ("Test of +,-,/,* over randomly generated test vectors.");
  new_line(1); 
  put ("The tests search for failures in +,-,/,*.  They don't attempt");
  new_line(1); 
  put ("to measure ultimate error in floating point operations.");
  new_line(1); 
  put ("The testing goes on for as many hours as you let it:");
  new_line(2);

  for repetitions in 1 .. 1_000_000 loop


  Big_Digits := Get_999;

  Exp_First :=  e_Real_Machine_Emin / 2 + e_Real_Machine_Mantissa;
  Exp_Last  :=  e_Real_Machine_Emax / 2 - 0;
  --  Above settings are for division/mult tests:
  
  --  The random num we make by scaling with Exp in [Exp_First, Exp_Last]
  --  should have the property that when squared it is still in range,
  --  since the function Random returns a rand in the range [0.0  .. 1.0).
  --  The smallest exp of Random is e_Real_Machine_Mantissa - 1.
  
   --*******************************************  
   --  NOTE: X = 1 + Rand or  Y = 1 + Rand; guarantees loss of 3 digits
   --  accuracy.  This method doesn't really measure error in * or /.  
   --  It does however show the importance of guard digits.
   --*******************************************
   Exp_Error  := e_Integer'First;
   Frac_Error := 1.0E-100;
  
   for I in Integer32 range 1 .. Mult_Limit loop
  
      X       := Random;
      X       := Scaling (X, Random_Exp (Exp_First, Exp_Last));
    
      Y       := One - Half_Digit*Random;
      Y       := Scaling (Y, Random_Exp (Exp_First, Exp_Last));
  
      Z1 := X;
      Mult (Z1, Y);  -- Z1 := X * Y;
      Z1 := Z1 / Y;  -- (so Z1 should equal X)
      Delta_X := Z1 - X;
    
      --Z1 := Y;
      --Square (Z1);
      --Z1 := Z1 / Y;
      --Delta_X := (Z1 - Y);
    
      Get_Normalized_Delta (X, Delta_X, Frac_Part, Exp_Part);
    
      IF Exp_Part  > Exp_Error  THEN Exp_Error := Exp_Part;   END IF;
      IF Frac_Part > Frac_Error THEN Frac_Error := Frac_Part; END IF;
         
   end loop;
   
   new_line; put ("Max Error =");
   Print_e_Real_data (Frac_Error, Exp_Error);

   --*******************************************
   Exp_Error  := e_Integer'First;
   Frac_Error := 1.0E-100;
  
   for I in Integer32 range 1 .. Mult_Limit loop
  
      X       := Random;
      X       := Scaling (X, Random_Exp (Exp_First, Exp_Last));
    
      Y       := One - Half_Digit*Random;
      Y       := Scaling (Y, Random_Exp (Exp_First, Exp_Last));
    
      Z1      := (X * Y) / Y;
      Delta_X := Z1 - X;
    
      Get_Normalized_Delta (X, Delta_X, Frac_Part, Exp_Part);
    
      IF Exp_Part  > Exp_Error  THEN Exp_Error := Exp_Part;   END IF;
      IF Frac_Part > Frac_Error THEN Frac_Error := Frac_Part; END IF;
         
   end loop;
   
   new_line; put ("Max Error =");
   Print_e_Real_data (Frac_Error, Exp_Error);
   
   --*******************************************  
   Exp_Error  := e_Integer'First;
   Frac_Error := 1.0E-100;
  
   for I in Integer32 range 1 ..  Mult_Limit loop
  
      X       := Random;
      X       := Scaling (X, Random_Exp (Exp_First, Exp_Last));
    
      Y       := Scaling (Big_Digits, Random_Exp (Exp_First, Exp_Last));
    
      Z1      := (X * Y) / Y;
      Delta_X := Z1 - X;
    
      Get_Normalized_Delta (X, Delta_X, Frac_Part, Exp_Part);
    
      IF Exp_Part  > Exp_Error  THEN Exp_Error := Exp_Part;   END IF;
      IF Frac_Part > Frac_Error THEN Frac_Error := Frac_Part; END IF;
         
   end loop;

   new_line; put ("Max Error =");
   Print_e_Real_data (Frac_Error, Exp_Error);
   
   --*******************************************  
    
   Exp_First := -e_Real_Machine_Emax + 2;
   Exp_Last  :=  e_Real_Machine_Emax - 2;
   --  Above settings are for +/- tests:
  
   --*******************************************  
   Exp_Error  := e_Integer'First;
   Frac_Error := 1.0E-100;
  
   for I in Integer32 range 1 .. Mult_Limit * 2**2 loop
  
      X       := Random;
      X       := Scaling (X, Random_Exp (Exp_First, Exp_Last));
    
      Y_Exp_First := Max (Exp_First, Exponent(X)+1);  
      Y_Exp_Last  := Min (Exp_Last,  Exponent(X) + e_Real_Machine_Mantissa + 1);  
      Y           := Scaling (Big_Digits, Random_Exp (Y_Exp_First, Y_Exp_Last));
    
      Z1      := (X + Y) - X;
      Delta_X := Z1 - Y;
    
      Get_Normalized_Delta (Y, Delta_X, Frac_Part, Exp_Part);
    
      IF Exp_Part  > Exp_Error  THEN Exp_Error := Exp_Part;   END IF;
      IF Frac_Part > Frac_Error THEN Frac_Error := Frac_Part; END IF;
         
   end loop;
 
   new_line; put ("Max Error =");
   Print_e_Real_data (Frac_Error, Exp_Error);
   
   
   --*******************************************  
   Exp_Error  := e_Integer'First;
   Frac_Error := 1.0E-100;
  
   for I in Integer32 range 1 .. Mult_Limit * 2**2 loop
  
      X       := Random;
      X       := Scaling (X, Random_Exp (Exp_First, Exp_Last));
    
      Y_Exp_First := Max (Exp_First, Exponent(X) - e_Real_Machine_Mantissa - 1);  
      Y_Exp_Last  := Min (Exp_Last,  Exponent(X) - 1);  
      Y           := Scaling (Big_Digits, Random_Exp (Y_Exp_First, Y_Exp_Last));
    
      Z1      := (X + Y) - Y;
      Delta_X := Z1 - X;
    
      Get_Normalized_Delta (X, Delta_X, Frac_Part, Exp_Part);
    
      IF Exp_Part  > Exp_Error  THEN Exp_Error := Exp_Part;   END IF;
      IF Frac_Part > Frac_Error THEN Frac_Error := Frac_Part; END IF;
         
   end loop;
   
   new_line; put ("Max Error =");
   Print_e_Real_data (Frac_Error, Exp_Error);
   
   
   --*******************************************  
   Exp_Error  := e_Integer'First;
   Frac_Error := 1.0E-100;
  
   for I in Integer32 range 1 .. Mult_Limit * 2**2 loop
  
      X       := Random;
      X       := Scaling (X, Random_Exp (Exp_First, Exp_Last));

      Y_Exp_First := Max (Exp_First, Exponent(X) - e_Real_Machine_Mantissa - 1);
      Y_Exp_Last  := Min (Exp_Last,  Exponent(X) - 1);  
      Y           := Random;
      Y           := Scaling (Y, Random_Exp (Y_Exp_First, Y_Exp_Last));
    
      Z1      := (X + Y) - Y;
      Delta_X := Z1 - X;
    
      Get_Normalized_Delta (X, Delta_X, Frac_Part, Exp_Part);
    
      IF Exp_Part  > Exp_Error  THEN Exp_Error := Exp_Part;   END IF;
      IF Frac_Part > Frac_Error THEN Frac_Error := Frac_Part; END IF;
         
   end loop;
   
   new_line; put ("Max Error =");
   Print_e_Real_data (Frac_Error, Exp_Error);
   
   end loop; -- repetitions
   
 
end;
   
