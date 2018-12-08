
-- Simple test of "*", "/", "+", and "-" using random arguments.
-- Doesn't test endpoints very well.

with Extended_Real;
with Extended_Real.e_Rand; 
with Extended_Real.IO; 
with Text_IO; use Text_IO;

procedure  e_real_io_tst_2 is

  type Real is digits 15; 
  
  package ext is new Extended_Real (Real);
  use ext;
  package eio is new ext.IO;
  use eio;
  package rnd is new ext.e_Rand;
  use rnd;
  package rio is new Text_IO.Float_IO(Real);
  use rio;
  package iio is new Text_IO.Integer_IO(e_Integer);
  use iio;
  
  X, Z1   : e_Real;
  Last : Integer;

  Max_Error, Err  : e_Real;

  Exp_First, Exp_Last : e_Integer; 

  type Integer32 is range -2**31+1..2**31-1;
  
  Mult_Limit : constant Integer32 := 8_000_000; -- usually > N * 10^6
  --  Number of iterations of multiplication div test.  The +/- tests
  --  are 8 times this numbers.
  
  Some_Seed : constant Integer := 7251738;
  --  Start at different rand stream.
  
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
     --X := Real (rnd.Next_Random_Int mod 2**12) * 2.0**(-12);

     --  returns random Real in range (0, 1]:
     X := Real (1 + rnd.Next_Random_Int mod 2**24) * 2.0**(-24);

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

     Next_Digit : e_Digit;
     Delta_Exp  : e_Integer;
     Result     : e_Real;     -- Init to 0 important
     
  begin
  
     for I in 0 .. e_Real_Machine_Mantissa-1 loop
     
        Delta_Exp  := -I - 1;
        Next_Digit := Scaling (Max_Digit, Delta_Exp);
        Result     := Next_Digit + Result;
        
     end loop;
      
     return Result;
     
  end;
  
  ----------------------------------
  -- Print_Extended_Real_Settings --
  ----------------------------------
 
   procedure Print_Extended_Real_Settings
   is
      Bits_In_Radix : constant := Desired_No_Of_Bits_In_Radix;
   begin
      new_line(1); 
      put ("    Desired_Decimal_Digit_Precision =");
      put (Integer'Image(Desired_Decimal_Digit_Precision)); 
      new_line(1); 
      new_line(1); 
      put ("Number of decimal digits of precision requested:    ");
      put (Integer'Image(Desired_Decimal_Digit_Precision));
      new_line(1); 
      put ("Number of digits in use (including 2 guard digits): "); 
      put (e_Integer'Image(e_Real_Machine_Mantissa));
      new_line(1); 
      put ("These digits are not decimal; they have Radix:       2**(");
      put (e_Integer'Image(Bits_In_Radix)); put(")");
      new_line(1); 
      put ("In other words, each of these digits is in range:    0 .. 2**(");
      put (e_Integer'Image(Bits_In_Radix)); put(")"); put (" - 1.");
      new_line(1); 
      put ("Number of decimal digits per actual digit is approx: 9"); 
      new_line(2);
      put("Guard digits (digits of extra precision) are appended to the end of");
      new_line(1); 
      put("each number.  There are always 2 guard digits. This adds up to 18"); 
      new_line(1); 
      put("decimal digits of extra precision. The arithmetic operators, (""*"",");
      new_line(1); 
      put("""/"", ""+"" etc) usually produce results that are correct to all");
      new_line(1); 
      put("digits except the final (guard) digit.");
      new_line(2);
      put("If a number is correct to all digits except the final (guard) digit,");
      new_line(1); 
      put("expect errors of the order:");
      new_line(2);
      put(e_Real_Image (e_Real_Model_Epsilon / (One+One)**Bits_In_Radix, aft => 10));
      new_line(2); 
      put("If you lose 2 digits of accuracy (i.e. both guard digits) instead");
      new_line(1); 
      put("of 1 (as in the above case) then you lose another 9 decimal digits");
      new_line(1); 
      put("of accuracy.  In this case expect errors of the order:");
      new_line(2);
      put(e_Real_Image (e_Real_Model_Epsilon, aft => 10));
      new_line(2); 
      put("The above number, by the way, is e_Real_Model_Epsilon.");
      new_line(3);
 
  end Print_Extended_Real_Settings;
  
begin

  rnd.Reset (Some_Seed);

  Print_Extended_Real_Settings;

  put ("The test translates binary to text, then back to binary, and prints");
  new_line(1); 
  put ("the difference: prints X - e_Real_Val (e_Real_Image (X)) over randomly");
  new_line(1); 
  put ("generated X's.  8_000_000 X's are used, and the max error is printed.");
  new_line(1); 
  put ("The testing goes on for as many hours as you let it:");
  new_line(2);

  for repetitions in 1 .. 1_000_000 loop


  Exp_First :=  e_Real_Machine_Emin+1; -- cause X has exp of -1 intitially.
  Exp_Last  :=  e_Real_Machine_Emax-1;
  
  --  The random num we make by scaling with Exp in [Exp_First, Exp_Last]
  --  The function Random returns a rand in the range [0.0  .. 1.0).
  --  The smallest exp of Random is Max_Available_Digits - 1.
  
   Max_Error := Zero;
  
   for I in Integer32 range 1 .. Mult_Limit loop
  
      X       := Random * Get_999;
      X       := Scaling (X, Random_Exp (Exp_First, Exp_Last));
    
      e_Real_Val (e_Real_Image (X), Z1, Last);
      Err := (Z1/X - One);
    
      if Err > Max_Error then
        Max_Error := Err;
      end if;

   end loop;
   
   new_line; 
   put ("Max Error =");
   put (e_Real_Image (Max_Error));

   end loop;
   
end;
