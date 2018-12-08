
with Text_IO; use Text_IO;
with Ada.Numerics.Generic_Elementary_Functions;
with Extended_Real;
with Extended_Real.Elementary_Functions;
with Extended_Real.IO;

procedure  e_real_demo_1 is

   type Real is digits 15; 
   
   package mth  is new Ada.Numerics.Generic_Elementary_Functions(Real);
   use mth;
   package ext is new Extended_Real (Real);
   use ext;
   package fnc is new Ext.Elementary_Functions (Sqrt, Log, Exp, Arcsin);
   use fnc;
   package rio is new Text_IO.Float_IO (Real);
   use rio;
   package iio is new Text_IO.Integer_IO (E_Integer);
   use iio;
   package eio is new Ext.IO;
   use eio;
   
   Ten : constant E_Real:= +10.0;
   
   Last : Natural;
   Y : E_Real;
   
   Blank_Str      : constant String(1..1024) := (others => ' ');
   Number_String  : String(1..1024) := (others => ' ');
   Seventy_Digits : constant String(1..74) := 
     "1.234567890123456789012345678901234567890123456789012345678901234567890E-1";

  Radix                : constant Real := E_Real_Machine_Radix;  
  
  Z1, Z2 : E_Real;
  
  --  Attempt to make Junk constants of order one with non-zero's all over.
  Junk1 : constant E_Real := One / (+3.0);
  Junk2 : constant E_Real := (+17.0) / (+19.0);
  
  Test_Vector_Seed : Real;
  
  No_Decimal_Digits : constant := Desired_Decimal_Digit_Precision;

   -----------
   -- Pause --
   ----------- 

   procedure Pause (s0,s1,s2,s3,s4,s5,s6,s7,s8,s9 : string := "") is
     Continue : Character := ' ';
   begin
     new_line;
     if S0 /= "" then put_line (S0); end if;
     if S1 /= "" then put_line (S1); end if;
     if S2 /= "" then put_line (S2); end if;
     if S3 /= "" then put_line (S3); end if;
     if S4 /= "" then put_line (S4); end if;
     if S5 /= "" then put_line (S5); end if;
     if S6 /= "" then put_line (S6); end if;
     if S7 /= "" then put_line (S7); end if;
     if S8 /= "" then put_line (S8); end if;
     if S9 /= "" then put_line (S9); end if;
     new_line;
     begin
	put ("Type a character to continue: ");
	get_immediate (Continue);
     exception
	when others => null;
     end;
   end Pause;

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
      put("In most large floating pt. calculations you lose both guard digits");
      new_line(1); 
      put("at the minimum.  The best precision you can expect is given by");
      new_line(1); 
      put("the 2nd number above.");
      new_line(1); 

      if e_Real_Machine_Mantissa > 7 then
         new_line(2); 
	 put ("PRESENTLY USING:");
	 put (e_Integer'Image(e_Real_Machine_Mantissa)); 
	 put (" DIGITS.");
         new_line(1); 
	 put ("Parts of this demo don't work great with more than 7 digits.");
      end if;

      Pause;

   end Print_Extended_Real_Settings;


  procedure Test_Rem (Test_Vector_Seed : Real) is 
    Difference, Max_Err : E_Real;
    Z1, Z2 : E_Real := Zero;
    Z4 : E_Real;
    No_of_Trials : constant Integer := 20_000;
  begin
    
    for I in 1 .. No_of_Trials loop
      Z1 := Z1 + Junk1;
      Z2 := Z2 + (+Test_Vector_Seed) * Junk2;
      
      -- try to reproduce Z1:
      
      Z4 := Remainder (Z1, Z2) + Z2 * Unbiased_Rounding (Z1 / Z2);  
      -- Z4 should equal Z1.
      
      Difference := Abs (One - Z1 / Z4);

      --put (Make_Real (Z3 / Z2)); put (Make_Real (Z4 / Z2)); New_Line;
      --  These quantities should always be in range [-0.5..0.5]

      if Difference > Max_Err then
          Max_Err := Difference; 
      end if;
       
   end loop;
   
   put ("Result = "); 
   put (e_Real_Image (Max_Err));
   new_line;
   
  end Test_Rem;

  -----------------------------
  -- Test_Digit_Mult_and_Div --
  -----------------------------
  
  procedure Test_Digit_Mult_and_Div 
   (Test_Vector_Seed : Real) 
  is 
    Difference, Max_Err : E_Real; -- init to 0 essential
    Delta_Digit : constant Real := 50_000.0;
    Real_Digit : Real := Radix + Delta_Digit - 1.0;
    Digit1 : e_Digit;
    Z1, Z2 : E_Real := Zero;
    No_of_Trials : constant Integer := 20_000;
  begin
    Z1   := (+Test_Vector_Seed) * Junk1;
    
    for I in 1 .. No_of_Trials loop

      Real_Digit := Real_Digit - Delta_Digit;
      Digit1     := Make_e_Digit (Real_Digit);

      Z2         := Digit1 * (Z1 / Digit1);
      Difference := Abs (One - Z2 / Z1);

      if Difference > Max_Err then
         Max_Err := Difference;
      end if;
       
    end loop;
   
    put ("Result = "); 
  --put (e_Real_Image (Max_Err);
    put (e_Real_Image (Max_Err, Aft => Integer'Min (No_Decimal_Digits, 50)));
    new_line;
   
  end Test_Digit_Mult_and_Div;
  
  -----------------------
  -- Test_Mult_and_Add --
  -----------------------

   -- uses Newton's method to calculate 1/(X).  Square and invert to
   -- compare with X; square and multiply with with X to compare w. 1.0.

  procedure Test_Mult_and_Add 
   (Test_Vector_Seed : Real) 
  is 
    
    Difference, Max_Err : E_Real;
    No_of_Trials : constant Integer := 20_000;

    function Inverse (X : E_Real) return E_Real is
      X_isqr     : E_Real;
      Iterations : constant Integer := 24;
      X_start    : constant Real    := 1.0 / Make_Real(X);
    begin
      X_isqr := Make_Extended (X_Start);
      for I in 1..Iterations loop
        X_isqr := X_isqr + (One - X * X_isqr) * X_isqr;
      end loop;
      return X_isqr;
    end Inverse;
    
  begin
    Z1  := Zero;
    
    for I in 1..No_of_Trials loop
      Z1 := Z1 + (+Test_Vector_Seed) * Junk1;
      
      Z2         := Inverse (Z1);
      Difference := Abs (One - Z1 * Z2);

      if Difference > Max_Err then
          Max_Err := Difference;
      end if;
       
   end loop;

   put ("Result = "); 
   put (e_Real_Image (Max_Err));
   new_line;
   
  end Test_Mult_and_Add;

  -----------------------
  -- Test_Mult_and_Div --
  -----------------------
  
  procedure Test_Mult_And_Div 
   (Test_Vector_Seed : Real) 
  is 
    Difference, Max_Err : E_Real;
    No_of_Trials : constant Integer := 20_000;
  begin
    for I in 1 .. No_of_Trials loop
      Z1  := Z1 + Junk2 * (+Test_Vector_Seed);
      Z2  := (Z1 * (One / Z1));
    --Z2  := (Z1 / Z1);
      Difference := Abs (One - Z2);

      if Difference > Max_Err then
          Max_Err := Difference;
      end if;
       
    end loop;

    put ("Result = "); 
    put (e_Real_Image (Max_Err));
    new_line;
  end Test_Mult_And_Div;
  

  ------------------------
  -- Test_Make_Extended --
  ------------------------
  
  procedure Test_Make_Extended (Test_Vector_Seed : Real) is 
    Max : Real;
    Difference   : Real;
    Junk1, Junk2 : Real;
    No_of_Trials : constant Integer := 20_000;
  begin
    Junk1 := 0.0;
    Max   := 0.0;
    for I in 1 .. No_of_Trials loop
      Junk1  := Junk1 + Test_Vector_Seed;
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
   

begin

   Print_Extended_Real_Settings;

   new_line(2); 
   put ("Next step: consider the following 70 digit number:");
   new_line(2); 
   put (" " & Seventy_Digits);

   new_line(1); 
   pause ( 
     "The 70 digit number printed above will be translated to binary,",
     "(e_Real) then translated back to text and printed just below.");

   E_Real_Val (Seventy_Digits, Y, Last);   
   new_line; 
   put (E_Real_Image (Y));   

   new_line(1); 
   pause ( 
     "The 1st number at the top is correctly reproduced (usually) a few digits",
     "beyond the advertized precision.  Beyond that a few more wrong digits",
     "are retained.  The point at which the correct digits stop and the wrong",
     "start varies, so it is usually best to keep them all for best accuracy." );
   new_line(1); 
  

 --pause ("You can round away a guard digit by calling function Machine(X):");
 --new_line; 
 --put (E_Real_Image (Machine (Y)));
 --new_line(1); 

   new_line(1); 
   new_line(1); 
   pause ( 
     "Usually you don't want to print all these digits, so use for example,",
     "E_Real_Image (X, Aft => 15) to print 15 digits after the decimal.",
     "This is not rounded to 15 digits.  It just truncates the string.",
     "Also, the min value for Aft is 12.  It will always print at least 12."
     );
   new_line(1); 
   E_Real_Val (Seventy_Digits, Y, Last);   
   new_line; 
   put (E_Real_Image (Y, Aft => 15));
   new_line(1); 

  
   --procedure E_Real_Val  (X    : in String; 
   --                     Result: out E_Real;
   --                      Last : out Natural);

   Pause ("More of the same.");

   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..50) := "12345678901234567890123456789012345678901234567890";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..52) := "  12345678901234567890123456789012345678901234567890";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..53) := "12345678901234567890123456789012345678901234567890   ";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..54) := "12345678901234567890123456789012345678901234567890.   ";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..54) := "123456789012345678901234.56789012345678901234567890   ";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..56) := "   .12345678901234567890123456789012345678901234567890  ";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..59) := "00000.12345678901234567890123456789012345678901234567890   ";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..55) := "+.12345678901234567890123456789012345678901234567890   ";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..55) := "-.12345678901234567890123456789012345678901234567890   ";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..58) := "-.12345678901234567890123456789012345678901234567890e-0111";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..57) := "-.12345678901234567890123456789012345678901234567890E0111";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line(1); 
   Number_String        := (others => ' ');
   Number_String(1..59) := ".1234567890123456789012345678901234567890123456789e99999999";
   new_line; put ("  "); put (Number_String(1..60));
   new_line; put ("   Translating this string to binary and back to text, I get:");
   E_Real_Val (Number_String, Y, Last);   
   new_line; put ("  " & E_Real_Image (Y));   
 
   new_line;
   pause ("Print a few numbers:");

   new_line; put ("One:");
   new_line; put (E_Real_Image (One));
   new_line; put ("10**400:");
   new_line; put (E_Real_Image (Ten**400));
   new_line; put ("10**-400:");
   new_line; put (E_Real_Image (Ten**(-400)));
   new_line; put ("10**-2899:");
   new_line; put (E_Real_Image (Zero + Ten**(-2899)));
   new_line; put ("E_Real_Model_Epsilon:");
   new_line; put (E_Real_Image (E_Real_Model_Epsilon));
   new_line; put ("E_Real_Model_Epsilon + 10**-600:");
   new_line; put (E_Real_Image (E_Real_Model_Epsilon + Ten**(-600)));
   new_line; put ("1 + E_Real_Model_Epsilon:");
   new_line; put (E_Real_Image (E_Real_Model_Epsilon + One));
   new_line; put ("1 + E_Real_Machine_Epsilon:");
   new_line; put (E_Real_Image (E_Real_Machine_Epsilon + One));
   new_line; put ("0.0:");
   new_line; put (E_Real_Image (Zero));
   new_line;


  new_line(1); 
  pause ("  The following tests the 2 functions:   Make_Extended and Make_Real.",
         "  These 2 routines translate between ordinary 15 decimal digit floats",
         "  (Real) and extended precision floats (e_Real). Numbers in the ordinary",
         "  floating point type (Real) are transformed into extended precision floats",
         "  (e_Real) by calls to Make_Extended. e_Real's are transformed back to",
         "  Real's (with a loss of precision) by calls to Make_Real. Below we print:",
         "         ",
         "       Result :=  X - Make_Real (Make_Extended(X))",
         "         ",
         "  2_000 values of X are used each test, and max Abs(Result) is printed.");
  new_line;

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

  Test_Vector_Seed := 1.234567891234567891E+31;
  Test_Make_Extended (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.234567891234567891E+8;
  Test_Make_Extended (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.234567891234567891E-8;
  Test_Make_Extended (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.234567891234567891E-31;
  Test_Make_Extended (Test_Vector_Seed);


  new_line(1); 
  pause ("  Some tests of ""+"" and ""*"".  The following is calculated:  ",
         "         ",
         "         Result :=  One - X * Reciprocal (X)",
         "         ",
         "  where Reciprocal (X) = 1/X is obtained from Newton's method.",
         "  2_000 values of X are used each test, and max Abs(Result) is printed.");
  new_line;
    
  Test_Vector_Seed := 1.2345678912345678E+31;
  Test_Mult_and_Add (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+8;
  Test_Mult_and_Add (Test_Vector_Seed);
    
  Test_Vector_Seed := 1.2345678912345678E-8;
  Test_Mult_and_Add (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-31;
  Test_Mult_and_Add (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-14;
  for I in 1..10 loop
     Test_Vector_Seed := Test_Vector_Seed + 1.2345678912345678E-4;
     Test_Mult_and_Add (Test_Vector_Seed);
  end loop;
  Test_Vector_Seed := 1.2345678912345678E+14;
  for I in 1..10 loop
     Test_Vector_Seed := Test_Vector_Seed * 1.2345678912345678E+4;
     Test_Mult_and_Add (Test_Vector_Seed);
  end loop;

   pause ( 
     "Notice that the error is usually much smaller than you might expect",
     "from the number of decimal digits requested.  For simple operations",
     "like ""*"", ""+"", and ""/"" the first of the 2 guard digits is usually,",
     "calculated correctly, (and that's an extra 9 decimal digits of precision.",
     "In fact 2 guard digits (18 decimal digits) are always used.  The 2nd",
     "guard digits is for safety, and so that Sqrt's etc are of similar accuracy",
     "to operations like ""*"".");
   new_line(1); 


  new_line(1); 
  pause ("  Some tests of ""*"" and ""/"".  The following is calculated:",
         "         ",
         "     Result :=  One - X * (One / X)",
         "         ",
         "  2_000 values of X are used each test, and max Abs(Result) is printed.");
  new_line;

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

  Test_Vector_Seed := 1.2345678912345678E+31;
  Test_Mult_And_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+8;
  Test_Mult_And_Div (Test_Vector_Seed);
    
  Test_Vector_Seed := 1.2345678912345678E-8;
  Test_Mult_And_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-31;
  Test_Mult_And_Div (Test_Vector_Seed); 
  

  new_line;
  pause ("  Some tests of operations between e_Digit types and e_Real types.",
         "  e_Digits are special e_Reals that are small and efficient.",
         "  X is an e_Real and D is an e_Digit.  The following is calculated:",
         "         ",
         "       Result := One - (D * (X / D)) / X)",
         "         ",
         "  2_000 values of D are used each test, and max Abs(Result) is printed.");
  new_line(2);
  
  Test_Vector_Seed := 1.2345678912345678E+61;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+28;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+27;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+26;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E+4;
  for I in 1..10 loop
     Test_Vector_Seed := Test_Vector_Seed * 1.2345678912345678E+4;
     Test_Digit_Mult_and_Div (Test_Vector_Seed);
  end loop;

  Test_Vector_Seed := 1.2345678912345678E-4;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-7;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-18;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-28;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-31;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  Test_Vector_Seed := 1.2345678912345678E-91;
  Test_Digit_Mult_and_Div (Test_Vector_Seed);
  
  
  new_line;
  pause ("  Some simple tests of the exponentiation operator ""**"":");
  new_line;

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
  pause ("  A test of Round_Away_Smallest_Guard_Digit and e_Real_Machine_Epsilon:");  
  new_line;

  Z1 := (One + e_Real_Machine_Epsilon) - One;
  Z1 := Z1 / e_Real_Machine_Epsilon;
  put(" should be 1.0:"); put (Make_Real(Z1)); new_line;

  Z1 := Machine (One + e_Real_Machine_Epsilon) - One;
  put(" should be 0.0:"); put (Make_Real(Z1)); new_line;

  Z2 := Make_Extended(0.99999999999999);
  Z1 := Machine (Z2 + e_Real_Machine_Epsilon) - Z2;
  Z1 := Z1 / e_Real_Machine_Epsilon;
  put(" should be 1.0:"); put (Make_Real(Z1)); new_line;

  loop
     begin
       new_line(2); put ("Enter a number: ");
       Number_String := Blank_Str;
       --get_line (Number_String, LengthStr);
       get_Line (Number_String, Last);
       E_Real_Val (Number_String, Y, Last);   
       new_line; put ("I read this as: ");
       new_line; put (E_Real_Image (Y));
       exit;
     exception
       when others =>
        put_line("Some error.  Try again.");
     end;
   end loop;

  
end;
