with Extended_Real;
with Extended_Real.Elementary_Functions;
with Extended_Real.IO;
with Ada.Numerics.Generic_Elementary_Functions;
with Text_IO; use Text_IO;

procedure  e_function_demo_1 is

   type Real is digits 15;

   package mth is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use mth;
   package ext is new Extended_Real (Real);
   use ext;
   package fnc is new Ext.Elementary_Functions (Sqrt, Log, Exp, Arcsin);
   use fnc;
   package eio is new Ext.IO;
   use eio;
   package rio is new Text_IO.Float_IO (Real);
   use rio;


   Delta_2 : e_Real;
   Z0, Z1, Z2, Z3, Z4, Difference1  : e_Real;
   Junk1 : constant e_Real := E_Quarter_Pi * Make_Extended (0.273);
   Junk2 : constant e_Real := E_Quarter_Pi * Make_Extended (0.233);
   Junk3 : constant e_Real := E_Inverse_Sqrt_2; -- 1/SQRT(2)


   Test_Vector_Seed : Real;

   N : Positive;

   Two_Digit  : constant e_Digit := Make_e_Digit (2.0);
   Half       : constant e_Real  := Make_Extended (0.5);

   Limit : constant Integer := 1200;
   --  Number of test runs for some procedures. Can't exceed 1200 because
   --  some tests use (10.0**0.25)**Limit.

   e_Real_Decimals : constant := Desired_Decimal_Digit_Precision;
   

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
     new_line;
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

      Pause;

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
      new_line(1);

      new_line(1);
      put ("The above number, by the way, is:   e_Real_Model_Epsilon.");
      new_line(2);
      put ("Most computationally intensive floating pt. calculations will");
      new_line(1);
      put ("lose 2 guard digits of accuracy at the minimum.");

      Pause;

   end Print_Extended_Real_Settings;

   ---------------------------
   -- Test_Cos_and_Arccos_4 --
   ---------------------------

   -- Test rel. error in large Z1 limit.

   procedure Test_Cos_and_Arccos_4 is
     I1 : constant Integer := -Limit;
     I2 :constant  Integer := 0;
     Max_Error : e_Real; -- init essential

   begin
     Z1 := Junk1;

     for I in I1..I2 loop
       Z1     := (Z1 / (Z1 + (+(1.77827941**I))));
       Z2     := Cos (Arccos (Z1));

       Difference1 := Abs ((Z2 - Z1) / Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Cos_and_Arccos_4;


    -- Test rel. error in small Z1 limit.

   procedure Test_Cos_and_Arccos_3 is
     I1 : constant Integer := -Limit;
     I2 :constant  Integer := 0;
     Max_Error : e_Real; -- init essential

   begin
     Z1 := Junk1;

     for I in I1..I2 loop
       Z1     := Abs(Z1);
       Z1     := -(Z1 / (Z1 + (+(1.77827941**I))**2));
       Z2     := Cos (Arccos (Z1));

       Difference1 := Abs ((Z2 - Z1) / Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Cos_and_Arccos_3;


    -- Can't test relative error in small Z1 limit, only absolute,
    -- because Cos (Arccos (Z1)) simply doesn't reproduce Z1 for finite
    -- numbers of digits.  Recall Arcos just calls Arcsin.

   procedure Test_Cos_and_Arccos_2 is
     I1 : constant Integer := -Limit/2;
     I2 :constant  Integer := Limit/2;
     Max_Error : e_Real; -- init essential

   begin
     Z1 := Junk1;

     for I in I1..I2 loop
       Z1     := (Z1 / (Z1 + (+(1.77827941**I))));
       Z2     := Cos (Arccos (Z1));

       Difference1 := Abs (Z2 - Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Cos_and_Arccos_2;

    -- Can't test relative error in small Z1 limit, only absolute,
    -- because Cos (Arccos (Z1)) simply doesn't reproduce Z1 for finite
    -- numbers of digits.  Recall Arcos just calls Arcsin.

   procedure Test_Cos_and_Arccos_1 is
     I1 : constant Integer := -Limit/2;
     I2 :constant  Integer := Limit/2;
     Max_Error : e_Real; -- init essential
   begin
     Z1 := Junk1;

     for I in I1..I2 loop
       Z1     := Abs(Z1);
       Z1     := -(Z1 / (Z1 + (+(1.77827941**I))**2));
       Z2     := Cos (Arccos (Z1));

       Difference1 := Abs (Z2 - Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Cos_and_Arccos_1;


   -- Superficial test of Sin(X) and Arcsin(X).

   procedure Test_Sin_and_Arcsin_1 is
     I1 : constant Integer := -Limit;
     I2 : constant Integer := 0;
     Max_Error : e_Real; -- init essential
   begin
     Z1 := Junk1;

     for I in I1..I2 loop
       Z1     := Z1 / (Z1 + (+(1.77827941**I)));
       Z2     := Sin (Arcsin (Z1));

       Difference1 := Abs ((Z2 - Z1) / Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Sin_and_Arcsin_1;


   procedure Test_Sin_and_Arcsin_2 is
     I1 : constant Integer := -Limit;
     I2 : constant Integer := Limit+5;
     Max_Error : e_Real; -- init essential
   begin
     Z1 := Junk1;

     for I in I1..I2 loop
       Z1     := Z1 / (Z1 + (+(1.77827941**I)));
       Z2     := Sin (Arcsin (Z1));

       Difference1 := Abs ((Z2 - Z1) / Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Sin_and_Arcsin_2;




   procedure Test_Sin_and_Arcsin_3 is
     Max_Error : e_Real; -- init essential
     I1 : constant Integer := -Limit/2;
     I2 : constant Integer := Limit+5;
   begin
     Z1 := Junk1;

     for I in I1..I2 loop
       Z1     := Abs(Z1);
       Z1     := -Z1 / (Z1 + (+(1.77827941**I)));
       Z2     := Sin (Arcsin (Z1));

       Difference1 := Abs (Z2/Z1 - One);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Sin_and_Arcsin_3;


    -- test: Sin(X)**2 + Cos(X)**2 - 1.

   procedure Test_Sin_and_Cos_0
    (Test_Vector_Seed : Real) 
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 10_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);

     for I in 1 .. No_of_Trials loop
       Z1  :=  Junk2 + Z1;

       Z2  := Sin(Z1)**2 + Cos(Z1)**2;

       Difference1 := Abs (One - Z2);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Sin_and_Cos_0;

    -- test: Sin(X)**2 + Cos(X)**2 - 1.

   procedure Test_Sin_and_Cos_2
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 10_000;
     Shift : constant e_Digit := Make_e_Digit (2.0**(-27));
   begin
     Z1 := e_Inverse_Sqrt_2;

     for I in 1 .. No_of_Trials loop
       Z1  :=  Z1 + Shift * (Shift * e_Real_Model_Epsilon);

       Z2  := Sin(Z1)**2 + Cos(Z1)**2;

       Difference1 := Abs (One - Z2);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Sin_and_Cos_2;


   procedure Test_Tan_and_Arctan_2
    (Test_Vector_Seed : Real) 
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 4_000;
   begin
     Z1 := Abs (Junk1 * (+Test_Vector_Seed));

     for I in 1 .. No_of_Trials loop
        Z1  :=  Z1 + Junk2;
        Z2  := Arctan (Z1);
      --Z4  := Sin(Z2) / Cos(Z2);  -- Z1 - Sin(Z2) / Cos(Z2) not ok if cos=0
        Z3  := Sin(Z2);
        Z4  := Cos(Z2) * Z1;
        Difference1 := Abs ((Z3 - Z4) / Z3);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Tan_and_Arctan_2;


   procedure Test_Tan_and_Arctan
    (Test_Vector_Seed : Real) 
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 1000;
     Digit_Trials : constant e_Digit := Make_e_Digit (Real (No_of_Trials));
     dZ : constant e_Real := Two_Digit * e_Quarter_pi / Digit_Trials;
   begin

     for I in 1 .. No_of_Trials loop
       Z1 := Make_e_Digit (Real(I)) * dZ;

     --Z2  := Sin(Z1) / Cos(Z1);
       Z2  := Divide (Sin(Z1), Cos(Z1));
       Z4  := Arctan (Z2);
       Difference1 := Abs ((Z4 - Z1) / (Z1));

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
            Max_Error := Difference1;
         end if;
       end if;

     end loop;

     Z1 := Abs (Junk1 * (+Test_Vector_Seed));
     for I in 0 .. No_of_Trials loop
       Z1 := Z1 / (Z1 + One);
       Z1 := Two_Digit * Z1 * e_Quarter_pi;

     --Z2  := Sin(Z1) / Cos(Z1);
       Z2  := Divide (Sin(Z1), Cos(Z1));
       Z4  := Arctan (Z2);
       Difference1 := Abs ((Z4 - Z1) / (Z1));

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
            Max_Error := Difference1;
         end if;
       end if;

     end loop;

     Z1 := One / Abs (Junk1 * (+Test_Vector_Seed));
     for I in 0 .. No_of_Trials loop
       Z1 := Z1 / (Z1 + One);
       Z1 := Two_Digit * Z1 * e_Quarter_pi;

     --Z2  := Sin(Z1) / Cos(Z1);
       Z2  := Divide (Sin(Z1), Cos(Z1));
       Z4  := Arctan (Z2);
       Difference1 := Abs ((Z4 - Z1) / (Z1));

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     Z1 := Abs (Junk1 * (+Test_Vector_Seed));
     for I in 0 .. No_of_Trials loop
       Z1 := Z1 / (Z1 + (+(1.77827941**I)));
       Z1 := Two_Digit * Z1 * e_Quarter_pi;
     --Z1 := Z1 / (Z1 + One);

     --Z2  := Sin(Z1) / Cos(Z1);
       Z2  := Divide (Sin(Z1), Cos(Z1));
       Z4  := Arctan (Z2);
       Difference1 := Abs ((Z4 - Z1) / (Z1));

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;


     Z1 := Abs (Junk1 * (+Test_Vector_Seed));
     for I in -No_of_Trials .. 0 loop
       Z1 := Z1 / (Z1 + (+(1.77827941**I)));
       Z1 := Two_Digit * Z1 * e_Quarter_pi;
     --Z1 := Z1 / (Z1 + One);

     --Z2  := Sin(Z1) / Cos(Z1);
       Z2  := Divide (Sin(Z1), Cos(Z1));
       Z4  := Arctan (Z2);
       Difference1 := Abs ((Z4 - Z1) / (Z1));

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Tan_and_Arctan;


   procedure Test_Sin_and_Cos_1
    (Test_Vector_Seed : Real) 
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 10_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);

     for I in 1 .. No_of_Trials loop
       Z1  :=  Z1 + Junk2;

       Z2  := Two_Digit * Sin(Z1) * Cos(Z1);
       Z3  := Sin (Two_Digit * Z1);

     --Difference1 := Abs (Z3 / Z2 - One);
       Difference1 := Abs ((Z3 - Z2) / Z2);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Sin_and_Cos_1;

   --------------
   -- Test_Log --
   --------------

   procedure Test_Log 
    (Test_Vector_Seed : Real) 
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 10_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk2 * (+Test_Vector_Seed);

     for I in 1 .. No_of_Trials loop
       Z1     := Z1 + Z0;
       Z2     := Log (Z1);
       Z2     := Exp (Z2);  -- should bring it back to Z1

       Difference1 := Abs (One - Z2 / Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Log;


   ----------------
   -- Test_Log_2 --
   ----------------

    -- Verified that Log(2) is correct by independent calc of Log_2.
    -- Now check other arguments using Log(A*B) = ...

   procedure Test_Log_2 (Test_Vector_Seed : Real) is
     Max_Error : e_Real; -- init essential
     Log_Z0    : e_Real;
     No_of_Trials : constant Integer := 10_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk2 * (+Test_Vector_Seed);

     Log_Z0 := Log (Z0);

     for I in 1 .. No_of_Trials loop
       Z1     := Z1 + Z0;
       Z2     := Log (Z1*Z0);
       Z3     := Log (Z1) + Log_Z0;

       Difference1 := Abs (One - Z3 / Z2);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);


   end Test_Log_2;


   -----------------------------
   -- Test_Exp_w_Integer_Args --
   -----------------------------

   -- Verified that Exp is correct by independent calc of Log_2.
   -- Now check other arguments using Exp(A+B) = ...

   procedure Test_Exp_w_Integer_Args 
    (Test_Vector_Seed : Integer := 0) 
   is
     Max_Error : e_Real; -- init essential
     Exp_Z0    : e_Real;
     No_of_Trials : constant Integer := 500_000;
   begin
     Z0 := (+Real (2 + Test_Vector_Seed));
     Z1 := (+0.0);

     Exp_Z0 := Exp (Z0);

     for I in 1 .. No_of_Trials loop
       Z1     := Z1 + Z0;
       Z2     := Exp (Z1 + Z0);
       Z3     := Exp (Z1) * Exp_Z0;

       Difference1 := Abs (Z3/Z2 - One);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Exp_w_Integer_Args;

   --------------
   -- Test_Exp --
   --------------

   procedure Test_Exp 
    (Test_Vector_Seed : Real) 
   is
     Max_Error, Exp_Z0 : e_Real; -- init essential
     No_of_Trials : constant Integer := 10_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk3 * (+Test_Vector_Seed);
     Exp_Z0 := Exp (Z0);

     for I in 1..No_of_Trials loop
       Z1     := Z1 + Z0;
       Z2     := Exp (Z1 + Z0);
       Z3     := Exp (Z1) * Exp_Z0;

       Difference1 := Abs (Z3/Z2 - One);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Exp;

   -----------------
   -- Test_Sqrt_2 --
   -----------------

   -- uses ** to calculate SQRT(X).  Square and
   -- compare with X; square and divide by X to compare w. 1.0.

   procedure Test_Sqrt_2 
    (Test_Vector_Seed : Real) 
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 10_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk2 * (+Test_Vector_Seed);

     for I in 1 .. No_of_Trials  loop
       Z1   := Z1 + Z0;
       Z2   := Z1**(Half);

       Difference1 := Abs (One - (Z2 * Z2) / Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Sqrt_2;

   ----------------------------
   -- Test_Reciprocal_Sqrt_2 --
   ----------------------------

   -- uses ** to calculate SQRT(X).  Square and
   -- compare with X; square and divide by X to compare w. 1.0.

   procedure Test_Reciprocal_Sqrt_2 
    (Test_Vector_Seed : Real) 
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 10_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk2 * (+Test_Vector_Seed);

     for I in 1 .. No_of_Trials  loop
       Z1   := Z1 + Z0;
       Z2   := Z1**(-Half);

     --Difference1 := Abs (One/Z1 - (Z2 * Z2)) * Z1;
       Difference1 := Abs (One - (Z2 * Z2) * Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Reciprocal_Sqrt_2;

   --------------------------
   -- Test_Reciprocal_Sqrt --
   --------------------------

   -- uses Newton's method to calculate SQRT(X).  Square and
   -- compare with X; square and divide by X to compare w. 1.0.

   procedure Test_Reciprocal_Sqrt 
    (Test_Vector_Seed : Real) 
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 50_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk2 * (+Test_Vector_Seed);

     for I in 1 .. No_of_Trials  loop
       Z1   := Z1 + Z0;
       Z2   := Reciprocal_Sqrt (Z1);

       Difference1 := Abs (One - (Z2 * Z2) * Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Reciprocal_Sqrt;

   ---------------
   -- Test_Sqrt --
   ---------------

   -- uses Newton's method to calculate SQRT(X).  Square and
   -- compare with X; square and divide by X to compare w. 1.0.

   procedure Test_Sqrt 
    (Test_Vector_Seed : Real) 
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 50_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk2 * (+Test_Vector_Seed);

     for I in 1..No_of_Trials loop
       Z1   := Z1 + Z0;
       Z2   := Sqrt (Z1);
       Z2   := Z2 * Z2;

     --Difference1 := Abs ((Z1 - Z2) / Z1);
       Difference1 := Abs (One - Z2 / Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Sqrt;

   ---------------------
   -- Test_Reciprocal --
   ---------------------

   -- uses Newton's method to calculate Inverse of (X).
   -- compare with X; mult. by X to compare w. 1.0.

   procedure Test_Reciprocal 
    (Test_Vector_Seed : Real)
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 50_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk2 * (+Test_Vector_Seed);

     for I in 1 .. No_of_Trials loop
       Z1     := Z1 + Z0;
       Z2     := Reciprocal(Z1);

       Difference1 := Abs (One - Z1 * Z2);
     --Difference1 := Abs (One / Z1 - Z2) * Z1;

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Reciprocal;

   ----------------------
   -- Test_Divide_stnd --
   ----------------------

   -- uses Newton's method to calculate X/Y
   -- compare with X; mult. by Y to compare w. X.

   procedure Test_Divide_stnd 
    (Test_Vector_Seed : Real)
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 50_000;
   begin
     Z2 := Junk1 * (+Test_Vector_Seed);
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk2 * (+Test_Vector_Seed);

     for I in 1 .. No_of_Trials loop
       Z1     := Z1 + Z0;
       Z2     := Z2 + Z2 + Junk1; 
       Z3     := Z2 / Z1;
     --Z3     := Divide (Z2, Z1);

       Difference1 := Abs (Z2 - Z3 * Z1) / Z2;

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;
     end loop;

     for I in 1 .. No_of_Trials loop
       Z1     := Z1 + Z0;
       Z2     := Z2 + Z2 + Junk1; 
       Z3     := Z1 / Z2;
     --Z3     := Divide (Z1, Z2);

       Difference1 := Abs (Z1 - Z3 * Z2) / Z1;

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;
     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Divide_stnd;


   -----------------
   -- Test_Divide --
   -----------------

   -- uses Newton's method to calculate X/Y
   -- compare with X; mult. by Y to compare w. X.

   procedure Test_Divide 
    (Test_Vector_Seed : Real)
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 50_000;
   begin
     Z2 := Junk1 * (+Test_Vector_Seed);
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk2 * (+Test_Vector_Seed);

     for I in 1 .. No_of_Trials loop
       Z1     := Z1 + Z0;
       Z2     := Z2 + Z2 + Junk1; 
     --Z3     := Z2 / Z1;
       Z3     := Divide (Z2, Z1);

     --Difference1 := Abs (Z2 - Z3 * Z1) / Z2;
       Difference1 := Divide (Abs (Z2 - Z3 * Z1), Z2);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;
     end loop;

     for I in 1 .. No_of_Trials loop
       Z1     := Z1 + Z0;
       Z2     := Z2 + Z2 + Junk1; 
     --Z3     := Z1 / Z2;
       Z3     := Divide (Z1, Z2);

     --Difference1 := Abs (Z1 - Z3 * Z2) / Z1;
       Difference1 := Divide (Abs (Z1 - Z3 * Z2), Z1);

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;
     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Divide;

   ---------------
   -- Test_Root --
   ---------------

   -- uses Newton's method to calculate Inverse of Nth root of (X).
   -- compare with X; ** and mult by X to compare w. 1.0.

   procedure Test_Root
    (Test_Vector_Seed : Real; N : Positive) 
   is
     Max_Error : e_Real; -- init essential
     No_of_Trials : constant Integer := 50_000;
   begin
     Z1 := Junk1 * (+Test_Vector_Seed);
     Z0 := Junk2 * (+Test_Vector_Seed);

     for I in 1..No_of_Trials loop
       Z1     := Z0 + Z1;

       Z2     := Reciprocal_Nth_Root (Z1, N);
       Z3     := Z2 ** N; -- so should be just Reciprocal of Z1

       Difference1 := One - Z1 * Z3;

       if Are_Not_Equal (Difference1, Zero) then
         if Difference1 > Max_Error then
             Max_Error := Difference1;
         end if;
       end if;

     end loop;

     put ("Estimated Max Error: "); 
     put (e_Real_Image (Max_Error, Aft => Integer'Min (40, e_Real_Decimals))); 
     new_line(1);

   end Test_Root;

begin

   Print_Extended_Real_Settings;

   new_line;
   Pause (
     "Some tests of Sin and Cos. ",
     "10_000 trials are performed for each number printed below.",
     "Test 1. calculate: Sin**2 + Cos**2 - 1."
     );
   new_line(2);

   Test_Sin_and_Cos_2;

   Test_Vector_Seed := 1.2345657891234E+24; -- near max arg for Cos
   Test_Sin_and_Cos_0 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+8;
   Test_Sin_and_Cos_0 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Sin_and_Cos_0 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+0;
   Test_Sin_and_Cos_0 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-3;
   Test_Sin_and_Cos_0 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-7;
   Test_Sin_and_Cos_0 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-31;
   Test_Sin_and_Cos_0 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-93;
   Test_Sin_and_Cos_0 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-233;
   Test_Sin_and_Cos_0 (Test_Vector_Seed);
 

   new_line;
   put ("Sin (2X) - 2*Sin (X)*Cos (X):");
   new_line(2);


   Test_Vector_Seed := 1.2345657891234E+24;
   Test_Sin_and_Cos_1 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+10;
   Test_Sin_and_Cos_1 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Sin_and_Cos_1 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+0;
   Test_Sin_and_Cos_1 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-3;
   Test_Sin_and_Cos_1 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-6;
   Test_Sin_and_Cos_1 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-8;
   Test_Sin_and_Cos_1 (Test_Vector_Seed);


   --  Arcsin_Stuff


   Pause ("Some tests of Sin and Arcsin. ");
   new_line(2);

   Test_Sin_and_Arcsin_1;
   Test_Sin_and_Arcsin_2;
   Test_Sin_and_Arcsin_3;

   new_line;
   put_Line ("Some tests of Cos and Arccos. ");
   new_line;

   Test_Cos_and_Arccos_1;
   Test_Cos_and_Arccos_2;
   Test_Cos_and_Arccos_3;
   Test_Cos_and_Arccos_4;

   Pause 
    ("Test of Sin (Arctan(X)) / Cos (Arctan(X)) = X",
     "4_000 trials are performed for each number printed below:"
     );

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Tan_and_Arctan_2 (Test_Vector_Seed);

   Test_Vector_Seed := 2.2345657891234E00;
   Test_Tan_and_Arctan_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E00;
   Test_Tan_and_Arctan_2 (Test_Vector_Seed);

   Test_Vector_Seed := 5.2345657891234E-1;
   Test_Tan_and_Arctan_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-3;
   Test_Tan_and_Arctan_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-10;
   Test_Tan_and_Arctan_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-99;
   Test_Tan_and_Arctan_2 (Test_Vector_Seed);

   Pause 
    ("Test of Arctan (Sin(X) / Cos(X)) = X.",
     "5_000 trials are performed for each number printed below:"
     );

   Test_Vector_Seed := 1.2345657891234E+54;
   Test_Tan_and_Arctan (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Tan_and_Arctan (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+10;
   Test_Tan_and_Arctan (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Tan_and_Arctan (Test_Vector_Seed);

   Test_Vector_Seed := 2.2345657891234E00;
   Test_Tan_and_Arctan (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E00;
   Test_Tan_and_Arctan (Test_Vector_Seed);

   Test_Vector_Seed := 5.2345657891234E-1;
   Test_Tan_and_Arctan (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-3;
   Test_Tan_and_Arctan (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-10;
   Test_Tan_and_Arctan (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-99;
   Test_Tan_and_Arctan (Test_Vector_Seed);


   Pause 
    ("Test of Exp (Log (X)) = X for small and large arguments.",
     "10_000 trials are performed for each number printed below:"
     );

   Test_Vector_Seed := 1.2345657891234E+99;
   Test_Log (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Log (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+10;
   Test_Log (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Log (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E00;
   Test_Log (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-3;
   Test_Log (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-10;
   Test_Log (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-34;
   Test_Log (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-99;
   Test_Log (Test_Vector_Seed);



   Pause ("Test of Log (X*Y) = Log (X) + Log (Y).",
     "10_000 trials are performed for each number printed below:"
     );

   Test_Vector_Seed := 1.2345657891234E+99;
   Test_Log_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Log_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+18;
   Test_Log_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+10;
   Test_Log_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Log_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+0;
   Test_Log_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-3;
   Test_Log_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-10;
   Test_Log_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-34;
   Test_Log_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-99;
   Test_Log_2 (Test_Vector_Seed);


   -- Exp_Stuff

   Pause (
     "Test Exp (X+Y) = Exp (X) * Exp (Y).",
     "10_000 trials are performed for each number printed below:"
     );


   Test_Vector_Seed := +0.2345657891234E+1;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := +1.2345657891234E+4;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := +1.2345657891234E+3;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := -1.2345657891234E+3;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := -1.2345657891234E-4;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := +1.2345657891234E-4;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := -1.2345657891234E+1;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := +1.2345657891234E+1;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := -1.2345657891234E+0;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := -1.2345657891234E-1;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := +1.2345657891234E-1;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := -1.2345657891234E-10;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := +1.2345657891234E-16;
   Test_Exp(Test_Vector_Seed);

   Test_Vector_Seed := -1.2345657891234E-16;
   Test_Exp(Test_Vector_Seed);


   new_line;
   pause ("Start simple test of Log and Exp routines. ");
   new_line;

   put_line("compare exp(0.25) with actual:");
   put (Exp(0.25)); new_line;
   put (Make_Real(Exp(Make_Extended(0.25)))); new_line;

   put_line("compare exp(0.2) with actual:");
   put (Exp(0.2)); new_line;
   put (Make_Real(Exp(Make_Extended(0.2)))); new_line;

   put_line("compare exp(0.28) with actual:");
   put (Exp(0.28)); new_line;
   put (Make_Real(Exp(Make_Extended(0.28)))); new_line;

   put_line("compare exp(1.0) with actual:");
   put (Exp(1.0)); new_line;
   put (Make_Real(Exp(Make_Extended(1.0)))); new_line;

   put_line("compare exp(20.0) with actual:");
   put (Exp(20.0)); new_line;
   put (Make_Real(Exp(Make_Extended(20.0)))); new_line;

   put_line("compare exp(-20.0) with actual:");
   put (Exp(-20.0)); new_line;
   put (Make_Real(Exp(Make_Extended(-20.0)))); new_line;

   put_line("compare log(2.0) with actual:");
   put (Log(2.0)); new_line;
   put (Make_Real(Log(Make_Extended(2.0)))); new_line;

   put_line("compare log(1.01) with actual:");
   put (Log(1.01)); new_line;
   put (Make_Real(Log(Make_Extended(1.01)))); new_line;

   put_line("compare log(1.0E34) with actual:");
   put (Log(1.0E34)); new_line;
   put (Make_Real(Log(Make_Extended(1.0E34)))); new_line;
   new_line(1);

   new_line(1);
   Pause ("Check Sin and Cos at a few specific arguments.");
   new_line(1);

   put ("Check Sin (0.25*Pi) = Sqrt(0.5). Estimated error = ");
   Delta_2 := Sin (E_Quarter_Pi) - Sqrt(Make_Extended(0.5));
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Sin (0.50*Pi) = 1.0.       Estimated error = ");
   Delta_2 := One - Sin (Two_Digit * E_Quarter_Pi);
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Sin (0.75*Pi) = Sqrt(0.5). Estimated error = ");
   Delta_2 := Sin (Make_e_Digit(5.0)*E_Quarter_Pi) + e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Sin (1.25*Pi) =-Sqrt(0.5). Estimated error = ");
   Delta_2 := Sin (Make_e_Digit(5.0)*E_Quarter_Pi) + e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Sin (1.5*Pi) = (-1.0).     Estimated error = ");
   Delta_2 := Sin (Make_e_Digit(6.0)*E_Quarter_Pi) + One;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Sin (1.75*Pi) =-Sqrt(0.5). Estimated error = ");
   Delta_2 := Sin (Make_e_Digit (7.0)*E_Quarter_Pi) + e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Sin (2.25*Pi) = Sqrt(0.5). Estimated error = ");
   Delta_2 := Sin (Make_e_Digit(9.0)*E_Quarter_Pi) - e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Sin (2.5*Pi) = 1.0.        Estimated error = ");
   Delta_2 := Sin (Make_e_Digit(10.0) * E_Quarter_Pi) - One;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Sin (2.75*Pi) = Sqrt(0.5). Estimated error = ");
   Delta_2 := Sin (Make_e_Digit(11.0)*E_Quarter_Pi) - e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Sin (3.0*Pi) = 0.0.        Estimated error = ");
   Delta_2 := Sin (Make_e_Digit(12.0)*E_Quarter_Pi);
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (3.0*Pi) = (-1.0).     Estimated error = ");
   Delta_2 := Cos (Make_e_Digit(12.0)*E_Quarter_Pi) + One;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (2.75*Pi) =-Sqrt(0.5). Estimated error = ");
   Delta_2 := Cos (Make_e_Digit(11.0)*E_Quarter_Pi) + e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (2.5*Pi) = 0.0.        Estimated error = ");
   Delta_2 := Cos (Make_e_Digit(10.0)*E_Quarter_Pi);
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (2.25*Pi) = Sqrt(0.5). Estimated error = ");
   Delta_2 := Cos (Make_e_Digit(9.0)*E_Quarter_Pi) - e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (2.0*Pi) = 1.0.        Estimated error = ");
   Delta_2 := Cos (Make_e_Digit(8.0)*E_Quarter_Pi) - One;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (1.75*Pi) = Sqrt(0.5). Estimated error = ");
   Delta_2 := Cos (Make_e_Digit(7.0)*E_Quarter_Pi) - e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (1.5*Pi) = 0.0.        Estimated error = ");
   Delta_2 := Cos (Make_e_Digit(6.0)*E_Quarter_Pi);
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (1.25*Pi) =-Sqrt(0.5). Estimated error = ");
   Delta_2 := Cos (Make_e_Digit(5.0)*E_Quarter_Pi) + e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (1.0*Pi) =-1.0.        Estimated error = ");
   Delta_2 := Cos (Make_e_Digit(4.0)*E_Quarter_Pi) + One;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (0.5*Pi) = 0.0.        Estimated error = ");
   Delta_2 := Cos (Two_Digit*E_Quarter_Pi);
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (0.25*Pi) = Sqrt(0.5). Estimated error = ");
   Delta_2 := Cos (E_Quarter_Pi) - Sqrt(Make_Extended(0.5));
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (-0.25*Pi) = Sqrt(0.5) Estimated error = ");
   Delta_2 := Cos (-E_Quarter_Pi) - Sqrt(Make_Extended(0.5));
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Cos (129.25*Pi) =-Sqrt(0.5). Estimated error = ");
   Delta_2 := Cos (Make_e_Digit(517.0)*E_Quarter_Pi) + e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check Sin (129.75*Pi) =-Sqrt(0.5). Estimated error = ");
   Delta_2 := Sin (Make_e_Digit (519.0)*E_Quarter_Pi) + e_Inverse_Sqrt_2;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   new_line(1);
   Pause ("Some more spot checks.");
   new_line(1);

   Z0 := E_Pi;
   for I in 1..10 loop
      Z0 := Z0 * Z0;
   end loop;
   Z1 := E_Pi**(2**10);

   put ("Check Pi**1024 - Pi*Pi*Pi*..*Pi*Pi*Pi: Estimated error = ");
   Delta_2 := One - Z1/Z0;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Check 2.0 * (1/SQRT(2.0)**2) - 1.0:    Estimated error = ");
   Delta_2 := e_Inverse_Sqrt_2**2 - (+0.5);
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Arccos(1/sqrt(2)) - Arcsin(1/sqrt(2)): Estimated error = ");
   Delta_2 := Arccos(E_Inverse_Sqrt_2) - Arcsin(E_Inverse_Sqrt_2);
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Arccos(-1/sqrt(2)) - 1.5 Pi/4:         Estimated error = ");
   Delta_2 := Arccos(-e_Inverse_Sqrt_2) - Two_Digit*e_Quarter_Pi - e_Quarter_Pi;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Arcsin(-1/sqrt(2)) + Pi/4:             Estimated error = ");
   Delta_2 := Arcsin(-e_Inverse_Sqrt_2) + e_Quarter_Pi;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Arccos(0) - Pi/2:                      Estimated error = ");
   Delta_2 := Arccos(Zero) - Two_Digit*e_Quarter_Pi;
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);

   put ("Arccos(-1) - Pi:                       Estimated error = ");
   Delta_2 := (Arccos(-One) - e_Pi);
   put (e_Real_Image (Delta_2, Aft => Integer'Min (10, e_Real_Decimals))); 
   new_line(1);


   Pause (
     "Test Exp (X+Y) = Exp (X) * Exp (Y) using integer valued Arguments.",
     "500_000 trials are performed for each number printed below:"
     );

   Test_Exp_w_Integer_Args;


   new_line;
   Pause (
      "Test of Sqrt routine.",
      "50_000 tests are performed for each number printed below:"
      );

   Test_Vector_Seed := 1.2345657891234E+134;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+74;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+27;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+8;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+1;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-1;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-4;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-8;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-14;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-34;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-70;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-173;
   Test_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-273;
   Test_Sqrt (Test_Vector_Seed);

   new_line;
   Pause (
      "Test of Sqrt routine using the ""**"" operator.",
      "10_000 tests are performed for each number printed below:"
      );

   Test_Vector_Seed := 1.2345657891234E+134;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+74;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+27;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+8;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+1;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-1;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-4;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-8;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-14;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-34;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-70;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-173;
   Test_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-273;
   Test_Sqrt_2 (Test_Vector_Seed);

   new_line;
   Pause (
      "Test of the standard ""/"" routine.",
      "Test calculates A - (A/B)*B.",
      "Use of both ""/"" and ""*"" usually means a loss of 2 guard digits accuracy.",
      "100_000 tests are performed for each number printed below:"
      );

   Test_Vector_Seed := 1.2345657891234E+273;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+123;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+10;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+08;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+1;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-1;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-3;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-8;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-18;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-28;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-39;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-78;
   Test_Divide_stnd(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-178;
   Test_Divide_stnd(Test_Vector_Seed);

   new_line;
   Pause (
      "Test of Newton-Raphson  Divide (A,B)  routine.",
      "Test calculates A - Divide (A,B)*B. The Newton-Raphson Divide is designed",
      "to minimize A - Divide (A,B)*B, and does so consistently better than the",
      "stnd ""/"" (schoolboy algorithm). It doesn't do the division more accurately.",
      "It is however quite a bit faster than the standard ""/"" operator.",
      "100_000 tests are performed for each number printed below:"
      );

   Test_Vector_Seed := 1.2345657891234E+273;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+123;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+10;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+08;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+1;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-1;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-3;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-8;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-18;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-28;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-39;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-78;
   Test_Divide(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-178;
   Test_Divide(Test_Vector_Seed);


   new_line;
   Pause (
      "Test of Reciprocal routine Reciprocal(A) = 1/A.",
      "Test calculates 1.0 - A*Reciprocal(A). The Newton-Raphson iteration is designed",
      "to minimize 1.0 - A*Reciprocal(A), and does so consistently better than using",
      "the standard ""/"" to get 1.0/A. It doesn't do the division more accurately.",
      "It is however quite a bit faster than the standard ""/"" operator.",
      "50_000 tests are performed for each number printed below:"
      );

   Test_Vector_Seed := 1.2345657891234E+273;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+123;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+10;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+08;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+1;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-1;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-3;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-8;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-18;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-28;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-39;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-78;
   Test_Reciprocal(Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-178;
   Test_Reciprocal(Test_Vector_Seed);
   Pause (
      "Test of Reciprocal_Sqrt routine.",
      "50_000 tests are performed for each number printed below:"
      );

   Test_Vector_Seed := 1.2345657891234E+273;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+173;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+27;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+8;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+1;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-1;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-4;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-8;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-14;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-34;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-70;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-173;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-273;
   Test_Reciprocal_Sqrt (Test_Vector_Seed);

   new_line;
   Pause (
      "Test of 1 / Sqrt routine using the ""**"" operator.",
      "10_000 tests are performed for each number printed below:"
      );

   Test_Vector_Seed := 1.2345657891234E+273;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+173;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+27;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+8;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E+1;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-1;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-4;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-8;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-14;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-34;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 9.9345657891234E-70;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-173;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);

   Test_Vector_Seed := 1.2345657891234E-273;
   Test_Reciprocal_Sqrt_2 (Test_Vector_Seed);


   Pause (
      "Test of Reciprocal_Nth_Root routine.",
      "5_000 tests are performed for each number printed below:"
      );

   for N in 1..4 loop

   new_line;
   put ("N = "); put(Integer'Image(N)); new_line;

   Test_Vector_Seed := 1.2345657891234E+123;
   Test_Root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E+27;
   Test_root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E+18;
   Test_root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E+8;
   Test_root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E+3;
   Test_root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E+1;
   Test_root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E-1;
   Test_root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E-3;
   Test_Root(Test_Vector_Seed, N);

   Test_Vector_Seed := 9.9345657891234E-14;
   Test_Root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E-17;
   Test_Root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E-33;
   Test_Root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E-70;
   Test_Root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E-112;
   Test_Root(Test_Vector_Seed, N);

   end loop;

   for M in 17..18 loop

   N := 123*M;

   new_line;
   put ("N = "); put(Integer'Image(N)); new_line;

   Test_Vector_Seed := 1.2345657891234E+34;
   Test_Root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E+27;
   Test_root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E+8;
   Test_root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E-8;
   Test_Root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E-17;
   Test_Root(Test_Vector_Seed, N);

   Test_Vector_Seed := 1.2345657891234E-33;
   Test_Root(Test_Vector_Seed, N);

   end loop;

   --  Reciprocal(Z2) is generally 2 or more times faster than "/"

   -- declare  Char : Character := ' '; begin
   --Z1     := Junk1 * (+1.23456789123E+12);
   --Z2     := Junk2 * (+1.234567891234E+00);
   --new_line;
   --put_Line ("Benchmark of division.  Enter a real 1.0 to start: ");
   --get (Char);
   --put_line("Start bench 3, 5_000_000 divisions: ");
   --for I in 1..5_000_000 loop
       --Z1 := Z1 / Z2;
   --end loop;
   --put_line("End of bench 3.");
   --new_line;
   --put_Line ("Benchmark of alt. division.  Enter a real 1.0 to start: ");
   --get (Char);
   --put_line("Start bench 4, 5_000_000 divisions: ");
   --for I in 1..5_000_000 loop
       --Z1 := Z1 * Reciprocal(Z2);
   --end loop;
   --put_line("End of bench 4.");
   --end;

end;

