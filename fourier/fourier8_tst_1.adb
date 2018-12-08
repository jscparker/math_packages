--***************************************************************
-- Procedure benchmarks FFT.  50000 calls to FFT and 50000 calls to
-- to the inverse FFT.
-- The FFT does 4N*ln(N) flops per call, where N=array length.
-- So 100000 FFT's do:
--  1840*10**6 flops where N=512.
--  4100*10**6 flops where N=1024.
--***************************************************************
with Fourier8;
with Text_IO;  use Text_IO;
with Ada.Numerics.Generic_Elementary_Functions;

procedure fourier8_tst_1 is

   type Real is digits 15;

   package math is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use math;
   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Integer);
   use iio;

 --No_of_Repetitions : constant := 25_000; -- std, as described above.
   No_of_Repetitions : constant := 1;

   Log_Of_Max_Data_Length : constant := 22;
   --  Means max data array size is 2**Log_Of_Max_Data_Length.
   --  The larger you make this array, the slower program might run.

   subtype Array_Index is Integer range 0..2**Log_Of_Max_Data_Length + 0;
   --  Sometimes its faster to add a bit to the end.

 --subtype Array_Index is Integer range 0..2**Log_Of_Max_Data_Length+1;
 --subtype Array_Index is Integer range 0..2**(Log_Of_Max_Data_Length+1)-2;

   type Data_Array is array(Array_Index) of Real;

   package fft8 is new Fourier8
      (Real, Array_Index, Data_Array, Log_Of_Max_Data_Length);
   use fft8;


 --Do_the_Bit_Reversal : constant Boolean := False;
   Do_the_Bit_Reversal : constant Boolean := True;


   D_Re : Data_Array;
   D_Im : Data_Array;

   Transformed_Data_Last, Data_Set_Last : Data_Index;
   Theta : Real;
   Inverse_No_Points, Inverse_No_Points2 : Real;
   Num                     : Integer;
   Error, D_Old, Max_Error : Real;
   Transformed_Data_Length : Real;
   
   Exp_Table : Exp_Storage;

 -- and add the following to the FFT parameter list:
 --
 --         Exp_Table  => Exp_Table,
 --

   ------------------
   -- Integer_Log2 --
   ------------------

   function Integer_Log2 (I : Array_Index)
      return Integer 
   is
      Logarithm : Integer := 0;
      Argument  : constant Integer := Integer(I);
      --  Rounds down, so on range (eg) 0..2**15-1 it
      --  returns 0 .. 14,  It returns 0 for I in 0..1
      --  returns 1 for I in 2..3, etc.
   begin
      for Exponent in 0..63 loop  
         exit when 2**Exponent > Argument; 
         Logarithm := Exponent;
      end loop;
      return Logarithm;
   end Integer_Log2;

begin

   new_line; 
   put ("Enter number of data points in the FFT.");
   new_line; 
   put ("Using a power of 2 like 1024 is good idea if you are doing benchmarks.");
   new_line; 
   put ("Enter the desired number: ");
   get (Num);
   new_line;
   
   Data_Set_Last := Array_Index (Num-1);
 
   --Inverse_No_Points := 1.0 / (Real (Data_Set_Last) + 1.0);
   -- But if Data_Set_Last+1 is not a power of 2, then FFT will pad data set
   -- out to the nearest power of two.  That is what we must divide by:
   
   Transformed_Data_Length := 2.0**(Integer_Log2(Data_Set_Last) + 1);
   Inverse_No_Points := 1.0 / Transformed_Data_Length;
   Inverse_No_Points2 := Inverse_No_Points**2;
 
   for I in Array_Index  loop
      D_Re(I)  := 0.0;
      D_Im(I)  := 0.0;
   end loop;
 
   for I in Array_Index range 0..Data_Set_Last loop
      Theta := Real (I);
      D_Re(I)  := 1.0 / (Theta + 1.0);
      D_Im(I)  := 1.0E-3 * Theta;
   end loop;
 
   for I in 1 .. No_of_Repetitions loop
     FFT (Data_Re                   => D_Re,
          Data_Im                   => D_Im,
          Transformed_Data_Last     => Transformed_Data_Last,
          Input_Data_Last           => Data_Set_Last,
          Exp_Table                 => Exp_Table,
          Inverse_FFT_Desired       => False,
          Normalized_Data_Desired   => False,
          Bit_Reversal_Desired      => Do_the_Bit_Reversal);
 
     FFT (Data_Re                   => D_Re,
          Data_Im                   => D_Im,
          Transformed_Data_Last     => Transformed_Data_Last,
          Input_Data_Last           => Transformed_Data_Last,
          Exp_Table                 => Exp_Table,
          Inverse_FFT_Desired       => True,
          Normalized_Data_Desired   => False,
          Bit_Reversal_Desired      => Do_the_Bit_Reversal);
 
     FFT (Data_Re                   => D_Re,
          Data_Im                   => D_Im,
          Transformed_Data_Last     => Transformed_Data_Last,
          Input_Data_Last           => Data_Set_Last,
          Exp_Table                 => Exp_Table,
          Inverse_FFT_Desired       => False,
          Normalized_Data_Desired   => False,
          Bit_Reversal_Desired      => Do_the_Bit_Reversal);
 
     FFT (Data_Re                   => D_Re,
          Data_Im                   => D_Im,
          Transformed_Data_Last     => Transformed_Data_Last,
          Input_Data_Last           => Transformed_Data_Last,
          Exp_Table                 => Exp_Table,
          Inverse_FFT_Desired       => True,
          Normalized_Data_Desired   => False,
          Bit_Reversal_Desired      => Do_the_Bit_Reversal);
 
     for J in 0..Transformed_Data_Last loop
        D_Re(J) := Inverse_No_Points2 * D_Re(J);
     end loop;
      
     for J in 0..Transformed_Data_Last loop
        D_Im(J) := Inverse_No_Points2 * D_Im(J);
     end loop;
      
   end loop;
 
   Max_Error := 0.0;  
   for J in 0..Data_Set_Last loop
      Theta := Real (J);
      D_Old  := 1.0 / (Theta + 1.0);
      Error := Abs (D_Re(J) - D_Old);
      if Max_Error < Error then Max_Error := Error; end if;
      
      D_Old  := 1.0E-3 * Theta;
      Error := Abs (D_Im(J) - D_Old);
      if Max_Error < Error then Max_Error := Error; end if;
   end loop;
   
   --Max_Error := 0.0;    
   if Data_Set_Last < Data_Index'Last then
      for J in Data_Set_Last+1 .. Transformed_Data_Last loop
         Error := Abs (D_Re(J) - 0.0);
         if Max_Error < Error then Max_Error := Error; end if;
         Error := Abs (D_Im(J) - 0.0);
         if Max_Error < Error then Max_Error := Error; end if;
      end loop;
    end if;
     
   if Do_the_Bit_Reversal then
      put("Maximum error = "); put(Max_Error); put (" ");
      new_line;
   end if;
   
end;
