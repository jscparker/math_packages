
-- Procedure tests DFT_Chirped.  Code fragments
-- demonstrate some features of DFT's.

with Chirped;
with Text_IO; use Text_IO;
with Ada.Numerics.Generic_Elementary_Functions;

procedure chirped_demo_1 is

   type Real is digits  15;

   package mth is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use mth;
   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Integer);
   use iio;

   Log_Of_Max_Data_Length : constant := 12;
   --  ***NOTICE this means max data array size is 4096.***
   --  The larger you make this array, the slower program
   --  might run.

   type Array_Index is range 0..2**(Log_Of_Max_Data_Length+1)-1;
 --subtype Array_Index is Integer range 0..2**(Log_Of_Max_Data_Length+1)-1;
   type Data_Array is array(Array_Index) of Real;

   package fft is new Chirped
     (Real, Array_Index, Data_Array, Log_Of_Max_Data_Length);
   use fft;

   D_Re, D_Im : Data_Array;

   Data_Set_Last : Data_Index;
   Basis_Function_Last      : Data_Index;
   Theta, Frequency, Norm   : Real;
   Two_Pi_Over_N            : Real;
   Num, Mode                : Integer;
   Area1, Area2             : Real;
   DeltaRE, DeltaIM         : Real;
   Dat_Re, Dat_Im           : Real;

   Max_Error, Del : Real;

   Pii : constant Real := 3.14159_26535_89793_23846;

   -----------
   -- Pause --
   -----------

   procedure Pause (s1,s2,s3,s4,s5,s6,s7,s8,s9,s10 : string := "") is
      Continue : Character := ' ';
   begin
      new_line(2);
      if S1 /= "" then put_line (S1); end if;
      if S2 /= "" then put_line (S2); end if;
      if S3 /= "" then put_line (S3); end if;
      if S4 /= "" then put_line (S4); end if;
      if S5 /= "" then put_line (S5); end if;
      if S6 /= "" then put_line (S6); end if;
      if S7 /= "" then put_line (S7); end if;
      if S8 /= "" then put_line (S8); end if;
      if S9 /= "" then put_line (S9); end if;
      if S10 /= "" then put_line (S10); end if;

      dialog: loop
	 begin
            new_line(1);
	    put ("Type a character to continue: ");
	    get_immediate (Continue);
	    exit dialog;
	 exception
	    when others => null;
	 end;
      end loop dialog;
      new_line(2);
   end pause;

begin
  new_line;
  put ("Input number of data points to use in the following tests.");
  new_line;
  put ("Maximum allowed number is: ");
  put (2**Log_Of_Max_Data_Length);
  new_line;
  put ("A nice prime number will demonstate the algorithm, for example 37.");
  new_line;
  put ("Enter the desired number: ");
  get (Num);

  Basis_Function_Last := Data_Index (Num - 1);

  Data_Set_Last := Basis_Function_Last;

  --*****************************************************************
  -- Test 1.
  -- Make a basis element of the set of Fourier basis states:
  -- Its DFT should zero everywhere except for the element = Mode.
  -- There, the value of the DFT should be (1.0, 0.0).  Below
  -- we normalize the basis element with the factor Norm = 1/Sqrt(N).
  --*****************************************************************

  Pause
   ("Test 1: Use the Chirped_FFT to take the",
    "Discrete Fourier Transform (DFT) the Fourier basis function:",
    "  ",
    "         Exp (i*T*Mode*2*Pi/N) / SQRT(N),",
    "  ",
    "where N is the number of points.  The result of the DFT",
    "should be a delta-function peaked at position Mode.");

  put_line ("Input Mode (0, 1, 2, ...) of the basis function you wish to DFT: ");
  put_line ("For example, if you are using 37 points, try entering mode 31. ");
  put ("Enter the desired number: ");
  Get (Mode);
  Frequency := Real (Mode);

  --*****************************************************************

  Two_Pi_Over_N := 2.0 * Pii / (Real (Basis_Function_Last) + 1.0);
  Norm          := 1.0 / SQRT ((Real (Basis_Function_Last) + 1.0));
  for I in 0..Basis_Function_Last loop
       Theta   := Two_Pi_Over_N * Frequency * Real (I);
       D_Re(I) := Norm * Cos (Theta);
       D_Im(I) := Norm * Sin (Theta);
  end loop;

  put_line ("Starting Test 1:");

  new_line; put ("Basis function was constructed using "); put(Num);
  put(" data points.");

  FFT_Chirped
      (Data_Re                  => D_Re,
       Data_Im                  => D_Im,
       Input_Data_Last          => Basis_Function_Last,
       Inverse_FFT_Desired      => False,
       Normalized_Data_Desired  => True);

  Pause
   ("Ending the Discrete Fourier Transform (DFT).",
    "We have taken the DFT of Exp (i*T*Mode*2*Pi/N) / SQRT(N).",
    "The result should be (0.0, 0.0) everywhere except for ",
    "a (1.0, 0.0) at data point Mode.");

  for I in Data_Index range 0..Basis_Function_Last loop
    new_line; put (Integer(I)); put(' ');
    put (D_Re(I)); put(' '); put(D_Im(I));
  end loop;


  --*****************************************************************
  -- Test 1b.
  -- We 1st make a basis element of the set of Fourier basis states:
  -- Its DFT should zero everywhere except for the element = Mode.
  -- There, the value of the DFT should be (1.0, 0.0).  Below
  -- we normalize the basis element with the factor Norm = 1/Sqrt(N).
  --*****************************************************************

  Pause
   ("Test 1b: Use the Chirped_FFT to take the inverse",
    "Discrete Fourier Transform (DFT) the Fourier basis function:",
    "  ",
    "         Exp (-i*T*Mode*2*Pi/N) / SQRT(N),",
    "  ",
    "where N is the number of points. The result of the DFT",
    "should be a delta-function peaked at position Mode.");

  put_line ("Input Mode (0, 1, 2, ...) of the basis function you wish to DFT.");
  put_line ("For example, if you are using 37 points, try entering mode 31.");
  put ("Enter the desired number: ");
  Get (Mode);
  Frequency := Real (Mode);



  Two_Pi_Over_N := 2.0 * Pii / (Real (Basis_Function_Last) + 1.0);
  Norm          := 1.0 / SQRT ((Real (Basis_Function_Last) + 1.0));
  for I in 0..Basis_Function_Last loop
     Theta   :=  Two_Pi_Over_N * Frequency * Real (I);
     D_Re(I) :=  Norm * Cos (Theta);
     D_Im(I) := -Norm * Sin (Theta);
  end loop;

  put_line ("Starting Test 1b:");

  new_line; put ("Basis function was constructed using "); put(Num);
  put(" data points.");

  FFT_Chirped
    (Data_Re                  => D_Re,
     Data_Im                  => D_Im,
     Input_Data_Last          => Basis_Function_Last,
     Inverse_FFT_Desired      => True,
     Normalized_Data_Desired  => True);

  Pause
   ("Ending the Discrete Fourier Transform (DFT).",
    "We have taken the inverse DFT of Exp (-i*T*Mode*2*Pi/N) / SQRT(N).",
    "The result should be (0.0, 0.0) everywhere except for ",
    "a (1.0, 0.0) at data point Mode.");

  for I in Data_Index range 0..Basis_Function_Last loop
     new_line; put (Integer(I)); put(' ');
     put (D_Re(I)); put(' '); put(D_Im(I));
  end loop;

  --*****************************************************************
  --  Test2.
  --  Verify Parceval's Theorem.  The area under the DFT coefficients
  --  should equal the area under the data.  Notice that was true in the
  --  above test, because Sum (|exp(i*T*Mode*2*Pi/N) / Sqrt(N)|**2) is
  --  one on the interval [0,N-1].
  --  First we create a new data set:
  --*****************************************************************

  for I in 0..Data_Set_Last loop
     Theta   := Real (I);
     D_Re(I) := 1.0 / (Theta + 1.0);
     D_Im(I) := Theta**2;
  end loop;

  Area1 := 0.0;
  for I in 0..Data_Set_Last loop
      Area1 := Area1 + D_Re(I)*D_Re(I) + D_Im(I)*D_Im(I);
  end loop;

  Pause
  ("Test 2: test of Parceval's theorem.  The area under the",
   "DFT curve (sum of coefficients modulus squared) should equal",
   "area under the original curve.  This requires the use of the",
   "normalized version of the DFT.");

  put_line ("Starting Test 2:");

  FFT_Chirped
    (Data_Re                  => D_Re,
     Data_Im                  => D_Im,
     Input_Data_Last          => Data_Set_Last,
     Inverse_FFT_Desired      => False,
     Normalized_Data_Desired  => True);

  put_line ("Ending DFT");

  Area2 := 0.0;
  for I in 0..Data_Set_Last loop
      Area2 := Area2 + D_Re(I)*D_Re(I) + D_Im(I)*D_Im(I);
  end loop;

  new_line; put ("Area under the original curve:        "); put(Area1);
  new_line; put ("Area under Fourier transformed curve: "); put(Area2);

  --  Test3.
  --  Inverse DFT of the DFT.

  Pause
   ("Test 3: take the inverse DFT of the DFT of some",
    "artificial data, and compare with the original data.",
    "The test calculates: Data - Inverse_DFT (DFT (Data),",
    "then prints the max error.");

  for I in 0..Data_Set_Last loop
     Theta   := Real (I);
     D_Re(I) := Sqrt (Theta) / (Theta + 1.0);
     D_Im(I) := Cos (0.03737*Theta**2/(Theta + 1.0)) * Theta / (Theta + 1.0);
  end loop;

  put_line("Starting Test 3:");

  FFT_Chirped
    (Data_Re                  => D_Re,
     Data_Im                  => D_Im,
     Input_Data_Last          => Data_Set_Last,
     Inverse_FFT_Desired      => False,
     Normalized_Data_Desired  => True);

  FFT_Chirped
    (Data_Re                  => D_Re,
     Data_Im                  => D_Im,
     Input_Data_Last          => Data_Set_Last,
     Inverse_FFT_Desired      => True,
     Normalized_Data_Desired  => True);


  Max_Error := 0.0;
  for I in 0..Data_Set_Last loop

     Theta  := Real (I);
     Dat_Re := Sqrt (Theta) / (Theta + 1.0);
     Dat_Im := Cos (0.03737*Theta**2/(Theta + 1.0)) * Theta / (Theta + 1.0);

     Del    := Abs (Dat_Re - D_Re(I));
     if Max_Error < Del then Max_Error := Del; end if;
     Del := Abs (Dat_Im - D_Im(I));
     if Max_Error < Del then Max_Error := Del; end if;

  end loop;

  new_line(2);
  put ("Max error in  Data - Inverse_DFT (DFT (Data)):"); put (Max_Error);
  new_line;


  Pause
   ("Test 4 is a long test..runs through entire range of Array_Index ",
    "to see if DFT_Inverse (DFT (Data)) = Data. Nothing is printed",
    "unless large errors are detected.  (Use 15 digit floating point.)",
    "Interrupt the program here if you do not want a long wait.");

  for N in Array_Index range 0..Data_Set_Last loop

     for I in Array_Index range 0..N loop
        Theta   := Real (I);
        D_Re(I) := 1.0 / (Theta + 1.0);
        D_Im(I) := Cos (0.03737 * Theta**2 / (Theta + 1.0));
     end loop;

     FFT_Chirped
       (Data_Re                  => D_Re,
        Data_Im                  => D_Im,
        Input_Data_Last          => N,
        Inverse_FFT_Desired      => False,
        Normalized_Data_Desired  => True);

     FFT_Chirped
       (Data_Re                  => D_Re,
        Data_Im                  => D_Im,
        Input_Data_Last          => N,
        Inverse_FFT_Desired      => True,
        Normalized_Data_Desired  => True);

     for I in Array_Index range 0..N loop
        Theta   := Real (I);
        Dat_Re  := 1.0 / (Theta + 1.0);
        Dat_Im  := Cos (0.03737 * Theta**2 / (Theta + 1.0));
        DeltaRE := Abs (Dat_Re - D_Re(I));
        DeltaIM := Abs (Dat_Im - D_Im(I));
        if DeltaRE > 1.0E-11 then
           new_line;
           put("FAILURE: "); put (DeltaRE); put(" at "); put(Real(I));
        end if;
        if DeltaIM > 1.0E-11 then
           new_line;
           put("FAILURE: "); put (DeltaIM); put(" at "); put(Real(I));
        end if;
     end loop;

  end loop;

end;
