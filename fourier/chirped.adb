
-------------------------------------------------------------------------------
-- package body Chirped, a chirped fast fourier transform.
-- Copyright (C) 1995-2018 Jonathan S. Parker
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
-------------------------------------------------------------------------------

with Fourier8;
with Text_IO; use Text_IO;
with Ada.Numerics.Generic_Elementary_Functions;

package body Chirped is

  package mth is new Ada.Numerics.Generic_Elementary_Functions (Real);
  use mth;

  -- Make an FFT that operates on an array with index Array_Index,
  -- (which has twice the range of Data_Index).

  package fft8 is new 
    Fourier8 (Real, Array_Index, Data_Array, Log_Of_Max_Data_Length+1);
  use fft8;


  Exp_Table : Exp_Storage;


  subtype Exponent_Of_Two_Type is Natural range 0..Log_Of_Max_Data_Length;


  -- Global tables of SIN's and COS's used by FFT_Chirped:

  type Exp_Function is record
    Re : Real;
    Im : Real;
  end record;
  type Sinusoid_Storage is Array(Data_Index) of Exp_Function;

  Zero : constant Exp_Function := (0.0, 0.0);

  W_Table : Sinusoid_Storage := (others => Zero);

  --  The following is the set of chirped exponentials Exp(i*(2Pi/N)*K*K/2).
  --  It is constructed on the first call to FFT_Chirped, and only rebuilt
  --  if the size of the data set is changed, resulting in a change to
  --  the global variable: Table_Status.Of_Current_Size_Of_W_Table.

  V_Re : Data_Array := (others => 0.0);
  V_Im : Data_Array := (others => 0.0);

  --  Following record holds info on the status of the arrays V and W_Table.
  --  Sometimes they have to updated, and whether or not they do depends
  --  on the settings held in the record Table_Status:

  type Table_Status is record
     Is_Not_Yet_Initialized                       : Boolean    := True;
     Inverse_FFT_Was_Desired_During_Previous_Call : Boolean    := False;
     Of_Current_Size_Of_W_Table                   : Data_Index := 0;
  end record;

  Memory : Table_Status;

  -----------------------
  -- Construct_W_Table --
  -----------------------

  -- The procedure calculates a table of W = Exp(-i(2Pi/N)*(K*K)/2)
  -- where K goes from 0..N-1, and N = LastIndexOfData + 1.
  --
  procedure Construct_W_Table 
    (Last_Data_Point : in Data_Index)
  is
     C, S, Theta : Real;

     pragma Assert(Real'Digits >= 15);

     Pii            : constant Real := 3.14159_26535_89793_23846_26433_83279_50288;
     TwoN           : constant Real := 2.0 * (Real(Last_Data_Point) + 1.0);
     Two_Pi_Over_2N : constant Real := 2.0 * Pii / TwoN;

     type Integer_64 is range 0 .. 2**63-1;
     B    : constant Integer_64 := 2 * (Integer_64 (Last_Data_Point) + 1); -- 2*N
     Base : Integer_64 := Integer_64 (Real'Ceiling (Sqrt (Real(B))));

     pragma Assert (B <= 2**40); -- K_Squared_Modulo_B needs this.

     ------------------------
     -- K_Squared_Modulo_B --
     ------------------------

     function K_Squared_Modulo_B (K : Data_Index)
        return Real
     is
        e, f, t1, t2, t3 : Integer_64;
        K64  : constant Integer_64 := Integer_64 (K);
     begin

        --  K is in 0..N-1, and B = 2*N

        if K64 < 2**31-1 then
           return Real ((k64 * k64) mod B);
        end if;

        if Base**2 < B then Base := Base + 1; end if;

        e  := K64 mod Base;
        f  := K64  /  Base;

        -- K = e + f*Base
        -- so K*K = e*e + 2*e*f*Base + f*f*Base**2

        t1 := (e*e) mod B;
        t2 := (e*f*2*Base) mod B;
        t3 := (f*f*(Base**2 mod B)) mod B;-- Base**2 mod B <= 2*Base

        return Real ((t1 + t2 + t3) mod B);

     end K_Squared_Modulo_B;

  begin

     W_Table (Data_Index'First) := (1.0, 0.0);

     -- For accuracy in the calculation of Cos ((2Pi)*(K*K)/(2N))
     -- use the equivalent quantity: Cos (2Pi * ((K*K) Mod (2N))/(2N))
     -- which equals Cos ((2Pi/(2N)) * (K*K) Mod (2N))
     -- (using A/B = (A Mod B) / B + Some_Integer.)

     for K in Data_Index range 1 .. Last_Data_Point loop
        Theta    := Two_Pi_Over_2N * K_Squared_Modulo_B (K); -- B = 2N
        C        := Cos (Theta);
        S        := Sin (Theta);
        W_Table (K) := (C, -S);
     end loop;

  end Construct_W_Table;

  -----------------
  -- FFT_Chirped --
  -----------------

  procedure FFT_Chirped
    (Data_Re, Data_Im        : in out Data_Array;
     Input_Data_Last         : in     Data_Index;
     Inverse_FFT_Desired     : in     Boolean := False;
     Normalized_Data_Desired : in     Boolean := True)
  is
     L_Minus_1, L_Temp, Half_Data  : Array_Index;
     Test                          : Array_Index;
     Table_Index, Starting_Index   : Array_Index;
     Last_Data_Point               : constant Array_Index := Input_Data_Last;
     Log_Of_Data_Size_Minus_1      : Exponent_Of_Two_Type;
     NormFactor                    : Real;
     W_Table_Was_Reconstructed     : Boolean;
     Choice_Of_Inverse_Has_Changed : Boolean;

     W_Re, W_Im, D_Re, D_Im, V1_Re, V1_Im : Real;

     ------------------
     -- Integer_Log2 --
     ------------------
 
     -- Rounds down, so on range (eg) 0..2**15-1 it returns 0 .. 14
     -- returns 0 for I in 0..1 and returns 1 for I in 2..3, etc.
 
     function Integer_Log2 (I : Array_Index) 
        return Exponent_Of_Two_Type 
     is
        Log2 : Exponent_Of_Two_Type := 0;
     begin
        for Exponent in Natural loop
           exit when 2**Exponent > I;
           Log2 := Exponent;
        end loop;
        return Log2;
     end Integer_Log2;

  begin

     -- Below use the variable names in the reference given above.

     if Last_Data_Point = 0 then return; end if;
  
     if Last_Data_Point > Data_Index'Last then
        put_line ("Can't FFT a data set that large. Input_Data_Last is too big.");
        raise Constraint_Error;
     end if;
     --  If checks are suppressed, then this check is very helpful.
  
     --**************************************************************
     -- Step 1.  We use a radix 2 FFT to calculate the convolution
     -- used in the chirped FFT.  Must find the first
     -- power of two (L) that is greater than or equal to 2*N-1 where
     -- N = Last_Data_Point+1.  So L is to equal the 1st
     -- power of 2 that is greater than or equal 2*Last_Data_Point + 1.
     -- What really want is L-1, which call L_Minus_1.  Notice
     -- that if Last_Data_Point = 255, then L_Minus_1 = 511.
     -- If Last_Data_Point = 256, then L_Minus_1 = 1023.
     --**************************************************************
  
     L_temp := 2*Last_Data_Point + 1;
     Log_Of_Data_Size_Minus_1 := Integer_Log2 (L_temp);
  
     --  Integer_Log2 rounds down so LogOf.. is in Exponent_Of_Two_Type range.
     --  This is LogOfDataSize - 1, not LogOf(DataSize-1).
     --  The length of the padded data is twice 2**Log_Of_Data_Size_Minus_1
  
     Half_Data := Array_Index (2**Integer(Log_Of_Data_Size_Minus_1));
  
     L_Minus_1 := Half_Data + (Half_Data - 1);
  
     --**************************************************************
     -- Step 2. Construct the W table, if a Table of the correct length
     -- has not already been constructed.  Perform Step 3 whenever we
     -- perform step 2. (And also if the setting of Inverse_FFT_Desired
     -- changes from what it was in the previous call.)
     --**************************************************************
  
     W_Table_Was_Reconstructed := False;
  
     if Memory.Of_Current_Size_Of_W_Table /= Last_Data_Point
          			or Memory.Is_Not_Yet_Initialized then
  
        Construct_W_Table (Last_Data_Point);
  
        W_Table_Was_Reconstructed         := True;
        Memory.Of_Current_Size_Of_W_Table := Last_Data_Point;
        Memory.Is_Not_Yet_Initialized     := False;
  
     end if;
  
     --**************************************************************
     -- Step 3. FFT the global chirp array V.  Don't remake it if it is already
     -- OK. V depends on size of Last_Data_Point. If size of Last_Data_Point
     -- not constant between calls to FFT_Chirped, then V and FFT(V) must be
     -- reconstructed.  Also, if the setting of Inverse_FFT_Desired changes
     -- from the previous call to FFT_Chirped then reconstruct V and FFT(V).
     -- V is global array that is often reused between calls.
     -- (W_Table is also a global that is saved between calls.)
     --**************************************************************
  
     Choice_Of_Inverse_Has_Changed :=
      (Inverse_FFT_Desired /= Memory.Inverse_FFT_Was_Desired_During_Previous_Call);
  
     if W_Table_Was_Reconstructed or Choice_Of_Inverse_Has_Changed then
  
        if Inverse_FFT_Desired then
           for I in Data_Index range 0..Last_Data_Point loop
              V_Re(I) := W_Table(I).Re;
              V_Im(I) := W_Table(I).Im;
           end loop;
        else
           for I in Data_Index range 0..Last_Data_Point loop
              V_Re(I) :=  W_Table(I).Re;
              V_Im(I) := -W_Table(I).Im;
           end loop;
        end if;
  
        Starting_Index := L_Minus_1 - Last_Data_Point + 1;
        for I in Array_Index range Starting_Index..L_Minus_1 loop
           Table_Index := L_Minus_1 - I;
           Table_Index := Table_Index + 1;
           if Inverse_FFT_Desired then
              V_Re(I) := W_Table(Table_Index).Re;
              V_Im(I) := W_Table(Table_Index).Im;
           else
              V_Re(I) :=  W_Table(Table_Index).Re;
              V_Im(I) := -W_Table(Table_Index).Im;
           end if;
        end loop;
  
        if Starting_Index > Last_Data_Point+1 then
           for I in Last_Data_Point+1..Starting_Index-1 loop
              V_Re(I) := 0.0; V_Im(I) :=  0.0;
           end loop;
        end if;
  
        FFT 
          (Data_Re                   => V_Re,
           Data_Im                   => V_Im,
           Transformed_Data_Last     => Test,
           Input_Data_Last           => L_Minus_1,
           Exp_Table                 => Exp_Table,
           Inverse_FFT_Desired       => False,
           Normalized_Data_Desired   => False);
  
     end if;
  
     -- Update the global memory of the status of the Inverse_FFT_Desired setting.
  
     Memory.Inverse_FFT_Was_Desired_During_Previous_Call := Inverse_FFT_Desired;
  
     --**************************************************************
     -- Step 4. Multiply Data by W, put the product into Data, and FFT it:
     --**************************************************************
  
     if Inverse_FFT_Desired then
        for I in Array_Index range 0..Last_Data_Point loop
           --  Data(I) := Data(I) * Conjugate (W_Table (I));
           W_Re := W_Table(I).Re; W_Im := -W_Table(I).Im;
           D_Re := Data_Re(I);    D_Im := Data_Im(I);
           Data_Re(I) := D_Re*W_Re - D_Im*W_Im;
           Data_Im(I) := D_Re*W_Im + D_Im*W_Re;
        end loop;
     else
        for I in Array_Index range 0..Last_Data_Point loop
           --  Data(I) := Data(I) * W_Table (I);
           W_Re := W_Table(I).Re; W_Im := W_Table(I).Im;
           D_Re := Data_Re(I);    D_Im := Data_Im(I);
           Data_Re(I) := D_Re*W_Re - D_Im*W_Im;
           Data_Im(I) := D_Re*W_Im + D_Im*W_Re;
        end loop;
     end if;
  
     for I in Array_Index range Last_Data_Point+1..L_Minus_1 loop
        Data_Re(I) := 0.0; Data_Im(I) := 0.0;
     end loop;
  
     FFT 
       (Data_Re                   => Data_Re,
        Data_Im                   => Data_Im,
        Transformed_Data_Last     => Test,
        Input_Data_Last           => L_Minus_1,
        Exp_Table                 => Exp_Table,
        Inverse_FFT_Desired       => False,
        Normalized_Data_Desired   => False);
  
     --**************************************************************
     -- Step 5.  Multiply Data by V.  Call the product Data,
     -- and then inverse FFT it.  (This completes the convolution.)
     --**************************************************************
  
     for I in Array_Index range 0..L_Minus_1 loop
        V1_Re:= V_Re(I);   V1_Im := V_Im(I);
        D_Re := Data_Re(I); D_Im := Data_Im(I);
        Data_Re(I) := D_Re*V1_Re - D_Im*V1_Im;
        Data_Im(I) := D_Re*V1_Im + D_Im*V1_Re;
     end loop;
  
     FFT 
       (Data_Re                   => Data_Re,
        Data_Im                   => Data_Im,
        Transformed_Data_Last     => Test,
        Input_Data_Last           => L_Minus_1,
        Exp_Table                 => Exp_Table,
        Inverse_FFT_Desired       => True,
        Normalized_Data_Desired   => False);
  
     --**************************************************************
     -- Step 6. Multiply result by the appropriate Chirped exp's, and
     -- simultaneously put the FFT'd data back into the output array Data.
     --**************************************************************
  
     if Inverse_FFT_Desired then
        for I in Data_Index range 0..Last_Data_Point loop
           --  Data(I) := Data(I) * Conjugate (W_Table(I));
           W_Re := W_Table(I).Re; W_Im := -W_Table(I).Im;
           D_Re := Data_Re(I);    D_Im := Data_Im(I);
           Data_Re(I) := D_Re*W_Re - D_Im*W_Im;
           Data_Im(I) := D_Re*W_Im + D_Im*W_Re;
        end loop;
     else
        for I in Data_Index range 0..Last_Data_Point loop
           --  Data(I) := Data(I) * W_Table (I);
           W_Re := W_Table(I).Re; W_Im := W_Table(I).Im;
           D_Re := Data_Re(I);    D_Im := Data_Im(I);
           Data_Re(I) := D_Re*W_Re - D_Im*W_Im;
           Data_Im(I) := D_Re*W_Im + D_Im*W_Re;
        end loop;
     end if;
  
     for I in Last_Data_Point+1..Array_Index'Last loop
        Data_Re(I) := 0.0; Data_Im(I) := 0.0;
     end loop;
  
     --**************************************************************
     -- Step 7. Normalize data if so desired.  We have saved some time
     -- by failing to normalize in the above FFT's.  Two of the three
     -- calls should have normalized, each dividing the data by
     -- Sqrt (L_Minus_1 + 1).  In addition, must divide by another
     -- factor: Sqrt (Last_Data_Point + 1).  We do them all together:
     --**************************************************************
  
     if Normalized_Data_Desired then
        NormFactor := (Real(L_Minus_1) + 1.0) * Sqrt (Real(Last_Data_Point) + 1.0);
        NormFactor := 1.0 / NormFactor;
        for I in Data_Index range 0..Last_Data_Point loop
           Data_Re(I) := NormFactor * Data_Re(I);
           Data_Im(I) := NormFactor * Data_Im(I);
        end loop;
     end if;
  
  end FFT_Chirped;
  
end Chirped;

