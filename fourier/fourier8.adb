
-------------------------------------------------------------------------------
-- package body Fourier8, a fast fourier transform
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

with Ada.Numerics.Generic_Elementary_Functions;

package body Fourier8 is

  package math is new Ada.Numerics.Generic_Elementary_Functions (Real);
  use math;

  subtype Exponent_Of_Two_Type is Natural range 0 .. Log_Of_Max_Data_Length;

  --------------------
  -- Make_Exp_Table --
  --------------------

  --  Make_Exp_Table calculates points on the unit circle
  --  starting at Theta = 0 and ending near Theta = -Pi. The
  --  points are separated by Theta = 2*Pi/N.  In other words, the
  --  bottom half of the unit circle is divided into N/2 equally spaced
  --  points.  These points are then stored in the array Exp_Table.
  --  For efficiency, only the 1st eighth of the unit circle
  --  is calculated.  The other points are inferred from that. Exp_Table
  --  is defined Exp (-i * Two_Pi_Over_N * Data_Index) where Data_Index
  --  is in the range 0 .. N/2 - 1 and N is the Data length, the power of
  --  2 that appears in Two_Pi_Over_N.

  procedure Make_Exp_Table
    (Padded_Data_Index_Last : in     Data_Index;
     Exp_Table              : in out Exp_Storage)
  is
     Pii        : constant Real := 3.14159_26535_89793_23846_26433_83279_50288;
     Sqrt_Half  : constant Real := 0.70710_67811_86547_52440_08443_62104_84904;
     Two_Pi_Over_N : constant Real := 2.0 * Pii / (Real (Padded_Data_Index_Last) + 1.0);
     C, S       : Real;
     Theta      : Real;

     -- Divide the unit circle into 2**M = Padded_Data_Index_Last+1 points.
     -- Divide the unit circle into eighths and quadrants:
     StartOf1stQuadrant : constant Data_Index := 0;
     StartOf2ndEighth   : constant Data_Index := Padded_Data_Index_Last / 8 + 1;
     StartOf2ndQuadrant : constant Data_Index := Padded_Data_Index_Last / 4 + 1;
     StartOf3rdQuadrant : constant Data_Index := Padded_Data_Index_Last / 2 + 1;
     StartOf4rthEighth  : constant Data_Index := StartOf2ndEighth * 3;
  begin

     if    Padded_Data_Index_Last = 0 then          -- 1 data point
	null;
     elsif Padded_Data_Index_Last = 1 then          -- 2 data points
        Exp_Table.Re(StartOf1stQuadrant) :=  1.0;
        Exp_Table.Im(StartOf1stQuadrant) :=  0.0;
     elsif Padded_Data_Index_Last = 3 then          -- 4 data points
        Exp_Table.Re(StartOf1stQuadrant) :=  1.0;
        Exp_Table.Im(StartOf1stQuadrant) :=  0.0;
        Exp_Table.Re(StartOf2ndQuadrant) :=  0.0;
        Exp_Table.Im(StartOf2ndQuadrant) := -1.0;
     elsif Padded_Data_Index_Last >= 7 then         -- 8 data points and higher:
        Exp_Table.Re(StartOf1stQuadrant) :=  1.0;
        Exp_Table.Im(StartOf1stQuadrant) :=  0.0;
        Exp_Table.Re(StartOf2ndEighth)   :=  Sqrt_Half;
        Exp_Table.Im(StartOf2ndEighth)   := -Sqrt_Half;
        Exp_Table.Re(StartOf2ndQuadrant) :=  0.0;
        Exp_Table.Im(StartOf2ndQuadrant) := -1.0;
        Exp_Table.Re(StartOf4rthEighth)  := -Sqrt_Half;
        Exp_Table.Im(StartOf4rthEighth)  := -Sqrt_Half;

        for Mode in Exp_Mode_Index range 1..StartOf2ndEighth-1 loop
           Theta :=  Two_Pi_Over_N * Real(Mode);
	   C     := Cos (Theta);
	   S     := Sin (Theta);
	   Exp_Table.Re(StartOf1stQuadrant + Mode) :=  C;
	   Exp_Table.Im(StartOf1stQuadrant + Mode) := -S;
	   Exp_Table.Re(StartOf2ndQuadrant + Mode) := -S;
	   Exp_Table.Im(StartOf2ndQuadrant + Mode) := -C;
	   Exp_Table.Re(StartOf2ndQuadrant - Mode) :=  S;
	   Exp_Table.Im(StartOf2ndQuadrant - Mode) := -C;
	   Exp_Table.Re(StartOf3rdQuadrant - Mode) := -C;
	   Exp_Table.Im(StartOf3rdQuadrant - Mode) := -S;
	end loop;

     end if;

     Exp_Table.Current_Size_Of_Exp_Table := Padded_Data_Index_Last;
     --  This global variable is initalized to 0. This tells the FFT routine
     --  whether to reconstruct Exp_Table when it's called again.  The table
     --  must be reconstructed whenever Padded_Data_Index_Last changes.

  end Make_Exp_Table;

  --------------------------
  --  Get_Bit_Reversal_Of --
  --------------------------

  -- Radix 4 Bit_Reversal using bits 0..Top_Bit.

  procedure Get_Bit_Reversal_Of
    (Data_Re, Data_Im            : in out Data_Array;
     Two_To_The_Top_Bit          : in     Data_Index;
     Half_Two_To_The_Top_Bit     : in     Data_Index;
     Quarter_Two_To_The_Top_Bit  : in     Data_Index;
     Top_Bit                     : in     Exponent_Of_Two_Type)
  is
     Temp  : Real;
     K, Bit_Reversed_K, Power_Of_Two : Data_Index;
     Two_To_Top_Bit_minus_Two_To_lesser_Bit : constant Data_Index
                              := Two_To_The_Top_Bit - Half_Two_To_The_Top_Bit;
     Two_To_Top_Bit_plus_Two_To_lesser_Bit : constant Data_Index
                              := Two_To_The_Top_Bit + Half_Two_To_The_Top_Bit;
  begin

     -- No_Of_Data_points = 2**(Top_Bit+1) <= 2:

     if Top_Bit = 0 then
        return;
     end if;

     -- No_Of_Data_points = 2**(Top_Bit+1) = 4:

     if Top_Bit = 1  then
        Temp       := Data_Re(2);
        Data_Re(2) := Data_Re(1);
        Data_Re(1) := Temp;
        Temp       := Data_Im(2);
        Data_Im(2) := Data_Im(1);
        Data_Im(1) := Temp;
        return;
     end if;

     -- No_Of_Data_points = 2**(Top_Bit+1) >= 8:

     K              := 0;
     Bit_Reversed_K := 0;

     for count in 0..Half_Two_To_The_Top_Bit-2 loop

        --  now K = 00xxxxxxxxxxx     and  Bit_Reversed_K = xxxxxxxxxxxx00

        K              := K + 1;
        Bit_Reversed_K := Bit_Reversed_K + Two_To_The_Top_Bit;

        if K < Bit_Reversed_K then
           Temp                    := Data_Re (Bit_Reversed_K);
           Data_Re(Bit_Reversed_K) := Data_Re(K);
           Data_Re(K)              := Temp;
           Temp                    := Data_Im (Bit_Reversed_K);
           Data_Im(Bit_Reversed_K) := Data_Im(K);
           Data_Im(K)              := Temp;
        end if;

       --  now K = 10xxxxxxxxxxx    and  Bit_Reversed_K =  xxxxxxxxxxxx01

        K              := K + 1;
        Bit_Reversed_K := Bit_Reversed_K - Two_To_Top_Bit_minus_Two_To_lesser_Bit;

        if K < Bit_Reversed_K then
           Temp                    := Data_Re (Bit_Reversed_K);
           Data_Re(Bit_Reversed_K) := Data_Re(K);
           Data_Re(K)              := Temp;
           Temp                    := Data_Im (Bit_Reversed_K);
           Data_Im(Bit_Reversed_K) := Data_Im(K);
           Data_Im(K)              := Temp;
        end if;

       --  now K = 01xxxxxxxxxxx         xxxxxxxxxxxx10

        K              := K + 1;
        Bit_Reversed_K := Bit_Reversed_K + Two_To_The_Top_Bit;

        if K < Bit_Reversed_K then
           Temp                    := Data_Re (Bit_Reversed_K);
           Data_Re(Bit_Reversed_K) := Data_Re(K);
           Data_Re(K)              := Temp;
           Temp                    := Data_Im (Bit_Reversed_K);
           Data_Im(Bit_Reversed_K) := Data_Im(K);
           Data_Im(K)              := Temp;
        end if;

       --  now K = 11xxxxxxxxxxx         xxxxxxxxxxxx11

        K              := K + 1;
        Bit_Reversed_K := Bit_Reversed_K - Two_To_Top_Bit_plus_Two_To_lesser_Bit;

        Power_Of_Two   := Quarter_Two_To_The_Top_Bit;          --  2**(Top_Bit-2)

        Bit_Reverse_K_Plus_1:
        for Exponent in reverse 0..Top_Bit-2 loop
           if Bit_Reversed_K < Power_Of_Two then
              -- B_R_K has a '0' at position Exponent.  Put a '1' there:
              Bit_Reversed_K := Bit_Reversed_K + Power_Of_Two;
              exit Bit_Reverse_K_Plus_1;
           else
              -- B_R_K has a '1' at position Exponent.  Put a '0' there:
              Bit_Reversed_K := Bit_Reversed_K - Power_Of_Two;
              Power_Of_Two   := Power_Of_Two / 2;
           end if;
        end loop Bit_Reverse_K_Plus_1;

        if K < Bit_Reversed_K then
           Temp                    := Data_Re (Bit_Reversed_K);
           Data_Re(Bit_Reversed_K) := Data_Re(K);
           Data_Re(K)              := Temp;
           Temp                    := Data_Im (Bit_Reversed_K);
           Data_Im(Bit_Reversed_K) := Data_Im(K);
           Data_Im(K)              := Temp;
        end if;

     end loop;

     --  now K = 00xxxxxxxxxxx         xxxxxxxxxxxx00

     K              := K + 1;
     Bit_Reversed_K := Bit_Reversed_K + Two_To_The_Top_Bit;

     if K < Bit_Reversed_K then
        Temp                    := Data_Re (Bit_Reversed_K);
        Data_Re(Bit_Reversed_K) := Data_Re(K);
        Data_Re(K)              := Temp;
        Temp                    := Data_Im (Bit_Reversed_K);
        Data_Im(Bit_Reversed_K) := Data_Im(K);
        Data_Im(K)              := Temp;
     end if;

     --  now K = 10xxxxxxxxxxx         xxxxxxxxxxxx01

     K              := K + 1;
     Bit_Reversed_K := Bit_Reversed_K - Two_To_Top_Bit_minus_Two_To_lesser_Bit;

     if K < Bit_Reversed_K then
        Temp                    := Data_Re (Bit_Reversed_K);
        Data_Re(Bit_Reversed_K) := Data_Re(K);
        Data_Re(K)              := Temp;
        Temp                    := Data_Im (Bit_Reversed_K);
        Data_Im(Bit_Reversed_K) := Data_Im(K);
        Data_Im(K)              := Temp;
     end if;


  end Get_Bit_Reversal_Of;

  pragma Inline (Get_Bit_Reversal_Of);

  ---------
  -- FFT --
  ---------

  procedure FFT
    (Data_Re, Data_Im        : in out Data_Array;
     Transformed_Data_Last   :    out Data_Index;
     Input_Data_Last         : in     Data_Index;
     Exp_Table               : in out Exp_Storage;
     Inverse_FFT_Desired     : in     Boolean     := False;
     Normalized_Data_Desired : in     Boolean     := False;
     Bit_Reversal_Desired    : in     Boolean     := True)
  is
    Sqrt_Half   : constant Real := 0.70710_67811_86547_52440_08443_62104_84904;

    -- Indexing for the FFT

    Flock_Width, Start_Of_Butterfly, Start_Of_Flock : Data_Index;
    Butterflies_Per_Flock     : Data_Index;
    Butterflies_Per_Flock2    : Data_Index;
    Butterflies_Per_Flock3    : Data_Index;
    Butterflies_Per_Flock4    : Data_Index;
    Butterflies_Per_Flock5    : Data_Index;
    Butterflies_Per_Flock6    : Data_Index;
    Butterflies_Per_Flock7    : Data_Index;

    No_Of_Butterfly_Flocks    : Data_Index;
    Next_Butterfly            : Data_Index;
    Padded_Data_Index_Last    : Data_Index;
    Half_Data, Quarter_Data, Eighth_Data : Data_Index;
    Log_Of_Data_Length_Minus_1 : Exponent_Of_Two_Type;
    Radix_2_Stage_Last : Exponent_Of_Two_Type;
    Radix_8_Stage_Last : Exponent_Of_Two_Type;

    No_Of_Radix_8_Stages     : Data_Index;
    No_Of_Radix_4_Stages     : Data_Index;
    No_Of_Radix_2_Stages     : Data_Index;
    Remaining_Radix_2_Stages : Data_Index;
    Log_Of_Data_Length       : Data_Index;
    No_Of_Preliminary_Radix_8_Stages : Data_Index;

    Final_Stage_Is_A_Radix_8_Stage : Boolean;
    Final_Stage_Is_A_Radix_4_Stage : Boolean;
    Final_Stage_Is_A_Radix_2_Stage : Boolean;

    NormalizationFactor :  Real   := 1.0;

    -- For Radix 2 butterflies

    DataTop_Re, DataBot_Re : Real;
    DataTop_Im, DataBot_Im : Real;

    -- For Radix 4 butterflies

    Data02sum_Re, Data02sum_Im, Data02del_Re, Data02del_Im : Real;
    Data13sum_Re, Data13sum_Im, Data13del_Re, Data13del_Im : Real;
	
    -- For Radix 8 butterflies

    D0_Re, D1_Re, D2_Re, D3_Re, D4_Re, D5_Re, D6_Re, D7_Re : Real;
    D0_Im, D1_Im, D2_Im, D3_Im, D4_Im, D5_Im, D6_Im, D7_Im : Real;

    D04sum_Re, D15sum_Re, D26sum_Re, D37sum_Re : Real;
    D04sum_Im, D15sum_Im, D26sum_Im, D37sum_Im : Real;

    D04del_Re, D15del_Re, D26del_Re, D37del_Re : Real;
    D04del_Im, D15del_Im, D26del_Im, D37del_Im : Real;

    D1537sum_Re, D1537sum_Im, D1537del_Re, D1537del_Im : Real;
    Sqrt_sum_Re, Sqrt_sum_Im, Sqrt_del_Re, Sqrt_del_Im : Real;

    D04del26sum_Im, D04del26sum_Re, D04del26del_Re, D04del26del_Im : Real;

    D04sum26sum_Re, D04sum26sum_Im, D04sum26del_Re, D04sum26del_Im : Real;
    D15sum37sum_Re, D15sum37sum_Im, D15sum37del_Re, D15sum37del_Im : Real;

    Temp_Re, Temp_Im  : Real;
    Exp1_Re, Exp2_Re, Exp3_Re           : Real;
    Exp1_Im, Exp2_Im, Exp3_Im           : Real;
    Exp4_Re, Exp5_Re, Exp6_Re, Exp7_Re  : Real;
    Exp4_Im, Exp5_Im, Exp6_Im, Exp7_Im  : Real;
    Index0, Index1, Index2, Index3 : Data_Index;
    Index4, Index5, Index6, Index7 : Data_Index;
    Exp_Id, E_low : Data_Index;


    ------------------
    -- Integer_Log2 --
    ------------------

    -- Rounds down, so on range (eg) 0..2**15-1 it returns 0 .. 14
    -- returns 0 for I in 0..1
    -- returns 1 for I in 2..3, etc.
    function Integer_Log2 (I : Data_Index)
       return Exponent_Of_Two_Type
    is
       Log2 : Exponent_Of_Two_Type := 0;
    begin
       for Exponent in Exponent_Of_Two_Type loop
          exit when 2**Exponent > I;
          log2 := Exponent;
       end loop;
    return Log2;
    end Integer_Log2;

 begin



   -- Step 1. The FFT of a data set of length 1 is itself.
   -- So leave the data unchanged, and return.

   if Input_Data_Last < 1 then
      Padded_Data_Index_Last := 0;
      Transformed_Data_Last  := 0;
      return;
   end if;


   -- Step 2.
   -- Pad the data with zeros
   -- out to the nearest power of 2.  First must find what
   -- that power of 2 is, so  use Integer_Log2 (Input_Data_Last).
   -- Integer_Log2 rounds down so LogOf.. is in Exponent_Of_Two_Type range.
   -- This is LogOfDataSize - 1, not LogOf(DataSize-1). The length of the
   -- padded data is twice 2**Log_Of_Data_Length_Minus_1

   Log_Of_Data_Length_Minus_1 := Integer_Log2 (Input_Data_Last);

   Half_Data    := 2 ** Integer (Log_Of_Data_Length_Minus_1);
   Quarter_Data := Half_Data / 2; -- Could be 0
   Eighth_Data  := Half_Data / 4; -- Could be 0

   Padded_Data_Index_Last := Half_Data - 1;
   Padded_Data_Index_Last := Padded_Data_Index_Last + Half_Data;

   --  Pad with zeros:
   if Input_Data_Last < Data_Index'Last then
      for I in Input_Data_Last+1..Padded_Data_Index_Last loop
         Data_Re (I) := 0.0;
      end loop;
      for I in Input_Data_Last+1..Padded_Data_Index_Last loop
          Data_Im (I) := 0.0;
      end loop;
   end if;


   -- Step 2b. Must must decide how many Radix 8, 4, and 2 stages to take.
   -- Radix 2 stage will be the final stage if it is required. For a
   -- Radix 2 FFT, Stages is in 0..Log_Of_Data_Length_Minus_1.

   Log_Of_Data_Length := Data_Index (Log_Of_Data_Length_Minus_1) + 1;

   Remaining_Radix_2_Stages := Log_Of_Data_Length;
   No_Of_Radix_8_Stages     := Remaining_Radix_2_Stages / 3;
   Remaining_Radix_2_Stages := Remaining_Radix_2_Stages-No_Of_Radix_8_Stages*3;
   No_Of_Radix_4_Stages     := Remaining_Radix_2_Stages / 2;
   Remaining_Radix_2_Stages := Remaining_Radix_2_Stages-No_Of_Radix_4_Stages*2;
   No_Of_Radix_2_Stages     := Remaining_Radix_2_Stages;

   Final_Stage_Is_A_Radix_8_Stage := False;
   Final_Stage_Is_A_Radix_4_Stage := False;
   Final_Stage_Is_A_Radix_2_Stage := False;

   if No_Of_Radix_2_Stages > 0 then
     Final_Stage_Is_A_Radix_2_Stage := True;
   elsif No_Of_Radix_4_Stages > 0 then
     Final_Stage_Is_A_Radix_4_Stage := True;
   elsif No_Of_Radix_8_Stages > 0 then
     Final_Stage_Is_A_Radix_8_Stage := True;
   end if;

   No_Of_Preliminary_Radix_8_Stages := No_Of_Radix_8_Stages;
   if Final_Stage_Is_A_Radix_8_Stage then  -- have  Radix_8_Stages > 0
      No_Of_Preliminary_Radix_8_Stages := No_Of_Radix_8_Stages - 1;
   end if;


  -- Step 3. If haven't done so already, construct the table
  -- of SIN and COS functions (actually exp(-i*theta)).
  -- Dont need to remake it if it has already been constructed
  -- for the right value of Padded_Data_Index_Last.

   if not (Exp_Table.Current_Size_Of_Exp_Table = Padded_Data_Index_Last) then
      Make_Exp_Table (Padded_Data_Index_Last, Exp_Table);
   end if;

   -- here inline the Inverse (backward) fft by hand:

   if Inverse_FFT_Desired then


   -- Step 5a. Perform the FFT radix 8 Stages on 0..Radix_8_Stage_Last
   -- Each "Butterfly" has 8 inputs on the left, and 8 outputs on the right.
   -- The 8 on the left (input) are at indices n, n+N/8, n+N/4, n+3N/8 etc.
   -- Index n is in the range 0..N/8-1, where N is the Data_Length, a power of 2.
   -- The output data points (on the right) are at the same points.
   -- The above is for stage 0.  At state 1, just replace the N by N/8.
   --
   -- The Radix 8 stage is really a series of 8 point DFT's.  The Radix
   -- 2 stage is a series of 2 point DFT's.  We want the output of a single
   -- Radix 8 stage to be identical to the output of 3 Radix 2 stages.
   -- Therefore have to make sure that the 8 x 8 DFT's produces
   -- bit-reversed output just like results of 3 radix 2 stages, so
   -- have to artifically bit-reverse the results of the 8 x 8.  To bit
   -- reverse an 8 pt DFT, just swap items 1 and 4, and items 3 and 6.
   --
   -- Stage 0 prepares Data so that the N pt DFT can be done by eight N/8
   -- pt DFT's in stage 1.  Instead of performing the 8 DFT's, Stage 1 prepares
   -- the Data so the eight N/8 pt DFT's can each be
   -- done by eight N/64 pt. DFT's in stage 2. The value N/8 in Stage=0 is called
   -- Butterflies_Per_Flock.  Since each butterfly has 8 outputs, the
   -- N/8 butterflies of stage 0 comprise the entire data set.  At stage 0, with
   -- N/8 butterflies, say Flock_Width = N and No_Of_Butterfly_Flocks = 1.
   -- At Stage 1, with N/64 butterflies (for each of the 8 DFT's), say
   -- Flock_Width = N/8, and No_Of_Butterfly_Flocks = 8.

   if No_Of_Preliminary_Radix_8_Stages > 0 then

   Radix_8_Stage_Last:=Exponent_Of_Two_Type(No_Of_Preliminary_Radix_8_Stages-1);

   for Stage in 0 .. Radix_8_Stage_Last loop

     No_Of_Butterfly_Flocks := 2 ** Integer(3*Stage);
     Butterflies_Per_Flock  := Eighth_Data / No_Of_Butterfly_Flocks;
     Butterflies_Per_Flock2 := Butterflies_Per_Flock  + Butterflies_Per_Flock;
     Butterflies_Per_Flock3 := Butterflies_Per_Flock2 + Butterflies_Per_Flock;
     Butterflies_Per_Flock4 := Butterflies_Per_Flock2 + Butterflies_Per_Flock2;
     Butterflies_Per_Flock5 := Butterflies_Per_Flock2 + Butterflies_Per_Flock3;
     Butterflies_Per_Flock6 := Butterflies_Per_Flock2 + Butterflies_Per_Flock4;
     Butterflies_Per_Flock7 := Butterflies_Per_Flock2 + Butterflies_Per_Flock5;
     --  Butterflies_Per_Flock   is in the range   Eighth_Data..1
     --  No_Of_Butterfly_Flocks  is in the range   1..Eighth_Data

     Flock_Width := 0; -- Flock width overflows if Stage = 0.
     if Stage > 0 then Flock_Width := 8 * Butterflies_Per_Flock; end if;

     for Butterfly_ID in 0 .. Butterflies_Per_Flock-1 loop

	Next_Butterfly := No_Of_Butterfly_Flocks * Butterfly_ID;

        -- table holds fft exps, not inverse fft,
        -- so must change sign of imag. parts:

        Exp_Id := Next_Butterfly;
        Exp1_Re := Exp_Table.Re(Exp_Id); Exp1_Im := -Exp_Table.Im(Exp_Id);

        Exp_Id := Exp_Id + Next_Butterfly;
        Exp2_Re := Exp_Table.Re(Exp_Id); Exp2_Im := -Exp_Table.Im(Exp_Id);

        Exp_Id := Exp_Id + Next_Butterfly;
        Exp3_Re := Exp_Table.Re(Exp_Id); Exp3_Im := -Exp_Table.Im(Exp_Id);

        Exp_Id := Exp_Id + Next_Butterfly;
        Exp4_Re := Exp_Table.Re(Exp_Id); Exp4_Im := -Exp_Table.Im(Exp_Id);

        Exp_Id := Exp_Id + Next_Butterfly;
        if Exp_Id < Half_Data then
           Exp5_Re := Exp_Table.Re(Exp_Id); Exp5_Im := -Exp_Table.Im(Exp_Id);
        else
           E_low := Exp_Id - Half_Data;
           Exp5_Re := -Exp_Table.Re(E_low); Exp5_Im := Exp_Table.Im(E_low);
        end if;

        Exp_Id := Exp_Id + Next_Butterfly;
        if Exp_Id < Half_Data then
           Exp6_Re := Exp_Table.Re(Exp_Id); Exp6_Im := -Exp_Table.Im(Exp_Id);
        else
           E_low := Exp_Id - Half_Data;
           Exp6_Re := -Exp_Table.Re(E_low); Exp6_Im := Exp_Table.Im(E_low);
        end if;

        Exp_Id := Exp_Id + Next_Butterfly;
        if Exp_Id < Half_Data then
           Exp7_Re := Exp_Table.Re(Exp_Id); Exp7_Im := -Exp_Table.Im(Exp_Id);
        else
           E_low := Exp_Id - Half_Data;
           Exp7_Re := -Exp_Table.Re(E_low); Exp7_Im := Exp_Table.Im(E_low);
        end if;

        --  Next want to perform the 8 pt. Fourier Trans. on butterfly
        --  with ID = Butterfly_ID, in the Flock with ID = Flock_ID.  The
        --  span of a flock is 8 times the number of butterflies in the
        --  flock.  So the start of the Flock is at Flock_ID * Span,
        --  and add on to that Butterfly_ID to get the start of the
        --  desired butterfly.  The index of the start of the butterfly
        --  is called K:

        --  Here is the abstract version of the inner loop, commented out.
        --  Below it is unrolled by hand to make it more efficient under
        --  gcc/gnat.
        --
        -- must change sign of imag. partsfor inverse fft:
--
-- D0new := ((D0+D4) + (D2+D6) +   ( (D1+D5) + (D3+D7))) * Exp0;
-- D2new := ((D0+D4) - (D2+D6) - i*( (D1+D5) - (D3+D7))) * Exp2;
-- D4new := ((D0+D4) + (D2+D6) -   ( (D1+D5) + (D3+D7))) * Exp4;
-- D6new := ((D0+D4) - (D2+D6) + i*( (D1+D5) - (D3+D7))) * Exp6;
--
-- D1new := ((D0-D4) - i(D2-D6) + a( 1-i)*(D1-D5) + a(-1-i)*(D3-D7))*Exp1;
-- D3new := ((D0-D4) + i(D2-D6) + a(-1-i)*(D1-D5) + a( 1-i)*(D3-D7))*Exp3;
-- D5new := ((D0-D4) - i(D2-D6) + a(-1+i)*(D1-D5) + a( 1+i)*(D3-D7))*Exp5;
-- D7new := ((D0-D4) + i(D2-D6) + a( 1+i)*(D1-D5) + a(-1+i)*(D3-D7))*Exp7;
--
-- The last set can be written:
--
--D1new = ((D0-D4) - i(D2-D6) - a(i((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7))))*Exp1;
--D3new = ((D0-D4) + i(D2-D6) - a(i((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7))))*Exp3;
--D5new = ((D0-D4) - i(D2-D6) + a(i((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7))))*Exp5;
--D7new = ((D0-D4) + i(D2-D6) + a(i((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7))))*Exp7;
--

        for Flock_ID in 0 .. No_Of_Butterfly_Flocks-1 loop

           Start_Of_Flock     := Flock_ID * Flock_Width;
           Start_Of_Butterfly := Start_Of_Flock + Butterfly_ID;


           Index0 := Start_Of_Butterfly;
           Index4 := Start_Of_Butterfly + Butterflies_Per_Flock4;

           D0_Re     := Data_Re (Index0); D0_Im  := Data_Im (Index0);
           D4_Re     := Data_Re (Index4); D4_Im  := Data_Im (Index4);

           D04sum_Re := D0_Re + D4_Re; D04sum_Im := D0_Im + D4_Im;
           D04del_Re := D0_Re - D4_Re; D04del_Im := D0_Im - D4_Im;

           Index2 := Start_Of_Butterfly + Butterflies_Per_Flock2;
           Index6 := Start_Of_Butterfly + Butterflies_Per_Flock6;

           D2_Re     := Data_Re (Index2); D2_Im  := Data_Im (Index2);
           D6_Re     := Data_Re (Index6); D6_Im  := Data_Im (Index6);

           D26sum_Re := D2_Re + D6_Re; D26sum_Im := D2_Im + D6_Im;
           D26del_Re := D2_Re - D6_Re; D26del_Im := D2_Im - D6_Im;

           D04sum26sum_Re := D04sum_Re + D26sum_Re;
           D04sum26sum_Im := D04sum_Im + D26sum_Im;

           D04sum26del_Re := D04sum_Re - D26sum_Re;
           D04sum26del_Im := D04sum_Im - D26sum_Im;

           -- D04del26sum  := D04del + Im*D26del;
           D04del26sum_Re := D04del_Re - D26del_Im;
           D04del26sum_Im := D04del_Im + D26del_Re;

           -- D04del26del  := D04del - Im*D26del;
           D04del26del_Re := D04del_Re + D26del_Im;
           D04del26del_Im := D04del_Im - D26del_Re;



           Index1 := Start_Of_Butterfly + Butterflies_Per_Flock;
           Index5 := Start_Of_Butterfly + Butterflies_Per_Flock5;

           D1_Re     := Data_Re (Index1); D1_Im  := Data_Im (Index1);
           D5_Re     := Data_Re (Index5); D5_Im  := Data_Im (Index5);

           D15sum_Re := D1_Re + D5_Re; D15sum_Im := D1_Im + D5_Im;
           D15del_Re := D1_Re - D5_Re; D15del_Im := D1_Im - D5_Im;

           Index3 := Start_Of_Butterfly + Butterflies_Per_Flock3;
           Index7 := Start_Of_Butterfly + Butterflies_Per_Flock7;

           D3_Re     := Data_Re (Index3); D3_Im  := Data_Im (Index3);
           D7_Re     := Data_Re (Index7); D7_Im  := Data_Im (Index7);

           D37sum_Re := D3_Re + D7_Re; D37sum_Im := D3_Im + D7_Im;
           D37del_Re := D3_Re - D7_Re; D37del_Im := D3_Im - D7_Im;

           D15sum37sum_Re := D15sum_Re + D37sum_Re;
           D15sum37sum_Im := D15sum_Im + D37sum_Im;

           D15sum37del_Re := D15sum_Re - D37sum_Re;
           D15sum37del_Im := D15sum_Im - D37sum_Im;

           D1537sum_Re := D15del_Re + D37del_Re;
           D1537sum_Im := D15del_Im + D37del_Im;
           D1537del_Re := D15del_Re - D37del_Re;
           D1537del_Im := D15del_Im - D37del_Im;


           -- D0new := ((D0+D4) + (D2+D6) +   ( (D1+D5) + (D3+D7))) * Exp0;
           -- Data(Index0) := D04sum + D26sum + D15sum + D37sum;
           Temp_Re          := D04sum26sum_Re + D15sum37sum_Re;
           Temp_Im          := D04sum26sum_Im + D15sum37sum_Im;
           Data_Re (Index0) := Temp_Re;
           Data_Im (Index0) := Temp_Im;

           -- Inverse FFT: change sign of i:
           -- D2new       := ((D0+D4) - (D2+D6) + i*((D1+D5) - (D3+D7))) * Exp2;
           -- Data(Index2) := (D04sum - D26sum + Im * (D15sum - D37sum)) * Exp2;
           Temp_Re          := D04sum26del_Re - D15sum37del_Im;
           Temp_Im          := D04sum26del_Im + D15sum37del_Re;
           Data_Re (Index2) := Temp_Re*Exp2_Re - Temp_Im*Exp2_Im;
           Data_Im (Index2) := Temp_Re*Exp2_Im + Temp_Im*Exp2_Re;

           -- D4new     := ((D0+D4) + (D2+D6) - ((D1+D5) + (D3+D7))) * Exp4;
           -- Data(Index1) := (D04sum + D26sum -  (D15sum + D37sum)) * Exp4;
           -- To bit reverse the 8 x 8 DFT, swap items 4 and 1:
           Temp_Re          := D04sum26sum_Re - D15sum37sum_Re;
           Temp_Im          := D04sum26sum_Im - D15sum37sum_Im;
           Data_Re (Index1) := Temp_Re*Exp4_Re - Temp_Im*Exp4_Im;
           Data_Im (Index1) := Temp_Re*Exp4_Im + Temp_Im*Exp4_Re;

           -- Inverse FFT: change sign of i:
           -- D6new     := ((D0+D4) - (D2+D6) - i*((D1+D5) - (D3+D7))) * Exp6;
           -- Data(Index3) := (D04sum - D26sum - Im*(D15sum - D37sum)) * Exp6;
           -- To bit reverse the 8 x 8 DFT, swap items 6 and 3:
           Temp_Re          := D04sum26del_Re + D15sum37del_Im;
           Temp_Im          := D04sum26del_Im - D15sum37del_Re;
           Data_Re (Index3) := Temp_Re*Exp6_Re - Temp_Im*Exp6_Im;
           Data_Im (Index3) := Temp_Re*Exp6_Im + Temp_Im*Exp6_Re;



           -- Inverse FFT: change sign of i:
           -- Sqrt_Del := a(-i*((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7));
           Sqrt_Del_Re     := Sqrt_Half * ( D1537sum_Im - D1537del_Re);
           Sqrt_Del_Im     := Sqrt_Half * (-D1537sum_Re - D1537del_Im);

           -- Inverse FFT: change sign of i:
           -- Sqrt_Sum := a(-i*((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7));
           Sqrt_Sum_Re     := Sqrt_Half * ( D1537sum_Im + D1537del_Re);
           Sqrt_Sum_Im     := Sqrt_Half * (-D1537sum_Re + D1537del_Im);

           -- Inverse FFT: change sign of i:
           -- D04del26sum  := D04del + Im*D26del;
           -- Sqrt_Del     := a(-i*((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7));
           -- D1new        := ((D0-D4) + i(D2-D6)) - Sqrt_Del) * Exp1
           -- Data(Index4) := (D04del  + Im*D26del - Sqrt_Del) * Exp1;
           -- To bit reverse the 8 x 8 DFT, swap items 4 and 1:
           Temp_Re          := D04del26sum_Re - Sqrt_del_Re;
           Temp_Im          := D04del26sum_Im - Sqrt_del_Im;
           Data_Re (Index4) := Temp_Re*Exp1_Re - Temp_Im*Exp1_Im;
           Data_Im (Index4) := Temp_Re*Exp1_Im + Temp_Im*Exp1_Re;

           -- Inverse FFT: change sign of i:
           -- D04del26del  := D04del - Im*D26del;
           -- Sqrt_Sum     := a(-i*((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7));
           -- D3new        := ((D0-D4) - i(D2-D6)) - Sqrt_Sum) * Exp3;
           -- Data(Index6) := (D04del  - Im*D26del - Sqrt_Sum) * Exp3;
           -- To bit reverse the 8 x 8 DFT, swap items 6 and 3:
           Temp_Re          := D04del26del_Re - Sqrt_sum_Re;
           Temp_Im          := D04del26del_Im - Sqrt_sum_Im;
           Data_Re (Index6) := Temp_Re*Exp3_Re - Temp_Im*Exp3_Im;
           Data_Im (Index6) := Temp_Re*Exp3_Im + Temp_Im*Exp3_Re;

           -- D04del26sum  := D04del + Im*D26del;
           -- D5new        := ((D0-D4) + i(D2-D6)) + Sqrt_Del) * Exp5;
           -- Data(Index5) := (D04del  + Im*D26del + Sqrt_Del) * Exp5;
           Temp_Re          := D04del26sum_Re + Sqrt_del_Re;
           Temp_Im          := D04del26sum_Im + Sqrt_del_Im;
           Data_Re (Index5) := Temp_Re*Exp5_Re - Temp_Im*Exp5_Im;
           Data_Im (Index5) := Temp_Re*Exp5_Im + Temp_Im*Exp5_Re;

           -- D04del26del  := D04del - Im*D26del;
           -- D7new        := ((D0-D4) - i(D2-D6)) + Sqrt_Sum) * Exp7;
           -- Data(Index7) :=  (D04del - Im*D26del + Sqrt_Sum) * Exp7;
           Temp_Re          := D04del26del_Re + Sqrt_sum_Re;
           Temp_Im          := D04del26del_Im + Sqrt_sum_Im;
           Data_Re (Index7) := Temp_Re*Exp7_Re - Temp_Im*Exp7_Im;
           Data_Im (Index7) := Temp_Re*Exp7_Im + Temp_Im*Exp7_Re;

	end loop;
     end loop;
   end loop;
   end if;


   -- Step 5b. Perform the last stage of the Radix 8 set if it exists.

   if Final_Stage_Is_A_Radix_8_Stage then

   Radix_8_Stage_Last := Exponent_Of_Two_Type(Log_Of_Data_Length / 3 - 1);

   for Stage in Radix_8_Stage_Last .. Radix_8_Stage_Last loop

     No_Of_Butterfly_Flocks := 2 ** Integer (3*Stage);
     Butterflies_Per_Flock  := Eighth_Data / No_Of_Butterfly_Flocks; -- 1
     Butterflies_Per_Flock2 := Butterflies_Per_Flock  + Butterflies_Per_Flock;
     Butterflies_Per_Flock3 := Butterflies_Per_Flock2 + Butterflies_Per_Flock;
     Butterflies_Per_Flock4 := Butterflies_Per_Flock2 + Butterflies_Per_Flock2;
     Butterflies_Per_Flock5 := Butterflies_Per_Flock2 + Butterflies_Per_Flock3;
     Butterflies_Per_Flock6 := Butterflies_Per_Flock2 + Butterflies_Per_Flock4;
     Butterflies_Per_Flock7 := Butterflies_Per_Flock2 + Butterflies_Per_Flock5;
     --  Butterflies_Per_Flock   is in the range   Eighth_Data..1
     --  No_Of_Butterfly_Flocks  is in the range   1..Eighth_Data

     Flock_Width := 0; -- Flock width overflows if Stage = 0.
     if Stage > 0 then Flock_Width := 8 * Butterflies_Per_Flock; end if;

     for Butterfly_ID in 0 .. Butterflies_Per_Flock-1 loop

	Next_Butterfly := No_Of_Butterfly_Flocks * Butterfly_ID;

        for Flock_ID in 0 .. No_Of_Butterfly_Flocks-1 loop

           Start_Of_Flock     := Flock_ID * Flock_Width;
           Start_Of_Butterfly := Start_Of_Flock + Butterfly_ID;


           Index0 := Start_Of_Butterfly;
           Index4 := Start_Of_Butterfly + Butterflies_Per_Flock4;

           D0_Re     := Data_Re (Index0); D0_Im  := Data_Im (Index0);
           D4_Re     := Data_Re (Index4); D4_Im  := Data_Im (Index4);

           D04sum_Re := D0_Re + D4_Re; D04sum_Im := D0_Im + D4_Im;
           D04del_Re := D0_Re - D4_Re; D04del_Im := D0_Im - D4_Im;

           Index2 := Start_Of_Butterfly + Butterflies_Per_Flock2;
           Index6 := Start_Of_Butterfly + Butterflies_Per_Flock6;

           D2_Re     := Data_Re (Index2); D2_Im  := Data_Im (Index2);
           D6_Re     := Data_Re (Index6); D6_Im  := Data_Im (Index6);

           D26sum_Re := D2_Re + D6_Re; D26sum_Im := D2_Im + D6_Im;
           D26del_Re := D2_Re - D6_Re; D26del_Im := D2_Im - D6_Im;

           D04sum26sum_Re := D04sum_Re + D26sum_Re;
           D04sum26sum_Im := D04sum_Im + D26sum_Im;

           D04sum26del_Re := D04sum_Re - D26sum_Re;
           D04sum26del_Im := D04sum_Im - D26sum_Im;

           -- D04del26sum  := D04del + Im*D26del;
           D04del26sum_Re := D04del_Re - D26del_Im;
           D04del26sum_Im := D04del_Im + D26del_Re;

           -- D04del26del  := D04del - Im*D26del;
           D04del26del_Re := D04del_Re + D26del_Im;
           D04del26del_Im := D04del_Im - D26del_Re;



           Index1 := Start_Of_Butterfly + Butterflies_Per_Flock;
           Index5 := Start_Of_Butterfly + Butterflies_Per_Flock5;

           D1_Re     := Data_Re (Index1); D1_Im  := Data_Im (Index1);
           D5_Re     := Data_Re (Index5); D5_Im  := Data_Im (Index5);

           D15sum_Re := D1_Re + D5_Re; D15sum_Im := D1_Im + D5_Im;
           D15del_Re := D1_Re - D5_Re; D15del_Im := D1_Im - D5_Im;

           Index3 := Start_Of_Butterfly + Butterflies_Per_Flock3;
           Index7 := Start_Of_Butterfly + Butterflies_Per_Flock7;

           D3_Re     := Data_Re (Index3); D3_Im  := Data_Im (Index3);
           D7_Re     := Data_Re (Index7); D7_Im  := Data_Im (Index7);

           D37sum_Re := D3_Re + D7_Re; D37sum_Im := D3_Im + D7_Im;
           D37del_Re := D3_Re - D7_Re; D37del_Im := D3_Im - D7_Im;

           D15sum37sum_Re := D15sum_Re + D37sum_Re;
           D15sum37sum_Im := D15sum_Im + D37sum_Im;

           D15sum37del_Re := D15sum_Re - D37sum_Re;
           D15sum37del_Im := D15sum_Im - D37sum_Im;

           D1537sum_Re := D15del_Re + D37del_Re;
           D1537sum_Im := D15del_Im + D37del_Im;
           D1537del_Re := D15del_Re - D37del_Re;
           D1537del_Im := D15del_Im - D37del_Im;


           -- D0new := ((D0+D4) + (D2+D6) +   ( (D1+D5) + (D3+D7))) * Exp0;
           -- Data(Index0) := D04sum + D26sum + D15sum + D37sum;
           Temp_Re          := D04sum26sum_Re + D15sum37sum_Re;
           Temp_Im          := D04sum26sum_Im + D15sum37sum_Im;
           Data_Re (Index0) := Temp_Re;
           Data_Im (Index0) := Temp_Im;

           -- Inverse FFT: change sign of i:
           -- D2new       := ((D0+D4) - (D2+D6) + i*((D1+D5) - (D3+D7))) * Exp2;
           -- Data(Index2) := (D04sum - D26sum + Im * (D15sum - D37sum)) * Exp2;
           Temp_Re          := D04sum26del_Re - D15sum37del_Im;
           Temp_Im          := D04sum26del_Im + D15sum37del_Re;
           Data_Re (Index2) := Temp_Re;
           Data_Im (Index2) := Temp_Im;

           -- D4new     := ((D0+D4) + (D2+D6) - ((D1+D5) + (D3+D7))) * Exp4;
           -- Data(Index1) := (D04sum + D26sum -  (D15sum + D37sum)) * Exp4;
           -- To bit reverse the 8 x 8 DFT, swap items 4 and 1:
           Temp_Re          := D04sum26sum_Re - D15sum37sum_Re;
           Temp_Im          := D04sum26sum_Im - D15sum37sum_Im;
           Data_Re (Index1) := Temp_Re;
           Data_Im (Index1) := Temp_Im;

           -- Inverse FFT: change sign of i:
           -- D6new     := ((D0+D4) - (D2+D6) - i*((D1+D5) - (D3+D7))) * Exp6;
           -- Data(Index3) := (D04sum - D26sum - Im*(D15sum - D37sum)) * Exp6;
           -- To bit reverse the 8 x 8 DFT, swap items 6 and 3:
           Temp_Re          := D04sum26del_Re + D15sum37del_Im;
           Temp_Im          := D04sum26del_Im - D15sum37del_Re;
           Data_Re (Index3) := Temp_Re;
           Data_Im (Index3) := Temp_Im;



           -- Inverse FFT: change sign of i:
           -- Sqrt_Del := a(-i*((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7));
           Sqrt_Del_Re     := Sqrt_Half * ( D1537sum_Im - D1537del_Re);
           Sqrt_Del_Im     := Sqrt_Half * (-D1537sum_Re - D1537del_Im);

           -- Inverse FFT: change sign of i:
           -- Sqrt_Sum := a(-i*((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7));
           Sqrt_Sum_Re     := Sqrt_Half * ( D1537sum_Im + D1537del_Re);
           Sqrt_Sum_Im     := Sqrt_Half * (-D1537sum_Re + D1537del_Im);

           -- Inverse FFT: change sign of i:
           -- D04del26sum  := D04del + Im*D26del;
           -- Sqrt_Del     := a(-i*((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7));
           -- D1new        := ((D0-D4) + i(D2-D6)) - Sqrt_Del) * Exp1
           -- Data(Index4) := (D04del  + Im*D26del - Sqrt_Del) * Exp1;
           -- To bit reverse the 8 x 8 DFT, swap items 4 and 1:
           Temp_Re          := D04del26sum_Re - Sqrt_del_Re;
           Temp_Im          := D04del26sum_Im - Sqrt_del_Im;
           Data_Re (Index4) := Temp_Re;
           Data_Im (Index4) := Temp_Im;

           -- Inverse FFT: change sign of i:
           -- D04del26del  := D04del - Im*D26del;
           -- Sqrt_Sum     := a(-i*((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7));
           -- D3new        := ((D0-D4) - i(D2-D6)) - Sqrt_Sum) * Exp3;
           -- Data(Index6) := (D04del  - Im*D26del - Sqrt_Sum) * Exp3;
           -- To bit reverse the 8 x 8 DFT, swap items 6 and 3:
           Temp_Re          := D04del26del_Re - Sqrt_sum_Re;
           Temp_Im          := D04del26del_Im - Sqrt_sum_Im;
           Data_Re (Index6) := Temp_Re;
           Data_Im (Index6) := Temp_Im;

           -- D04del26sum  := D04del + Im*D26del;
           -- D5new        := ((D0-D4) + i(D2-D6)) + Sqrt_Del) * Exp5;
           -- Data(Index5) := (D04del  + Im*D26del + Sqrt_Del) * Exp5;
           Temp_Re          := D04del26sum_Re + Sqrt_del_Re;
           Temp_Im          := D04del26sum_Im + Sqrt_del_Im;
           Data_Re (Index5) := Temp_Re;
           Data_Im (Index5) := Temp_Im;

           -- D04del26del  := D04del - Im*D26del;
           -- D7new        := ((D0-D4) - i(D2-D6)) + Sqrt_Sum) * Exp7;
           -- Data(Index7) :=  (D04del - Im*D26del + Sqrt_Sum) * Exp7;
           Temp_Re          := D04del26del_Re + Sqrt_sum_Re;
           Temp_Im          := D04del26del_Im + Sqrt_sum_Im;
           Data_Re (Index7) := Temp_Re;
           Data_Im (Index7) := Temp_Im;

	end loop;
     end loop;
   end loop;
   end if;


   -- Step 5b. Perform an optimized final stage: a Radix 4 stage if it
   -- exists. All the Exps are 1.0, optimizations reflect that.

   if Final_Stage_Is_A_Radix_4_Stage then

 --Radix_4_Stage_Last := Exponent_Of_Two_Type(Log_Of_Data_Length / 2 - 1);

 --for Stage in 0 .. Radix_4_Stage_Last loop

   --No_Of_Butterfly_Flocks := 2 ** Integer (2*Stage);
   --Butterflies_Per_Flock  := Quarter_Data / No_Of_Butterfly_Flocks;


     No_Of_Butterfly_Flocks := Quarter_Data;
     Butterflies_Per_Flock  := 1;
     Butterflies_Per_Flock2 := 2*Butterflies_Per_Flock;
     Butterflies_Per_Flock3 := 3*Butterflies_Per_Flock;
     --  here: Butterflies_Per_Flock  := 1;
     --  here: No_Of_Butterfly_Flocks := Quarter_Data;

     Flock_Width := 0; -- Flock width overflows if Stage = 0.
     --if Stage > 0 then Flock_Width := 4 * Butterflies_Per_Flock; end if;
     Flock_Width := 4 * Butterflies_Per_Flock;

     for Butterfly_ID in 0 .. Butterflies_Per_Flock-1 loop

        --  All the Exp's equal One.


        for Flock_ID in 0 .. No_Of_Butterfly_Flocks-1 loop
	   Start_Of_Flock     := Flock_ID * Flock_Width;
	   Start_Of_Butterfly := Start_Of_Flock + Butterfly_ID;

	   Index0 := Start_Of_Butterfly;
	   Index1 := Start_Of_Butterfly + Butterflies_Per_Flock;
	   Index2 := Start_Of_Butterfly + Butterflies_Per_Flock2;
	   Index3 := Start_Of_Butterfly + Butterflies_Per_Flock3;
	
	   Data02Sum_Re := Data_Re(Index0) + Data_Re(Index2);
	   Data02Sum_Im := Data_Im(Index0) + Data_Im(Index2);
	   Data02Del_Re := Data_Re(Index0) - Data_Re(Index2);
	   Data02Del_Im := Data_Im(Index0) - Data_Im(Index2);
	   Data13Sum_Re := Data_Re(Index1) + Data_Re(Index3);
	   Data13Sum_Im := Data_Im(Index1) + Data_Im(Index3);
	   Data13Del_Re := Data_Re(Index1) - Data_Re(Index3);
	   Data13Del_Im := Data_Im(Index1) - Data_Im(Index3);

	   Data_Re (Index0) := Data02sum_Re + Data13sum_Re;
	   Data_Im (Index0) := Data02sum_Im + Data13sum_Im;
	   Data_Re (Index1) := Data02sum_Re - Data13sum_Re;
	   Data_Im (Index1) := Data02sum_Im - Data13sum_Im;
	   Data_Re (Index2) := Data02del_Re - Data13del_Im;
	   Data_Im (Index2) := Data02del_Im + Data13del_Re;
	   Data_Re (Index3) := Data02del_Re + Data13del_Im;
	   Data_Im (Index3) := Data02del_Im - Data13del_Re;
	end loop;
     end loop;
 --end loop;

   end if;


   -- Step 5c. Perform an optimized final stage, a radix 2 stage
   -- if it exists.  In this stage Stage := Log_Of_Size_Minus_1,
   -- all of the twiddle factors are 1 (Exp_Fctn = (1.0, 0.0)).

   if Final_Stage_Is_A_Radix_2_Stage then

   Radix_2_Stage_Last := Exponent_Of_Two_Type(Log_Of_Data_Length - 1);

   for Stage in Radix_2_Stage_Last .. Radix_2_Stage_Last loop

      No_Of_Butterfly_Flocks := 2 ** Stage;
      Butterflies_Per_Flock  := Half_Data / No_Of_Butterfly_Flocks;
      --  Butterflies_Per_Flock = 1, and No_Of_Butterfly_Flocks = Half_Data

      Flock_Width := 0;
      if Stage > 0 then Flock_Width := 2*Butterflies_Per_Flock; end if;

      for Butterfly_ID in 0 .. Butterflies_Per_Flock-1 loop  -- in 0..0

	--  Next_Butterfly := No_Of_Butterfly_Flocks * Butterfly_ID;
        --  Exp_Id := Next_Butterfly;	
	--  Get_Exp_At (Exp_Id, Inverse_FFT_Desired, Exp1_Re, Exp1_Im);

	for Flock_ID in 0 .. No_Of_Butterfly_Flocks-1 loop
	
	   Start_Of_Flock     := Flock_ID * Flock_Width;
	   Start_Of_Butterfly := Start_Of_Flock + Butterfly_ID;
	   Index0 := Start_Of_Butterfly;
	   Index1 := Start_Of_Butterfly + Butterflies_Per_Flock;
	   DataTop_Re := Data_Re (Index0);
	   DataTop_Im := Data_Im (Index0);
	   DataBot_Re := Data_Re (Index1);
	   DataBot_Im := Data_Im (Index1);

	   --  Here is the 2 pt. FFT.  Exp = (1,0) here.
	   --  Data0 :=  DataTop + DataBot;
	   --  Data1 := (DataTop - DataBot)*(Exp1_Re, Exp1_Im);
	
	   Data_Re (Index0)  := DataTop_Re + DataBot_Re;
	   Data_Im (Index0)  := DataTop_Im + DataBot_Im;
	   Data_Re (Index1)  := DataTop_Re - DataBot_Re;
	   Data_Im (Index1)  := DataTop_Im - DataBot_Im;
	end loop;
     end loop;
    end loop;

   end if;

   else -- not Inverse FFT; here inline the forward fft by hand:


   -- Step 5a. Perform the FFT radix 8 Stages on 0..Radix_8_Stage_Last
   -- Each "Butterfly" has 8 inputs on the left, and 8 outputs on the right.
   -- The 8 on the left (input) are at indices n, n+N/8, n+N/4, n+3N/8 etc.
   -- Index n is in the range 0..N/8-1, where N is the Data_Length, a power
   -- of 2. The output data points (on the right) are at the same points.
   -- The above is for stage 0.  At state 1, just replace the N by N/8.
   --
   -- The Radix 8 stage is really a series of 8 point DFT's.  The Radix
   -- 2 stage is a series of 2 point DFT's.  Want the output of a single
   -- Radix 8 stage to be identical to the output of 3 Radix 2 stages.
   -- Therefore have to make sure that the 8 x 8 DFT's produces
   -- bit-reversed output just like results of 3 radix 2 stages, so
   -- have to artifically bit-reverse the results of the 8 x 8.  To bit
   -- reverse an 8 pt DFT, just swap items 1 and 4, and items 3 and 6.
   --
   -- Stage 0 prepares Data so that the N pt DFT can be done by eight N/8
   -- pt DFT's in stage 1.  Instead of performing the 8 DFT's, Stage 1
   -- prepares the Data so the eight N/8 pt DFT's can each be
   -- done by eight N/64 pt. DFT's in stage 2. The value N/8 in Stage=0 is
   -- called Butterflies_Per_Flock.  Since each butterfly has 8 outputs, the
   -- N/8 butterflies of stage 0 comprise the entire data set.  At stage 0,
   -- with N/8 butterflies, say Flock_Width = N and No_Of_Butterfly_Flocks = 1.
   -- At Stage 1, with N/64 butterflies (for each of the 8 DFT's), say
   -- Flock_Width = N/8, and No_Of_Butterfly_Flocks = 8.

   if No_Of_Preliminary_Radix_8_Stages > 0 then

   Radix_8_Stage_Last :=Exponent_Of_Two_Type(No_Of_Preliminary_Radix_8_Stages-1);

   for Stage in 0 .. Radix_8_Stage_Last loop

     No_Of_Butterfly_Flocks := 2 ** Integer (3*Stage);
     Butterflies_Per_Flock  := Eighth_Data / No_Of_Butterfly_Flocks;
     Butterflies_Per_Flock2 := Butterflies_Per_Flock  + Butterflies_Per_Flock;
     Butterflies_Per_Flock3 := Butterflies_Per_Flock2 + Butterflies_Per_Flock;
     Butterflies_Per_Flock4 := Butterflies_Per_Flock2 + Butterflies_Per_Flock2;
     Butterflies_Per_Flock5 := Butterflies_Per_Flock2 + Butterflies_Per_Flock3;
     Butterflies_Per_Flock6 := Butterflies_Per_Flock2 + Butterflies_Per_Flock4;
     Butterflies_Per_Flock7 := Butterflies_Per_Flock2 + Butterflies_Per_Flock5;
     --  Butterflies_Per_Flock   is in the range   Eighth_Data..1
     --  No_Of_Butterfly_Flocks  is in the range   1..Eighth_Data

     Flock_Width := 0; -- Flock width overflows if Stage = 0.
     if Stage > 0 then Flock_Width := 8 * Butterflies_Per_Flock; end if;

     for Butterfly_ID in 0 .. Butterflies_Per_Flock-1 loop

	Next_Butterfly := No_Of_Butterfly_Flocks * Butterfly_ID;

        Exp_Id := Next_Butterfly;
        Exp1_Re := Exp_Table.Re(Exp_Id); Exp1_Im := Exp_Table.Im(Exp_Id);

        Exp_Id := Exp_Id + Next_Butterfly;
        Exp2_Re := Exp_Table.Re(Exp_Id); Exp2_Im := Exp_Table.Im(Exp_Id);

        Exp_Id := Exp_Id + Next_Butterfly;
        Exp3_Re := Exp_Table.Re(Exp_Id); Exp3_Im := Exp_Table.Im(Exp_Id);

        Exp_Id := Exp_Id + Next_Butterfly;
        Exp4_Re := Exp_Table.Re(Exp_Id); Exp4_Im := Exp_Table.Im(Exp_Id);

        Exp_Id := Exp_Id + Next_Butterfly;
        if Exp_Id < Half_Data then
           Exp5_Re := Exp_Table.Re(Exp_Id); Exp5_Im := Exp_Table.Im(Exp_Id);
        else
           E_low := Exp_Id - Half_Data;
           Exp5_Re := -Exp_Table.Re(E_low); Exp5_Im := -Exp_Table.Im(E_low);
        end if;

        Exp_Id := Exp_Id + Next_Butterfly;
        if Exp_Id < Half_Data then
           Exp6_Re := Exp_Table.Re(Exp_Id); Exp6_Im := Exp_Table.Im(Exp_Id);
        else
           E_low := Exp_Id - Half_Data;
           Exp6_Re := -Exp_Table.Re(E_low); Exp6_Im := -Exp_Table.Im(E_low);
        end if;

        Exp_Id := Exp_Id + Next_Butterfly;
        if Exp_Id < Half_Data then
           Exp7_Re := Exp_Table.Re(Exp_Id); Exp7_Im := Exp_Table.Im(Exp_Id);
        else
           E_low := Exp_Id - Half_Data;
           Exp7_Re := -Exp_Table.Re(E_low); Exp7_Im := -Exp_Table.Im(E_low);
        end if;

        --  Next perform the 8 pt. Fourier Trans. on butterfly
        --  with ID = Butterfly_ID, in the Flock with ID = Flock_ID.  The
        --  span of a flock is 8 times the number of butterflies in the
        --  flock.  So the start of the Flock is at Flock_ID * Span,
        --  and add on to that Butterfly_ID to get the start of the
        --  desired butterfly.  The index of the start of the butterfly
        --  is called K:

        --  Here is the abstract version of the inner loop, commented out.
        --  Below it is unrolled by hand to make it more efficient under
        --  gcc/gnat.
        --
--
-- D0new := ((D0+D4) + (D2+D6) +   ( (D1+D5) + (D3+D7))) * Exp0;
-- D2new := ((D0+D4) - (D2+D6) - i*( (D1+D5) - (D3+D7))) * Exp2;
-- D4new := ((D0+D4) + (D2+D6) -   ( (D1+D5) + (D3+D7))) * Exp4;
-- D6new := ((D0+D4) - (D2+D6) + i*( (D1+D5) - (D3+D7))) * Exp6;
--
-- D1new := ((D0-D4) - i(D2-D6) + a( 1-i)*(D1-D5) + a(-1-i)*(D3-D7))*Exp1;
-- D3new := ((D0-D4) + i(D2-D6) + a(-1-i)*(D1-D5) + a( 1-i)*(D3-D7))*Exp3;
-- D5new := ((D0-D4) - i(D2-D6) + a(-1+i)*(D1-D5) + a( 1+i)*(D3-D7))*Exp5;
-- D7new := ((D0-D4) + i(D2-D6) + a( 1+i)*(D1-D5) + a(-1+i)*(D3-D7))*Exp7;
--
-- The last set can be written:
--
--D1new = ((D0-D4) - i(D2-D6) - a(i((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7))))*Exp1;
--D3new = ((D0-D4) + i(D2-D6) - a(i((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7))))*Exp3;
--D5new = ((D0-D4) - i(D2-D6) + a(i((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7))))*Exp5;
--D7new = ((D0-D4) + i(D2-D6) + a(i((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7))))*Exp7;
--


        for Flock_ID in 0 .. No_Of_Butterfly_Flocks-1 loop


           Start_Of_Flock     := Flock_ID * Flock_Width;
           Start_Of_Butterfly := Start_Of_Flock + Butterfly_ID;


           Index0 := Start_Of_Butterfly;
           Index4 := Start_Of_Butterfly + Butterflies_Per_Flock4;

           D0_Re     := Data_Re (Index0); D0_Im  := Data_Im (Index0);
           D4_Re     := Data_Re (Index4); D4_Im  := Data_Im (Index4);

           D04sum_Re := D0_Re + D4_Re; D04sum_Im := D0_Im + D4_Im;
           D04del_Re := D0_Re - D4_Re; D04del_Im := D0_Im - D4_Im;

           Index2 := Start_Of_Butterfly + Butterflies_Per_Flock2;
           Index6 := Start_Of_Butterfly + Butterflies_Per_Flock6;

           D2_Re     := Data_Re (Index2); D2_Im  := Data_Im (Index2);
           D6_Re     := Data_Re (Index6); D6_Im  := Data_Im (Index6);

           D26sum_Re := D2_Re + D6_Re; D26sum_Im := D2_Im + D6_Im;
           D26del_Re := D2_Re - D6_Re; D26del_Im := D2_Im - D6_Im;

           D04sum26sum_Re := D04sum_Re + D26sum_Re;
           D04sum26sum_Im := D04sum_Im + D26sum_Im;

           D04sum26del_Re := D04sum_Re - D26sum_Re;
           D04sum26del_Im := D04sum_Im - D26sum_Im;

           -- D04del26sum  := D04del + Im*D26del;
           D04del26sum_Re := D04del_Re - D26del_Im;
           D04del26sum_Im := D04del_Im + D26del_Re;

           -- D04del26del  := D04del - Im*D26del;
           D04del26del_Re := D04del_Re + D26del_Im;
           D04del26del_Im := D04del_Im - D26del_Re;



           Index1 := Start_Of_Butterfly + Butterflies_Per_Flock;
           Index5 := Start_Of_Butterfly + Butterflies_Per_Flock5;

           D1_Re     := Data_Re (Index1); D1_Im  := Data_Im (Index1);
           D5_Re     := Data_Re (Index5); D5_Im  := Data_Im (Index5);

           D15sum_Re := D1_Re + D5_Re; D15sum_Im := D1_Im + D5_Im;
           D15del_Re := D1_Re - D5_Re; D15del_Im := D1_Im - D5_Im;

           Index3 := Start_Of_Butterfly + Butterflies_Per_Flock3;
           Index7 := Start_Of_Butterfly + Butterflies_Per_Flock7;

           D3_Re     := Data_Re (Index3); D3_Im  := Data_Im (Index3);
           D7_Re     := Data_Re (Index7); D7_Im  := Data_Im (Index7);

           D37sum_Re := D3_Re + D7_Re; D37sum_Im := D3_Im + D7_Im;
           D37del_Re := D3_Re - D7_Re; D37del_Im := D3_Im - D7_Im;

           D15sum37sum_Re := D15sum_Re + D37sum_Re;
           D15sum37sum_Im := D15sum_Im + D37sum_Im;

           D15sum37del_Re := D15sum_Re - D37sum_Re;
           D15sum37del_Im := D15sum_Im - D37sum_Im;

           D1537sum_Re := D15del_Re + D37del_Re;
           D1537sum_Im := D15del_Im + D37del_Im;
           D1537del_Re := D15del_Re - D37del_Re;
           D1537del_Im := D15del_Im - D37del_Im;



           -- D0new := ((D0+D4) + (D2+D6) + ((D1+D5) + (D3+D7))) * Exp0;
           -- Data(Index0) := D04sum + D26sum + D15sum + D37sum;
           Temp_Re          := D04sum26sum_Re + D15sum37sum_Re;
           Temp_Im          := D04sum26sum_Im + D15sum37sum_Im;
           Data_Re (Index0) := Temp_Re;
           Data_Im (Index0) := Temp_Im;

           -- D2new       := ((D0+D4) - (D2+D6) - i*((D1+D5) - (D3+D7))) * Exp2;
           -- Data(Index2) := (D04sum - D26sum - Im * (D15sum - D37sum)) * Exp2;
           Temp_Re          := D04sum26del_Re + D15sum37del_Im;
           Temp_Im          := D04sum26del_Im - D15sum37del_Re;
           Data_Re (Index2) := Temp_Re*Exp2_Re - Temp_Im*Exp2_Im;
           Data_Im (Index2) := Temp_Re*Exp2_Im + Temp_Im*Exp2_Re;

           -- D4new     := ((D0+D4) + (D2+D6) - ((D1+D5) + (D3+D7))) * Exp4;
           -- Data(Index1) := (D04sum + D26sum -  (D15sum + D37sum)) * Exp4;
           -- To bit reverse the 8 x 8 DFT, swap items 4 and 1:
           Temp_Re          := D04sum26sum_Re - D15sum37sum_Re;
           Temp_Im          := D04sum26sum_Im - D15sum37sum_Im;
           Data_Re (Index1) := Temp_Re*Exp4_Re - Temp_Im*Exp4_Im;
           Data_Im (Index1) := Temp_Re*Exp4_Im + Temp_Im*Exp4_Re;

           -- D6new     := ((D0+D4) - (D2+D6) + i*((D1+D5) - (D3+D7))) * Exp6;
           -- Data(Index3) := (D04sum - D26sum + Im*(D15sum - D37sum)) * Exp6;
           -- To bit reverse the 8 x 8 DFT, swap items 6 and 3:
           Temp_Re          := D04sum26del_Re - D15sum37del_Im;
           Temp_Im          := D04sum26del_Im + D15sum37del_Re;
           Data_Re (Index3) := Temp_Re*Exp6_Re - Temp_Im*Exp6_Im;
           Data_Im (Index3) := Temp_Re*Exp6_Im + Temp_Im*Exp6_Re;



           -- Sqrt_Del := a(i*((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7));
           Sqrt_Del_Re   := Sqrt_Half * (-D1537sum_Im - D1537del_Re);
           Sqrt_Del_Im   := Sqrt_Half * ( D1537sum_Re - D1537del_Im);

           -- Sqrt_Sum := a(i*((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7));
           Sqrt_Sum_Re   := Sqrt_Half * (-D1537sum_Im + D1537del_Re);
           Sqrt_Sum_Im   := Sqrt_Half * ( D1537sum_Re + D1537del_Im);

           -- D04del26del  := D04del - Im*D26del;
           -- Sqrt_Del     := a(i*((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7));
           -- D1new        := ((D0-D4) - i(D2-D6) - Sqrt_Del) * Exp1
           -- Data(Index4) :=  (D04del - Im*D26del - Sqrt_Del) * Exp1;
           -- To bit reverse the 8 x 8 DFT, swap items 4 and 1:
           Temp_Re           := D04del26del_Re - Sqrt_del_Re;
           Temp_Im          := D04del26del_Im - Sqrt_del_Im;
           Data_Re (Index4) := Temp_Re*Exp1_Re - Temp_Im*Exp1_Im;
           Data_Im (Index4) := Temp_Re*Exp1_Im + Temp_Im*Exp1_Re;

           -- D04del26sum  := D04del + Im*D26del;
           -- Sqrt_Sum     := a(i*((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7));
           -- D3new        := ((D0-D4) + i(D2-D6) - Sqrt_Sum) * Exp3;
           -- Data(Index6) :=  (D04del + Im*D26del - Sqrt_Sum) * Exp3;
           -- To bit reverse the 8 x 8 DFT, swap items 6 and 3:
           Temp_Re          := D04del26sum_Re - Sqrt_sum_Re;
           Temp_Im          := D04del26sum_Im - Sqrt_sum_Im;
           Data_Re (Index6) := Temp_Re*Exp3_Re - Temp_Im*Exp3_Im;
           Data_Im (Index6) := Temp_Re*Exp3_Im + Temp_Im*Exp3_Re;

           -- D04del26del  := D04del - Im*D26del;
           -- D5new        := ((D0-D4) - i(D2-D6) + Sqrt_Del) * Exp5;
           -- Data(Index5) :=  (D04del - Im*D26del + Sqrt_Del) * Exp5;
           Temp_Re          := D04del26del_Re + Sqrt_del_Re;
           Temp_Im          := D04del26del_Im + Sqrt_del_Im;
           Data_Re (Index5) := Temp_Re*Exp5_Re - Temp_Im*Exp5_Im;
           Data_Im (Index5) := Temp_Re*Exp5_Im + Temp_Im*Exp5_Re;

           -- D04del26sum  := D04del + Im*D26del;
           -- D7new        := ((D0-D4) + i(D2-D6) + Sqrt_Sum) * Exp7;
           -- Data(Index7) :=  (D04del + Im*D26del + Sqrt_Sum) * Exp7;
           Temp_Re          := D04del26sum_Re + Sqrt_sum_Re;
           Temp_Im          := D04del26sum_Im + Sqrt_sum_Im;
           Data_Re (Index7) := Temp_Re*Exp7_Re - Temp_Im*Exp7_Im;
           Data_Im (Index7) := Temp_Re*Exp7_Im + Temp_Im*Exp7_Re;

       end loop;
     end loop;
   end loop;
   end if;


  -- Step 5b. Perform the last stage of the Radix 8 set if it exists.

   if Final_Stage_Is_A_Radix_8_Stage then

   Radix_8_Stage_Last := Exponent_Of_Two_Type(Log_Of_Data_Length / 3 - 1);

   for Stage in Radix_8_Stage_Last .. Radix_8_Stage_Last loop

     No_Of_Butterfly_Flocks := 2 ** Integer (3*Stage);
     Butterflies_Per_Flock  := Eighth_Data / No_Of_Butterfly_Flocks;
     Butterflies_Per_Flock2 := Butterflies_Per_Flock  + Butterflies_Per_Flock;
     Butterflies_Per_Flock3 := Butterflies_Per_Flock2 + Butterflies_Per_Flock;
     Butterflies_Per_Flock4 := Butterflies_Per_Flock2 + Butterflies_Per_Flock2;
     Butterflies_Per_Flock5 := Butterflies_Per_Flock2 + Butterflies_Per_Flock3;
     Butterflies_Per_Flock6 := Butterflies_Per_Flock2 + Butterflies_Per_Flock4;
     Butterflies_Per_Flock7 := Butterflies_Per_Flock2 + Butterflies_Per_Flock5;
     --  Butterflies_Per_Flock   is in the range   Eighth_Data..1
     --  No_Of_Butterfly_Flocks  is in the range   1..Eighth_Data

     Flock_Width := 0; -- Flock width overflows if Stage = 0.
     if Stage > 0 then Flock_Width := 8 * Butterflies_Per_Flock; end if;

     for Butterfly_ID in 0 .. Butterflies_Per_Flock-1 loop

	Next_Butterfly := No_Of_Butterfly_Flocks * Butterfly_ID;

        for Flock_ID in 0 .. No_Of_Butterfly_Flocks-1 loop


           Start_Of_Flock     := Flock_ID * Flock_Width;
           Start_Of_Butterfly := Start_Of_Flock + Butterfly_ID;


           Index0 := Start_Of_Butterfly;
           Index4 := Start_Of_Butterfly + Butterflies_Per_Flock4;

           D0_Re     := Data_Re (Index0); D0_Im  := Data_Im (Index0);
           D4_Re     := Data_Re (Index4); D4_Im  := Data_Im (Index4);

           D04sum_Re := D0_Re + D4_Re; D04sum_Im := D0_Im + D4_Im;
           D04del_Re := D0_Re - D4_Re; D04del_Im := D0_Im - D4_Im;

           Index2 := Start_Of_Butterfly + Butterflies_Per_Flock2;
           Index6 := Start_Of_Butterfly + Butterflies_Per_Flock6;

           D2_Re     := Data_Re (Index2); D2_Im  := Data_Im (Index2);
           D6_Re     := Data_Re (Index6); D6_Im  := Data_Im (Index6);

           D26sum_Re := D2_Re + D6_Re; D26sum_Im := D2_Im + D6_Im;
           D26del_Re := D2_Re - D6_Re; D26del_Im := D2_Im - D6_Im;

           D04sum26sum_Re := D04sum_Re + D26sum_Re;
           D04sum26sum_Im := D04sum_Im + D26sum_Im;

           D04sum26del_Re := D04sum_Re - D26sum_Re;
           D04sum26del_Im := D04sum_Im - D26sum_Im;

           -- D04del26sum  := D04del + Im*D26del;
           D04del26sum_Re := D04del_Re - D26del_Im;
           D04del26sum_Im := D04del_Im + D26del_Re;

           -- D04del26del  := D04del - Im*D26del;
           D04del26del_Re := D04del_Re + D26del_Im;
           D04del26del_Im := D04del_Im - D26del_Re;



           Index1 := Start_Of_Butterfly + Butterflies_Per_Flock;
           Index5 := Start_Of_Butterfly + Butterflies_Per_Flock5;

           D1_Re     := Data_Re (Index1); D1_Im  := Data_Im (Index1);
           D5_Re     := Data_Re (Index5); D5_Im  := Data_Im (Index5);

           D15sum_Re := D1_Re + D5_Re; D15sum_Im := D1_Im + D5_Im;
           D15del_Re := D1_Re - D5_Re; D15del_Im := D1_Im - D5_Im;

           Index3 := Start_Of_Butterfly + Butterflies_Per_Flock3;
           Index7 := Start_Of_Butterfly + Butterflies_Per_Flock7;

           D3_Re     := Data_Re (Index3); D3_Im  := Data_Im (Index3);
           D7_Re     := Data_Re (Index7); D7_Im  := Data_Im (Index7);

           D37sum_Re := D3_Re + D7_Re; D37sum_Im := D3_Im + D7_Im;
           D37del_Re := D3_Re - D7_Re; D37del_Im := D3_Im - D7_Im;

           D15sum37sum_Re := D15sum_Re + D37sum_Re;
           D15sum37sum_Im := D15sum_Im + D37sum_Im;

           D15sum37del_Re := D15sum_Re - D37sum_Re;
           D15sum37del_Im := D15sum_Im - D37sum_Im;

           D1537sum_Re := D15del_Re + D37del_Re;
           D1537sum_Im := D15del_Im + D37del_Im;
           D1537del_Re := D15del_Re - D37del_Re;
           D1537del_Im := D15del_Im - D37del_Im;



           -- D0new := ((D0+D4) + (D2+D6) + ((D1+D5) + (D3+D7))) * Exp0;
           -- Data(Index0) := D04sum + D26sum + D15sum + D37sum;
           Temp_Re          := D04sum26sum_Re + D15sum37sum_Re;
           Temp_Im          := D04sum26sum_Im + D15sum37sum_Im;
           Data_Re (Index0) := Temp_Re;
           Data_Im (Index0) := Temp_Im;

           -- D2new       := ((D0+D4) - (D2+D6) - i*((D1+D5) - (D3+D7))) * Exp2;
           -- Data(Index2) := (D04sum - D26sum - Im * (D15sum - D37sum)) * Exp2;
           Temp_Re          := D04sum26del_Re + D15sum37del_Im;
           Temp_Im          := D04sum26del_Im - D15sum37del_Re;
           Data_Re (Index2) := Temp_Re;
           Data_Im (Index2) := Temp_Im;

           -- D4new     := ((D0+D4) + (D2+D6) - ((D1+D5) + (D3+D7))) * Exp4;
           -- Data(Index1) := (D04sum + D26sum -  (D15sum + D37sum)) * Exp4;
           -- To bit reverse the 8 x 8 DFT, swap items 4 and 1:
           Temp_Re          := D04sum26sum_Re - D15sum37sum_Re;
           Temp_Im          := D04sum26sum_Im - D15sum37sum_Im;
           Data_Re (Index1) := Temp_Re;
           Data_Im (Index1) := Temp_Im;

           -- D6new     := ((D0+D4) - (D2+D6) + i*((D1+D5) - (D3+D7))) * Exp6;
           -- Data(Index3) := (D04sum - D26sum + Im*(D15sum - D37sum)) * Exp6;
           -- To bit reverse the 8 x 8 DFT, swap items 6 and 3:
           Temp_Re          := D04sum26del_Re - D15sum37del_Im;
           Temp_Im          := D04sum26del_Im + D15sum37del_Re;
           Data_Re (Index3) := Temp_Re;
           Data_Im (Index3) := Temp_Im;



           -- Sqrt_Del := a(i*((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7));
           Sqrt_Del_Re   := Sqrt_Half * (-D1537sum_Im - D1537del_Re);
           Sqrt_Del_Im   := Sqrt_Half * ( D1537sum_Re - D1537del_Im);

           -- Sqrt_Sum := a(i*((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7));
           Sqrt_Sum_Re   := Sqrt_Half * (-D1537sum_Im + D1537del_Re);
           Sqrt_Sum_Im   := Sqrt_Half * ( D1537sum_Re + D1537del_Im);

           -- D04del26del  := D04del - Im*D26del;
           -- Sqrt_Del     := a(i*((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7));
           -- D1new        := ((D0-D4) - i(D2-D6) - Sqrt_Del) * Exp1
           -- Data(Index4) :=  (D04del - Im*D26del - Sqrt_Del) * Exp1;
           -- To bit reverse the 8 x 8 DFT, swap items 4 and 1:
           Temp_Re           := D04del26del_Re - Sqrt_del_Re;
           Temp_Im          := D04del26del_Im - Sqrt_del_Im;
           Data_Re (Index4) := Temp_Re;
           Data_Im (Index4) := Temp_Im;

           -- D04del26sum  := D04del + Im*D26del;
           -- Sqrt_Sum     := a(i*((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7));
           -- D3new        := ((D0-D4) + i(D2-D6) - Sqrt_Sum) * Exp3;
           -- Data(Index6) :=  (D04del + Im*D26del - Sqrt_Sum) * Exp3;
           -- To bit reverse the 8 x 8 DFT, swap items 6 and 3:
           Temp_Re          := D04del26sum_Re - Sqrt_sum_Re;
           Temp_Im          := D04del26sum_Im - Sqrt_sum_Im;
           Data_Re (Index6) := Temp_Re;
           Data_Im (Index6) := Temp_Im;

           -- D04del26del  := D04del - Im*D26del;
           -- D5new        := ((D0-D4) - i(D2-D6) + Sqrt_Del) * Exp5;
           -- Data(Index5) :=  (D04del - Im*D26del + Sqrt_Del) * Exp5;
           Temp_Re          := D04del26del_Re + Sqrt_del_Re;
           Temp_Im          := D04del26del_Im + Sqrt_del_Im;
           Data_Re (Index5) := Temp_Re;
           Data_Im (Index5) := Temp_Im;

           -- D04del26sum  := D04del + Im*D26del;
           -- D7new        := ((D0-D4) + i(D2-D6) + Sqrt_Sum) * Exp7;
           -- Data(Index7) :=  (D04del + Im*D26del + Sqrt_Sum) * Exp7;
           Temp_Re          := D04del26sum_Re + Sqrt_sum_Re;
           Temp_Im          := D04del26sum_Im + Sqrt_sum_Im;
           Data_Re (Index7) := Temp_Re;
           Data_Im (Index7) := Temp_Im;

	end loop;
     end loop;
   end loop;
   end if;


   -- Step 5b. Perform an optimized final stage: a Radix 4 stage if it
   -- exists. All the Exps are 1.0, optimizations reflect that.

   if Final_Stage_Is_A_Radix_4_Stage then

 --Radix_4_Stage_Last := Exponent_Of_Two_Type(Log_Of_Data_Length / 2 - 1);

 --for Stage in 0 .. Radix_4_Stage_Last loop

   --No_Of_Butterfly_Flocks := 2 ** Integer (2*Stage);
   --Butterflies_Per_Flock  := Quarter_Data / No_Of_Butterfly_Flocks;


     No_Of_Butterfly_Flocks := Quarter_Data;
     Butterflies_Per_Flock  := 1;
     Butterflies_Per_Flock2 := 2*Butterflies_Per_Flock;
     Butterflies_Per_Flock3 := 3*Butterflies_Per_Flock;
     --  here: Butterflies_Per_Flock  := 1;
     --  here: No_Of_Butterfly_Flocks := Quarter_Data;

     Flock_Width := 0; -- Flock width overflows if Stage = 0.
     --if Stage > 0 then Flock_Width := 4 * Butterflies_Per_Flock; end if;
     Flock_Width := 4 * Butterflies_Per_Flock;

     for Butterfly_ID in 0 .. Butterflies_Per_Flock-1 loop

        --  All the Exp's equal One.

        for Flock_ID in 0 .. No_Of_Butterfly_Flocks-1 loop
	   Start_Of_Flock     := Flock_ID * Flock_Width;
	   Start_Of_Butterfly := Start_Of_Flock + Butterfly_ID;

	   Index0 := Start_Of_Butterfly;
	   Index1 := Start_Of_Butterfly + Butterflies_Per_Flock;
	   Index2 := Start_Of_Butterfly + Butterflies_Per_Flock2;
	   Index3 := Start_Of_Butterfly + Butterflies_Per_Flock3;
	
	   Data02Sum_Re := Data_Re(Index0) + Data_Re(Index2);
	   Data02Sum_Im := Data_Im(Index0) + Data_Im(Index2);
	   Data02Del_Re := Data_Re(Index0) - Data_Re(Index2);
	   Data02Del_Im := Data_Im(Index0) - Data_Im(Index2);
	   Data13Sum_Re := Data_Re(Index1) + Data_Re(Index3);
	   Data13Sum_Im := Data_Im(Index1) + Data_Im(Index3);
	   Data13Del_Re := Data_Re(Index1) - Data_Re(Index3);
	   Data13Del_Im := Data_Im(Index1) - Data_Im(Index3);

	   Data_Re (Index0) := Data02sum_Re + Data13sum_Re;
	   Data_Im (Index0) := Data02sum_Im + Data13sum_Im;
	   Data_Re (Index1) := Data02sum_Re - Data13sum_Re;
	   Data_Im (Index1) := Data02sum_Im - Data13sum_Im;
	   Data_Re (Index2) := Data02del_Re + Data13del_Im;
	   Data_Im (Index2) := Data02del_Im - Data13del_Re;
	   Data_Re (Index3) := Data02del_Re - Data13del_Im;
	   Data_Im (Index3) := Data02del_Im + Data13del_Re;
	end loop;
	
     end loop;
 --end loop;
   end if;


   -- Step 5c. Perform an optimized final stage, a radix 2 stage
   -- if it exists.  In this stage Stage := Log_Of_Size_Minus_1,
   -- all of the twiddle factors are 1 (Exp_Fctn = (1.0, 0.0)).

   if Final_Stage_Is_A_Radix_2_Stage then

   Radix_2_Stage_Last := Exponent_Of_Two_Type(Log_Of_Data_Length - 1);

   for Stage in Radix_2_Stage_Last .. Radix_2_Stage_Last loop

      No_Of_Butterfly_Flocks := 2 ** Stage;
      Butterflies_Per_Flock  := Half_Data / No_Of_Butterfly_Flocks;
      --  Butterflies_Per_Flock = 1, and No_Of_Butterfly_Flocks = Half_Data

      Flock_Width := 0;
      if Stage > 0 then Flock_Width := 2*Butterflies_Per_Flock; end if;

      for Butterfly_ID in 0 .. Butterflies_Per_Flock-1 loop  -- in 0..0

	--  Next_Butterfly := No_Of_Butterfly_Flocks * Butterfly_ID;
        --  Exp_Id := Next_Butterfly;	
	--  Get_Exp_At (Exp_Id, Inverse_FFT_Desired, Exp1_Re, Exp1_Im);

	for Flock_ID in 0 .. No_Of_Butterfly_Flocks-1 loop
	
	   Start_Of_Flock     := Flock_ID * Flock_Width;
	   Start_Of_Butterfly := Start_Of_Flock + Butterfly_ID;
	   Index0 := Start_Of_Butterfly;
	   Index1 := Start_Of_Butterfly + Butterflies_Per_Flock;
	   DataTop_Re := Data_Re (Index0);
	   DataTop_Im := Data_Im (Index0);
	   DataBot_Re := Data_Re (Index1);
	   DataBot_Im := Data_Im (Index1);

	   --  Here is the 2 pt. FFT.  Exp = (1,0) here.
	   --  Data0 :=  DataTop + DataBot;
	   --  Data1 := (DataTop - DataBot)*(Exp1_Re, Exp1_Im);
	
	   Data_Re (Index0)  := DataTop_Re + DataBot_Re;
	   Data_Im (Index0)  := DataTop_Im + DataBot_Im;
	   Data_Re (Index1)  := DataTop_Re - DataBot_Re;
	   Data_Im (Index1)  := DataTop_Im - DataBot_Im;
	end loop;
     end loop;
   end loop;
   end if;

   end if;  -- not inverse fft


   if Normalized_Data_Desired then
      NormalizationFactor := 1.0 / SQRT (Real(Padded_Data_Index_Last)+1.0);
      for K in 0..Padded_Data_Index_Last loop
	  Data_Re (K) := Data_Re (K) * NormalizationFactor;
      end loop;
      for K in 0..Padded_Data_Index_Last loop
	  Data_Im (K) := Data_Im (K) * NormalizationFactor;
      end loop;
   end if;


   if Bit_Reversal_Desired then
   Get_Bit_Reversal_Of
     (Data_Re                     => Data_Re,
      Data_Im                     => Data_Im,
      Two_To_The_Top_Bit          => Half_Data,
      Half_Two_To_The_Top_Bit     => Quarter_Data,
      Quarter_Two_To_The_Top_Bit  => Eighth_Data,
      Top_Bit                     => Log_Of_Data_Length_Minus_1);
   end if;

   Transformed_Data_Last := Padded_Data_Index_Last;

end FFT;

end Fourier8;
