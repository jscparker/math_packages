
-- Warning: puts large array on stack.  You get segmentation fault if stack
-- space is insufficient. To get more stack space in bash shell, type:
--    ulimit -s 1024M
--    ulimit -s unlimited
-- at the prompt. Type ulimit -a for summary of limits.  
-- In c shell type
--      limit stacksize 1024M
--      limit stacksize unlimited
--
-- Translated from Marsaglia, Tsang  diehard suite.
--
-- Tests:  Disorderly.Basic_Rand

with Chi_Gaussian_CDF; --  cumulative distribution functions for getting p_vals
with Disorderly.Basic_Rand; use Disorderly.Basic_Rand;
with Disorderly.Basic_Rand.Clock_Entropy;
with Ada.Numerics.Generic_Elementary_Functions;
with Text_IO; use Text_IO;
with Ada.Numerics.Discrete_Random; -- for testing compiler's Generator

procedure gorilla_tst_demo_2 is

  pragma Assert (Random_Int'Last >= 2**32-1);
  --  Random_Int should be uniform on 0..2**n-1 in order to pass, with n >= 32. 

  type Real is digits 15;

  package math is new  Ada.Numerics.Generic_Elementary_Functions (Real); 
  use math;
  package Normal_p_vals is new Chi_Gaussian_CDF (Real); 
  use Normal_p_vals;

  Stream_1 : State;

  Random_x : Random_Int;
  k_th_bit : Random_Int;
  
  Log_of_32 : constant  := 5;
  type Gorilla_Data_Index is mod 2**(2**Log_of_32);

  Bytes_Bit_Count : array(Gorilla_Data_Index range 0..255) of 
                                 Gorilla_Data_Index := (others => 0);
  Count_of_Ones : Gorilla_Data_Index := 0;

  subtype Gorilla_Test_Bit_Range is Gorilla_Data_Index range 0 .. 31;
  Two_to_the : array(Gorilla_Test_Bit_Range) of Gorilla_Data_Index;

  --  can't change from 26 and 2**26:

  No_of_Bits_per_Gorilla_Word : constant := 26;
  No_of_Gorilla_Words_in_Test : constant Gorilla_Data_Index := 2**26;

  Gorilla_Word : Gorilla_Data_Index;
  j, q : Gorilla_Data_Index;
  No_of_Missing_Gorilla_Words : Gorilla_Data_Index;

  Log_Size_of_t_array : constant := No_of_Bits_per_Gorilla_Word - Log_of_32; -- 21

  t : array(Gorilla_Data_Index range 0..2**Log_Size_of_t_array-1) of Gorilla_Data_Index;

  ave_p_val, ad_ks_p_val, variance, ave_variance : Real;
  --  Probably won't converge to theoretical vals to high accuracy, because
  --  the calc of AD-KS p-vals is very approximate.

  z, v : Real;
  u : array (Gorilla_Test_Bit_Range) of Real;

  ------------------
  -- Get_Random_g --
  ------------------

  type Unsigned_48 is mod 2**48;

  package Discrete_48_bit is new Ada.Numerics.Discrete_Random (Unsigned_48);

  Stream_gnat : Discrete_48_bit.Generator;

  procedure Get_Random_g(X : out Random_Int; S : in Discrete_48_bit.Generator) 
  is
  begin
    X := Random_Int (Discrete_48_bit.Random (S) / 2**16);
  end Get_Random_g;


  ---------------------------
  -- AD_KS_test_p_value_32 --
  ---------------------------

  -- converts AD (anderson-darling) KS statistic into p-value, n=32

  function AD_KS_test_p_value_32 (z : Real) return Real is
    y, Result : Real := 0.0;
  begin

    if z < 0.0 then
      raise Constraint_Error;
    elsif (z = 0.0) then
      Result := Exp (0.266);
    elsif (z < 0.5) then
      y := Exp (-1.27 * Log (z));
      Result := Exp (0.266-(0.70025-(0.009804-0.000213*y)*y)*y);
    elsif (z<1.0) then
      Result :=  (0.53383+(1.65285-(1.988-0.6634*z)*z)*z)*z-0.21862;
    elsif (z<2.0) then
      Result :=  0.99987-0.6616*Exp (-1.0896*z)-0.953*Exp (-2.005*z);
    else
      Result :=  1.0-0.52686*Exp (-1.05276*z)-0.68053*Exp (-1.62034*z);
    end if;

    return Result;

  end AD_KS_test_p_value_32;

begin 

  declare 
    x, error : Real; All_OK : Boolean := true; 
  begin
  for i in 1..160 loop
    x := -30.0 + real(i)*1.0;
    Test_Normal_cdf(x, error);
    if abs error > 256.0*Real'Epsilon then 
      put ("Failure in Normal distribution CDF."); put(real'image(error));
      All_OK := False;
    end if;
  end loop;
  if All_OK then
    put ("Calculation of Normal distribution's cumulative distribution function OK.");
  end if;
  end;
  new_line(2);

  put ("Wait a few minutes. The p-vals should be uniformly distributed in (0,1):");
  new_line(2);

  Clock_Entropy.Reset (Stream_1);
  --  Call Reset only once; then repeat Test many times using  Stream_1.

  -- Gorilla test for 2^26 bits, positions 0 to 31:

  -- Bytes_Bit_Count(i) gives you number of 1's in byte i:

  for i in Gorilla_Data_Index range 1 .. 255  loop
    q := i;
    for m in 0 .. 7 loop
      if q rem 2 > 0 then 
        Bytes_Bit_Count(i) := Bytes_Bit_Count(i) + 1;
      end if;
      q := q / 2;
    end loop;
  end loop;

  Two_to_the(0) := 1;
  for i in 1 .. Gorilla_Test_Bit_Range'Last loop
    Two_to_the(i) := 2*Two_to_the(i-1); 
  end loop; -- set k_th_bits

  ave_p_val    := 0.0;
  ave_variance := 0.0;

  -- get  p_val  for each repetition of the test and print average p_val.

  for Test_id in 1 .. 2_000_000_000 loop

    Do_Each_Bit: for bit_id in Gorilla_Test_Bit_Range loop

        k_th_bit := Random_Int (Two_to_the(bit_id));     -- k_th_bit = 2**bit_id
	
	t := (others => 0);

        -- Push 26 bits onto the Gorilla_Word as tho' it were a stack. 
	-- Take the  k_th bit of 26 consecutive Random_X's, and fill
	-- bits 0..25 of the initial Gorilla_Word with these 26 bits.

        Gorilla_Word := 0;
        for i in 1 .. No_of_Bits_per_Gorilla_Word  loop
            Gorilla_Word := Gorilla_Word / 2;
            Get_Random (Random_X, Stream_1);
            if (Random_X AND k_th_bit) > 0 then
              Gorilla_Word := Gorilla_Word + 2**(No_of_Bits_per_Gorilla_Word-1);
            end if;
        end loop;

        --  Loop gets 2^26 (No_of_Gorilla_Words_in_Test) overlapping Gorilla_Word's; 
	--  Sets a unique bit in t() to store evidence of its observation. To do this
	--  we could make t(Gorilla_Word) either true or false, and pack the bits;
	--  but here we just translate the original code: does the equivalent with 
	--  1's  and 0's.

        for Word_id in 1 .. No_of_Gorilla_Words_in_Test  loop

          -- push k_th bit of Random_X  onto top of bit_stack = Gorilla_Word.

          Gorilla_Word := Gorilla_Word / 2;
          Get_Random (Random_X, Stream_1);
          if (Random_X AND k_th_bit) > 0 then  
            Gorilla_Word := Gorilla_Word + 2**(No_of_Bits_per_Gorilla_Word-1);
          end if;

          j    := Gorilla_Word AND (2**Log_Size_of_t_array-1);  
          --  The lower 21 bits of Gorilla_Word is the index j of array t.

          t(j) := t(j) OR Two_to_the (Gorilla_Word / 2**Log_Size_of_t_array);
          --  The upper 5 bits of Gorilla_Word give a unique val in 0..31, 
	  --  which bit is set to 1 in t(j). Use OR operator. A 1 at that bit
          --  means there were 1 or more occurances of the number Gorilla_Word.
	  --  This btw is how t got its size: 
	  --    Log_Size_of_t_array = No_of_Bits_per_Gorilla_Word - Log_of_32
	  --    Log_Size_of_t_array = 26 - 5

        end loop; 

        -- count 1's in  t(i) 

        Count_of_Ones := 0;
        for i in t'Range loop
          j := t(i);
          Count_of_Ones  := Count_of_Ones  
              + Bytes_Bit_Count( j AND 255)
              + Bytes_Bit_Count((j / 2**8 ) AND 255)
              + Bytes_Bit_Count((j / 2**16) AND 255)
              + Bytes_Bit_Count((j / 2**24) AND 255);
        end loop;

	No_of_Missing_Gorilla_Words := No_of_Gorilla_Words_in_Test - Count_of_Ones;
        --  put(Random_Int'Image(No_of_Missing_Gorilla_Words));

        z := (Real (No_of_Missing_Gorilla_Words) - 24687971.0) / 4170.0; 
	--  converts gaussian to standard normal distribution; z-score.
	--  ASSUMES 26 bit Gorilla_Words, and No_of_Gorilla_Words_in_Test = 2**26.

        u(bit_id) := Normal_CDF (z); 
	--  z-score converted to p-value; area under normal curve

       new_Line;
       put ("Gorilla Test number "); put (Integer'Image(Test_id));put(",  ");
       put ("Bit number "); put (Gorilla_Test_Bit_Range'Image(bit_id));put(":");
       new_Line;
       put ("p-val ="); 
       put (Real'Image(u(bit_id)));

    end loop Do_Each_Bit; --end bit_id loop

    -- Now do KS test on 32 u's.  first, sort u's:

    for i in Gorilla_Test_Bit_Range'First .. Gorilla_Test_Bit_Range'Last-1  loop
      for j in Gorilla_Test_Bit_Range'First .. Gorilla_Test_Bit_Range'Last-1-i loop
         if  u(j+1) < u(j)  then
            v      := u(j);
            u(j)   := u(j+1);
            u(j+1) := v;
         end if;
       end loop;
    end loop;

    -- calculate Anderson-Darling statistic:

    z := -1024.0;
    for i in  Gorilla_Test_Bit_Range  loop
      v := u(i) * (1.0 - u(Gorilla_Test_Bit_Range'Last - i));
      if (v < 1.0e-30) then   v := 1.0e-30;   end if;
      z := z - Real (i+i+1) * Log (v);
    end loop;

    ad_ks_p_val  := AD_KS_test_p_value_32 (z / 32.0);

    ave_p_val    := ave_p_val + ad_ks_p_val;
    variance     := (ad_ks_p_val - 0.5)**2;
    ave_variance := ave_variance + variance;

    new_Line(2);
    put ("Gorilla ks test number "); put (Integer'Image(Test_id));put(".");
    new_Line;
    put ("ks p-val average  (should approach 0.5 after 100's of ks tests): "); 
    put (Real'Image(ave_p_val / Real (Test_id)));
    new_Line;
    put ("ks p-val variance (should approach 1.0 after 100's of ks tests): "); 
    put (Real'Image(ave_variance / (Real (Test_id)*(0.25/3.0)))); -- should -> 1.0
    new_Line;

  end loop;  -- for Test_id

end;

