
--  Stnd Binary Rank Test.

with text_io; use text_io;
with Binary_Rank;
with Disorderly.Random; use Disorderly.Random;
with Disorderly.Random.Clock_Entropy;
with Ada.Numerics.Discrete_Random;
with Chi_Gaussian_CDF;

procedure rank_tst_demo_1 is

   type Real is digits 15;

   Bits_per_Random_Word : constant := 32;
   pragma Assert (Bits_per_Random_Word > 31 and Bits_per_Random_Word < 62);
   --  32 bits for general tests. (Test parameters below were designed for > 30 bits.)

   package Rank_Matrix is new Binary_Rank 
     (No_of_Bits_per_Segment => Bits_per_Random_Word,
      Segments_per_Vector    => 1);                    -- can't change from 1

   use Rank_Matrix;

   r : Binary_Matrix;

   Max_Row, Max_Col : constant Matrix_Index := Matrix_Index'Last;
   -- Matrix_Index goes from 1 ..  Bits_per_Vector.
   -- all array start at 1, since that's the way Marsaglia's diehard does it.
   -- Notice that we can do sub blocks of the full matrix by changing above.

   -- Variables for chi-square test:

   package chi_cdf is new Chi_Gaussian_CDF (Real);

   Sample_Size : constant := 2**12; -- for quick demo, 2**12 is min. acceptable.
 --Sample_Size : constant := 2**27; -- 2**27 is good for strong test.
   --  Usually, the best way to make the Chi-Square as strong as possible is
   --  to make the Sample_Size as large as possible. 

   x : Random_Int;
   Stream_1 : State;
   Rank : Integer;
   Outcome_id : Matrix_Index;

   ave_chi, ave_p, ave_normalized_var_of_p : Real := 0.0;


   True_Degrees_of_Freedom : constant := 3;

   subtype Range_of_Possible_Outcomes is 
                Matrix_Index range Max_Col-True_Degrees_of_Freedom .. Max_Col;

   type Statistical_Data is array (Range_of_Possible_Outcomes) of Real;

   Observed_Count : Statistical_Data := (others => 0.0);
   Expected_Count : Statistical_Data := (others => 0.0);

   Expected_Frequency : constant Statistical_Data :=

   -- FOR 31 bits and greater:

   --  True_Degrees_of_Freedom : constant := 2;
  --(0.13363571467, 0.57757619017, 0.28878809515);

   --  True_Degrees_of_Freedom : constant := 3;
     (0.00528545025, 0.12835026442, 0.57757619017, 0.28878809515);

   --  True_Degrees_of_Freedom : constant := 4;
   --(4.66639518324E-05, 0.0052387863054, 0.1283502644829, 
   -- 0.577576190173205, 0.2887880950866);

  --  True_Degrees_of_Freedom : constant := 5;
   -- (9.6962450869E-08, 4.65669893816E-05, 0.0052387863054,
   --  0.12835026448293, 0.577576190173205, 0.2887880950866);


   -- from known analytical distribution (true asymptotically also):

   -- expected rank_freq
   --   rank    probability
   --   <=29    0.00528545025
   --     30    0.12835026442
   --     31    0.57757619017
   --     32    0.28878809515
   --
   --     58 4.66639518324322E-05
   --     59 5.23878630542589E-03
   --     60 1.28350264482934E-01
   --     61 5.77576190173205E-01
   --     62 2.88788095086602E-01

   --     57 9.69624508687517E-08
   --     58 4.65669893815635E-05
   --     59 5.23878630542589E-03
   --     60 1.28350264482934E-01
   --     61 5.77576190173205E-01
   --     62 2.88788095086602E-01
   
   -- 61 bits:
   --
   --      56 9.69624508687517E-08
   --      57 4.65669893815635E-05
   --      58 5.23878630542589E-03
   --      59 1.28350264482934E-01
   --      60 5.77576190173205E-01
   --      61 2.88788095086602E-01
   --
   --      57 4.66639518324322E-05   and all < 57
   --      58 5.23878630542589E-03
   --      59 1.28350264482934E-01
   --      60 5.77576190173205E-01
   --      61 2.88788095086602E-01
   
   --25 6.05558642701658E-15 
   --26 4.88352774683776E-11
   --27 9.69136088580580E-08
   --28 4.65669892297724E-05
   --29 5.23878629810739E-03
   --30 1.28350264423167E-01
   --31 5.77576190173205E-01
   --32 2.88788095153841E-01

  --------------------
  -- Get_Random_Stnd -
  --------------------

  --  Compiler's native random number generator.

  type Unsigned_Stnd is mod 2**Bits_per_Random_Word;

  package rnd is new Ada.Numerics.Discrete_Random (Unsigned_Stnd);

  Stream_Stnd : rnd.Generator;

  procedure Get_Random_Stnd(X : out Random_Int; S : in rnd.Generator) 
  is
  begin
    X := Random_Int (rnd.Random (S));
  end Get_Random_Stnd;

  pragma Inline (Get_Random_Stnd);

begin 

  for i in Range_of_Possible_Outcomes loop
     Expected_Count (i) := Expected_Frequency (i) * Real(Sample_Size);
  end loop;
  --  Init array of expected counts: same really for 31 bit to 64 bit nums.
  --  True_Degrees_of_Freedom : constant := 3; These 4 sum to 1.0

  rnd.Reset (Stream_Stnd);  
  -- compiler's PRNG; Initialize the stream even if you don't use it.

  Clock_Entropy.Reset (Stream_1); 
  --  Init Stream_1 for calls to   Get_Random(x, Stream_1);

  --  do 1 chi-sqr test per Chi_Test_id value:

  for Chi_Test_id in Long_Integer range 1..2**28 loop  -- forever really
 
    Observed_Count := (others => 0.0);
 
    for Draw_id in Long_Integer range 1 .. Sample_Size loop
   
      for Col_id in Matrix_Index range  Matrix_Index'First .. Max_Col loop
      for Seg_id in Segments loop
       --Get_Random_Stnd (X, Stream_Stnd);
         Get_Random(X, Stream_1);
         r(Col_id)(Seg_id) := Unsigned_Segment(X mod 2**Bits_per_Random_Word);
      end loop;
      end loop;
   
      Get_Rank(r, Max_Row, Max_Col, Rank);
      --text_io.put(integer'Image(Max_Col - Rank));
   
      if Matrix_Index(Rank) > Range_of_Possible_Outcomes'First then
        Outcome_id := Matrix_Index(Rank);
      else
        Outcome_id := Range_of_Possible_Outcomes'First;
      end if;
      Observed_Count(Outcome_id) := Observed_Count(Outcome_id) + 1.0;
   
    end loop;
 
    declare
      chi, e, s, normalized_variance_of_p, p_val : Real;
      Degrees_of_Freedom : constant Real := Real (True_Degrees_of_Freedom);
    begin
      chi := 0.0;
      for Outcome in Range_of_Possible_Outcomes loop
        e :=  Expected_Count(Outcome);
        s := (Observed_Count(Outcome) - e)**2 / e;
        chi := chi + s;
      end loop;
      p_val := chi_cdf.Chi_Squared_CDF (Degrees_of_Freedom, chi);
      normalized_variance_of_p := (p_val - 0.5)**2 / (0.25/3.0);
  
      ave_chi := ave_chi + chi;
      ave_p   := ave_p + p_val;
      ave_normalized_var_of_p := ave_normalized_var_of_p + normalized_variance_of_p;
    end;
 
    new_line(1); 
    Put ("p-val average  (should be 0.5, after 1000's of iterations):");
    put (Real'Image (ave_p / Real (Chi_Test_id)));
    new_line; 
    Put ("p-val variance (should be 1.0, after 1000's of iterations):");
    put (Real'Image (ave_normalized_var_of_p / (Real (Chi_Test_id)))); 
    new_line; 
 
  end loop;

end;


