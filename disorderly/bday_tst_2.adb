

-- Set up at present to test the Basic_Rand.
--
-- Warning: puts large array on stack.  You get segmentation fault if stack
-- space is insufficient. To get more stack space in bash shell, type:
--    ulimit -s 1024M
--    ulimit -s unlimited
-- at the prompt. Type ulimit -a for summary of limits.  
-- In c shell type
--      limit stacksize 1024M
--      limit stacksize unlimited
--
-- after 30738012 trials have 50% probability of 2 or more collisions (if
-- using 48 bit uniformly distributed rands).

with Chi_Gaussian_CDF; --  cumulative distribution functions for getting p_vals
with Disorderly.Basic_Rand; use Disorderly.Basic_Rand;
with Disorderly.Basic_Rand.Clock_Entropy;
with Text_IO; use Text_IO;
with Ada.Numerics.Discrete_Random;
with Sorted_Array;

-- Test bday_tst_1 generates N random numbers (usually N = 2**25 + 2**24),
-- and then counts how many of these random numbers occurred
-- more than once in the set of N. 
--
-- Uses package Sorted_Array to sort the N random numbers. The sorted
-- data is then inspected to find the number of rands that occurred
-- more than once in the sample (collisions). The number of collisions
-- should follow a poisson distribution, which is tested with a Chi-square  
-- goodness-of-fit test. This routine does not count the number of collisions 
-- in the spacings of the N random numbers (Marsaglia's birthday-spacings-test).
--
-- You have a choice of 2 random num generators: 
--
--     Disorderly.Random    or    Ada.Numerics.Discrete_Random. 
--

procedure bday_tst_2 is

  --  Choose your Random Number Generator:

  Use_Compilers_Stnd_Generator : constant Boolean := False;
  --  IF False then use    Disorderly.Basic_Rand
  --  IF True  then use    Ada.Numerics.Discrete_Random

  Bits_per_Random_Word : constant := 48; -- can't change

  type Real is digits 15;

  package Chi_CDF is new Chi_Gaussian_CDF (Real);

  -- prob of no collisions (of an event of probability 2**(-48) = e) in n trials,
  -- (ie prob of no repeats of any of the 2**48 = 1/e possibilities):
  --  p(0) =  1*(1-e)*(1-2e)* ... *(1-(n-1)e)
  --  
  -- prob of exactly 1 collision:
  --  p(1) = (n(n+1)/2) * 1*e*(1-e)*(1-2e)* ... *(1-(n-2)e)
  --       = p(0) * (n(n+1)/2) * e * (1+(n-1)e)
  --  
  -- Poissonian p(k) = lambda^k exp(-lambda) / k! 
  --  
  -- p(0) =  exp(-lambda)
  -- p(1) =  exp(-lambda)* lambda = p(0) * lambda
  --  
  --  say  e = 2**(-48)  and  n = 2**25
  --  
  -- so by the 1st Poissonian p(0) above: 
  --    lambda = -log (1*(1-e)*(1-2e)* ... *(1-(n-1)e)) = 2.00000001987
  -- by the 2nd Poissonian p(1) above: 
  --    lambda = (n(n+1)/2) * e * (1+(n-1)e) =  2.000000298
  --  
  -- Agree to 7 significant figures, so we'll assume poissonian statistics
  -- with lambda obtained as above.
  --  

  Min_Recorded_Cnt : constant := 2;   -- counts 0..2   goto 1st  bin
  Max_Recorded_Cnt : constant := 7;   -- counts 7..inf goto last bin

  N : constant Parent_Random_Int := 2**25 + 2**24; 
  --  When N = 2**25 + 2**24, collision counts of 0..2 go into 1st bin,
  --  and then counts of 7 and higher go into bin 6th bin.

  Max_Val_of_Random_Ints : constant Parent_Random_Int := 2**Bits_per_Random_Word-1;

  type Table_Index is mod 2**26;

  package Sorted_Table is
     new Sorted_Array 
       (Item                    => Parent_Random_Int, 
        Max_Size_of_Item        => Max_Val_of_Random_Ints, 
        Max_Allowed_No_of_Items => N, 
        Table_Index             => Table_Index); 

  subtype Range_of_Recorded_Outcomes is Table_Index  
                           range Min_Recorded_Cnt..Max_Recorded_Cnt;

  type Birthday_Count_Statistics is array (Range_of_Recorded_Outcomes) of Real;

  Prob_of_k_Collisions : constant Birthday_Count_Statistics :=  
    (
     1.73578060853667E-01, 1.68717879896370E-01, 1.89807618654555E-01,
     1.70826860183124E-01, 1.28120147682862E-01, 1.68949432729422E-01
    );

  True_Degrees_of_Freedom : constant := 5.0; --  Prob  sums to 1.0.


  -- The following is stronger, (N=2**26) but requires a larger array.
  --
  --N : constant := 2**26; 
  --  When N = 2**26, collision counts of 0..4 go into 1st bin,
  --  and then counts of 12 and higher go into bin 9th bin.

  --Min_Recorded_Cnt : constant := 4;  -- counts 0..4    goto 1st  bin
  --Max_Recorded_Cnt : constant := 12; -- counts 12..inf goto last bin
  --
  --Prob_of_k_Collisions : Birthday_Count_Statistics :=  
  --  (
  --   9.96323936620416E-02, 9.16036574975765E-02, 1.22138211816770E-01,
  --   1.39586529870595E-01, 1.39586531950597E-01, 1.24076919138310E-01,
  --   9.92615367897598E-02, 7.21902096500887E-02, 1.11924009624262E-01
  --  );
  --
  --True_Degrees_of_Freedom : constant := 8.0;
  --  Prob  sums to 1.0.


  Bday_Year_Length    : constant := N;
  Sample_Size         : constant := 1600;    -- can use 1000 for N=2**26 case above.
  No_of_Chi_Tests     : constant := 2**31-1; -- just keep doing chi tests

  Cnt : Table_Index;
  Observed_Count : Birthday_Count_Statistics := (others => 0.0);
  Expected_Count : Birthday_Count_Statistics := (others => 0.0);

  p_val : Real;

  Stream_1 : State;

  X : Random_Int;

  --------------------
  -- Get_Random_Stnd -
  --------------------

  --  Compiler's built in random number generator.
  --  GNAT compiler complains if you ask for more than 48 bits.

  type Unsigned_Stnd is mod 2**Bits_per_Random_Word;

  package rnd is new Ada.Numerics.Discrete_Random (Unsigned_Stnd);

  Stream_stnd : rnd.Generator;

  procedure Get_Random_Stnd(X : out Random_Int; S : in rnd.Generator) 
  is
  begin
    X := Random_Int (rnd.Random (S));
  end Get_Random_Stnd;

  pragma Inline (Get_Random_Stnd);

begin
  --  Initialize states of the random number generators.
  --  Both generators use the clock to choose an initial seed.
  --  Call Reset only once; then repeat Test many times using  Stream_1.

  rnd.Reset (Stream_stnd);

  Clock_Entropy.Reset (Stream_1);

  --  Use Probability to init Expection values:

  for k in Range_of_Recorded_Outcomes loop
     Expected_Count(k) := Real (Sample_Size) * Prob_of_k_Collisions(k); 
  end loop;

  -- get  p_val  for each repetition of the chi test

  new_line;  
  if Use_Compilers_Stnd_Generator then
     put ("Using    Ada.Numerics.Discrete_Random   to generate random numbers");
  else
     put ("Using    Disorderly.Random    to generate random numbers");
  end if;
  new_line;  
  put ("Doing a Chi-squared Goodness-of-fit test with sample size ="); 
  put (Integer'Image (Sample_Size));
  new_line; 
  put ("Usually takes several hours on a 64-bit PC."); 
  new_line; 

  for Chi_Test_id in 1 .. No_of_Chi_Tests loop

     Observed_Count := (others => 0.0);

     Fill_the_Bins:
     for Trial_id in 1 .. Sample_Size loop 

        Sorted_Table.Initialize_Table_for_Restart;

        for i in 1 .. Bday_Year_Length loop
          if Use_Compilers_Stnd_Generator then
             Get_Random_Stnd (X, Stream_stnd);
          else
             Get_Random (X, Stream_1);
          end if;
          Sorted_Table.Insert_and_Sort (X mod 2**Bits_per_Random_Word);
        end loop; -- in i

        Cnt := Sorted_Table.No_of_Collisions_Detected;

        put (Table_Index'Image (Cnt));

        if Cnt >= Max_Recorded_Cnt then
          Observed_Count(Max_Recorded_Cnt) := Observed_Count(Max_Recorded_Cnt) + 1.0;
        elsif Cnt <= Min_Recorded_Cnt then
          Observed_Count(Min_Recorded_Cnt) := Observed_Count(Min_Recorded_Cnt) + 1.0;
        else
          Observed_Count(Cnt) := Observed_Count(Cnt) + 1.0;
        end if;

        if not Sorted_Table.Array_Sort_Successful then
           put_line ("Failure in array sort. Should never happen.");
           return;
        end if;

     end loop Fill_the_Bins;  -- Trial_id in 1 .. Sample_Size

     Get_Chi_p_val:
     declare
       chi, e, s, p : Real;
     begin
       chi := 0.0;
       for i in Range_of_Recorded_Outcomes loop
         e := Expected_Count(i);
         s := (Observed_Count(i) - e)**2 / e;
         chi := chi + s;
       end loop;
       p := chi_cdf . Chi_Squared_CDF (True_Degrees_of_Freedom, chi);
       p_val := p;
     end Get_Chi_p_val;

     new_line; 
     put ("Chi-squared Goodness-of-fit test, number"); put (Integer'Image (Chi_Test_id));
     put (" with sample size "); put (Integer'Image (Sample_Size));
     new_line; 
     put ("p-val (should be uniformly distributed in [0, 1) )");
     new_line; 
     put (Real'Image (p_val));
     new_line; 

  end loop;  -- for Chi_Test_id in 1..No_of_Chi_Tests

end;

