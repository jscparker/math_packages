
with Disorderly.Basic_Rand; use Disorderly.Basic_Rand;
with Disorderly.Basic_Rand.Clock_Entropy;
with Ada.Numerics.Discrete_Random;
with Chi_Gaussian_CDF;
with text_io; use text_io;

--  translated from Marsaglia, Tsang  diehard suite.

procedure gcd_6bytes_2 is

  Bits_per_Random_Word : constant := 48;
  --  Must set this correctly here.  There's no way to check this.

  Stream_1 : Disorderly.Basic_Rand.State;
  --  Create a stream of Random numbers.
  --  Initializing this after the begin, w/ a call to Reset_with_Calendar

  subtype gnat_Random_Int is Random_Int range 0 .. 2**Bits_per_Random_Word-1;
  package rnd is new ada.numerics.Discrete_Random (gnat_Random_Int);
  g : rnd.generator;
  --  Call rnd.Reset (g); to initialize.

  type Real is digits 15;

  package Chi_Analysis is new Chi_Gaussian_CDF (Real);  use Chi_Analysis;

  type Unsigned_64 is mod 2**64;

  type Statistical_Data is array (Unsigned_64 range <>) of Real;


  -- Greatest Common Divisor Count test.
  --
  -- 1st test counts No of occurances of GCD's calculated for pairs of Rands:

  Span_of_GCD_Count_Test : constant := 100;

  subtype GCD_Count_Test_Range is Unsigned_64 range 1 .. Span_of_GCD_Count_Test;

  subtype GCD_Counts is Statistical_Data (GCD_Count_Test_Range);

  True_DOF_for_GCD_Count_Test : constant := Span_of_GCD_Count_Test - 1;


  -- Greatest Common Divisor Iterations test.
  --
  -- 2nd test counts No of Iterations required to find GCD of a pair of Rands:

  subtype GCD_Iterations_Test_Range is Unsigned_64 range 5..53;

  subtype GCD_Iterations_Stats is Statistical_Data (GCD_Iterations_Test_Range);

  True_DOF_for_Iterations_Test : constant := 48;

  Probability_of_GCD_Iterations : constant GCD_Iterations_Stats := 

  (
   1.82133747963399E-09,
   9.58898969998596E-09, 4.95900373562431E-08, 2.27481731154183E-07,
   9.28626032480111E-07, 3.40040609676470E-06, 1.12699881420596E-05,
   3.40712181596310E-05, 9.44518676549460E-05, 2.41209404024540E-04,
   5.69742894566334E-04, 1.24896390954733E-03, 2.54813052386917E-03,
   4.85060876606904E-03, 8.63352069814511E-03, 1.43949377324915E-02,
   2.25197150861574E-02, 3.31014969046822E-02, 4.57718000608556E-02,
   5.96024363457740E-02, 7.31521241333012E-02, 8.46879485357183E-02,
   9.25374818484444E-02, 9.54841741832165E-02, 9.30779824121723E-02,
   8.57414430641431E-02, 7.46531672639515E-02, 6.14429557739744E-02,
   4.78042557873558E-02, 3.51559583569137E-02, 2.44342786072529E-02,
   1.60456052696471E-02, 9.95230933086590E-03, 5.82790566428320E-03,
   3.22011217073676E-03, 1.67784031979238E-03, 8.23785205774590E-04,
   3.80796274659059E-04, 1.65557171340292E-04, 6.76216404052808E-05,
   2.59157200839291E-05, 9.30882508186161E-06, 3.12834803067835E-06,
   9.80931601705069E-07, 2.86878314503330E-07, 7.78360272500109E-08,
   1.97336355105359E-08, 4.56083856575878E-09, 1.23887369078947E-09
  );

   --   distr based on 141 * 2^34 sample size

  ---------------------------------
  -- Get_Chi_Statistic_and_P_val --
  ---------------------------------

  procedure Get_Chi_Statistic_and_P_val 
    (Probability_Distribution : in Statistical_Data;
     Observed_Count           : in Statistical_Data;
     True_Degrees_of_Freedom  : in Positive;
     Sample_Size              : in Unsigned_64;
     Chi_squared              : out Real;
     P_val, P_val_Variance    : out Real)
  is
     Expected_Count, Sum : Real;
  begin
     Sum := 0.0;
     for i in Probability_Distribution'Range loop
       Expected_Count := Probability_Distribution(i) * Real (Sample_Size);
       Sum := Sum + (Observed_Count(i) - Expected_Count)**2 / Expected_Count;
     end loop;
     Chi_squared    := Sum;
     P_val          := Chi_Squared_CDF (Real (True_Degrees_of_Freedom), Chi_squared);
     P_val_Variance := (P_val-0.5)**2;
  end Get_Chi_Statistic_and_P_val;
 
  ------------------------------
  -- Greatest_Common_Divisors --
  ------------------------------

  -- translated from Marsaglia's diehard suite.
  -- GCD Test, uses pairs of Rand's: u and v
  -- where pairs = Sample_Size.
  -- ***Requires uniform rands on 0..2**Bits_per_Random_Word-1.***

  procedure Greatest_Common_Divisors 
   (Sample_Size : in Unsigned_64;
    Count_of_GCD_Iterations : out GCD_Iterations_Stats) is

    Observed_Count_of_GCDs :  GCD_Counts;

    s, e : Real;
    p99, chi99,variance_p99 : Real;
    ave_chi99, ave_p99, ave_variance_p99 : Real := 0.0;

    p, chi, variance_p : Real;
    ave_p, ave_chi, ave_variance_p : Real := 0.0;

    k : Unsigned_64;
    u, v, w : Unsigned_64;
    u0, v0 : Random_Int;

    No_of_Samples : constant Integer := 2**16;

  begin

    Observed_Count_of_GCDs  := (others => 0.0);
    Count_of_GCD_Iterations := (others => 0.0);

    Outer: for j in 1..No_of_Samples loop

    Observed_Count_of_GCDs  := (others => 0.0);
    Count_of_GCD_Iterations := (others => 0.0);

    for i in Unsigned_64 range 1 .. Sample_Size loop

      Get_Pair: loop
         --Get_Random(u0, Stream_1);
         --Get_Random(v0, Stream_1);
         --u := Unsigned_64 (u0 mod 2**Bits_per_Random_Word);
         --v := Unsigned_64 (v0 mod 2**Bits_per_Random_Word);
         u := Unsigned_64 (rnd.Random (g));
         v := Unsigned_64 (rnd.Random (g));
         exit Get_Pair when (u > 0 and then v > 0);
      end loop Get_Pair;
 
      k := 0;
 
      Euclid: loop
        w := u mod v;
        u := v;
        v := w;
        k := k + 1;
        exit Euclid when v = 0;
      end loop Euclid;
 
      --  k is Observed number of Iterations to obtain greatest common divisor (GCD).
      --  u is the greatest common divisor (GCD).
 
      if k < Count_of_GCD_Iterations'First then
         k := Count_of_GCD_Iterations'First;
      end if;
      if k > Count_of_GCD_Iterations'Last then
         k := Count_of_GCD_Iterations'Last;
      end if;
 
      Count_of_GCD_Iterations(k) := Count_of_GCD_Iterations(k)+1.0;
 
      if u > Observed_Count_of_GCDs'Last then
         u := Observed_Count_of_GCDs'Last;
      end if;
      if u < Observed_Count_of_GCDs'First then
         u := Observed_Count_of_GCDs'First;
      end if;
 
      Observed_Count_of_GCDs(u) := Observed_Count_of_GCDs(u) + 1.0;

   end loop;

   Get_Chi_Statistic_and_P_val 
     (Probability_Distribution => Probability_of_GCD_Iterations,
      Observed_Count           => Count_of_GCD_Iterations,
      True_Degrees_of_Freedom  => True_DOF_for_Iterations_Test,
      Sample_Size              => Sample_Size,
      Chi_squared              => chi,
      P_val                    => p,
      P_val_Variance           => variance_p);
 
    ave_chi := ave_chi + chi;
    ave_p   := ave_p + p;
    ave_variance_p := ave_variance_p + variance_p;
 
 
    -- on range 1..99 distribution seems to be:   (0.607926 + 6.0e-8 * i) / i^2
    -- theoretical value, with inf number of bits: 0.60792710 / i^2
    --
    -- e := Real (Sample_Size) * 0.6081842 / Real (i)**2;--asymptotically, i = 5410
 
    p99 := 0.0;
    variance_p99 := 0.0;
  --e := Real (Sample_Size) * 0.61097691e-2;  -- in theory, p >> 2**32
    e := Real (Sample_Size) * 0.61097e-2;     -- I get 0.61097e-2
    chi99 := (Observed_Count_of_GCDs(GCD_Count_Test_Range'Last) - e)**2 / e;
    for i in GCD_Count_Test_Range'First .. GCD_Count_Test_Range'Last-1 loop
       e := Real (Sample_Size) * (0.607926 + 6.0E-8 * Real (i)) / Real (i)**2;
       s := (Observed_Count_of_GCDs(i) - e)**2 / e;
       chi99 := chi99 + s;
    end loop;
    p99 := Chi_Squared_CDF (Real(True_DOF_for_GCD_Count_Test), chi99);
    variance_p99   := (p99-0.5)**2;
 
    ave_chi99 := ave_chi99 + chi99;
    ave_p99   := ave_p99 + p99;
    ave_variance_p99   := ave_variance_p99 + variance_p99;
 
 
    new_line(1); 
    put("Test"); put (Integer'Image(j));
    put(". Chi^2 (48 dof), ave p-val, and ave normalized variance of GCD iterations:");
    new_line; 
    put("          ");
    put (Real'Image (chi));
    put (Real'Image (ave_p / Real(j)));
    put (Real'Image (ave_variance_p / (Real(j)*(0.25/3.0)))); -- should -> 1.0
 
    new_line(1); 
    put("         Chi^2 (99 dof), ave p-val, and ave normalized variance of GCD's:");
    new_line; 
    put("          ");
    put (Real'Image (chi99));
    put (Real'Image (ave_p99 / Real(j)));
    put (Real'Image (ave_variance_p99 / (Real(j)*(0.25/3.0))));
 
  end loop Outer;


  end Greatest_Common_Divisors;

begin

  rnd.Reset (g);

  Disorderly.Basic_Rand.Clock_Entropy.Reset (Stream_1); 
  --  The state of the generator is Stream_1. (Starts up a random stream.)

  test: declare

    Sample_Size : constant Unsigned_64 := 2**40;  -- turn way up to best see failure
    -- 2**32 Sample_Size is OK, chi squared wise.
    -- 2**37 Sample_Size is gd stnd tst. Tks a day at least (2+ for gnat etc)!

    Full_Sample_Size : Real;
    Sample_Iteration_Stats : GCD_Iterations_Stats;
    Full_Iteration_Stats : GCD_Iterations_Stats := (others => 0.0);

  begin

    for i in 1..2**16 loop

      Greatest_Common_Divisors (Sample_Size, Sample_Iteration_Stats);

      Full_Sample_Size := Real(i)*Real(Sample_Size);

      for k in Full_Iteration_Stats'Range loop
        Full_Iteration_Stats(k) := Full_Iteration_Stats(k) + Sample_Iteration_Stats(k);
      end loop;

      new_line;
      put (Integer'Image (i)); put ("Total Sample_Size:"); 
      put (Real'Image (Real(i)*Real(Sample_Size)));
      new_line;
      for k in Full_Iteration_Stats'Range loop
        if (Integer(k)-Integer(Full_Iteration_Stats'First)) mod 3 = 0 then
          new_line;
        end if;
        put (Real'Image (Full_Iteration_Stats(k) / Full_Sample_Size)); put (",");
      end loop;  
      new_line; 
    end loop;  

  end test;

end;

