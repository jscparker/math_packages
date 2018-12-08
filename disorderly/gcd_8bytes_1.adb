
with Disorderly.Random; use Disorderly.Random;
with Disorderly.Random.Clock_Entropy;
with Chi_Gaussian_CDF;
with text_io; use text_io;

--  translated from Marsaglia, Tsang  diehard suite.

procedure gcd_8bytes_1 is

  Bits_per_Random_Word : constant := 61;
  --  Must set this correctly here.  There's no way to check this.

  Stream_1 : Disorderly.Random.State;
  --  Create a stream of Random numbers
  --  by initializing this after the begin, w/ a call to Reset_with_Calendar
  
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

  subtype GCD_Iterations_Test_Range is Unsigned_64 range 7..57;

  subtype GCD_Iterations_Stats is Statistical_Data (GCD_Iterations_Test_Range);

  True_DOF_for_Iterations_Test : constant := 50;

  type Permissable_Range_of_Bits_per_Word is range 61..64;

  Probability_of_GCD_Iterations : constant array(Permissable_Range_of_Bits_per_Word) 
                                      of GCD_Iterations_Stats := 

  --61bit,  1.5 * 2^40 sample_size
  ((
     3.60634420876918E-11, 1.84113362237163E-10, 9.27526422816774E-10,
     4.18462466610514E-09, 1.79804729945634E-08, 7.02996697762738E-08,
     2.50987637707073E-07, 8.29602156401328E-07, 2.54719381463592E-06,
     7.29226651500263E-06, 1.95308024083953E-05, 4.90890664008001E-05,
     1.16038067220791E-04, 2.58472948031419E-04, 5.43582520786794E-04,
     1.08110812383846E-03, 2.03641299373450E-03, 3.63780116996445E-03,
     6.16967723370813E-03, 9.94544514548540E-03, 1.52515013444075E-02,
     2.22684802100285E-02, 3.09799653111006E-02, 4.10936488194199E-02,
     5.20026693911464E-02, 6.28131812491515E-02, 7.24521456393869E-02,
     7.98344028653051E-02, 8.40652918002730E-02, 8.46145726650518E-02,
     8.14282151139038E-02, 7.49344837759449E-02, 6.59512901418012E-02,
     5.55169851045914E-02, 4.47006642858943E-02, 3.44248149624985E-02,
     2.53562087328527E-02, 1.78608364057373E-02, 1.20299972414100E-02,
     7.74637900821804E-03, 4.76738610228969E-03, 2.80356488466445E-03,
     1.57480260911499E-03, 8.44672639884622E-04, 4.32374513247455E-04,
     2.11148911267114E-04, 9.83079048835308E-05, 4.36203933887831E-05,
     1.84285277298287E-05, 7.40598087482478E-06, 4.38173035787338E-06
   ),

  --62bit,   1.25 * 2^40 sample_size 
   (
    1.0e-11, 1.0e-10, 1.0e-09,
    3.12356860376894E-09, 1.07153027784079E-08, 4.23540768679231E-08,
    1.53834844240919E-07, 5.18621527589858E-07, 1.62184005603194E-06,
    4.72031460958533E-06, 1.28798790683504E-05, 3.29731628880836E-05,
    7.93799634266179E-05, 1.80152335087769E-04, 3.86169017292559E-04,
    7.83033167681424E-04, 1.50444347600569E-03, 2.74219085113145E-03,
    4.74746785839670E-03, 7.81447832123376E-03, 1.22424124856480E-02,
    1.82686513180670E-02, 2.59868027780612E-02, 3.52614222552802E-02,
    4.56673437729478E-02, 5.64818627550267E-02, 6.67421547936101E-02,
    7.53817606542725E-02, 8.14047969237436E-02, 8.40791420108871E-02,
    8.30774562251463E-02, 7.85453180062177E-02, 7.10665999366029E-02,
    6.15405520489730E-02, 5.10077870050736E-02, 4.04669900148292E-02,
    3.07284356771561E-02, 2.23318681637466E-02, 1.55313623552502E-02,
    1.03353681937733E-02, 6.57944859049166E-03, 4.00580784917111E-03,
    2.33198539935984E-03, 1.29759441624628E-03, 6.89835629600566E-04,
    3.50286864704685E-04, 1.69781225849874E-04, 7.85187381552532E-05,
    3.46313805493992E-05, 1.45411046105437E-05, 9.24259074963629E-06
    ),

  --63bit 
   (
     1.0e-11, 1.0e-10, 1.0e-09,
      1.84718373930082E-09, 6.40193320577964E-09, 2.55613485933282E-08,
      9.43073246162385E-08, 3.22844243783039E-07, 1.02684680314269E-06,
      3.04464447253850E-06, 8.45538488647435E-06, 2.20202909986256E-05,
      5.39864840902737E-05, 1.24805929772265E-04, 2.72558390861377E-04,
      5.63293603590864E-04, 1.10343044616456E-03, 2.05119602287595E-03,
      3.62303316342150E-03, 6.08702397312300E-03, 9.73696121127432E-03,
      1.48417563932526E-02, 2.15752666608751E-02, 2.99300470460366E-02,
      3.96473091677763E-02, 5.01783137979146E-02, 6.07052047053003E-02,
      7.02302970366873E-02, 7.77271460929115E-02, 8.23192137850128E-02,
      8.34508585840013E-02, 8.09940879425994E-02, 7.52742528147792E-02,
      6.69985561098656E-02, 5.71142716544273E-02, 4.66340440261774E-02,
      3.64706130903869E-02, 2.73180850017525E-02, 1.95966771325402E-02,
      1.34618205693187E-02, 8.85359909625550E-03, 5.57415191451582E-03,
      3.35856920628430E-03, 1.93597860015871E-03, 1.06730554853129E-03,
      5.62642482691444E-04, 2.83431539173762E-04, 1.36380568619643E-04,
      6.26575392743689E-05, 2.74562980848714E-05, 1.87182404260966E-05
     ),

  --64bit 

    (
     1.0e-11, 1.0e-10, 1.0e-09,
     1.04345316584739E-09, 3.79978802003380E-09, 1.54057084324045E-08,
     5.78617295508997E-08, 2.00539862918150E-07, 6.48949814300674E-07,
     1.95518324390933E-06, 5.52275103271111E-06, 1.46452957172490E-05,
     3.65230218449142E-05, 8.59657468633183E-05, 1.91198751215577E-04,
     4.02582168842653E-04, 8.03587449586808E-04, 1.52299020523464E-03,
     2.74339848013672E-03, 4.70193459930467E-03, 7.67594423026215E-03,
     1.19457277723591E-02, 1.77364065986088E-02, 2.51409167733881E-02,
     3.40445156152176E-02, 4.40656369518870E-02, 5.45457989380035E-02,
     6.45975542736172E-02, 7.32208271168525E-02, 7.94619972674910E-02,
     8.25876607805387E-02, 8.22245617538202E-02, 7.84341962288130E-02,
     7.16953611735031E-02, 6.28069448083503E-02, 5.27326066778591E-02,
     4.24344945974250E-02, 3.27284291460981E-02, 2.41922337119225E-02,
     1.71370716094518E-02, 1.16315170244181E-02, 7.56359276545279E-03,
     4.71090568017745E-03, 2.80972215655816E-03, 1.60432155716990E-03,
     8.76614772803603E-04, 4.58251069641038E-04, 2.29045271690766E-04,
     1.09453786394119E-04, 4.99610013018052E-05, 3.64976355437345E-05)
     );

   --  2^64 distr based on 1.25 * 2^40 sample size

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

  -- from diehard.  
  -- GCD Test, uses pairs of Rand's: u and v
  -- where pairs = Sample_Size.
  -- ***Requires uniform rands on 0..2**No_of_Bits_in_Test-1.***

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

    No_of_Samples : constant Integer := 2**20;

  begin

    Observed_Count_of_GCDs  := (others => 0.0);
    Count_of_GCD_Iterations := (others => 0.0);

    Outer: for j in 1..No_of_Samples loop

    Observed_Count_of_GCDs  := (others => 0.0);
    Count_of_GCD_Iterations := (others => 0.0);

    for i in Unsigned_64 range 1 .. Sample_Size loop

      Get_Pair: loop
        Get_Random (u0, Stream_1);
        Get_Random (v0, Stream_1);
        u := Unsigned_64 (u0);
        v := Unsigned_64 (v0);
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
     (Probability_Distribution => Probability_of_GCD_Iterations (Bits_per_Random_Word),
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
    put(". Chi^2 (47 dof), ave p-val, and ave normalized variance of GCD iterations:");
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

  Disorderly.Random.Clock_Entropy.Reset (Stream_1); 
  --  The state of the generator is Stream_1. (Starts up a random stream.)

  test: declare

    Sample_Size : constant Unsigned_64 := 2**35;  -- turn way up to best see failure
    -- 2**35 Sample_Size is OK, chi squared wise.
    -- 2**36 Sample_Size is unimpeachable.
    -- 2**37 Sample_Size is gd stnd tst. Tks mny hrs!

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

