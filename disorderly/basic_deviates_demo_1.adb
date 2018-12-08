
--  Uses a giant array, so may need to set stacksize large at (eg) Linux 
--  prompt, as in: ulimit -s unlimited

with Text_io; Use Text_io;
with Ada.Numerics.Generic_Elementary_Functions;
with Disorderly.Basic_Rand; use Disorderly.Basic_Rand;
with Disorderly.Basic_Rand.Clock_Entropy;
with Disorderly.Basic_Rand.Deviates;

-- To use Random.Deviates, you need to "with" Disorderly.Basic_Rand.
-- which provides the random num generator (which you just ignore), 
-- and the essential    Reset (Stream_x, Seeds ...) 
-- which allows you to create multiple independent streams
-- of Random nums.
--
-- Below we use the Clock version of Reset to get the initial state (Stream).
-- The clock provides the Seeds. You can get all the streams you want:
--
--   Disorderly.Basic_Rand.Clock_Entropy.Reset (Stream_1);  
--   Disorderly.Basic_Rand.Clock_Entropy.Reset (Stream_2);  
--   Disorderly.Basic_Rand.Clock_Entropy.Reset (Stream_3);  


procedure basic_deviates_demo_1 is

   type Integer64 is range -2**63+1..2**63-1;

   type Real is digits 15;

   package math is new Ada.Numerics.Generic_Elementary_Functions (Real); 
   use math;
   package rio is new Float_IO(Real); 
   use rio;
   package dev is new Disorderly.Basic_Rand.Deviates (Real); 
   use dev;

   Stream_1 : State; -- exported by Disorderly.Basic_Rand
 
   Chosen_Distribution : Integer;
 
   Delta_X_stnd : constant Real  := 0.004; 
   -- If you change parameters in the distributions, may have make this much
   -- smaller. (For example, the Beta distribution becomes very sharply 
   -- peaked if aa<1, or if aa>>1; same with bb. Same with Gamma distribution
   -- if you make s<1; etc. For example, with Beta distribution b = .9, needed
   -- Delta_X_stnd = 0.00004, and HalfBins = 4.0 / Delta_X_stnd)

   HalfBins     : constant Integer64 := Integer64 (400.0 / Delta_X_stnd); 
   --  need 400+ for Cauchy; much less for others, 
   --  so use 400. Means that results are stored in a giant array.

   Delta_X : Real := Delta_X_stnd;
   --  Delta_X must be 1 for Poisson, Binomial, Neg_Binomial!  

   --  For Chi-Squared: 

   Degrees_of_Freedom : constant Real := 3.0;
   Chi_init    : Chi_Initializer;

   --  For Cauchy:

   S_Cauchy : constant Real := 0.25; -- Narrows up Lorentzian (tails fall off slow)
   Pi       : constant Real := 3.141_592_652_589_793_238;

   --  For Gamma:

   s : constant Real := 3.85;
   Gamma_init    : Gamma_Initializer;

   --  For Beta:

   aa : constant Real := 3.3;
   bb : constant Real := 2.2;
   Beta_init    : Beta_Initializer;

   --  For Student_t:

   Student_t_m : constant Positive := 8; -- degrees of freedom
   S_t_init    : Student_t_Initializer;

   --  For Neg_Binomial:

   Neg_Binomial_r : constant Real := 10.4;
   Neg_Binomial_p : constant Real := 0.21;
   NB_Init        : Neg_Binomial_Initializer;

   --  For Binomial:

   Binomial_n : constant Positive := 2000;
   Binomial_p : constant Real     := 0.1;
   B_Init     : Binomial_Initializer;

   --  For Poisson:

   Poisson_Mean : constant Real := 150.0;
   P_Init : Poisson_Initializer;

   --  For Exponential:

   Exponential_Mean : constant Real := 9.0;

   --  For Normal (Gaussian):

   Standard_Dev   : constant Real := 6.0;
   Gaussian_Mean  : constant Real := 24.0;
   N_Init : Normal_Initializer;
 
   --  For Log_Normal:

   Lognormal_Sigma : constant Real := 0.5;
   Lognormal_Mean  : constant Real := 0.0;
   LN_Init : Log_Normal_Initializer;

 
   Sample_Size : Integer64;
 
   subtype Bin_Range is Integer64 range -HalfBins+1 .. HalfBins;
   type Data is array(Bin_Range) of Real;

   Histogram    : Data := (others => 0.0);
   Distribution : Data := (others => 0.0);
 
   Choice    : Real := 0.0;
   Norm, Sum : Real := 0.0;
   X, D, a   : Real := 0.0;
 
   Observed_Bin_id : Integer64;
   IX : Integer;

begin

   --  use Calendar to get initial state: Stream_1:
   Disorderly.Basic_Rand.Clock_Entropy.Reset (Stream_1);  

   new_line;
   put ("Choose A Distribution. Enter a number:");
   new_line(2);
   put ("Uniform  = 0, Normal    = 1,  Lorentzian  = 2");
   new_Line;
   put ("Rayleigh = 3, Student_t = 4,  Exponential = 5");
   new_Line;
   put ("Weibull  = 6, Poisson   = 7,  Binomial    = 8");
   new_Line;
   put ("Beta     = 9, Gamma     = 10, Chi_Squared = 11");
   new_line;
   put ("Neg_Binomial = 12, Log_Normal = 13");
   new_line;
   get (Choice);
   Chosen_Distribution := Integer(Choice);

   new_line(2);
   put ("Presently calculating variance of distance between the observed");
   new_line(1);
   put ("distribution of a sample of N random deviates, and the exact"); 
   new_line(1);
   put ("distribution they are meant to obey:");
   new_line(1);


   if Choice = 7.0 or Choice = 8.0 or Choice = 12.0 then
      Delta_X := 1.0; -- Poisson, Binomial
   else
      Delta_X := Delta_X_stnd;
   end if;


   Sample_Size := 2_000; -- initial sample size is 10x this.


   for Resized_Sample_Size in 1 .. 16 loop

      Histogram    := (others => 0.0);
      Distribution := (others => 0.0);
   
      Sample_Size :=  Sample_Size * 10;

      new_line(2); 
      put ("Using Sample size N                   =  ");
      put (Integer64'Image (Sample_Size));


      for i in Integer64 range 1 .. Sample_Size loop
   
         if Chosen_Distribution = 0 then -- Uniform
            Get_Random_Real (X, Stream_1);
         elsif Chosen_Distribution = 1 then -- Normal
            Get_Normal (Gaussian_Mean, Standard_Dev, N_Init, Stream_1, X);
         elsif Chosen_Distribution = 2 then -- Cauchy (Lorentzian)
            Get_Cauchy (S_Cauchy, Stream_1, X);
         elsif Chosen_Distribution = 3 then -- Rayleigh
            Get_Rayleigh (Stream_1, X);
         elsif Chosen_Distribution = 4 then -- Student_t
            Get_Student_t (Student_t_m, S_t_init, Stream_1, X);
         elsif Chosen_Distribution = 5 then -- Exponential
            Get_Exponential (Exponential_Mean, Stream_1, X);
         elsif Chosen_Distribution = 6 then -- Rayleigh really
            Get_Weibull (2.0, Stream_1, X);
         elsif Chosen_Distribution = 7 then -- Poisson
            Get_Poisson (Poisson_Mean, P_Init, Stream_1, X);
         elsif Chosen_Distribution = 8 then -- Binomial
            Get_Binomial (Binomial_n, Binomial_p, B_Init, Stream_1, X);
         elsif Chosen_Distribution = 9 then -- Beta
            Get_Beta (aa, bb, Beta_Init, Stream_1, X);
         elsif Chosen_Distribution = 10 then -- Gamma
            Get_Gamma (s, Gamma_Init, Stream_1, X);
         elsif Chosen_Distribution = 11 then -- Chi_Squared
            Get_Chi_Squared (Degrees_of_Freedom, Chi_Init, Stream_1, X);
         elsif Chosen_Distribution = 12 then -- Neg_Binomial
            Get_Neg_Binomial (Neg_Binomial_r, Neg_Binomial_p, NB_Init, Stream_1, X);
         elsif Chosen_Distribution = 13 then -- Log_Normal
            Get_Log_Normal (Lognormal_Mean, Lognormal_Sigma, LN_Init, Stream_1, X);
         end if;
     
         Observed_Bin_id := Integer64 (X  / Delta_X + 0.5);
         if Observed_Bin_id in Bin_Range then
            Histogram(Observed_Bin_id) := Histogram(Observed_Bin_id) + 1.0;
         end if;
   
      end loop;
   
      -- Normalize the curves.  (Normalize the curves
      -- the same way that the distribution curves generated below
      -- are normalized: integrate over X with dX.) Here
      -- dX is Delta_X, so multiply by Delta_X at the end.
   
      Norm := 0.0;
      for Bin_id in Bin_Range loop
         Norm := Norm + Histogram(Bin_id);
      end loop;
      Norm := Norm*Delta_X;
   
      for Bin_id in Bin_Range loop
         Histogram(Bin_id) := Histogram(Bin_id) / Norm;
      end loop;
   
       -- Next make normalized distribution:
  
      for Bin_id in Bin_Range loop
   
         X := Delta_X * (Real (Bin_id) - 0.5);
         D := 0.0;
    
         if Chosen_Distribution = 0 then -- Uniform of [0,1)
	    if X >= 0.0  and then X < 1.0 then
               D := 1.0;
	    else
	       D := 0.0;
	    end if;
         elsif Chosen_Distribution = 1 then -- Normal
	    D   := Normal_Probability (Gaussian_Mean, Standard_Dev, X);
         elsif Chosen_Distribution = 2 then -- Cauchy (Lorentzian)
            a := S_Cauchy;
	    D := a / ((a*a + x*x) * Pi);
         elsif Chosen_Distribution = 3 then -- Rayleigh 
            if X > 0.0 then
	       D := 2.0 * x * Exp (-x*x);
            else
               D := 0.0;
            end if;
         elsif Chosen_Distribution = 4 then -- student_t
	    D := Student_t_Probability (Student_t_m, X);
         elsif Chosen_Distribution = 5 then -- Exponential
	    if X > 0.0 then
	       D := Exp (-X / Exponential_Mean) / Exponential_Mean;
	    else
	       D := 0.0;
	    end if;
         elsif Chosen_Distribution = 6 then -- Weibull 
            if X > 0.0 then
	       D := 2.0 * x * Exp (-x*x);
            else
               D := 0.0;
            end if;
         elsif Chosen_Distribution = 7 then -- Poisson
	    IX := Integer (Bin_id) - 1;
            D  := Poisson_Probability (Poisson_Mean, IX);
         elsif Chosen_Distribution = 8 then -- Binomial
	    IX := Integer (Bin_id) - 1;
	    D  := Binomial_Probability (Binomial_n, IX, Binomial_p);
         elsif Chosen_Distribution = 9 then -- Beta
            if X > 0.0  and then X < 1.0 then
	       D := Beta_Probability (aa, bb, X);
            else
               D := 0.0;
            end if;
         elsif Chosen_Distribution = 10 then -- Gamma
            if X > 0.0 then
	       D := Gamma_Probability (s, X);
            else
               D := 0.0;
            end if;
         elsif Chosen_Distribution = 11 then -- Chi_Squared
            if X > 0.0 then
	       D := Chi_Squared_Probability (Degrees_of_Freedom, X);
            else
               D := 0.0;
            end if;
         elsif Chosen_Distribution = 12 then -- Neg_Binomial
	    IX := Integer (X - 0.5); -- Integer (Bin_id) - 1
            D  := Neg_Binomial_Probability (Neg_Binomial_r, IX, Neg_Binomial_p);
         elsif Chosen_Distribution = 13 then -- Log_Normal
            D  := Log_Normal_Probability (Lognormal_Mean, Lognormal_Sigma, X);
         end if;
    
         Distribution (Bin_id) := D;

       --if bin_id >= 0 then
          --t := Abs (Poissonian(Bin_id-1) - D);
          --if t>0.0 then
             --put(Integer64'Image (bin_id-1));  put (Real'Image (t));
          --end if;
       --end if;

      end loop;
   
      Sum := 0.0;
      for j in Bin_Range loop
         Sum := Sum + (Distribution(j) - Histogram(j))**2;
      end loop;
      Sum := Sum*Delta_X;
      
      new_line; 
      put ("Variance (standard deviation squared) =");
      put ("  "); put (Real'Image (Sum));
      new_line; 
      put ("Stnd_Deviation**2 should go approximately as 1/N.");
      new_line; 
      put ("The actual value of the variance will fluctuate statistically.");
      new_line; 
      put ("Fluctuations are very large for Poisson, Binomial - run repeatedly.");
      new_line; 
      put ("This will run as long as you let it. Control-c to escape");
      new_line; 

   end loop;   
   
end;
