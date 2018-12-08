
--  Procedure multi_var_deviates_demo_1, test and demonstrate random deviates.
--
--  Uses a giant array, so may need to set stacksize large at (eg) Linux 
--  prompt, as in: ulimit -s unlimited

with Text_io; Use Text_io;
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


procedure multi_var_deviates_demo_1 is

   type Integer64 is range -2**63+1..2**63-1;

   type Real is digits 15;

   package dev is new Disorderly.Basic_Rand.Deviates (Real); 
   use dev;

   Stream_1 : State; -- exported by Disorderly.Basic_Rand
 
   Delta_X_stnd : constant Real  := 0.01; 
   -- If you change parameters in the distributions, may have make this much
   -- smaller. (For example, the Beta distribution becomes very sharply 
   -- peaked if aa<1, or if aa>>1; same with bb. Same with Gamma distribution
   -- if you make s<1; etc. For example, with Beta distribution b = .9, needed
   -- Delta_X_stnd = 0.00004, and HalfBins = 4.0 / Delta_X_stnd)

   HalfBins     : constant Integer64 := Integer64 (25.0 / Delta_X_stnd); 

   Delta_X : Real := Delta_X_stnd;
   --  Delta_X must be 1 for Poisson, Binomial, Neg_Binomial!  


   --  For Multivariate_Normal:

   Index_First : constant Positive := 1;
   Index_Last  : constant Positive := 2;
   Mean        : constant Vector(Index_First..Index_Last) := (1.0, 2.0);
   X_vec       : Vector(Index_First..Index_Last);
   Covariance       : Matrix(Mean'Range, Mean'Range) := (others => (others => 1.0));
   LU_of_Covariance : Matrix(Mean'Range, Mean'Range) := (others => (others => 0.0));
   MVN_Init         : MV_Normal_Initializer;

 
   Sample_Size : Integer64;
 
   subtype Bin_Range is Integer64 range -HalfBins+1 .. HalfBins;
   type Data is array(Bin_Range, Bin_Range) of Real;

   Histogram    : Data;
   Distribution : Data;
 
   Norm, Sum : Real := 0.0;
   X, D   : Real := 0.0;
 
begin

   --  Init covariance matrix for multivariate normal.
   --  Try to make a positive definite matrix:
   X := 0.12345;
   for i in Index_First .. Index_Last loop
   for j in Index_First .. Index_Last loop
      X := X + 0.12345;
      Covariance(i,j) := X;
   end loop;
   end loop;

   for i in Index_First .. Index_Last loop
   for j in i .. Index_Last loop
      Covariance(i,j) := Covariance(j,i);
   end loop;
   end loop;

   for i in Index_First .. Index_Last loop
      Covariance(i,i) := Covariance(i,i) + 3.0 * Real (i);
   end loop;

   -- For this we need the Cholesky (LU) Decomposition of the Covariance matrix:

   Choleski_Decompose 
     (Covariance,
      LU_of_Covariance); -- Choleski Decomp of Covariance matrix.
   
   --  use Calendar to get initial state: Stream_1:

   Disorderly.Basic_Rand.Clock_Entropy.Reset (Stream_1);  

   -- Next make normalized distribution:
  
   Distribution := (others => (others => 0.0));

   for Bin_id_2 in Bin_Range loop
   for Bin_id_1 in Bin_Range loop

      X_vec(1) := Delta_X * (Real (Bin_id_1) - 0.5);
      X_vec(2) := Delta_X * (Real (Bin_id_2) - 0.5);

      D := Multivariate_Normal_Probability (Mean, LU_of_Covariance, X_Vec);
 
      Distribution (Bin_id_1, Bin_id_2) := D;

   end loop;
   end loop;
  
 
   new_line(2);
   put ("Presently calculating variance of distance between the observed");
   new_line(1);
   put ("distribution of a sample of N random deviates, and the exact"); 
   new_line(1);
   put ("distribution they are meant to obey:");
   new_line(1);


   Delta_X := Delta_X_stnd;

   Sample_Size := 2_000; -- initial sample size is 10x this.


   for Resized_Sample_Size in 1 .. 16 loop

      Histogram    := (others => (others => 0.0));
   
      Sample_Size :=  Sample_Size * 10;

      new_line(2); 
      put ("Using Sample size N                   =  ");
      put (Integer64'Image (Sample_Size));

      for i in Integer64 range 1 .. Sample_Size loop
   
         Get_Multivariate_Normal (Mean, LU_of_Covariance, MVN_Init, Stream_1, X_vec);
         declare 
            Bin_id_1 : constant Integer64 := Integer64 (X_vec(1) / Delta_X + 0.5);
            Bin_id_2 : constant Integer64 := Integer64 (X_vec(2) / Delta_X + 0.5);
         begin
            if Bin_id_1 in Bin_Range then
            if Bin_id_2 in Bin_Range then
               Histogram(Bin_id_1, Bin_id_2) := Histogram(Bin_id_1, Bin_id_2) + 1.0;
            end if;
            end if;
         end;
   
      end loop;
   
      -- Normalize the curves.  (Normalize the curves
      -- the same way that the distribution curves generated below
      -- are normalized: integrate over X with dX.) Here
      -- dX is Delta_X, so multiply by Delta_X at the end.
   
      Norm := 0.0;
      for Bin_id_2 in Bin_Range loop
      for Bin_id_1 in Bin_Range loop
         Norm := Norm + Histogram(Bin_id_1, Bin_id_2);
      end loop;
      end loop;
      Norm := Norm*Delta_X**2;
   
      for Bin_id_2 in Bin_Range loop
      for Bin_id_1 in Bin_Range loop
         Histogram(Bin_id_1, Bin_id_2) := Histogram(Bin_id_1, Bin_id_2) / Norm;
      end loop;
      end loop;
   
      Sum := 0.0;
      for Bin_id_2 in Bin_Range loop
      for Bin_id_1 in Bin_Range loop
       Sum := Sum + (Distribution(Bin_id_1, Bin_id_2)-Histogram(Bin_id_1, Bin_id_2))**2;
      end loop;
      end loop;
      Sum := Sum*Delta_X**2;
      
      new_line; 
      put ("Variance (standard deviation squared) =");
      put ("  "); put (Real'Image (Sum));
      new_line; 
      put ("Stnd_Deviation**2 should go approximately as 1/N.");
      new_line; 
      put ("The actual value of the variance will fluctuate statistically.");
      new_line; 
      put ("It will run as long as you let it. Control-c to escape");
      new_line; 

   end loop;   
   
end;
