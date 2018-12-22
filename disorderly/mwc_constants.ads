
-- package MWC_Constants
--
-- Multiply_with_Carry_Constants 
--
-- Multiply with carry is an algorithm for calculating
--
--     a*x mod (a*b-1)
--
-- efficiently. Main benefit: even if a*x overflows 64 bit arithmetic,
-- the multiply with carry algorithm works entirely in 64 bit arithmetic.
-- No extended precision arithmetic is required.
--
-- use g*c mod(c-1) = g mod(c-1)
--     a*x mod(a*b) = a*(x mod b)  (using x = (x/b)*b + x mod b)
--
--     z = z mod c + (z/c)*c
--     z mod (c-1) = (z mod c + (z/c)) mod (c-1)
--
--  let c = a*b, z = a*x
--     a*x mod (a*b-1) = [a*(x mod b)  + x/b] mod (a*b-1)
--
--  or say:     x = x_hi*b + x_lo
--
--    a*x mod (a*b-1) = [(a*b)*x_hi + a*x_lo] mod (a*b-1).
--
--                    = x_hi mod (a*b-1) +  a*x_lo mod (a*b-1) 
--                    = [x_hi + a*x_lo] mod (a*b-1) 
--
--    x_hi = x mod (b) and x_lo = x / b.
--
--  x = a*b-1, and x=0 produce streams of period = 1.
--  x_hi + a*x_lo reaches max of a*b-1 when x = a*b-1. 
--
--  x is in [1 .. ab-2].   NOTE (ab-2) / 2 is also prime.
--  x/b = (ab-2)/b = ((a-1)b + b-2)/b = a-1.
--  x mod b = (ab-2) mod b = ((a-1)b + b-2) mod b = b-2
--  So if x = ab-2, then x_hi + a*x_lo = a-1  + a*(b-2) 
--
--
-- EVALUATING  x mod (2**n - k).  
--
--     (for example in an LCG you want:    a*y mod (2**n - k) 
--
--      b = 2**n = base
--
--      write x in base b:  x = (x/b)*b + x mod b
--
--      x mod m = ((x/b)*b + x mod b) mod (b - k)
--              = ((x/b)*k)  + (x/b)*(b - k) + (x mod b)) mod  (b - k)
--              = ((x/b)*k) + (x mod b)) mod  (b - k)
--
-- e.g.  b = 2**32 and k = 5.
--
--
-- VERIFIED SAFE PRIMES:
--
--
--  Do[ If[ PrimeQ[a*(2^32)  - 1] && 
--          PrimeQ[(a*(2^32) - 1 - 1)/2]  , 
--          Print[(a*(2^32)  - 1), "  ", a, "  ",
--                  ((a - 1594719654)*1.81338)], Null ], 
--          {a, 1594719654 - 30000, 1594719654 + 0000}]
--
--
--    ratio of periods   = (1+sqrt(5))/2       (1.618033..  golden mean)
--    product of periods =  2^124 * 0.9999999960368
--
--    (7250911033471008767*11732400383950061567) / 2^126   = 0.9999999960368
--
--  b=2**32
--
--  m0 = 11732400383950061567  a0 = 2731662333
--  m1 = 7250911033471008767   a1 = 1688234283
--

package MWC_Constants
   with Spark_Mode => On
is
   pragma Pure (MWC_Constants);

   --  If a single parameter mutates, the generator fails, so reassert.
   --
   b0 : constant := 2**32; 
   a0 : constant := 2731662333; -- a0 must be the larger a 
   pragma Assert (a0 = 2731662333);  -- 2731662333, 2731662333
   pragma Assert (b0 = 2**32);       -- 2**32, 2**32
   -- m0 = 11732400383950061567  a0 = 2731662333

   b1 : constant := 2**32; 
   a1 : constant := 1688234283;
   pragma Assert (a1 = 1688234283);  -- 1688234283, 1688234283
   pragma Assert (b1 = 2**32);       -- 2**32, 2**32
   --  m1 = a1*2^32-1 = 7250911033471008767   a1 = 1688234283

   --  Ratio of periods = golden mean, just for fun.
   --  (7250911033471008767*11732400383950061567)/4 = 2^124 * 0.99999999604

   m0 : constant :=  a0 * b0 - 1;  -- LCG modulus (safe prime)
   m1 : constant :=  a1 * b1 - 1;  -- LCG modulus (safe prime)

   p0 : constant := (m0 - 1) / 2; -- period (prime)
   p1 : constant := (m1 - 1) / 2; -- period (prime)

   pragma Assert (a0 < b0);
   pragma Assert (a1 < b1);
   pragma Assert (a1 < a0); 
   pragma Assert ((a1/2)*2 < a1);  -- but only for b1 = 2^32

   --pragma Assert (m0 < Parent_Random_Int'Last); 
   --pragma Assert (m1 < Parent_Random_Int'Last); 

   --pragma Assert (m0 > Random_Int'Last + 1);  -- for Reset to work.
   --pragma Assert (m1 > Random_Int'Last + 1);  -- ie m1 > 2**61

   -- b0 : constant := 2**32-5;   -- too slow usually
   -- a0 : constant := 2731663104; -- a0 must be the larger a 
   -- m0 = a0*(2^32-5)-1 =  11732403681711531263 
   --
   -- b0 : constant := 2**32-5;
   -- a0 : constant := 2147476584; -- FOR TESTING
   -- m0 = a0*(2^32-5)-1 =  9223341686468413943  = 2^63 * .999999978
   --
   -- b0 : constant := 2**16+1;
   -- a0 : constant := 57684; -- FOR TESTING
   -- m0 = a0*b0-1 = 3780436307; p0 = 1890218153
   --
   -- b0 : constant := 2**16+1;
   -- a0 : constant := 59592; -- FOR TESTING
   -- m0 = a0*b0-1 = 3905480903; p0 = 1952740451
   --
   -- b0 : constant := 2**16+1;
   -- a0 : constant := 59964; -- FOR TESTING
   -- m0 = a0*b0-1 = 3929860667; p0 = 1964930333
   --
   -- b0 : constant := 2**16+1;
   -- a0 : constant := 49380; -- FOR TESTING
   -- m0 = a0*b0-1 = 3236217059; p0 = 1618108529
   --
   -- m0 : constant := a0 * b0 - 1; -- For the testing parameters above
   -- p0 : constant := (m0-1) / 2;  -- For the testing parameters above

end MWC_Constants;

