
-------------------------------------------------------------------------------
-- package body Disorderly.Random, Non-linear Random Number Generator
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

package body Disorderly.Random is

  --  MWC stands for Multiply_with_Carry. There are 2 MWC generators,
  --  which are combined to make a uniform full-period generator.
  --  We call this combined generator one of the 3 component generators.
  --  The state (composed of 256 bits) is contained in an array of four
  --  64 bit ints indexed by the following 4 ints:

  MWC_id_0     : constant State_Index := 0;
  MWC_id_1     : constant State_Index := 1;
  Xor_Shift_id : constant State_Index := 2;
  Quadratic_id : constant State_Index := 3;

  --------------------------------
  -- Get_Random_X_squared_Mod_P --
  --------------------------------

  --  X_new = (X_old * X_old) mod P,  followed by P - X_new if X_new > (P-1)/2. 
  --
  -- The final step makes the generator (with full_period_desired = true)
  -- full period if P mod 16 = 7 and if P is a doubly safe prime; half period
  -- if P mod 16 = 15, (which seems to be ok even as a stand-alone 
  -- generator). If full_period_desired = false then they are half and quarter
  -- period respectively.
  --
  -- The real cost here is the 128 bit modulo operation: X*X MOD P  
  -- for prime P, and 128 bit X*X in:
  --
  --            (X*X = 128bit) mod P.
  --
  -- The  X*X  in 64 bit modular (unsigned) arithmetic is now cheap on pc's,
  -- but not the rest. The trick is (1) write X = (X/b)*b + X mod b for some 
  -- b<=2^32.  b is called the base. Then (2) find a special (safe-safe) prime 
  -- such that P = a*b-1, and (P-1)/2 is prime and (P-3)/4 is prime. If
  -- P mod 16 is 7 then generator can be made full period, w/ period = 2*Large-
  -- prime. If P mod 16 = 15, then generator can be at most half period, 
  -- with period = Large prime.
  -- The rest is straightforward but not necessarily pleasant.
  --
  -- use g*c mod(c-1) = g mod(c-1)
  --     a*x mod(a*b) = a*(x mod b)  (using x = (x/b)*b + x mod b)
  --     z = z mod c + (z/c)*c
  --     z mod (c-1) = (z mod c + (z/c)) mod (c-1)
  --  let c = a*b, z = a*x
  --     a*x mod (a*b-1) = [a*(x mod b)  + x/b] mod (a*b-1)
  --
  -- z*z mod (c-1) ?   use z = z mod a + (z/a)*a = alfa + beta*a;  c= ab, a<b
  --                       alfa = z mod a < a so alfa**2 < ab-1
  --
  -- we now have terms like a*x mod (a*b-1) thanks to the  (z/a)*a term above.
  -- Want ((alfa + beta*a)^2)  mod  (a*b-1) 
  --           = [ alfa^2 + a*[(2*alfa*beta) mod b] + (2*alfa*beta)/b 
  --                + a*[(beta^2*a) mod b] + (beta^2*a)/b  ] mod  (a*b-1)
  --
  -- This can all be done in 64 bit arithmetic even if a, b, and z are near 2**32,
  -- and for example (beta^2*a) is of the order 96 bits.  We use doubly safe
  -- primes of the form a*b-1.
  --
  --  IMPORTANT: must have P = a*b-1 < 2^63
  --
  --  Mathematica script for b = (2^32-5): 
  -- 
  --    Do[If[ PrimeQ[a*(2^32-5) - 1] && 
  --           PrimeQ[(a*(2^32-5) - 1 - 1)/2]  && 
  --           PrimeQ[(a*(2^32-5) - 1 - 3)/4] , 
  --           Print[(a)], Null], 
  --           {a, 2^31 - 5000000, 2^31 + 0 }]
  --
  --  for b = (2^31-1): 
  -- 
  --    Do[If[ PrimeQ[a*(2^31-1) - 1] && 
  --           PrimeQ[(a*(2^31-1) - 1 - 1)/2]  && 
  --           PrimeQ[(a*(2^31-1) - 1 - 3)/4] , 
  --     Print[(a)], Null], {a, 2^32 - 50000, 2^32 + 0 }]
  --
  -- the hard terms are (beta^2*a) mod b and beta^2*a/b. 
  -- Use (ab) mod c = (a mod c)(b mod c) mod c  for the 1st term 
  -- (just all in mod b arith).
  -- For the (beta^2*a/b) use 
  --     beta^2 = beta^2 mod b + (beta^2/b)b =  (lo32 bits + hi32)
  -- and so (beta^2*a/b) = (a*(beta^2 mod b))/b + a* (beta^2/b)
  --
  procedure Get_Random_X_squared_Mod_P
    (Random_x :    out Parent_Random_Int;
     S        : in out State)
  is
     Result : Parent_Random_Int;

     alfa, beta, beta_a, alfa_sq, two_alfa_beta : Parent_Random_Int;
     beta_a_div_b, beta_a_mod_b                 : Parent_Random_Int;
     beta_sq_a_div_b, beta_sq_a_mod_b           : Parent_Random_Int;
     two_alfa_beta_div_b, two_alfa_beta_mod_b   : Parent_Random_Int;
     beta_a_mod_b_times_beta                    : Parent_Random_Int;
     beta_a_mod_b_times_beta_div_b              : Parent_Random_Int;

     --  Option X: with near power-of-2 period (for tests)
     --
     --  b : constant := 2**32;
     --  a : constant := 2147483085;    -- a = 2147483085    
     --  P = a*b-1 = 9223369618788188159 = 0.99999974 * 2^63 
     --  P mod 16 = 15

     --  Option A: with near power-of-2 period
     --
     --  b : constant := 2**32-5;       -- 4294967291
     --  a : constant := 2146582440;    -- a = 2146582440    
     --  P = a*b-1 = 9219501367234970039 = 0.99958 * 2^63 
     --  P mod 16 = 7
     --  if the parameters mutate, then generator fails, so write them thrice:
     --  pragma Assert (a=2146582440 and (a=2146582440) and (a=2146582440)); 
     --  pragma Assert ((b=2**32-5) and (b=2**32-5) and (b=2**32-5)); 

     --
     --  Option B: power-of-2 base b = 2**log_of_b
     --  Bit faster, but can never be full period.
     --  On 32 bit machines, gnat doesn't give the expected speed-up, unless
     --  you write shifts (eg xx / 2**log_of_b)  as  Shift_Right (xx, log_of_b);
     --
     --  b : constant := 2**32; 
     --  a : constant := 1947103035;  
     --  P = a*b-1 = 8362743857267343359 = 1.813381 * 2^62 = 0.9067 * 2^63
     --  period = (P-3)/4 = 2090685964316835839 = 1.813381 * 2^60
     --  P mod 16 = 15

     -- Option C: preferred provided an additional step is used to make results 
     --           uniform on 0..2**n-1, (since then the implicit mod 2**61 does
     --           significant additional bit mixing if P is far from power-of-2).

     b : constant := 2**32-5;       -- b = 4294967291
     a : constant := 1952409288;    -- a = 1952409288    
     --  P =  a*b-1 = 8385534030604598807 = 0.90916 * 2^63 
     --  P mod 16 = P mod 8 = 7
     --  if the parameters mutate, then generator fails, so write them thrice:
     pragma Assert (a=1952409288 and (a=1952409288) and (a=1952409288)); 
     pragma Assert ((b=2**32-5) and (b=2**32-5) and (b=2**32-5)); 

     P  : constant := a*b-1;      -- prime
     P2 : constant := (P-1)/2;    -- prime
     P_minus_1 : constant := P-1;
 
     pragma Assert (P > Random_Int'Last + 1); -- for Reset to work
     pragma Assert (P < Parent_Random_Int'Last / 2);
     pragma Assert (a < b); 

     --  (not a) or b  is same as a implies b
     --  want full_period_desired implies P mod 16 = 7
     --  ie when full_period_desired is true, it must be true that P mod 16 = 7.
     pragma Assert (not Full_Period_Quadratic_Desired or (P mod 16 = 7)); 
 
  begin
 
     beta   := S.X(Quadratic_id) / a; 
 
     beta_a := beta*a;
     alfa   := S.X(Quadratic_id) - beta_a;  -- alfa = X mod a < a 
 
     alfa_sq       := alfa*alfa;
     two_alfa_beta := 2*alfa*beta;
 
     beta_a_div_b := beta_a  /  b;
     beta_a_mod_b := beta_a - b*beta_a_div_b;
 
     beta_a_mod_b_times_beta       := beta_a_mod_b * beta;
     beta_a_mod_b_times_beta_div_b := beta_a_mod_b_times_beta / b;
 
     beta_sq_a_div_b := beta_a_mod_b_times_beta_div_b + beta * beta_a_div_b; 
     -- tricky part. to avoid overflow, we used (beta_a*beta) / b,
     -- and expanded beta_a = b*beta_a_div_b + beta_a_mod_b
 
     beta_sq_a_mod_b := beta_a_mod_b_times_beta - b*beta_a_mod_b_times_beta_div_b;
     -- should be   (beta*beta*a) mod b = ((beta*a) mod b) * beta) mod b 
     -- (just all in mod b arithmetic; can do that for arbitrary b).
     -- Could  alfa*beta*2 = 2*a*b-eps  overflow 64 bit?
     -- Not as long as we have a<2^31 if b is 2^32: so again a*b-1< 2^63
 
     two_alfa_beta_div_b := two_alfa_beta  /  b; 
     two_alfa_beta_mod_b := two_alfa_beta - b*two_alfa_beta_div_b;
 
     Result := alfa_sq + two_alfa_beta_div_b + beta_sq_a_div_b; 
     if Result > P then Result := Result - P; end if; 
 
     Result := Result + a * two_alfa_beta_mod_b; 
     if Result > P then Result := Result - P; end if; 
    
     Result := Result + a * beta_sq_a_mod_b; 
     if Result > P then Result := Result - P; end if; 

     -- during testing:
     --if Result > P then raise Constraint_Error; end if; 

     --if Result > P then Result := Result - P; end if; 
     --Note we use P (not P-1) since if  Result = P then period = 0 anyway.

     S.X(Quadratic_id) := Result;
     Random_x          := Result;

     --  Random_x is now in range 2 .. p-1, but only about half the numbers in
     --  that range are returned. p2 = (p-1)/2 of them to be precise. Which half 
     --  depends on the initial condition. All it takes to put them in range
     --  1..p2 is p - Random_x. Surprizingly, with that transformation, 
     --  each number in the range  1..p2 is returned exactly once during the
     --  period of length p2.

     if Full_Period_Quadratic_Desired then
        if Random_x > P2 then 
           Random_x := P_minus_1 - Random_x;
        else
           Random_x := Random_x - 1;
        end if;
     end if;
     --  Now Random_x is in  0 .. P2-1. 
     --  (x=p-1 is transformed into 0.)
     --  (x=p-2 is transformed into 1.)
     --  (x=2   is transformed into 1.)
     --  (x=p2+0 is transformed into p2-1.)
     --  (x=p2+1 is transformed into p2-1.)
     --  (x=p2-1 is transformed into p2-2.)
     --  (x=p2+2 is transformed into p2-2.)
     --  Optional, since its not clear if this is a benefit in a
     --  combination generator in which the other components of the generator
     --  are full period (guaranteeing overall full-period generation).

  end Get_Random_X_squared_Mod_P;

  pragma Inline (Get_Random_X_squared_Mod_P);


  -----------------------------------
  -- Multiply_with_Carry_Constants --
  -----------------------------------

  -- MULTIPLY WITH CARRY:
  --
  -- want a*x mod (a*b-1)
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

  package Multiply_with_Carry_Constants is

     --  If a single parameter mutates, the generator fails, so assert thrice.
     --
     b0 : constant := 2**32; 
     a0 : constant := 2731662333; -- a0 must be the larger a 
     pragma Assert ((a0 = 2731662333) and (a0 = 2731662333) and (a0 = 2731662333));
     pragma Assert ((b0 = 2**32) and (b0 = 2**32) and (b0 = 2**32));
     -- m0 = 11732400383950061567  a0 = 2731662333

     b1 : constant := 2**32; 
     a1 : constant := 1688234283;
     pragma Assert ((a1 = 1688234283) and (a1 = 1688234283) and (a1 = 1688234283));
     pragma Assert ((b1 = 2**32) and (b1 = 2**32) and (b1 = 2**32));
     --  m1 = a1*2^32-1 = 7250911033471008767   a1 = 1688234283

     --b0 : constant := 2**32-5;
     --a0 : constant := 2731663104; -- a0 must be the larger a 
     --pragma Assert ((a0 = 2731663104) and (a0 = 2731663104) and (a0 = 2731663104));
     --pragma Assert ((b0 = 2**32-5) and (b0 = 2**32-5) and (b0 = 2**32-5));
     -- m0 = a0*(2^32-5)-1 =  11732403681711531263 

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

     pragma Assert (m0 < Parent_Random_Int'Last); 
     pragma Assert (m1 < Parent_Random_Int'Last); 

     pragma Assert (m0 > Random_Int'Last + 1);  -- for Reset to work.
     pragma Assert (m1 > Random_Int'Last + 1);  -- ie m1 > 2**61

  end Multiply_with_Carry_Constants;


  -------------------------
  -- Get_Random_MWC_64_2 --
  -------------------------

  -- Marsaglia-Zaman multiply-with-carry algorithm is used to construct
  -- two LCG's (linear congruential generators) with prime number moduli.
  -- The 2 prime number moduli p1 and p2 satisfy p1*p2 ~ 2^p. They are
  -- combined to make a generator uniform on 0..2**p-1. 
  -- They're not really LCG's.
  -- A slight modification at the end makes each of the 2 LCG's full-period,
  -- so that they can be combined in way that mixes bits rather more
  -- strongly than simple addition (modulo some prime), and ensures very 
  -- high degree of uniformity on 0 .. 2**p-1. 
  -- Reminder: The full 63+ bits produced by a single MWC generator 
  -- should never be used alone. (The lower 32 bits are pretty good, but
  -- its the combination generator below that meets minimal standards.) 

  procedure Get_Random_MWC_64_2
    (Random_x :    out Parent_Random_Int;
     S        : in out State)
  is
     S_lo, S_hi : Parent_Random_Int;
     S0, S1     : Parent_Random_Int;

     use Multiply_with_Carry_Constants;

  begin

     S_hi   := S.X(MWC_id_0)  /  b0;
     S_lo   := S.X(MWC_id_0) mod b0;
     S0     := a0 * S_lo + S_hi;

     S.X(MWC_id_0) := S0;


     --  S0 is in  1 .. m0-1.   m0-1 equals 2*p0, where p0 is prime.

     if S0 > p0 then 
        S0 := m0 - S0 - 1;
     else
        S0 := S0 - 1;
     end if;
     --  Each number in range  0 .. p0-1  appears exactly once during period p0.
     --  So S0 is full-period with prime number period p0.
     --  (If x=p+1,  m - x = 2p+1  - (p+1)  = p, so m-x-1 = p-1.)


     S_hi   := S.X(MWC_id_1)  /  b1;
     S_lo   := S.X(MWC_id_1) mod b1;
     S1     := a1 * S_lo + S_hi;

     S.X(MWC_id_1) := S1;

     if S1 > p1 then 
        S1 := m1 - S1 - 1;
     else
        S1 := S1 - 1;
     end if;
     --  Now S1 is in  0 .. p1-1 

     Random_x := S0 + p0 * S1;

     -- Requires modular arithmetic, 2^n = 2^64. Requires p0 > p1.
     -- The above is equivalent to doing the arithmetic in full precision, 
     -- followed by mod 2^64, because (with for example m = 2^64):
     --   (a*b + c) mod m = ((a mod m) * (b mod m) + c mod m) mod m.
     -- This is ok because we are going to  mod 2^61  it as a final step.
     -- In full precision math, Random_x here would return all nums in range
     -- 0 .. (p0*p1-1) exactly once during  period = p0*p1. For example,
     -- if x0 and x1 are produced by full period rearrangement generators:
     --    x0 in 0..9, with period p0 = 10, and
     --    x1 in 0..6, with period p1 = 7, then
     -- you'ld get z = x0 + p0*x1 in range 0 .. 69, no missing nums during
     -- period of 70. So z is produced with period p0*p1 = 70, and the new 
     -- generator is a full period rearrangement generator: each number in 
     -- the range 0..69 appears exactly once during each period of length 70.

  end Get_Random_MWC_64_2;

  pragma Inline (Get_Random_MWC_64_2);

  -----------------------------
  -- Get_Random_XOR_Shift_61 --
  -----------------------------

  -- Marsaglia's XOR-SHIFT generator mixes bits by the following recipe:
  --
  --  S2  := (S2 XOR (S2 * 2**20))mod 2**61; (a shift-left-by-20  step)
  --  S2  :=  S2 XOR (S2 / 2**32);           (a shift-right-by-32 step)
  --
  -- Marsaglia's Algorithm is used to find shifts that yield generators with
  -- full (prime) period of 2**61-1.  The following are full period shifts that
  -- get best scores on a stringent avalanche test.  (But you really need about
  -- 12 of these shifts to get best resemblance of true randomness. There are
  -- other deficiencies in the XOR-SHIFT generator that can never be removed
  -- by adding more xor-shifts, so we stop at 8 and use it as part of a 
  -- combination generator.)
  --
  -- best 8's:  1 3 6 13 21 29 31 12             57 59 2
  --            1 3 7 15 30 32 20 6              56 58 2 (preferred)
  --            1 3 7 16 21 30 26 6              55 55 0
  --            1 3 9 16 21 23 23 1              43 54 11 (looks promising)
  --            1 3 12 18 23 25 29 15            61 65 4
  --            1 3 13 24 31 35 8 1              63 53 10 
  -- 8 or less is probably appropriate for combination generator.


  procedure Get_Random_XOR_Shift_61
    (Random_x :    out Random_Int;
     S        : in out State)
  is
     X2 : Random_Int;

     S8 : constant := 6;  S7 : constant := 20; S6 : constant := 32; 
     S5 : constant := 30; S4 : constant := 15; S3 : constant := 7; 
     S2 : constant := 3;  S1 : constant := 1; 

     --  Error detection by assertion:

     pragma Assert 
      (S8=6 and S7=20 and S6=32 and S5=30 and S4=15 and S3=7 and S2=3 and S1=1);

     --  Error correction by inspection:

     --  (if mutated parameters are detected, use data given below to correct).
     --   1 3 7 15 30 32 20 6
     --   1 3 7 15 30 32 20 6
     --   1 3 7 15 30 32 20 6

  begin

     X2  := S.X(Xor_Shift_id);
     X2  :=  X2 XOR (X2 / 2**S8);
     X2  := (X2 XOR (X2 * 2**S7))mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S6); 
     X2  := (X2 XOR (X2 * 2**S5))mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S4);
     X2  := (X2 XOR (X2 * 2**S3))mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S2);
     X2  := (X2 XOR (X2 * 2**S1))mod 2**61;
     S.X(Xor_Shift_id) := X2;

     Random_x := X2;

  end Get_Random_XOR_Shift_61;

  pragma Inline (Get_Random_XOR_Shift_61);

  ----------------
  -- Get_Random --
  ----------------

  --  Get_Random is a combination generator: 3 component generators w, x, and y
  --  are combined to produce a generator z whose period is the product of
  --  the periods of w, x, and y.  z is generated by z = (w + x + y) mod p,
  --  where y is full period in range 0..p-1, and the periods of the 
  --  generators for w, x and y have no prime factors in common.
  --  Only one of the generators w, x, and y
  --  needs to be full-period to guarantee a full-period z. In the present 
  --  case, 2 of the 3 component generators are full-period.
  --  The third (the nonlinear generator) is optionally full-period.
  --  z is therefore uniform and full period (in the sense that every number
  --  in the range 0..p-1 appears the same number of times during z's period). 
  --  In the present case p = 2^61. (Actually y is in the range 1..p1, missing
  --  the 0, but the missing 0 has negligible affect.)

  procedure Get_Random
    (Random_x :    out Random_Int;
     S        : in out State) 
  is
     X1, X2, X3 : Parent_Random_Int;
  begin
     Get_Random_MWC_64_2 (X1, S);
     Get_Random_XOR_Shift_61 (X2, S);
     Get_Random_X_squared_Mod_P (X3, S);
   
     Random_x := Random_Int ((X1 + X2 + X3) mod 2**61);

  end Get_Random;


  -----------
  -- Reset --
  -----------

  -- Step 1. Guarantee that each unique and valid choice 
  -- of the Initiators (Initiator1, Initiator2, .. )
  -- produces a unique and valid initial state S: S.X(0), S.X(1), ...
  --
  -- Step 2. Guarantee that a change of 1 bit or more in any one Initiator
  -- will cause changes in all elements of the initial state S: S.X(0), S.X(1), ... 
  --
  -- Harder than I expected...
  --
  procedure Reset
    (S          :    out State;
     Initiator1 : in     Seed_Random_Int := 1111;
     Initiator2 : in     Seed_Random_Int := 2222;
     Initiator3 : in     Seed_Random_Int := 3333;
     Initiator4 : in     Seed_Random_Int := 4444)
  is
     
     Seeds : Vals := (Initiator1, Initiator2, Initiator3, Initiator4);

     --  variables for doing a transpose:

     Chomp, Seed, X2 : Parent_Random_Int;
     Chomp_Shift_Length : constant Integer := (1 + Bits_per_Random_Number / S.X'Length);

     --  the transpose requires:

     pragma Assert (State_Index'First = 0);
     pragma Assert (S.X'Length < 9);
     pragma Assert (Random_Int'Last = 2**Bits_per_Random_Number-1);
     pragma Assert (Bits_per_Random_Number = 61);

     use Multiply_with_Carry_Constants;

     --  constants for an 8 stage full-period XOR-Shift function:

     S8 : constant := 15; S7 : constant := 29; S6 : constant := 25;   
     S5 : constant := 23; S4 : constant := 18; S3 : constant := 12;
     S2 : constant := 3;  S1 : constant := 1; 

     pragma Assert 
      (S8=15 and S7=29 and S6=25 and S5=23 and S4=18 and S3=12 and S2=3 and S1=1);
     --  if mutated parameters are detected, use data given in Step 1 to correct.

  begin

     -- Want the Seeds in range, even if checks are off:
 
     for i in State_Index loop
        if not Seeds(i)'Valid then raise Constraint_Error; end if;
     end loop;
 
     -- Keep scrambling until S.X(2) /= 0:
 
     Get_non_Zero_State: 
     for Keep_Trying in 1..17 loop  -- astronomically unlikely this will take >2 tries
 
        -- Step 1:
        -- Start by transforming each seed with a function 
        -- that is a 1-1 mapping between all elements of 0..2**61-1.
        -- (0 maps to 0, but is weeded out).
        -- 
        -- period = 2**61-1: maps 0 to 0 but thats ok
        --                   1 3 12 18 23 25 29 15 
        --                   1 3 12 18 23 25 29 15 
    
        for i in State_Index loop   
           X2 := Seeds(i);
           X2  :=  X2 XOR (X2 / 2**S8);
           X2  := (X2 XOR (X2 * 2**S7))mod 2**61;
           X2  :=  X2 XOR (X2 / 2**S6); 
           X2  := (X2 XOR (X2 * 2**S5))mod 2**61;
           X2  :=  X2 XOR (X2 / 2**S4);
           X2  := (X2 XOR (X2 * 2**S3))mod 2**61;
           X2  :=  X2 XOR (X2 / 2**S2);
           X2  := (X2 XOR (X2 * 2**S1))mod 2**61;
           Seeds(i) := X2;
        end loop;
 
 
        -- Step 2: The Transpose.
        -- Make an N x N matrix out of the N elements of array Seed
        -- (by breaking each element into N Parts), and transpose it;
        -- write it to state S.X. Here N = 4, and Chomp_Shift_Length = 16.
 
        S.X := (others => 0);
 
        for j in State_Index loop 
 
           Seed := Seeds(j);

           --  Each of the X's gets a 16-bit bite of Seed, placed at 2**(j*16): 

           for i in State_Index loop 
              Chomp  := Seed MOD 2**Chomp_Shift_Length;
              S.X(i) := S.X(i) + Chomp * 2**(j*Chomp_Shift_Length);
              Seed   := Seed / 2**Chomp_Shift_Length;
           end loop;
 
        end loop;
 
        -- Scramble the state S.X again:
 
        for i in State_Index loop   
           X2 := S.X(i);
           X2  :=  X2 XOR (X2 / 2**S8);
           X2  := (X2 XOR (X2 * 2**S7))mod 2**61;
           X2  :=  X2 XOR (X2 / 2**S6); 
           X2  := (X2 XOR (X2 * 2**S5))mod 2**61;
           X2  :=  X2 XOR (X2 / 2**S4);
           X2  := (X2 XOR (X2 * 2**S3))mod 2**61;
           X2  :=  X2 XOR (X2 / 2**S2);
           X2  := (X2 XOR (X2 * 2**S1))mod 2**61;
           S.X(i) := X2;
        end loop;
 
        if Keep_Trying > 1 and then S.X(Xor_Shift_id) /= 0 then   
           exit Get_non_Zero_State;
        else
           Seeds := S.X; -- make new seeds and scramble again.
        end if;

     end loop Get_non_Zero_State;


     -- Step 3: Error correction.
     --
     -- Weed out Seed values that give a period of 1 to component generators.

     for i in State_Index range MWC_id_0 .. MWC_id_1 loop
     if S.X(i) = 0 then                 -- Period=1
        S.X(i) := Random_Int'Last + 1;  -- Is valid init as asserted above
     end if;
     end loop;

     if S.X(MWC_id_0) = (a0 * b0 - 1) then -- Period=1
        raise Constraint_Error;  -- can't happen, if assertion governing a0 is enforced
     end if;

     if S.X(MWC_id_1) = (a1 * b1 - 1) then -- Period=1
        raise Constraint_Error;  -- can't happen, if assertion governing a1 is enforced
     end if;

     --  xor/shift: already done. Weeded out the 0's above.
     --if S.X(Xor_Shift_id) = 0 then
        --S.X(Xor_Shift_id) := Random_Int'Last;
     --end if;

     declare 
        q : constant State_Index := Quadratic_id;
     begin
        if S.X(q) = 0  or  S.X(q) = 1 then           -- Period=1 for x*x mod p
           S.X(q) := Random_Int'Last + 1 + S.X(q);   -- Is valid init as asserted above
        end if;
     end;

  end Reset;

  -----------
  -- Value --
  -----------

  function Value (Coded_State : State_String) return State is
     S  : State;
     Start_Of_Seed : Integer := 1;

     use Multiply_with_Carry_Constants; -- for a0, a1, b

  begin

     for i in State_Index loop
        S.X(i) := Random_Int'Value (
            Coded_State(Start_Of_Seed .. Start_Of_Seed+Seeds_Image_Width-1));
        Start_Of_Seed := Start_Of_Seed + Seeds_Image_Width;
     end loop;

     -- Finally, invalid state vals ought to raise Constraint_Error (should be 
     -- impossible and therefore something went unacceptibly wrong.):

     if (S.X(MWC_id_0) = 0) or (S.X(MWC_id_0) = (a0*b0-1)) then -- Period=1
        raise Constraint_Error;
     end if;
     if (S.X(MWC_id_1) = 0) or (S.X(MWC_id_1) = (a1*b1-1)) then -- Period=1
        raise Constraint_Error;
     end if;

     --  Period=1 of XOR_shift component generator:
     if S.X(Xor_Shift_id) = 0  or not (S.X(Xor_Shift_id) in Random_Int) then 
       raise Constraint_Error; 
     end if;

     --  Period=1 of X*X mod N component generator:
     if S.X(Quadratic_id) = 0  or S.X(Quadratic_id) = 1 then 
       raise Constraint_Error; 
     end if;

     return S;

  end Value;

  -----------
  -- Image --
  -----------

  function Image (Of_State : State) return State_String is
     Result : State_String := (others => ' ');
     Start_Of_Seed : Integer := 1;
  begin
     for i in State_Index loop
        declare
           Y : constant String := Random_Int'Image (Of_State.X (i));
        begin
           Result(Start_Of_Seed .. Start_Of_Seed+Y'length-1) := Y;
           Start_Of_Seed := Start_Of_Seed + Seeds_Image_Width;
        end;
     end loop;
     return Result;
  end Image;

end Disorderly.Random;

