
-------------------------------------------------------------------------------
-- package body Disorderly.Basic_Rand, Linear Random Number Generator
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

with text_utilities;

package body Disorderly.Basic_Rand
   with Spark_Mode => On
is

  package text_handling is new text_utilities (No_Of_Seeds, Parent_Random_Int);
  use text_handling;

  --  The state (composed of 3 x 64 bits) is contained in an array of 3
  --  64 bit ints indexed by the following 3 ints:

  MWC_id_0     : constant State_Index := 0;
  MWC_id_1     : constant State_Index := 1;
  Xor_Shift_id : constant State_Index := 2;

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

     --  If a single parameter mutates, the generator fails, so assert twice.
     --
     b0 : constant := 2**32;
     a0 : constant := 2731662333; -- a0 must be the larger a
     pragma Assert ((a0 = 2731662333) and (b0 = 2**32));
     -- modulus not period:  m0 = 11732400383950061567  a0 = 2731662333

     b1 : constant := 2**32;
     a1 : constant := 1688234283;
     pragma Assert ((a1 = 1688234283) and (b1 = 2**32));
     -- modulus not period:  m1 = a1*2^32-1 = 7250911033471008767   a1 = 1688234283

     --  Ratio of periods = golden mean, just for fun.
     --  Period: (7250911033471008767*11732400383950061567)/4 = 2^124 * 0.99999999604
     --
     --b0 : constant := 2**32-5;
     --a0 : constant := 2731663104; -- a0 must be the larger a
     --pragma Assert ((a0 = 2731663104) and (a0 = 2731663104) and (a0 = 2731663104));
     --pragma Assert ((b0 = 2**32-5) and (b0 = 2**32-5) and (b0 = 2**32-5));
     -- m0 = a0*(2^32-5)-1 =  11732403681711531263

     m0 : constant :=  a0 * b0 - 1;  -- LCG modulus (safe prime)
     m1 : constant :=  a1 * b1 - 1;  -- LCG modulus (safe prime)

     --p0 : constant := (m0 - 1) / 2; -- period (prime)
     --p1 : constant := (m1 - 1) / 2; -- period (prime)

     pragma Assert (a0 < b0);
     pragma Assert (a1 < b1);
     pragma Assert (a1 < a0);
     pragma Assert ((a1/2)*2 < a1);  -- but only for b1 = 2^32

     pragma Assert (m0 < Parent_Random_Int'Last);
     pragma Assert (m1 < Parent_Random_Int'Last);

     pragma Assert (m0 > Random_Int'Last + 1);  -- for Reset to work.
     pragma Assert (m1 > Random_Int'Last + 1);  -- ie m1 > 2**61

  end Multiply_with_Carry_Constants;

  use Multiply_with_Carry_Constants; -- for m0, m1

  subtype Valid_MWC_0_Range is Parent_Random_Int range 1 .. m0-1;
  subtype Valid_MWC_1_Range is Parent_Random_Int range 1 .. m1-1;
  subtype Valid_Xor_Shift_Range is Random_Int;

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

--  procedure Get_Random_MWC_64_2
--    (Random_x :    out Parent_Random_Int;
--     S        : in out State)
--  is
--     S_lo, S_hi : Parent_Random_Int;
--     S0, S1     : Parent_Random_Int;
--
--     use Multiply_with_Carry_Constants;
--
--  begin
--
--     S_hi   := S.X(MWC_id_0)  /  b0;
--     S_lo   := S.X(MWC_id_0) mod b0;
--     S0     := a0 * S_lo + S_hi;
--
--     S.X(MWC_id_0) := S0;
--
--
--     --  S0 is in  1 .. m0-1.   m0-1 equals 2*p0, where p0 is prime.
--
--     if S0 > p0 then
--       S0 := m0 - S0 - 1;
--     else
--       S0 := S0 - 1;
--     end if;
--     --  Each number in range  0 .. p0-1  appears exactly once during period p0.
--     --  So S0 is full-period with prime number period p0.
--     --  (If x=p+1,  m - x = 2p+1  - (p+1)  = p, so m-x-1 = p-1.)
--
--
--     S_hi   := S.X(MWC_id_1)  /  b1;
--     S_lo   := S.X(MWC_id_1) mod b1;
--     S1     := a1 * S_lo + S_hi;
--
--     S.X(MWC_id_1) := S1;
--
--     if S1 > p1 then
--       S1 := m1 - S1 - 1;
--     else
--       S1 := S1 - 1;
--     end if;
--     --  Now S1 is in  0 .. p1-1
--
--     Random_x := S0 + p0 * S1;
--
--     -- Requires modular arithmetic, 2^n = 2^64. Requires p0 > p1.
--     -- The above is equivalent to doing the arithmetic in full precision,
--     -- followed by mod 2^64, because (with for example m = 2^64):
--     --   (a*b + c) mod m = ((a mod m) * (b mod m) + c mod m) mod m.
--     -- This is ok because we are going to  mod 2^61  it as a final step.
--     -- In full precision math, Random_x here would return all nums in range
--     -- 0 .. (p0*p1-1) exactly once during  period = p0*p1. For example,
--     -- if x0 and x1 are produced by full period rearrangement generators:
--     --    x0 in 0..9, with period p0 = 10, and
--     --    x1 in 0..6, with period p1 = 7, then
--     -- you'ld get z = x0 + p0*x1 in range 0 .. 69, no missing nums during
--     -- period of 70. So z is produced with period p0*p1 = 70, and the new
--     -- generator is a full period rearrangement generator: each number in
--     -- the range 0..69 appears exactly once during each period of length 70.
--
--  end Get_Random_MWC_64_2;
--
--  pragma Inline (Get_Random_MWC_64_2);

  -------------------------
  -- Get_Random_MWC_64_x --
  -------------------------

  -- Stnd Zaman Marsalglia multiply-with-carry LCG: just use the
  -- lower 30 or 31 bits of 2 generators and concatenate to get greater
  -- than 61 bits; there's nothing wrong with this - passes all the
  -- usual tests. The bits are nearly uniform because we only use the
  -- lower 31 or 30 bits of 63 or 62 bit nums.

  procedure Get_Random_MWC_64_x
    (Random_x :    out Parent_Random_Int;
     S        : in out State)
  is
     S_lo, S_hi : Parent_Random_Int;
     S0, S1     : Parent_Random_Int;
  begin

     S_hi   := S.X(MWC_id_0)  /  b0; -- b0, b1 from Multiply_with_Carry_Constants
     S_lo   := S.X(MWC_id_0) mod b0;
     S0     := a0 * S_lo + S_hi;

     S.X(MWC_id_0) := S0;

     S_hi   := S.X(MWC_id_1)  /  b1;
     S_lo   := S.X(MWC_id_1) mod b1;
     S1     := a1 * S_lo + S_hi;

     S.X(MWC_id_1) := S1;

     -- concatenate 31 bits from S0, + 30 from S1 on range 0..x (shifted 31 bits):

     Random_x := (S0 mod 2**31) + (S1 - 1) * 2**31;    -- type Parent_Random_Int

  end Get_Random_MWC_64_x;

  pragma Inline (Get_Random_MWC_64_x);

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
     X2 : Parent_Random_Int;

     S8 : constant := 6;  S7 : constant := 20; S6 : constant := 32;
     S5 : constant := 30; S4 : constant := 15; S3 : constant := 7;
     S2 : constant := 3;  S1 : constant := 1;

     --  Error detection is by assertion:

     pragma Assert
      (S8=6 and S7=20 and S6=32 and S5=30 and S4=15 and S3=7 and S2=3 and S1=1);

     --  Error correction is by inspection:
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

     Random_x := X2 mod 2**61;

     S.X(Xor_Shift_id) := Random_x;

  end Get_Random_XOR_Shift_61;

  ----------------
  -- Get_Random --
  ----------------

  --  Get_Random is a combination generator: 2 component generators x, and y
  --  are combined to produce a generator z whose period is the product of
  --  the periods of x, and y.  z is generated by z = (x + y) mod p,
  --  where y is full period in range 0..p-1, and the periods of the
  --  generators for x and y have no prime factors in common.
  --  Only one of the two generators x and y
  --  needs to be full-period to guarantee a full-period z.

  procedure Get_Random
    (Random_x :    out Random_Int;
     S        : in out State)
  is
     X1, X2 : Parent_Random_Int;
  begin
   --Get_Random_MWC_64_2 (X1, S); 
   -- Slower but higher degree of uniformity than below; all the bits in the final
   -- result will be full-period. The alternative below results in a generator
   -- in which the bits have 1/2 the period of the full generator.

     Get_Random_MWC_64_x (X1, S); -- faster.
     Get_Random_XOR_Shift_61 (X2, S);

     Random_x := Random_Int ((X1 + X2) mod 2**61);

  end Get_Random;


  function Valid_State (S : in State) return Boolean is
     Result : constant Boolean :=
       (S.X(MWC_id_0)     in Valid_MWC_0_Range) and
       (S.X(MWC_id_1)     in Valid_MWC_1_Range) and
       (S.X(Xor_Shift_id) in Valid_Xor_Shift_Range);
  begin
     return Result;
  end Valid_State;

  procedure Correct (S : in out State) is
  begin
     -- order important
     if S.X(MWC_id_0) > Valid_MWC_0_Range'Last then -- Period=1 if equality holds
        S.X(MWC_id_0) := S.X(MWC_id_0) mod (Valid_MWC_0_Range'Last+1);
     end if;

     if S.X(MWC_id_0) = 0 then -- Period=1
        S.X(MWC_id_0) := Valid_MWC_0_Range'Last;
     end if;

     if S.X(MWC_id_1) > Valid_MWC_1_Range'Last then -- Period=1 if equality holds 
        S.X(MWC_id_1) := S.X(MWC_id_1) mod (Valid_MWC_1_Range'Last+1);
     end if;

     if S.X(MWC_id_1) = 0 then -- Period=1
        S.X(MWC_id_1) := Valid_MWC_1_Range'Last;
     end if;

     if S.X(Xor_Shift_id) > Valid_Xor_Shift_Range'Last then
        S.X(Xor_Shift_id) := S.X(Xor_Shift_id) mod (Valid_Xor_Shift_Range'Last+1);
     end if;

     if S.X(Xor_Shift_id) = 0 then -- Period=1
        S.X(Xor_Shift_id) := Valid_Xor_Shift_Range'Last;
     end if;

  end Correct;

  procedure Get_Random_XOR_Shift_61_b (Random_x : in out Parent_Random_Int) is
     X2 : Parent_Random_Int := Random_x;

     S8 : constant := 15; S7 : constant := 29; S6 : constant := 25;
     S5 : constant := 23; S4 : constant := 18; S3 : constant := 12;
     S2 : constant := 3;  S1 : constant := 1;
     pragma Assert
      (S8=15 and S7=29 and S6=25 and S5=23 and S4=18 and S3=12 and S2=3 and S1=1);
     --  if mutated parameters are detected, use these for error correction.
     --     1 3 12 18 23 25 29 15
     --     1 3 12 18 23 25 29 15
  begin

     X2  :=  X2 XOR (X2 / 2**S8);
     X2  := (X2 XOR (X2 * 2**S7))mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S6);
     X2  := (X2 XOR (X2 * 2**S5))mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S4);
     X2  := (X2 XOR (X2 * 2**S3))mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S2);
     X2  := (X2 XOR (X2 * 2**S1))mod 2**61;
     
     Random_x := X2;

  end Get_Random_XOR_Shift_61_b;

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
  procedure Reset
    (S          :    out State;
     Initiator1 : in     Seed_Random_Int := 1111;
     Initiator2 : in     Seed_Random_Int := 2222;
     Initiator3 : in     Seed_Random_Int := 3333)
  is
     Seeds : Vals := (Initiator1, Initiator2, Initiator3);

     Cut, Seed : Parent_Random_Int;
     Cut_Shift_Length : constant Integer := (1 + Bits_per_Random_Number / S.X'Length);

     --  the transpose requires:

     pragma Assert (State_Index'First = 0);
     pragma Assert (S.X'Length < 9);
     pragma Assert (Bits_per_Random_Number = 61);
     pragma Assert (Random_Int'Last = 2**Bits_per_Random_Number-1);
  begin

     for Keep_Trying in 1 .. 17 loop  -- arbitrary loop limit. More is good.

        -- Step 1:
        -- Start by transforming each seed with a function
        -- that is a 1-1 mapping between all elements of 0..2**61-1.
        -- (0 maps to 0, but is weeded out).

        for i in State_Index loop
           Get_Random_XOR_Shift_61_b (Seeds(i));
        end loop;

        -- Step 2:
        -- Make an N x N matrix out of the N elements of array Seed
        -- (by breaking each element into N Parts), and transpose it;
        -- write it to to the state S.X:

        S.X := (others => 0);

        for j in State_Index loop

           Seed := Seeds(j);

           for i in State_Index loop
              Cut    := Seed mod 2**Cut_Shift_Length;
              S.X(i) := S.X(i) + Cut * 2**(j*Cut_Shift_Length);
              Seed   := Seed / 2**Cut_Shift_Length;
           end loop;

        end loop;

        -- Scramble the state S.X again:

        for i in State_Index loop
           Get_Random_XOR_Shift_61_b (S.X(i));
        end loop;

        Seeds := S.X; -- use the new seeds and scramble again.

     end loop;

     -- Notice each S.X is in 0..2^61-1, so have S.X < a0 * b0 - 1 etc.
     -- Still may have S.X = 0, which gives period of 1 to component generators.

     -- Step 3: Error correction.

     if not Valid_State (S) then
        Correct (S);
     end if;

  end Reset;

  -----------
  -- Value --
  -----------

  function Value (Coded_State : in State_String) return State is
     S  : State;
     Seed_1_1st : constant Positive := Coded_State'First;
     Seed_j_1st : Positive;
  begin

     for j in State_Index loop
        Seed_j_1st  := Seed_1_1st  + (j-State_Index'First) * Rand_Image_Width;
        S.X(j) := Value (Coded_State(Seed_j_1st .. Seed_j_1st + Rand_Image_Width - 1));
     end loop;

     return S;

  end Value;

  -----------
  -- Image --
  -----------

  function Image (Of_State : State) return State_String is
     Result : State_String := (others => '0');
     Y : Random_Int_String;
     Seed_1_1st : constant Positive := Result'First;
     Seed_j_1st : Positive;
  begin

     for j in State_Index loop
        Y := Image (Of_State.X (j));
        Seed_j_1st  := Seed_1_1st  + (j-State_Index'First) * Rand_Image_Width;
        Result(Seed_j_1st .. Seed_j_1st+Rand_Image_Width-1) := Y;
     end loop;

     return Result;
  end Image;

  -- Puts a ' ' in front of each 20 digit number.

  function Formatted_Image_Of (Of_State : State_String) return Formatted_State_String is
     Result : Formatted_State_String := (others => ' ');
     Seed_1_1st : constant Positive := Result'First;
     Seed_y_1st, Seed_j_1st : Positive;
     Y : Random_Int_String;
  begin
     for j in State_Index loop
        Seed_j_1st := Seed_1_1st + (j-State_Index'First) * (Rand_Image_Width + Leading_Spaces);
        Seed_y_1st := Seed_1_1st + (j-State_Index'First) *  Rand_Image_Width;
        Y := Of_State(Seed_y_1st .. Seed_y_1st+Rand_Image_Width-1);
        Result(Seed_j_1st .. Seed_j_1st+Rand_Image_Width-1) := Y;
     end loop;

     return Result;
  end Formatted_Image_Of;

  function Formatted_Image (Of_State : State) return Formatted_State_String is
  begin
     return Formatted_Image_Of (Image (Of_State));
  end Formatted_Image;

end Disorderly.Basic_Rand;

