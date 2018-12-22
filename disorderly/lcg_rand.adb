
-------------------------------------------------------------------------------
-- package body LCG_Rand, Linear Congruential Generator
-- Copyright (C) 2018 Jonathan S. Parker
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

package body LCG_Rand
   with Spark_Mode => On
is

  -- Marsaglia-Zaman multiply-with-carry algorithm is used to construct
  -- two LCG's (linear congruential generators) with prime number moduli
  -- m0 and m1, and prime number periods (m0-1)/2 and (m1-1)/2. 

  -------------------------
  -- Get_Random_LCG_64_0 --
  -------------------------

  -- a0 = 2731662333
  -- b0 = 2^32
  -- m0 = a0 * b0 - 1  
  -- m0 = 11732400383950061567
  -- result should be
  --    x  = (x*a0) mod m0 
  -- but (x0*a0) overflows 64 bit computation, so tests need 128 bit ints.

  procedure Get_Random_LCG_64_0 (X0 : in out Parent_Random_Int) is
     S_lo, S_hi : Parent_Random_Int;
  begin

     S_lo   := X0 mod b0;
     S_hi   := X0  /  b0;
     X0     := a0 * S_lo + S_hi;
     --  X0_final is in 1 .. m0-1 if x0_initial is in 1 .. m0-1.
     --  m0 is prime.
     --  Period = p0 = (m0-1) / 2. Period is prime.

  end Get_Random_LCG_64_0;

  -- a1 = 1688234283
  -- b1 = 2^32
  -- m1 = a1 * b1 - 1  
  -- m1 = 7250911033471008767
  -- result should be 
  --    x1 := (x1*a1) mod m1,
  -- but (x1*a1) overflows 64 bit computation, so tests need 128 bit ints.

  procedure Get_Random_LCG_64_1 (X1 : in out Parent_Random_Int) is
     S_lo, S_hi : Parent_Random_Int;
  begin
     S_lo   := X1 mod b1;
     S_hi   := X1  /  b1;
     X1     := a1 * S_lo + S_hi;
     --  X1 is in  1 .. m1-1. 
     --  m1 is prime.
     --  Period = p1 = (m1-1) / 2. Period is prime.

  end Get_Random_LCG_64_1;


  procedure Get_Random_LCG_64_Combined
    (S0       : in out Parent_Random_Int;
     S1       : in out Parent_Random_Int;
     Random_x :    out Parent_Random_Int) -- result
  is
     x0 : Parent_Random_Int; 
     -- verify that summing x0+s1 doesn't go out-of-range of 2^64: p0+m1<2^64
     pragma Assert (p0 + m1 <= Parent_Random_Int'Last);
  begin

     Get_Random_LCG_64_0 (S0); -- this updates state S0
     --  x0 = S0 is in 1 .. m0-1, with prime number period p0 = (m0-1)/2

     x0 := S0;
     if x0 > p0 then 
        x0 := (m0 - 1) - x0;
     else
        x0 := x0 - 1;
     end if;
     --  Now x0 is in 0 .. p0-1, with prime number period p0 = (m0-1)/2

     Get_Random_LCG_64_1 (S1); -- this updates state S1
     --  S1 is in 1 .. m1-1, with prime number period p1 = (m1-1)/2

     Random_x := (x0 + S1) mod p0;
     -- Stnd combination generator. Only 1 of the 2 component generators (the 0th) 
     -- needs to be full period. If that's the case, then the combination
     -- generator returns numbers uniformly distributed in 1 .. p0-1, with no
     -- subcycles shorter than p0*p1; in other words it's a full period generator.
     --
     -- Random_LCG_64_0, as modified above, is a rearrangement generator. It takes
     -- each of the numbers in the range 0 .. p0-1, and returns them out of order.
     -- During a cycle of length p0, it returns each of the numbers exactly
     -- once. (Random_LCG_64_1 has a similar property: it takes 1/2 the numbers
     -- in 0 .. m1-1, and returns each of them exactly once during a period of
     -- length (m1-1)/2.)
     --
     -- Let's say the 1st generator has period p0, and returns each number in 0..p0-1
     -- exactly once during that period. If you then add the results of the 1st generator
     -- to those of another generator modulo p0, and if the 2nd generator's period, p1,
     -- has no prime factors in common with p0, then the resulting combination
     -- generator is full period, with period p0*p1. The combination generator returns
     -- each number in 0..p0-1 exactly p1 times during its period of length p0*p1.

  end Get_Random_LCG_64_Combined;
   
end LCG_Rand;

