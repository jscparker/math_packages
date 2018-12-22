
-------------------------------------------------------------------------------
-- package Disorderly.Random, Non-linear Random Number Generator
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

-- PACKAGE Disorderly.Random
--
-- Procedure  Disorderly.Random.Get_Random  is a non-linear pseudo random number 
-- generator designed to be task-safe, and optimized for high statistical quality. 
--
--  1. Uniform output in the range 0 .. 2**61-1.  
--     (ie, 61 random bits per call.)  
--     (2**61-1 is a Mersenne Prime; explanation below.) 
--  2. Period > 2^246. (Actually, period is about 1.813381 * 2^246.)
--  3. Period is the product of 4 large (19 decimal digit) prime numbers:
--     4192767015302299403 * 3625455516735504383 * 5866200191975030783 * (2^61-1).
--  4. Generator is non-linear.
--  5. Generator is full-period, and 2 of the 3 component generators are
--     full-period. (The 3rd, the non-linear one, is optionally full-period.)
--  6. Generator is pure (the package is stateless) for convenient use
--     in multi-tasking simulations.
--  7. CPU time per call is constant (again for multi-tasking simulations).
--  8. Size of state per generator is 4 x 64 bits.
--  9. Speed:  Typical benchmarks (on a 64-bit intel PC) measure CPU time per
--     call at around 1/2 to 2/3 that of a call to a 64-bit Sin(x).
--  10. The State can be initialized by a version of procedure Reset that calls 
--      Ada.Calendar.Clock. This version of Reset is confined to a child package,
--      since its not pure. The child package is called Random.Clock_Entropy.
--
-- Items 1-5 are general characteristics of good generators.
-- Items 6-8 are attributes that are desirable in a language with built-in 
-- concurrency.
--
-- For a general RATIONALE, see the parent package:   Disorderly.ads
--

package Disorderly.Random
   with Spark_Mode => On
is

   pragma Pure (Disorderly.Random);

   type Parent_Random_Int is mod 2**64;
   --  Internally, all arithmetic is done on this type.  The generator 
   --  is designed for machines with efficient 64-bit integer arithmetic.

   Bits_per_Random_Number : constant := 61;

   subtype Random_Int is Parent_Random_Int range 0 .. 2**Bits_per_Random_Number-1;
   --  The random number generator returns Ints of this type.  Its
   --  best to think of these Ints as 61 bins, each containing a random bit.
   --  Notice that it is *not* a mod 2**61 type.
   --
   --  The number p = 61 is a special kind of prime that makes 2**p-1 a 
   --  prime number. Primality of 2**p-1 is a basic requirement of the 
   --  algorithm that generates the random numbers.
   --
   --  The number of bits returned by a RNG is rarely the actual
   --  number required by a particular application. Additional
   --  work is necessary to put the random stream into the desired range
   --  and numeric type.

   subtype Seed_Random_Int is Random_Int range 1 .. Random_Int'Last;
   --  Seeds (Initiators) must never be 0.

   type State is private;

   procedure Get_Random
     (Random_x :    out Random_Int;
      S        : in out State);

   procedure Reset 
     (S          :    out State;
      Initiator1 : in     Seed_Random_Int := 1111;
      Initiator2 : in     Seed_Random_Int := 2222;
      Initiator3 : in     Seed_Random_Int := 3333;
      Initiator4 : in     Seed_Random_Int := 4444);
   --  procedure Reset initializes State S. 
   --  There are no calls to Calendar, so it is reproducible: use the same 
   --  seeds each time and you get the same state S each time.


   --  The following routines translate state S into text, and back again. 
   --  Useful for saving state S on HD in text format.

   No_Of_Seeds       : constant := 4;
   Rand_Image_Width  : constant := 20; -- 2^64-1 = 1.8446744073709551615 * 10^19
   Max_Image_Width   : constant := Rand_Image_Width * No_Of_Seeds;

   subtype State_String is String(1 .. Max_Image_Width);

   function Value (Coded_State : State_String) return State with
      Pre => (for all i in State_String'Range => Coded_State(i) in '0' .. '9');

   function Image (Of_State : State) return State_String;

   --  Detect invalid States, for example after input from HDD.

   function Valid_State (S : in State) return Boolean;

   --  Make an easier to read version of State_String by putting a space in front
   --  of each of the 20 digit numbers: Formatted_Image().

   Leading_Spaces               : constant := 1;
   Formatted_State_String_Width : constant := (Rand_Image_Width + Leading_Spaces)*No_Of_Seeds;
   subtype Formatted_State_String is String(1 .. Formatted_State_String_Width);

   function Formatted_Image (Of_State : State) return Formatted_State_String;

   function Are_Equal (State_1, State_2 : State) return Boolean;

private

   subtype State_Index is Integer range 0..No_Of_Seeds-1;

   type Vals is array(State_Index) of Parent_Random_Int;

   type State is record
      X : Vals := (others => 701);
   end Record;

end Disorderly.Random;

