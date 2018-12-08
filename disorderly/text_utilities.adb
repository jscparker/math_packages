
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

package body Text_Utilities
   with Spark_Mode => On
is

  subtype Integer_Digit is Parent_Random_Int range 0 .. 9;

  Digit_Image : constant array(Integer_Digit) of Character_Digit :=
     ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9');

  Digit_Value : constant array(Character_Digit) of Integer_Digit :=
     (0, 1, 2, 3, 4, 5, 6, 7, 8, 9);

  subtype Log_Range is Natural range 0 .. 19;

  Ten_to_the : constant array(Log_Range) of Parent_Random_Int :=
    (10**00, 10**01, 10**02, 10**03, 10**04, 10**05, 10**06, 10**07, 10**08, 10**09,
     10**10, 10**11, 10**12, 10**13, 10**14, 10**15, 10**16, 10**17, 10**18, 10**19);

  -----------
  -- Value --
  -----------

  -- Random_Int_String contains a 20 digit unsigned integer.
  -- Every Parent_Random_Int (0 .. 2**64-1) can be represented in 20 decimal digits.
  --
  -- All state values satisfy State.X(i) < 2^64.
  --
  -- The calculation uses the modular arithmetic of modular type Parent_Random_Int. 
  -- If the image string encodes a number > 2^64-1, then of course the calculation
  -- won't return that value, but all State.X(i) values are in a range that 
  -- is correctly evaluated by this function.
  --
  -- For more general use, one may want to check if the string is out-of-range
  -- of Parent_Random_Int (i.e. > 2^64-1). To do that, sum the least significant 19
  -- digits first. Then check if the most significant digit is '0', '1', or greater.
  -- If greater the string is out of range, if '1' it might be out of range, and if
  -- '0' then it's in range.

  function Value (Random_Int_Image : in Random_Int_String) return Parent_Random_Int
  is
     Val : Parent_Random_Int := 0;
  begin

     for j in Random_Int_String_Index loop
       Val := Val +
          Digit_Value (Random_Int_Image(j))*Ten_to_the(Random_Int_String_Index'Last - j);
     end loop;
     return Val;

  end Value;

  -----------
  -- Image --
  -----------

  -- Every 64 bit Parent_Random_Int can be represented by a 20 decimal digit 
  -- number.

  function Image (X : in Parent_Random_Int) return Random_Int_String
  is
     Ten : constant Parent_Random_Int := 10;
     Y : Parent_Random_Int := X;
     Result : Random_Int_String;
  begin

     for j in reverse Random_Int_String_Index loop
        Result(j) := Digit_Image (Y mod Ten); 
        Y := Y / Ten;
     end loop;

     return Result;

  end Image;

end Text_Utilities;
