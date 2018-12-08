
-----------------------------------------------------------------------
-- package body Extended_Real.E_Rand, extended precision random numbers.
-- Copyright (C) 2008-2018 Jonathan S. Parker
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
---------------------------------------------------------------------------


package body Extended_Real.E_Rand is

  State : Random_Int := 701;

  ----------------
  -- Get_Random --
  ----------------

  -- 61 bit rands:

  function Next_Random_Int 
     return Random_Int
  is
     X2 : Random_Int;

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

     X2  := State;
     X2  :=  X2 XOR (X2 / 2**S8);
     X2  := (X2 XOR (X2 * 2**S7))mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S6); 
     X2  := (X2 XOR (X2 * 2**S5))mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S4);
     X2  := (X2 XOR (X2 * 2**S3))mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S2);
     X2  := (X2 XOR (X2 * 2**S1))mod 2**61;
     State := X2;

     return X2;

  end Next_Random_Int;

  -----------
  -- Reset --
  -----------
  
  -- Initiator and state must never be negative or 0!

  procedure Reset 
    (Initiator : in Positive := 7777777)
  is
     X : Integer := Initiator mod 2**30; 
     -- if Ints are 64 bit, keep it under 61 bit, while still portable to 32 bit int.
  begin
     if X = 0 then X := 1; end if;
     State := Random_Int (X);
  end Reset;
  
  ------------
  -- Random --
  ------------

  function Random 
     return E_Real  
  is
     Result : E_Real;
  begin

     for I in Digit_Index loop
        Result.Digit(I) := Digit_Type (Next_Random_Int mod 2**No_Of_Bits_in_Radix);
     end loop;

     if Result.Digit(Digit_Index'First) = Digit_Zero then
        Result.Digit(Digit_Index'First) := Digit_One;
     end if;
   
     Result.Is_Zero     := False;
     Result.Is_Positive := True;
     Result.Is_Infinite := False;
     Result.Exp         := -1;

     return Result;  

  end Random;
 
end Extended_Real.E_Rand;
