
-------------------------------------------------------------------------------
-- package body Xor_Shift, xor shift generator
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

package body Xor_Shift
   with Spark_Mode => On
is

  -- Marsaglia's XOR-SHIFT generator mixes bits by the following recipe:
  --
  --  S2  := (S2 XOR (S2 * 2**20))mod 2**61; (a shift-left-by-20  step)
  --  S2  :=  S2 XOR (S2 / 2**32);           (a shift-right-by-32 step)
  --
  -- Below are some shifts that produce full period generators on 1 .. 2**61-1.
  -- 2**61-1 is a Mersenne prime. The following are the best full period
  -- shifts I could find (based on scores from the avalanche test).
  --
  -- best 8's:  1 3 6 13 21 29 31 12             57 59 2
  --            1 3 7 15 30 32 20 6              56 58 2 (preferred)
  --            1 3 7 16 21 30 26 6              55 55 0
  --            1 3 9 16 21 23 23 1              43 54 11 (looks promising)
  --            1 3 12 18 23 25 29 15            61 65 4
  --            1 3 13 24 31 35 8 1              63 53 10 
  --

  -----------------------------
  -- Get_Random_XOR_Shift_61 --
  -----------------------------

  procedure Get_Random_XOR_Shift_61 (Random_x : in out Parent_Random_Int) is
     X2 : Parent_Random_Int := Random_x;

     S8 : constant := 6;  S7 : constant := 20; S6 : constant := 32; 
     S5 : constant := 30; S4 : constant := 15; S3 : constant := 7; 
     S2 : constant := 3;  S1 : constant := 1; 

     pragma Assert 
      (S8=6 and S7=20 and S6=32 and S5=30 and S4=15 and S3=7 and S2=3 and S1=1);
     --  if mutated parameters are detected, use the following for error correction.
     --   1 3 7 15 30 32 20 6
     --   1 3 7 15 30 32 20 6

  begin
     X2  :=  X2 XOR (X2 / 2**S8);
     X2  := (X2 XOR (X2 * 2**S7)) mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S6); 
     X2  := (X2 XOR (X2 * 2**S5)) mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S4);
     X2  := (X2 XOR (X2 * 2**S3)) mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S2);
     X2  := (X2 XOR (X2 * 2**S1)) mod 2**61;

     Random_x := X2;

  end Get_Random_XOR_Shift_61;


  procedure Get_Random_XOR_Shift_61_b (Random_x : in out Parent_Random_Int) is
     X2 : Parent_Random_Int := Random_x;

     S8 : constant := 15; S7 : constant := 29; S6 : constant := 25;
     S5 : constant := 23; S4 : constant := 18; S3 : constant := 12;
     S2 : constant := 3;  S1 : constant := 1;

     pragma Assert
      (S8=15 and S7=29 and S6=25 and S5=23 and S4=18 and S3=12 and S2=3 and S1=1);
     --  if mutated parameters are detected, use the following for error correction.
     --     1 3 12 18 23 25 29 15
     --     1 3 12 18 23 25 29 15

  begin

     X2  :=  X2 XOR (X2 / 2**S8);
     X2  := (X2 XOR (X2 * 2**S7)) mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S6);
     X2  := (X2 XOR (X2 * 2**S5)) mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S4);
     X2  := (X2 XOR (X2 * 2**S3)) mod 2**61;
     X2  :=  X2 XOR (X2 / 2**S2);
     X2  := (X2 XOR (X2 * 2**S1)) mod 2**61;
     
     Random_x := X2;

  end Get_Random_XOR_Shift_61_b;

end Xor_Shift;
