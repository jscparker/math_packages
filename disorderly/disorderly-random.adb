
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

with Text_Utilities;
with XOR_Shift;
with LCG_Rand;

package body Disorderly.Random
   with Spark_Mode => On
is

  package text_handling is new text_utilities (No_Of_Seeds, Parent_Random_Int);
  use text_handling;

  package MWC is new LCG_Rand (Parent_Random_Int);
  use MWC;

  package Xor_Rand is new XOR_Shift (Parent_Random_Int);
  use Xor_Rand;

  --  The 256 bit state of the full generator is contained in an array of four
  --  64 bit ints indexed by the following 4 ints:

  LCG_id_0     : constant State_Index := 0;
  LCG_id_1     : constant State_Index := 1;
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
  package X_squared_Mod_P_Constants is

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
     pragma Assert (a=1952409288); -- 1952409288, 1952409288
     pragma Assert (b=2**32-5);    -- 2**32-5, 2**32-5

     P  : constant := a*b-1;      -- prime
   --P2 : constant := (P-1)/2;    -- prime
   --P_minus_1 : constant := P-1;
 
     pragma Assert (P > Random_Int'Last + 2);
     pragma Assert (P < Parent_Random_Int'Last / 2);
     pragma Assert (a < b); 

  end X_squared_Mod_P_Constants;


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

     use X_squared_Mod_P_Constants;

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
     --  depends on the initial condition.

  end Get_Random_X_squared_Mod_P;

  pragma Inline (Get_Random_X_squared_Mod_P);


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
     S1, S0 : Parent_Random_Int;
  begin
     S0 := S.X(LCG_id_0);
     S1 := S.X(LCG_id_1);
     Get_Random_LCG_64_Combined (S0, S1, X1);
     S.X(LCG_id_0) := S0;
     S.X(LCG_id_1) := S1;

     X2 := S.X(XOR_Shift_Id);
     Get_Random_XOR_Shift_61 (X2);
     S.X(XOR_Shift_Id) := X2;

     Get_Random_X_squared_Mod_P (X3, S);
   
     Random_x := Random_Int ((X1 + X2 + X3) mod 2**61);

  end Get_Random;

  -- Code for initializing states:

--subtype Valid_LCG_0_Range     is Parent_Random_Int range 1 .. m0-1;
--subtype Valid_LCG_1_Range     is Parent_Random_Int range 1 .. m1-1;
--subtype Valid_XOR_Shift_Range is Parent_Random_Int range 1 .. 2**61-1;
  subtype Valid_Quadratic_Range is Parent_Random_Int range 2 .. X_squared_Mod_P_Constants.P-1;

  function Valid_State (S : in State) return Boolean is
     Result : constant Boolean :=
       (S.X(LCG_id_0)     in Valid_LCG_0_Range)     and
       (S.X(LCG_id_1)     in Valid_LCG_1_Range)     and
       (S.X(Xor_Shift_id) in Valid_Xor_Shift_Range) and
       (S.X(Quadratic_id) in Valid_Quadratic_Range);
  begin
     return Result;
  end Valid_State;

  procedure Make_Correct (S : in out State) is
  begin

     if S.X(LCG_id_0)  > Valid_LCG_0_Range'Last then
        S.X(LCG_id_0) := Valid_LCG_0_Range'Last;
     end if;

     if S.X(LCG_id_0)  < Valid_LCG_0_Range'First then -- Period=1 if = 0
        S.X(LCG_id_0) := Valid_LCG_0_Range'First;
     end if;

     if S.X(LCG_id_1)  > Valid_LCG_1_Range'Last then
        S.X(LCG_id_1) := Valid_LCG_1_Range'Last;
     end if;

     if S.X(LCG_id_1)  < Valid_LCG_1_Range'First then -- Period=1 if = 0
        S.X(LCG_id_1) := Valid_LCG_1_Range'First;
     end if;

     if S.X(Xor_Shift_id)  > Valid_Xor_Shift_Range'Last then
        S.X(Xor_Shift_id) := Valid_Xor_Shift_Range'Last;
     end if;

     if S.X(Xor_Shift_id)  < Valid_Xor_Shift_Range'First then
        S.X(Xor_Shift_id) := Valid_Xor_Shift_Range'First;
     end if;

     if S.X(Quadratic_id)  > Valid_Quadratic_Range'Last then
        S.X(Quadratic_id) := Valid_Quadratic_Range'Last;
     end if;

     if S.X(Quadratic_id)  < Valid_Quadratic_Range'First then -- Period=1 for 0,1
        S.X(Quadratic_id) := Valid_Quadratic_Range'First;
     end if;

  end Make_Correct;

  -----------
  -- Reset --
  -----------

  -- Step 1. Want each unique and valid choice
  -- of the Initiators (Initiator1, Initiator2, .. )
  -- to produce a unique and valid initial state S: S.X(0), S.X(1), ...
  --
  -- Step 2. Want a change of 1 bit or more in any one Initiator to
  -- cause changes in all elements of the initial state S: S.X(0), S.X(1), ...
  --
  procedure Reset
    (S          :    out State;
     Initiator1 : in     Seed_Random_Int := 1111;
     Initiator2 : in     Seed_Random_Int := 2222;
     Initiator3 : in     Seed_Random_Int := 3333;
     Initiator4 : in     Seed_Random_Int := 4444)
  is
     Seeds : Vals := (Initiator1, Initiator2, Initiator3, Initiator4);

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
        Make_Correct (S);
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

  function Are_Equal (State_1, State_2 : State) return Boolean is
     Result : Boolean := True;
  begin
     for i in State_1.X'Range loop
        if State_1.X(i) /= State_2.X(i) then Result := False; end if;
     end loop;
     return Result;
  end Are_Equal;
       
end Disorderly.Random;

