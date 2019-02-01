
-----------------------------------------------------------------------
-- package body e_Derivs, high-order automatic differentiation of functions
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

with Ada.Numerics;

package body e_Derivs is

   pragma Assert (Deriv_Index'Last > 1);

   Real_Val, Factorial : Derivatives;
   --  global constant arrays, initialized after the begin statement.

   Zero : constant Real := (+0.0);
   One  : constant Real := (+1.0);
   Two  : constant Real := (+2.0);
   Half : constant Real := (+0.5);

  --------------------------
  -- Make_Constant_Arrays --
  --------------------------

  -- Try to minimize inexact arithmetic.

  procedure Make_Constant_Arrays is
  begin

     Real_Val(0)   := Zero;
     for I in 1..Deriv_Index'Last loop
        Real_Val(I)   := Real_Val(I-1) + One;
     end loop;

     Factorial(0)  := One;
     for I in 1..Deriv_Index'Last loop
        Factorial(I)  := Real_Val(I) * Factorial(I-1);
     end loop;

  end Make_Constant_Arrays;

  ---------------
  -- Un_Reduce --
  ---------------

  -- Multiply by (Order of Derivative) factorial to make reduced derivative
  -- into ordinary Derivative.

  procedure Un_Reduce (F : in out Derivatives) is
  begin
     for I in 2..Deriv_Index'Last loop
        F(I) := F(I) * Factorial(I);
     end loop;
  end Un_Reduce;

  ------------------
  -- Make_Reduced --
  ------------------

  -- Divide by (Order of Derivative) factorial to make true derivative into
  -- Leibnitz reduced Derivative.

  procedure Make_Reduced (F : in out Derivatives) is
  begin
    for I in 2..Deriv_Index'Last loop
       F(I) := F(I) / Factorial(I);
    end loop;
  end Make_Reduced;

  ---------
  -- "+" --
  ---------

  function "+" (F, G : in  Derivatives) return Derivatives is
     Result : Derivatives;
  begin
    for I in Deriv_Index loop
       Result(I) := F(I) + G(I);
    end loop;
    return Result;
  end "+";

  ---------
  -- "-" --
  ---------

  function "-" (F, G : in  Derivatives) return Derivatives is
     Result : Derivatives;
  begin
    for I in Deriv_Index loop
       Result(I) := F(I) - G(I);
    end loop;
    return Result;
  end "-";

  ---------
  -- "*" --
  ---------

  -- Derivatives, up to n, of H(t) :=  F(t) * G(t).  Notation: (d/dt)**n = d^n.
  -- so write (d/dt)**n (F(t) * G(t)) == d^n (F*G) = d^n (H).
  -- Returns the Leibnitz Reduced derivative of H.  This routine assumes
  -- that F and G are already in reduced form (i.e. (d/dt)**k F(t) / k! is
  -- the k-th reduced derivative of F, and is stored in the k-th place of
  -- array F).  Notice that this means that the j-th reduced derivative of
  -- H is BackwardDotProduct (F, G) := SUM(k in 0..j){F(k) * G(j-k)}.  The
  -- pascal coeffients are cancelled out by the factorials in the reduced
  -- derivatives.  If you need all derivatives from 0 to n, then this
  -- routine is the most efficient way.
  -- If you just need the n-th derivative it might be faster to avoid
  -- Leibnitz reduced derivatives.
  --
  function "*" (F, G : in  Derivatives) return Derivatives is
     H   : Derivatives := (others => Zero);
     Sum : Real        := Zero;
  begin
     for Order in Deriv_Index loop
        Sum := Zero;
        for I in 0..Order loop
           Sum := Sum + F(I) * G(Order - I);
        end loop;
        H(Order) := Sum;
     end loop;

     return H;

  end "*";

  -------------
  -- Compose --
  -------------

  -- Routine gets derivatives of composition of F and G: H(t) := (F(G(t)).
  -- Want all derivatives up to n.  Notation: the j-th derivative of H is
  -- written d^j (H) := (d/dt)**j H(t), and the j-th derivative of F is written
  -- Fj.  To get recursive formula for d^n (H), use
  --                    d^n (H) := d^(n-1) (G1 * F1(G(t))).
  -- Next use the formula for the (n-1)th derivative of the product of
  -- G1 and F1(G(t)).  then have a formula for the n-th derivative
  -- of F(G(t)) in terms of the (n-1)th derivative of the composition of 2
  -- functions, F1 and G: F1(G(t)).  Can do it with recursive calls, but
  -- here calculate the intermediate quantities explicitly, so that inner
  -- inner loops are optimized.
  -- Suppose want d^n H(t) := d^n (F(G(t)).  Start with
  --
  -- d^0 Fn-0(G(t))                                                     (j=n-0)
  -- d^0 Fn-1(G(t))   d^1 Fn-1(G(t))                                    (j=n-1)
  -- d^0 Fn-2(G(t))   d^1 Fn-2(G(t))  d^2 Fn-2(G(t))                    (j=n-2)
  --
  -- ....
  --
  -- d^0 F1(G(t))     d^1 F1(G(t))    d^2 F1(G(t))  ...  d^n-1 F1(G(t))   (j=1)
  -- d^0 F0(G(t))     d^1 F0(G(t))    d^2 F0(G(t))  ...  d^n-1 F0(G(t))   (j=0)
  --
  -- with a missing d^n F0(G(t)) in the final row and column.  The final row
  -- is of course the desired result.
  --
  -- Each row above is calculated from the previous row by the recursive formula
  -- given above. For example, once have the j=1 row can get the
  -- desired result using:
  --
  --   d^n (H)  =  d^(n-1) (G1 * F1(G(t))) :=
  --
  --              n-1
  --            = SUM { d^(n-1-k) G1 * d^k (F1(G(t)) * Coeff(n-1,k) }
  --              k=0
  --
  -- Above, Coeff(n-1,j) is the Pascal's triangle coefficient (n-1)!/k!(n-1-k)!.
  -- Can remove this entirely from the sum by using Leibnitz reduced form:
  -- divide each argument of d^j by j!.  Then multiply the final sum by (n-1)!
  -- to get the true derivative.  Avoid the final multiplication by (n-1)! if
  -- you want the reduced derivative of H.  (Still need to multiply the result
  -- by a factor of 1/n to make d^n (H) the reduced derivative of H.)
  -- This saves much time if you only rarely convert from the reduced form to the
  -- true derivatives.  To do it this way, must put G1 into reduced
  -- form.  are given F and G in reduced form.  Must shift G array
  -- down one place, and multiply the components by (k+1) (k is the index) to
  -- adjust the 1/(k+1)! term to 1/k!.  If d^k (F1(G(t)) is the reduced
  -- derivative of the argument F1(G(t), then the formula for the reduced
  -- derivatives of H := F(G(t)) becomes,
  --
  --           n-1
  -- d^n (H) = SUM { d^(n-1-k) G1 * d^k (F1(G(t)) } / n
  --           k=0
  --
  --         == SUM { G1(n-1-k) * DF(1,k) } / n
  --
  -- Below the array DF(j,k) will hold d^k Fj(G(t)).  For example, the final row
  -- above will be placed into DF(1,k), where k is in 0..n-1.  Row j of DF(j,k)
  -- is calculated from row j+1 by setting H := Fj(G(t)) using
  --
  --           k-1
  -- d^k (H) = SUM { d^(k-1-i) G1 * d^i (Fj+1(G(t)) } * (j+1) / k
  --           i=0
  --
  -- The j+1 term comes from the fact that Fj is also reduced, and when we
  -- took the true derivative of it, had to divide and multiply by j+1 to
  -- keep it a reduced derivative.  Notice that the affect of this repeated
  -- multiplication of Fj+1 by j+1 is simply to unreduce F.  If F had been
  -- introduced as unreduced to begin with, then the formula would be identical
  -- but without the j+1 factor.  So it's both more efficient and more accurate
  -- to unreduce F at the beginning rather than along the way.  If do that
  -- then the final formula for the reduced derivative d^k Fj(G(t)) becomes:
  --
  --   d^k (H)    :=    SUM { G1(k-1-i) * DF(j+1,i) } / k    :=    DF(j,k)
  --

  function Compose
    (F, G : in Derivatives)
     return Derivatives
  is
     G1, F0, H : Derivatives;
     DF : array (Deriv_Index, Deriv_Index) of Real;
     Col, J, K : Deriv_Index;
     Sum : Real;
  begin
      -- The reduced derivatives of F and G are placed in arrays F and G.  Array
      -- F holds the reduced derivatives of function F evaluated at G(t), not
      -- the derivatives of F(G(t)), which want to calculate.  now need the
      -- array of reduced derivatives of the 1st derivative of G, called G1, and
      -- need to unreduce F, as explained above.

      for I in 1..Deriv_Index'Last loop
         G1(I-1) := G(I) * Real_Val(I);
      end loop;
      for I in Deriv_Index loop
         F0(I) := F(I) * Factorial(I);
      end loop;

      --************************************************************************
      -- Get the first Column of DF, the ones with the 0-th derivatives of Fk(G(t)).
      -- The col j is constant at 0, and k varies from Max_Order_Of_Deriv to 0.
      --************************************************************************
      Col := 0;
      for Row in Deriv_Index loop
         k := Deriv_Index'Last - Row;
         DF(k, Col) := F0(k);
      end loop;

      --************************************************************************
      -- Get the next rows of DF, 1..Max_Order_Of_Deriv.  Row and column refer to the
      -- Rows and Cols of the table above.  k and j are the indices of DF, which
      -- match the k and j in d^k Fj(G(t)).
      --************************************************************************
      for Row in 1..Deriv_Index'Last loop
           j := Deriv_Index'Last - Row;
           for Col in 1..Row loop
                k := Col;

      --************************************************************************
      -- Calculate d^k Fj(G(t)) := DF(j,k) given DF(j+1, k).
      --   SUM { G1(k-1-i) * DF(j+1,k) } / k.
      --************************************************************************
            Sum := Zero;
            for i in 0..k-1 loop
               Sum := Sum + G1(k-1-i) * DF(j+1,i);
            end loop;
            DF(j,k) := Sum / Real_Val (k);  -- "/" always more accurate, and k>0
          --DF(j,k) := Sum * Inverse_Of(k);

          end loop;
      end loop;

      --************************************************************************
      -- The desired derivatives are the final row := Max_Order_Of_Deriv, (j=0)
      -- of DF.  This is a reduced form.
      --************************************************************************
      for I in Deriv_Index loop
         H(I) := DF(0, I);
      end loop;

      return H;

  end Compose;

  ----------------
  -- Reciprocal --
  ----------------

  --  Returns reduced derivatives of 1 / F(t), given reduced derivs of F.
  --  Just use composition of A(t) = 1/t with F(t): H(t) = A(F(t)) = 1/F(t),
  --  so H(t) = Compose (Reciprocal, F).
  --
  function Reciprocal (F : in  Derivatives) return Derivatives is
  begin
     if F(0) = Zero then
       raise d_Argument_Error;
     end if;

     return Compose (Reciprocal (F(0)), F);

  end Reciprocal;

  ---------
  -- "/" --
  ---------

  --  Returns reduced derivatives of F(t) / G(t), given reduced derivs of F, G.
  --  Just use: H(t) = F(t) * Reciprocal (G(t)).
  --
  function "/" (F, G : in  Derivatives) return Derivatives is
  begin
     if G(0) = Zero then
       raise d_Argument_Error;
     end if;

     return F * Reciprocal (G);

  end "/";

  ----------
  -- "**" --
  ----------

  -- Get reduced derivatives of F(t) to an integer Power: F(t)**n.
  -- H(t) = En(F(t)) where En(t') = (t')**n, so that H(t) = F(t)**n.  So use
  -- Compose (En, F), where En = F(0)**n.
  -- (or H(t) = Pow_n(F(t)).

  function "**"
    (F        : in Derivatives;
     Exponent : in Natural)
     return Derivatives
  is
  begin

    return Compose (F(0)**Exponent, F);

  end "**";

  ----------
  -- "**" --
  ----------

  -- Get reduced derivatives of Time to a non-negative integer Power.
  -- Follow the Ada standard: Zero**0 is One by definition.
  -- Anything to the 0 is defined to be 1.
  -- (Reduced derivs are just pascals triangle * t**n.)

  function "**"
    (Time     : in Real;
     Exponent : in Natural)
     return Derivatives
  is
     Derivs        : Derivatives := (others => Zero);
     Power_Of_t    : Derivatives := (others => Zero);
     Coeff         : Derivatives;
     Exp_Of_t, Exp : Integer     := Exponent;
     Highest_Order : Deriv_Index := Deriv_Index'Last;
  begin

     if Exponent = 0 then
        Derivs (0) := One;  -- Even Zero**0 is 1. All derivs are Zero.
        return Derivs;
     end if;

     -- Exponent > 0, so Time=0 implies Time**N = Zero, but the N-th deriv is N!
     -- The N-th deriv in reduced form is One.  All other derivs are Zero.

     if Time = Zero then
        if Exponent <= Integer (Deriv_Index'Last) then
           Derivs (Deriv_Index (Exponent)) := One;  -- All other init to Zero.
        end if;
        return Derivs;  -- already reduced.
     end if;

     -- Now know that Exponent > 0, and Time /= Zero.
     -- Of course the high powers of time can under/overflow.

     Exp_Of_t      := Exponent;
     Power_Of_t(0) := Time ** Exp_Of_t;
     Highest_Order := Deriv_Index'Last;

     for I in 1..Deriv_Index'Last loop
        Exp_Of_t := Exp_Of_t - 1;
        if Exp_Of_t > 0 then
           Power_Of_t(I) := Time ** Exp_Of_t; -- slow but better accuracy.
        else
           Power_Of_t(I) := One;
           Highest_Order := I;
           exit;
        end if;
     end loop;


     -- Coeff has been initialized to Zero:

     Coeff (0) := One;
     Exp       := Exponent;

     for I in 1..Highest_Order loop
       if Exp > 0 then
          Coeff(I) := Coeff(I-1) * (+Real_8 (Exp));
       else
          exit;
       end if;
       Exp := Exp - 1;
     end loop;

     Derivs(0) := Power_Of_t(0);
     for I in 1..Highest_Order loop
        Derivs(I) := Power_Of_t(I) * Coeff(I);
     end loop;

     Make_Reduced (Derivs);

    return Derivs;

  end "**";

  ----------------
  -- Reciprocal --
  ----------------

  -- Get reduced derivatives of One / Time.  Remember, high order derivs go as
  -- (1/Time)**Order.  Easy get underflow/overflow if Time too large or small.
  --
  function Reciprocal 
    (Time : in  Real)
     return Derivatives
  is
     Result : Derivatives;
     Coeff : Real := One;
  begin

     if Time = Zero then
        raise d_Argument_Error;
     end if;

     Result (0) := One / Time;
     for I in 1..Deriv_Index'Last loop
        Coeff  := -Coeff;
        Result (I) := Coeff / Time ** (Integer(I+1)); --detectable accuracy imprv.
      --Result (I) := -Result (I-1) / Time;
     end loop;
     --  Notice that neglect the I! coefficient of Result, which is
     --  lost when Result(I) is divided by I! to get reduced deriv.

     return Result;

   end Reciprocal;

   ------------
   -- Sqrt_d --
   ------------

   function Sqrt_d
     (Time      : in Real;
      Frequency : in Real := One_d;
      Phase     : in Real := Zero_d)
      return Derivatives
   is
      Result : Derivatives;
      Arg  : Real;
   begin

      Arg := Time * Frequency + Phase;

     --  Can't take first deriv at 0.0:

      if Arg = Zero then
        raise d_Argument_Error;
      end if;

      Result (0) := Sqrt (Arg);

      Result(1)  := Result (0) * Half * Frequency / Arg;

      for I in Deriv_Index'First+2 .. Deriv_Index'Last loop
         Result(I) := -Result (I-1) * (Real_Val (I-1) - Half) * Frequency;
         Result(I) :=  Result (I) / (Arg * Real_Val (I));
      end loop;

      -- reduced by the  / Real_Val (I) above.

      return Result;

   end Sqrt_d;

   --------------
   -- Arctan_d --
   --------------

   --  Reduced derivs of Arctan.
   --
   function Arctan_d
     (Time : in Real) 
      return Derivatives
   is
      Result : Derivatives;
      One_Plus_Arg_Squared  : Derivatives := (others => Zero);
      One_Over : Derivatives;
      Arg : constant Real := Time;
   begin

      --  Reduced derivatives (divide by order of derivative factorial):
      One_Plus_Arg_Squared (Deriv_Index'First)   := One + Arg*Arg;
      One_Plus_Arg_Squared (Deriv_Index'First+1) := Two * Arg;
      One_Plus_Arg_Squared (Deriv_Index'First+2) := One; -- reduced by / 2!

      One_Over := Reciprocal (One_Plus_Arg_Squared (0));

      Result := Compose (One_Over, One_Plus_Arg_Squared);

      for I in reverse Deriv_Index'First+1 .. Deriv_Index'Last loop
         Result (I) := Result (I-1) / Real_Val (I);
      end loop;

      Result (Deriv_Index'First) := Arctan (Arg);

      return Result;

   end Arctan_d;

   --------------
   -- Arcsin_d --
   --------------

   --  Reduced derivs of Arcsin.  Arg must be in (-1.0, 1.0).  Result
   --  of Arcsin is in range -pi/2 .. pi/2.
   --
   function Arcsin_d
     (Time : in Real) 
      return Derivatives
   is
      Result : Derivatives;
      One_Minus_Arg_Squared  : Derivatives := (others => Zero);
      Reciprocal_Square_Root : Derivatives;
      Arg : constant Real := Time;
   begin

      --  Can't take derivs at -1.0 or 1.0:

      if (Arg <= -One) or (Arg >= One) then
         raise d_Argument_Error;
      end if;

      --  Reduced derivatives (divide by order of derivative factorial):
      One_Minus_Arg_Squared (Deriv_Index'First)   :=  One - Arg*Arg;
      One_Minus_Arg_Squared (Deriv_Index'First+1) := -Two * Arg;
      One_Minus_Arg_Squared (Deriv_Index'First+2) := -One;

      Reciprocal_Square_Root := Reciprocal (Sqrt_d (One_Minus_Arg_Squared (0)));

      Result := Compose (Reciprocal_Square_Root, One_Minus_Arg_Squared);

      for I in reverse Deriv_Index'First+1 .. Deriv_Index'Last loop
         Result (I) := Result (I-1) / Real_Val (I);
      end loop;

      Result (Deriv_Index'First) := ArcSin (Arg);

      return Result;

   end Arcsin_d;

   --------------
   -- Arccos_d --
   --------------

   --  Reduced derivs of Arccos.  Arg must be in (-1.0, 1.0).  Result
   --  of Arccos is in range 0.0 .. pi.  Arccos = Pi/2 - Arcsin, so call
   --  Arcsin.

   function Arccos_d
     (Time : in Real) return Derivatives
   is
      Result : Derivatives;
      Arg : constant Real := Time;
   begin

      --  Can't take derivs at -1.0 or 1.0:

      if (Arg <= -One) or (Arg >= One) then
          raise d_Argument_Error;
      end if;

      Result := Arcsin_d (Arg);

      for I in Deriv_Index'First+1 .. Deriv_Index'Last loop
         Result (I) := -Result (I);
      end loop;

      Result (Deriv_Index'First) := ArcCos (Arg);

      return Result;

   end Arccos_d;

  -----------
  -- Sin_d --
  -----------

  -- Get reduced derivatives of Sin
  --
  function Sin_d
    (Time : in Real;
     Frequency : in Real := One_d;
     Phase     : in Real := Zero_d)
     return Derivatives
  is
     SI, CO        : Real;
     FreqPower     : Real    := One;
     Index_Is_Even : Boolean := True;
     Sinusoid      : Derivatives;
  begin
     SI    :=  Sin (Frequency*Time + Phase);
     CO    :=  Cos (Frequency*Time + Phase);

     Index_Is_Even := True;
     FreqPower     := One;
     for I in Deriv_Index loop
        if Index_Is_Even then
           Sinusoid(I) := FreqPower * SI;
           FreqPower   := FreqPower * Frequency;
        else
           Sinusoid(I) := FreqPower * CO;
           FreqPower   := -FreqPower * Frequency;
        end if;
        Index_Is_Even := not Index_Is_Even;
     end loop;

     Make_Reduced (Sinusoid);

     return Sinusoid;

  end Sin_d;

  -----------
  -- Cos_d --
  -----------

  -- Get reduced derivatives of Cos
  --
  function Cos_d
    (Time : in Real;
     Frequency : in Real := One_d;
     Phase     : in Real := Zero_d)
     return Derivatives
  is
     SI, CO        : Real;
     FreqPower     : Real    := One;
     Index_Is_Even : Boolean := True;
     CoSinusoid : Derivatives;
  begin
     SI    :=  Sin (Frequency*Time + Phase);
     CO    :=  Cos (Frequency*Time + Phase);

     Index_Is_Even := True;
     FreqPower     := One;
     for I in Deriv_Index loop
        if Index_Is_Even then
           CoSinusoid(I) := FreqPower * CO;
           FreqPower     := -FreqPower * Frequency;
        else
           CoSinusoid(I) := FreqPower * SI;
           FreqPower     := FreqPower * Frequency;
        end if;
        Index_Is_Even := not Index_Is_Even;
     end loop;

     Make_Reduced (CoSinusoid);

     return CoSinusoid;

  end Cos_d;

  -----------
  -- Exp_d --
  -----------

  -- Get reduced derivatives Exp, the envelope.
  --
  function Exp_d
    (Time : in Real;
     Frequency : in Real := One_d;
     Phase     : in Real := Zero_d)
     return Derivatives
  is
     Val, Arg  : Real;
     FreqPower : Real := One;
     Expon     : Derivatives;
  begin
     Arg := Frequency*Time + Phase;
     --  May want to do some special casing with this.  If large negative
     --  for example, set Val to Zero, etc.

     Val := Exp (Arg);

     FreqPower := One;
     for I in Deriv_Index loop
        Expon(I)  := FreqPower * Val;
        FreqPower := FreqPower * Frequency;
     end loop;

     Make_Reduced (Expon);

     return Expon;

  end Exp_d;

  -----------
  -- Log_d --
  -----------

  -- Get reduced derivatives of natural Log (base e).
  --
  function Log_d
    (Time : in Real;
     Frequency : in Real := One_d;
     Phase     : in Real := Zero_d)
     return Derivatives
  is
     Arg, FreqPower : Real;
     Result : Derivatives;
  begin
     Arg := Frequency*Time + Phase;

     if Arg <= Zero then
       raise d_Argument_Error;
     end if;

     Result := Reciprocal (Arg);
     --  Derivatives of 1/t evaluated at t=arg

     for I in reverse 1..Deriv_Index'Last loop
        Result(I) := Result(I-1) / Real_Val(I);
     end loop;
     --  Derivatives of 1/t are higher derivatives of Log(t).
     --  Need to properly reduce tho: turn (I-1)! into I!.

     Result(0) := Log (Arg);

     FreqPower := One;
     for I in 1..Deriv_Index'Last loop
       FreqPower := FreqPower * Frequency;
       Result(I) := FreqPower * Result(I);
     end loop;

     return Result;

  end Log_d;

  ----------------
  -- Taylor_Sum --
  ----------------

  --  Want   Sum  =  B_0 + B_1*t + B_2*t**2 + ... + B_n*t**n.
  --  or in Horner's form: Sum  =  B_0 + t*(B_1 + ... + t*(B_n-1 + t*B_n)))))).
  --  This is easily written as matrix equation, with Sum = S_0:
  --
  --  S_n = B_n;  S_n-1 = B_n-1 + t*S_n;  S_1 = B_1 + t*S_2; S_0 = B_0 + t*S_1;
  --
  --  In matrix form, vector S is the solution to matrix equation M*S = B,
  --  where B = (B_0,...,B_n), S = (S_0,...,S_n) and matrix M is equal to
  --  the unit matrix I minus t*O1, where O1 is all 1's on the 1st upper off-
  --  diagonal.
  --
  --  M  is not diag. domin. if t>>1, so:
  --  This form is chosen because the solution vector S can
  --  be improved numerically by iterative refinement with Newton's
  --  method:
  --             S{k+1} = S{k} + M_inverse * (B - M*S{k})
  --
  --  where S = M_inverse * B is the calculation of S given above.  If the
  --  said calculation of S is numerically imperfect, then the iteration above
  --  may produce improved values of S.  If the Coefficients of
  --  the polynomial B are numerically poor, then this effort will be wasted.
  --  It's often the case with
  --
  --  Iterative refinement
  --
  --  If  No_Of_Iterations=0  then usual solution is returned.
  --  If  No_Of_Iterations=1  then solution is refined iteratively once.
  --
  --  This version refines residual rather than actual solution, (so that
  --  can subtract small quantities from small, rather than large.)
  --  If y = exact error in 1st iteration: y = S{1} - S{inf}, then y is the
  --  solution of M*y = d_1 where d_1 = B - M*S{1}.
  --  Let M_inv denote approximate inverse of M.  Iterate for y using
  --
  --      Delta_Y_{k+1} == Y_{k+1} - Y_k = M_inv*(d_1 - M*Y_k).
  --
  --  Remember Y = SUM (Delta_Y_k's).  Here's the actual method:
  --
  --   Delta_Y_1 = M_inv*d_1
  --   Let d_2 == d_1 - M*Delta_Y_1
  --   Delta_Y_2 = M_inv*(d_1 - M*Delta_Y_1) = M_inv*d_2
  --   Let d_3 == d_2 - M*Delta_Y_2
  --   Delta_Y_3 = M_inv*(d_1 - M*Delta_Y_1 - M*Delta_Y_2) = M_inv*d_3
  --
  --  so:    d_k = d_{k-1} - M*Delta_Y_{k-1}; Delta_Y_k = M_inv*d_k
  --
  --  Sum the Delta_Y_k's to get the correction to S(1): Y = SUM (Delta_Y_k's).
  --  notice are iterating for the error in the error in the error...
  --  First method doesn't seem to work.
  --  Second method works best in simple tests, so use it.
  --
  function Taylor_Sum
    (t                : in Real;
     F                : in Derivatives;
     Max_Index        : in Deriv_Index   := Deriv_Index'Last;
     No_Of_Iterations : in Natural       := 0)
     return Real
  is
     S_first : Derivatives;
     D_k, Delta_Y_k : Derivatives;
     Y_f : Derivatives := (others => Zero); --init impor.
     B : Derivatives renames F;

     --  Get Product = M*S(k):

     function M_times
       (S                : in Derivatives;
        t                : in Real;
        Coeff_Last       : in Deriv_Index   := Deriv_Index'Last)
        return Derivatives
     is
        Product : Derivatives;
     begin
        for n in Coeff_Last .. Deriv_Index'Last loop
           Product(n) := Zero;
        end loop;
        Product(Coeff_Last) := S(Coeff_Last);
        if Deriv_Index'First < Coeff_Last then
           for n in reverse Deriv_Index'First..Coeff_Last-1 loop
              Product(n) := S(n) - t * S(n+1);
           end loop;
        end if;
        return Product;
     end;

     --  Solve for S in the matrix equation M*S = B

     function M_inv_times
       (B                : in Derivatives;
        t                : in Real;
        Coeff_Last       : in Deriv_Index   := Deriv_Index'Last)
        return Derivatives
     is
        S : Derivatives;
     begin
        for n in Coeff_Last .. Deriv_Index'Last loop
           S(n) := Zero;
        end loop;
        S(Coeff_Last) := B(Coeff_Last);
        if Deriv_Index'First < Coeff_Last then
           for n in reverse Deriv_Index'First..Coeff_Last-1 loop
              S(n) := B(n) + t * S(n+1);
           end loop;
        end if;
        return S;
     end;

     function Difference
       (A                : in Derivatives;
        B                : in Derivatives;
        Coeff_Last       : in Deriv_Index   := Deriv_Index'Last)
        return Derivatives
     is
        S : Derivatives;
     begin
        for n in Coeff_Last .. Deriv_Index'Last loop
           S(n) := Zero;
        end loop;
        for n in reverse Deriv_Index'First..Coeff_Last loop
           S(n) := A(n) - B(n);
        end loop;
        return S;
     end;

  begin
     --  Get S{1} = S_first, 1st soln, as defined above.

     S_first := M_inv_times (B, t, Max_Index);   --does init ok

     If No_Of_Iterations > 0 then

        --  get d_k = B - M*S_first,  the 1st residual:

        D_k := Difference (B, M_times (S_first, t, Max_Index), Max_Index); -- D_1

        --  get Delta_Y = M_inv*D_k:

        Delta_Y_k := M_inv_times (D_k, t, Max_Index);
        Y_f := Delta_Y_k;

        for Iteration in 1..No_Of_Iterations loop

           --  get d_k = d_k - M*Delta_Y_k

           D_k := Difference(D_k, M_times (Delta_Y_k, t, Max_Index), Max_Index);

           --  get Delta_Y = M_inv*D_k:

           Delta_Y_k := M_inv_times (D_k, t, Max_Index);

           --  Increment solution for full correction to S_first, Y_f:

           for I in Deriv_Index'First..Max_Index loop
              Y_f(I) := Y_f(I) + Delta_Y_k(I);
           end loop;

        end loop;

     end if;

     --for I in Deriv_Index'First..Max_Index loop
     --   S(I) := Y_f(I) + S_first(I);
     --end loop;

     --if Max_Index < Deriv_Index'Last then
     --for I in Max_Index+1 ..Deriv_Index'Last loop
     --   S(I) := Zero;
     --end loop;
     --end if:

     return Y_f(Deriv_Index'First) + S_first(Deriv_Index'First);

  end Taylor_Sum;

  ------------------
  -- Taylor_Sum_2 --
  ------------------

  --  Want   Sum  =  B_0 + B_1*t + B_2*t**2 + ... + B_n*t**n.
  --  or in Horner's form: Sum  =  B_0 + t*(B_1 + ... + t*(B_n-1 + t*B_n)))))).
  --  This is easily written as matrix equation, with Sum = S_0:
  --
  --  S_n = B_n;  S_n-1 = B_n-1 + t*S_n;  S_1 = B_1 + t*S_2; S_0 = B_0 + t*S_1;
  --
  --  In matrix form, vector S is the solution to matrix equation M*S = B,
  --  where B = (B_0,...,B_n), S = (S_0,...,S_n) and matrix M is equal to
  --  the unit matrix I minus t*O1, where O1 is all 1's on the 1st upper off-
  --  diagonal.
  --
  function Taylor_Sum_2
    (t                : in Real;
     F                : in Derivatives;
     Max_Index        : in Deriv_Index   := Deriv_Index'Last) return Real is

     S_first, Bu : Derivatives;
     B : Derivatives renames F;

     --  Solve for S in the matrix equation M*S = B

     function M_inv_times
       (B                : in Derivatives;
        t                : in Real;
        Coeff_Last       : in Deriv_Index   := Deriv_Index'Last)
        return Derivatives
     is
        S : Derivatives;
     begin
        for n in Coeff_Last .. Deriv_Index'Last loop
           S(n) := Zero;
        end loop;
        S(Coeff_Last) := B(Coeff_Last);
        if Deriv_Index'First < Coeff_Last then
           for n in reverse Deriv_Index'First..Coeff_Last-1 loop
              S(n) := B(n) + t * S(n+1) / Real_val(n+1);
           end loop;
        end if;
        return S;
     end;

  begin

     Bu := B;
     un_reduce(Bu);

     --  Get S{1} = S_first, 1st soln, as defined above.

     S_first := M_inv_times (Bu, t, Max_Index);   --does init ok

     return S_first(Deriv_Index'First);

  end Taylor_Sum_2;

  function Derivative_Of
    (F              : in  Derivatives;
     Order_Of_Deriv : in Deriv_Index)
     return Derivatives is

     Result : Derivatives := F;

  begin

     if Order_Of_Deriv = 0 then return Result; end if;

     for Order in 1..Order_Of_Deriv loop

     -- shft rgt
     for n in Deriv_Index'First..Deriv_Index'Last-1 loop
        Result(n) := Result(n+1);
     end loop;
     Result(Deriv_Index'Last) := Zero;


     -- reduce
     for n in Deriv_Index'First+1..Deriv_Index'Last-Order loop
        Result(n) := Result(n) * Real_Val(n+1);
     end loop;

     end loop;

     return Result;

  end Derivative_Of;

  --  Integration const is set to 0.

  function Integral_Of
    (F : in  Derivatives)
     return Derivatives is

     Result : Derivatives := F;

  begin

     -- shft left
     for n in reverse Deriv_Index'First+1..Deriv_Index'Last loop
        Result(n) := Result(n-1);
     end loop;
     Result(Deriv_Index'First) := Zero;

     -- reduce
     for n in Deriv_Index'First+2..Deriv_Index'Last loop
        Result(n) := Result(n) / Real_Val (n);
     end loop;

     return Result;

  end Integral_Of;

begin

  Make_Constant_Arrays;

end e_Derivs;
