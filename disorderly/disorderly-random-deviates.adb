
-------------------------------------------------------------------------------
-- package body Disorderly.Random.Deviates, Floating point random deviates.
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

with Ada.Numerics.Generic_elementary_functions;
with Gamma;
package body Disorderly.Random.Deviates is

   Sqrt_of_Min_Real : constant Real := Two**(Real'Machine_Emin / 2);
   Sqrt_of_Two_Pi   : constant Real := 2.506628274631000502415765284811;
   Natural_Log_of_2 : constant Real := 0.693147180559945309417232121458;

   package Math is new Ada.Numerics.Generic_elementary_functions (Real);
   use Math;
   package Gam is new Gamma (Real);
   use Gam;

   -------------------
   -- Log_of_1_plus --
   -------------------

   --  More accurate than using Log (1 + x) for Abs x << 1.

   function Log_of_1_plus
      (x : Real)
      return Real
   is
      u   : constant Real := One + x;
      Eps : constant Real := Two ** (-Real'Machine_Mantissa / 2 - 6);
   begin
      if u <= Zero then
         raise Constraint_Error;
      end if;

      -- Use Log(1+x) = x - x^2/2 + x^3/3 - ... = x*(1 - x/2 + x^2/3 ...)
      -- So if x is somewhat smaller than Sqrt (Real'Epsilon) then 1 + x^2/3 = 1:

      if Abs (u - One) < Eps then
         return x - Half * x * x;
      end if;

      return Log(u) * x / (u - One);    -- u /= One; see above.

   end Log_of_1_plus;

   ----------
   -- "**" --
   ----------

   function "**"
      (Left, Right : Real)
      return Real
   is
      Arg : Real;
   begin
      if Left < Zero then
         raise Constraint_Error;
      end if;

      if Right = Zero then    -- notice this makes 0**0 = 1
         return One;
      end if;

      if Left = Zero then
         return Zero;
      end if;

      Arg := Right * Log (Left);

      if Arg < -Two_to_the_Ninth then
         return Zero;
      end if;

      if Arg > Two_to_the_Ninth then   -- or let Exp decide about this?
         raise Constraint_Error;
      end if;

      return Exp (Arg);

   end "**";

   -----------------
   -- Random_Real --
   -----------------

   --  Returns random reals in interval [0, 1).

   pragma Assert (Real'Digits >= 15);
   -- use 53+ bit Reals to get the most out of our 61 bit Generator.

   No_of_Bits_in_Mantissa : constant := 53;
   -- leave it here, even if using 18 digit reals.

   Scale : constant Real := Two**(-No_of_Bits_in_Mantissa);
   -- the Asserts below help optimize this. Its best to just leave
   -- No_of_Bits_in_Mantissa = 53 bits, and use Real'Digits = 15.

   procedure Get_Random_Real
     (Random_Real : out    Real;
      Stream      : in out State)
   is
      Test1 : Real := One - Scale;
      pragma Volatile (Test1); --if this is disabled, then assert test stops working
      pragma Assert (Test1 < One);
      -- Prevents:  Random_Real = 1.
      -- If Real'Digits = 15, and you get the usual IEEE floats, then
      -- No_of_Bits_in_Mantissa = 53 bits passes this; 54 doesn't.

      X : Random_Int;
   begin
      Get_Random (X, Stream);
      Random_Real := Real (X mod 2**No_of_Bits_in_Mantissa) * Scale;
   end Get_Random_Real;

   ---------------------
   -- Get_Exponential --
   ---------------------

   procedure Get_Exponential
     (Mean   : in     Real;
      Stream : in out State;
      Result : out    Real)
   is
      X : Real;
   begin
      if Mean <= Zero then raise Constraint_Error; end if;
      Get_non_Zero_Rand: loop
         Get_Random_Real (X, Stream);
	 exit Get_non_Zero_Rand when X > Zero;
      end loop Get_non_Zero_Rand;
      Result := -Mean * Log (X);
   end Get_Exponential;

   ------------------
   -- Get_Rayleigh --
   ------------------

   --  The Rayleigh distribution probability density:
   --
   --   f(x) = 2 * x * Exp (-x**2)
   --

   procedure Get_Rayleigh
     (Stream : in out State;
      Result : out    Real)
   is
      X : Real;
   begin

      Get_Exponential (Mean => One, Stream => Stream, Result => X);
      Result := Sqrt (X);

   end Get_Rayleigh;

   -----------------
   -- Get_Weibull --
   -----------------

   --  The Weibull distribution probability density:
   --
   --   f(x) = a * x**(a-1) * Exp (-x**a)
   --

   procedure Get_Weibull
     (a      : in     Real;
      Stream : in out State;
      Result : out    Real)
   is
      X : Real;
   begin
      if a = Zero then  raise Constraint_Error;  end if;

      Get_Exponential (Mean => One, Stream => Stream, Result => X);

      Result := X**(One / a);
      --  Get_Exponential  should never return X=0.

   end Get_Weibull;

   ------------------------------
   -- Neg_Binomial_Probability --
   ------------------------------

   --  [Gamma(r + k)/[(Gamma(r)*k!] * p^r * (1-p)^k
   --
   -- Probability of r successes and k failures in n = k+r Bernoulli trials.
   -- Assumes that the final trial is a success, and p = success probability.

   function Neg_Binomial_Probability
     (r : in Real;
      k : in Integer;
      p : in Real)
      return Real
   is
      Arg : Real;
   begin
      if  k < 0   then  return Zero;  end if;

      if  not p'Valid then  raise Constraint_Error;  end if;
      if  p <= Zero   then  raise Constraint_Error;  end if;
      if  p >= One    then  raise Constraint_Error;  end if;

      if  not r'Valid then  raise Constraint_Error;  end if;
      if  r <= Zero   then  raise Constraint_Error;  end if;

      Arg :=  (Log_Gamma (r + Real(k))   -
               Log_Gamma (r)             -
               Log_Gamma (Real(k+1))     +
               Log (p)  * r              +
               Real (k) * Log_of_1_plus (-p));

      if Arg < -Two_to_the_Ninth then
         return Zero;
      end if;

      return Exp (Arg);

   end Neg_Binomial_Probability;

   ----------------------
   -- Get_Neg_Binomial --
   ----------------------

   -- Adapted from Fortran 77 code from the book:
   --     Dagpunar, J. 'Principles of random variate generation'
   --     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
   -- Translated from  Allan Miller's  Fortran 90.
   --
   --    Must have:  r > 0.0
   --       r = the `power' parameter of the negative binomial
   --    Must have:  0 < p < 1
   --       p = Bernoulli success probability
   --
   -- To make this routine conform to distribution returned by function
   -- Neg_Binomial_Probability above, (the formula given above is from the
   -- Wikipedia article on the negative binomial distribution), need to
   -- swap p and q = 1-p. ie, want r to be number of successes, not failures,
   -- with last trial being a success.

   procedure Get_Neg_Binomial
     (r       : in     Real;
      p       : in     Real;
      NB_Init : in out Neg_Binomial_Initializer;
      Stream  : in out State;
      Result  :    out Real)
   is
      -- uln = -ln (2**Machine_Emin) = -Real'Machine_Emin * ln (2).
      uln : constant Real := -Real (Real'Machine_Emin) * Natural_Log_of_2;

      h       : constant Real    := 0.85;
      p1      : constant Real    := One - p;            -- the old q
      q1      : constant Real    := p;                  -- the old p
      Trunc_r : constant Real    := Real'Truncation (r);
      mr      : constant Integer := Integer (Trunc_r);

      r1     : Real := r;
      k1, k2 : Real;
      Reciprocal_Log_q1 : Real;
      Reciprocal_Log_p1 : Real;
   begin
      Result := Zero;

      if r <= Zero then
         raise Constraint_Error;
      end if;

      if NB_Init.Uninitialized or else p /= NB_Init.p then
         if p <= Zero or p >= One then
            raise Constraint_Error;
         end if;
         Reciprocal_Log_q1 := One / Log (q1 + Min_Allowed_Real);
         Reciprocal_Log_p1 := One / Log (p1 + Min_Allowed_Real);

         NB_Init.p                 := p;
         NB_Init.Uninitialized     := False;
         NB_Init.Reciprocal_Log_q1 := Reciprocal_Log_q1;
         NB_Init.Reciprocal_Log_p1 := Reciprocal_Log_p1;
      end if;

      Reciprocal_Log_q1 := NB_Init.Reciprocal_Log_q1;
      Reciprocal_Log_p1 := NB_Init.Reciprocal_Log_p1;

      -- The 2nd part (k2 calculation) is too slow if p1 is too near 1.
      -- The following (Get_k1) is used as an optimization if p1 > 0.85 = h.
      -- Increment k1, Decrement r1:

      k1 := Zero;

      if p1 > h and then mr > 0 then

      Get_k1: declare
         Z, del_k : Real;
      begin
         for i in 1 .. mr loop
            Get_Non_Zero: loop
               Get_Random_Real (Z, Stream);
               exit Get_Non_Zero when Z > Zero;
            end loop Get_Non_Zero;
            del_k := Real'Truncation (Reciprocal_Log_p1 * Log (Z));
            k1    := k1 + del_k;
          --r1    := r1 - One;
          -- Moved outside loop as r1 := r1 - Trunc_r; assumes mr>0.
         end loop;
         r1 := r1 - Trunc_r;
	 -- if Trunc_r = Real (mr) = Real'Truncation (r), then r1 is 0 always.
	 -- When r1=0, the following calculation of k2 does nothing: k2=0.
      end Get_k1;
      end if;

      if r1 > Zero and then r1 > -uln * Reciprocal_Log_q1 then
         raise Constraint_Error; -- p1 too large (q1 near 0) for this value of r
      end if;


      k2 := Zero;

      Get_k2: declare
         y : Real := q1**r1;
         g : Real := r1;
         Z : Real;
      begin
         Get_Random_Real (Z, Stream);
         loop
            exit when Z < y;
            Z  := Z - y;
            k2 := k2 + One;
            y  := y * p1 * g / k2;
            g  := g + One;
         end loop;
      end Get_k2;

      Result := Real'Truncation (k1 + k2 + Half);

   end Get_Neg_Binomial;

   --------------------------
   -- Binomial_Probability --
   --------------------------

   --  [n!/(k!(n-k)!)] * p^k * (1-p)^(n-k)     if 0 <= k <= n;
   --  returns 0 otherwise.

   function Binomial_Probability
     (n : in Positive;
      k : in Integer;
      p : in Real)
      return Real
   is
      Arg : Real;
   begin
      if  not p'Valid then  raise Constraint_Error;  end if;
      if  p <= Zero   then  raise Constraint_Error;  end if;
      if  p >= One    then  raise Constraint_Error;  end if;

      if  k > n then  return Zero; end if;
      if  k < 0 then  return Zero; end if;

      Arg :=  (Log_Gamma (Real(n+1))   -
               Log_Gamma (Real(k+1))   -
               Log_Gamma (Real(n-k+1)) +
               Real (k)   * Log (p)    +
               Real (n-k) * Log_of_1_plus (-p));

      if Arg < -Two_to_the_Ninth then
         return Zero;
      end if;

      return Exp (Arg);

   end Binomial_Probability;

   --------------------
   -- Reset_Binomial --
   --------------------

   -- procedure Get_Binomial  must call this each time n or p is changed.
   --    p = Bernoulli success probability
   --           0 < p < 1
   --    n = number of bernoulli trials : Positive
   --           n >= 1

   procedure Reset_Binomial
     (n : in Positive;
      p : in Real;
      B_Init : out Binomial_Initializer)
   is
   begin
      if p >= One  then   raise Constraint_Error;   end if;
      if p <= Zero then   raise Constraint_Error;   end if;

      B_Init.n := n;
      B_Init.p := p;

      B_Init.r0 := Integer (Real'Truncation (Real (n+1) * p));
      if B_Init.r0 > n then B_Init.r0 := n; end if;

      B_Init.p_r           := Binomial_Probability (n, B_Init.r0, p);
      B_Init.Odds_Ratio    := p / (One - p);
      B_Init.Uninitialized := False;

   end Reset_Binomial;

   ------------------
   -- Get_Binomial --
   ------------------

   -- function generates a random binomial deviate using C.D.Kemp's method.
   --    p = Bernoulli success probability
   --    Must have 0 < p < 1.
   --    n = number of Bernoulli trials
   --    Must have n >= 1, (n integer valued).
   -- Reference: Kemp, C.D. (1986). `A modal method for generating binomial
   --            variables', Commun. Statist. - Theor. Meth. 15(3), 805-813.
   -- Author: Allan Miller (Fortran90 version)

   procedure Get_Binomial
     (n      : in     Positive;
      p      : in     Real;
      B_Init : in out Binomial_Initializer;
      Stream : in out State;
      Result : out    Real)
   is
      pragma Assert(Integer'Last >= 2**31-1); -- really should make n 64 bit.
      Odds_Ratio : Real;
      p_r        : Real;
      ru, rd, r0 : Integer;
      u, pd, pu  : Real;
   begin

      if B_Init.Uninitialized or p /= B_Init.p or n /= B_Init.n then
         Reset_Binomial (n, p, B_Init);
      end if;

      r0         := B_Init.r0;         -- r0 almost always = np = mean
      Odds_Ratio := B_Init.Odds_Ratio;
      p_r        := B_Init.p_r;

      Result := Real (r0); -- Init

      Get_Random_Real (u, Stream);
      u := u - p_r;
      if u < Zero then
        Result := Real (r0); return;
      end if;

      pu := p_r;
      ru := r0;
      pd := p_r;
      rd := r0;
      loop
         rd := rd - 1;
         if rd >= 0 then
           pd := pd * Real (rd+1) / (Odds_Ratio * Real (n-rd));
           u := u - pd;
           if u < Zero then
             Result := Real (rd); return;
           end if;
         end if;

         ru := ru + 1;
         if (ru <= n) then
           pu := pu * Real (n-ru+1) * Odds_Ratio / Real (ru);
           u := u - pu;
           if u < Zero then
             Result := Real (ru); return;
           end if;
         end if;
      end loop;

   end Get_Binomial;

   -------------------------
   -- Poisson_Probability --
   -------------------------

   --  The Poisson distribution (probability mass function):
   --
   --      f(k) = Mean^k * Exp (-Mean) / k!
   --
   --  Input includes 0.

   function Poisson_Probability
     (Mean : in Real;
      k    : in Integer)
      return Real
   is
      Arg : Real;
   begin
      if  not Mean'Valid then  raise Constraint_Error;  end if;
      if  Mean <= Zero   then  raise Constraint_Error;  end if;

      if  k < 0 then  return Zero; end if;

      Arg :=  ( Log (Mean) * Real (k)
               -Mean
               -Log_Gamma (Real(k+1)) );

      if Arg < -Two_to_the_Ninth then
         return Zero;
      end if;

      return Exp (Arg);

   end Poisson_Probability;

   ------------------------
   -- Get_Poisson_Tocher --
   ------------------------

   -- Tocher's algorithm for generating Poisson deviates.
   -- Slow for Mean ~ 8 or greater.
   -- Notice should think of it as single-precision, 30 bit.

   procedure Get_Poisson_Tocher
     (Mean   : in     Real;
      Stream : in out State;
      Result : out    Real)
   is
      S : Real := One;
      A, Random_Real_1, Random_Real_2 : Real;
      X, X1, X2 : Random_Int;
      Bits_Per_X : constant      := 30;
      Top_Shift  : constant      := Bits_per_Random_Number - Bits_Per_X;
      Scale_30   : constant Real := Two**(-Bits_Per_X);
   begin
      Result := Zero;
      if Mean <= Zero then raise Constraint_Error; end if;
      A := Exp (-Mean);
      loop
         Get_Random (X, Stream);
         X1 := (X mod 2**Bits_Per_X);            -- bottom 30 bits
         Random_Real_1 := Real (X1) * Scale_30;

         S := S * Random_Real_1;
         exit when S < A;
         Result := Result + One;

         X2 := (X / 2**Top_Shift);               -- top 30 bits
         Random_Real_2 := Real (X2) * Scale_30;

         S := S * Random_Real_2;
         exit when S < A;
         Result := Result + One;
      end loop;

   end Get_Poisson_Tocher;

   pragma Inline (Get_Poisson_Tocher);

   -----------------
   -- Get_Poisson --
   -----------------

   --  Use Binomial with small p if mean > 8, else use tocher's algorithm.

   procedure Get_Poisson
     (Mean   : in     Real;
      P_Init : in out Poisson_Initializer;
      Stream : in out State;
      Result : out    Real)
   is
      Reciprocal_p : constant := 2.0**(-p_Shift);
      p : constant := 2.0**(p_Shift);
      n  : Integer;
   begin

      --if Mean <= Zero   then   raise Constraint_Error;   end if;
      -- Checked by Get_Poisson_Tocher

      if Mean < Eight then

         Get_Poisson_Tocher (Mean, Stream, Result);

      else

         --  if Mean > 2.0**11 then   raise Constraint_Error;   end if;
         --  if n0 >= Real (Integer'Last) then  raise Constraint_Error;  end if;
         --  Let Integer () raise these in next step:

         n := Integer (Mean * Reciprocal_p);

         --  p := Mean / Real (n);
         --  waste of time: error is about the same if p off by this amount.
         --  In other words can use n-1, n, n+1, doesn't matter.

         Get_Binomial
           (n      => n,
            p      => p,
            B_Init => P_Init,
            Stream => Stream,
            Result => Result);

      end if;

   end Get_Poisson;

   ------------------------
   -- Normal_Probability --
   ------------------------

   --  The Normal probability density:
   --
   --   f(X) = A * Exp (-(X - Mean)**2 / (Two*Sigma**2))
   --
   --     A = 1.0 / (Sigma * Sqrt (2*Pi))

   function Normal_Probability
     (Mean  : in Real;   -- Mean    of random variable X
      Sigma : in Real;   -- Std Dev of random variable X
      X     : in Real)
      return Real
   is
      A, B, Arg : Real;
   begin
      if Sigma < Sqrt_of_Min_Real then  raise Constraint_Error;  end if;

      Arg := (X - Mean)**2 / (Two*Sigma**2);

      if Arg > Two_to_the_Ninth then 
        return Zero;
      end if;

      B := Exp (-Arg);

      A := Sigma * Sqrt_of_Two_Pi;
      return  B / A;

   end Normal_Probability;

   ------------------
   -- Get_Normal_2 --
   ------------------

   --  Marsaglia-Bray method.
   --
   --  Sigma = Standard_Deviation

   procedure Get_Normal_2
     (Mean   : in     Real;
      Sigma  : in     Real;
      Stream : in out State;
      X1, X2 :    out Real)
   is
      S, Z1, Z2, Coeff_S : Real;
   begin
      loop
         Get_Random_Real (Z1, Stream);
         Get_Random_Real (Z2, Stream);
         Z1 := Two * Z1 - One;
         Z2 := Two * Z2 - One;
         S  := Z1*Z1 + Z2*Z2;
         exit when S < One and then S > Zero;
      end loop;

      Coeff_S := Sigma * Sqrt (-Two * Log (S) / S);
      X1      := Mean + Z1 * Coeff_S;
      X2      := Mean + Z2 * Coeff_S;

   end Get_Normal_2;

   pragma Inline (Get_Normal_2);

   ----------------
   -- Get_Normal --
   ----------------

   --  Marsaglia-Bray method.
   --
   --  Sigma = Standard_Deviation

   procedure Get_Normal
     (Mean   : in     Real;
      Sigma  : in     Real;
      N_Init : in out Normal_Initializer;
      Stream : in out State;
      Result :    out Real)
   is
      X1, X2 : Real;
   begin
      if N_Init.Uninitialized or else
         Mean  /= N_Init.Mean or else
         Sigma /= N_Init.Sigma
      then
         N_Init.Mean  := Mean;
         N_Init.Sigma := Sigma;

         Get_Normal_2 (Mean, Sigma, Stream, X1, X2);

         Result    := X1;               -- return X1
         N_Init.X2 := X2;               -- store X2
         N_Init.Uninitialized := False; -- Store of X's not empty.
      else
         Result    := N_Init.X2;        -- return X2
         N_Init.Uninitialized := True;  -- have depleted store of X's
      end if;

   end Get_Normal;

   ----------------------------
   -- Log_Normal_Probability --
   ----------------------------

   --  The Log_Normal probability density:
   --
   --   f(X) = Exp (-(Log(X) - Mean_Z)**2 / (2*Sigma_Z**2)) / (A * X)
   --
   --   f(X) = 0 (for X <= 0)
   --
   --     A = Sigma_Z * Sqrt (2*Pi)
   --     Z = Log(X)

   function Log_Normal_Probability
     (Mean_Z  : in Real;   -- Mean    of random variable Z = Log (X)
      Sigma_Z : in Real;   -- Std Dev of random variable Z = Log (X)
      X       : in Real)
      return Real
   is
      A, B : Real;
   begin
      if       X < Sqrt_of_Min_Real then  return Zero;             end if;
      if Sigma_Z < Sqrt_of_Min_Real then  raise Constraint_Error;  end if;

      B := Exp (-(Log(X) - Mean_Z)**2 / (Two*Sigma_Z**2));
      A := Sigma_Z * Sqrt_of_Two_Pi;

      return  B / (A * X);

   end Log_Normal_Probability;

   --------------------
   -- Get_Log_Normal --
   --------------------

   procedure Get_Log_Normal    -- outputs random variable X
     (Mean_Z  : in     Real;   -- Mean    of random variable Z = Log (X)
      Sigma_Z : in     Real;   -- Std Dev of random variable Z = Log (X)
      LN_Init : in out Log_Normal_Initializer;
      Stream  : in out State;
      Result  :    out Real)   -- X
   is
      Z1, Z2 : Real;
   begin
      if LN_Init.Uninitialized    or else
         Mean_Z  /= LN_Init.Mean  or else
         Sigma_Z /= LN_Init.Sigma
      then
         LN_Init.Mean  := Mean_Z;
         LN_Init.Sigma := Sigma_Z;

         Get_Normal_2 (Mean_Z, Sigma_Z, Stream, Z1, Z2);

         Result     := Exp (Z1);         -- return X1
         LN_Init.X2 := Exp (Z2);         -- store X2
         LN_Init.Uninitialized := False; -- Store of X's not empty.
      else
         Result     := LN_Init.X2;       -- return X2
         LN_Init.Uninitialized := True;  -- have depleted store of X's
      end if;
   end Get_Log_Normal;

   ----------------
   -- Get_Cauchy --
   ----------------

   --  The Cauchy (Lorentzian) distribution probability density:
   --
   --    f(X) = A / [(A*A + X*X) * Pi] -- normalized and scaled
   --
   --  Generates random deviates X in range (-inf, inf).

   procedure Get_Cauchy
     (A      : in     Real;
      Stream : in out State;
      Result :    out Real)
   is
      S, X1, X2 : Real;
   begin
      Result := Zero;
      loop
         Get_Random_Real (X1, Stream);
         Get_Random_Real (X2, Stream);
         X1 := Two * X1 - One;
         X2 := Two * X2 - One;
         S  := X1*X1 + X2*X2;
         exit when S < One and then X2 > Min_Allowed_Real;
	 -- X2=0 is extremely rare.
      end loop;

      Result := A * X1 / X2;

   end Get_Cauchy;

   ---------------------------
   -- Student_t_Probability --
   ---------------------------

   --  Gamma((m+1)/2) * (1 + x*x/m)^(-(m+1)/2) / [Sqrt (m*Pi) * Gamma(m/2)]

   function Student_t_Probability
     (m : in Positive;
      x : in Real)
      return Real
   is
      r    : constant Real := Half * Real (m + 1);
      s    : constant Real := Half * Real (m);
      m_Pi : constant Real := Real (m) * 3.141_592_653_589_793_238_463;
      Arg, Result : Real := Zero;
   begin
      Arg := (Log_Gamma (r)
             -Log_Gamma (s)
             -Log (m_Pi) * Half
         -r * Log_of_1_plus (x*x / Real(m)));

      if Arg < -Two_to_the_Ninth then
         Result := Zero;
      else
         Result := Exp (Arg);
      end if;

      return Result;

   end Student_t_Probability;

   -------------------
   -- Get_Student_t --
   -------------------

   -- generates a random deviate from a
   -- t distribution using Kinderman and Monahan's ratio method.
   --
   -- Must have m >= 1.
   -- m = degrees of freedom of distribution;
   --
   -- Adapted from Fortran 77 code from the book:
   -- Dagpunar, J. 'Principles of random variate generation'
   -- Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
   -- Translated from Allan Miller's Fortran 90 version.

   procedure Get_Student_t
     (m                  : in     Positive;
      Student_t_Init     : in out Student_t_Initializer;
      Stream             : in out State;
      Result             : out    Real)
   is
      s, c, a, f, g : Real;
      r, x, v : Real;
   begin

      s := Real (m);
      c := -Quarter * (s + One);

      if Student_t_Init.Uninitialized or Student_t_Init.m /= m then
         a := Four / (One + One / s)**c;
         f := Sixteen / a;
         if m > 1 then
            g := s - One;
            g := ((s + One) / g)**c * Sqrt ((s + s) / g);
         else
            g := One;
         end if;
         Student_t_Init.a := a;
         Student_t_Init.f := f;
         Student_t_Init.g := g;
         Student_t_Init.m := m;
      else
         a := Student_t_Init.a;
         f := Student_t_Init.f;
         g := Student_t_Init.g;
      end if;

      Student_t_Init.Uninitialized := False;

      Keep_Trying: loop
         This_Attempt: loop
            Get_Random_Real (r, Stream);
            exit This_Attempt when r <= Zero;
            Get_Random_Real (v, Stream);
            x := (Two*v - One) * g / r;
            v := x*x;
            if v > (Five - a*r) then
               exit This_Attempt when f < r*(v + Three);
               exit This_Attempt when r > (One + v/s)**c;
            end if;
            exit Keep_Trying; -- if we got this far we're done.
         end loop This_Attempt;
      end loop Keep_Trying;

      Result := x;

   end Get_Student_t;

   ----------------------
   -- Beta_Probability --
   ----------------------

   --  f(x) = x**(aa-1) * (1-x)**(bb-1) *
   --                       Gamma (aa + bb) / (Gamma (bb)*Gamma (aa))
   --  must have 1 > x > 0.

   function Beta_Probability
     (aa, bb : in Real;
      x : in Real)
      return Real
   is
      Arg, Result : Real := Zero;
   begin

      if x <= Zero then   raise Constraint_Error;   end if;
      if x >= One  then   raise Constraint_Error;   end if;

      Arg :=  Log_Gamma (aa + bb)
             -Log_Gamma (aa)
             -Log_Gamma (bb)
             +Log (x) * (aa - One)
             +Log_of_1_plus (-x) * (bb - One);

      if Arg < -Two_to_the_Ninth then
         Result := Zero;
      else
         Result := Exp (Arg);
      end if;

      return Result;

   end Beta_Probability;

   --------------
   -- Get_Beta --
   --------------

   -- Adapted from Fortran 77 code from the book:
   --     Dagpunar, J. 'Principles of random variate generation'
   --     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
   -- translated from Alan Miller's Fortran 90 version.
   --
   -- Get_Beta generates a random deviate in [0,1] from a beta distribution
   -- with density proportional to x**(aa-1) * (1-x)**(bb-1), using
   -- cheng's log logistic method.
   --
   -- Must have aa > 0, bb > 0.
   --
   --  f(x) = x**(aa-1) * (1-x)**(bb-1) *
   --                       Gamma (aa + bb) / (Gamma (bb)*Gamma (aa))
   --

   procedure Get_Beta
     (aa        : in     Real;
      bb        : in     Real;
      Beta_Init : in out Beta_Initializer;
      Stream    : in out State;
      Result    :    out Real)
   is
      log_of_4 : constant Real := 1.38629436111989061883; --1.38629436112 is gd
      a, b, g, r, s, x, y, z, test2 : Real;
      d, f, h, t, c : Real;
      Swap          : Boolean;
   begin
      if (aa <= Zero) or (bb <= Zero) then
         raise Constraint_Error;
      end if;

      if Beta_Init.Alpha /= aa then   Beta_Init.Uninitialized := True;   end if;
      if Beta_Init.Beta  /= bb then   Beta_Init.Uninitialized := True;   end if;

      if Beta_Init.Uninitialized then
         a := aa;
         b := bb;
         Swap := b > a;
         if Swap then
            g := b;
            b := a;
            a := g;
         end if;
         d := a / b;
         f := a + b;
         if (b > One) then
            h := Sqrt ((Two*a*b - f) / (f - Two));
            t := One;
         else
            h := b;
            t := One / (One + (a / (Max_Allowed_Real*b))**b);
         end if;
         c := a + h;
         Beta_Init.Alpha := aa;
         Beta_Init.Beta  := bb;
	 Beta_Init.d := d;
	 Beta_Init.f := f;
	 Beta_Init.h := h;
	 Beta_Init.t := t;
	 Beta_Init.c := c;
	 Beta_Init.Swap := Swap;
         Beta_Init.Uninitialized := False;
      else
         d := Beta_Init.d;
         f := Beta_Init.f;
         h := Beta_Init.h;
         t := Beta_Init.t;
         c := Beta_Init.c;
         Swap := Beta_Init.Swap;
      end if;

      Keep_Trying :loop
         This_Attempt: loop
            Get_Random_Real (r, Stream);
            Get_Random_Real (x, Stream);
            s := r*r*x;
            exit This_Attempt when (r < Min_Allowed_Real or else s <= Zero);
            if r < t then
               x := Log (r / (One - r)) / h;
               y := d*Exp (x);
               z := c*x + f*Log ((One + d) / (One + y)) - log_of_4;
               if (s - One) > z then
                  exit This_Attempt when (s - s*z) > One;
                  exit This_Attempt when  Log (s) > z;
               end if;
               Result := y / (One + y);
            else
               Test2 := (One + One / d)**f;
               exit This_Attempt when (Four*s > Test2);
               Result := One;
            end if;
            exit Keep_Trying;   -- We are done. Only way out.
         end loop This_Attempt; -- Make another attempt.
      end loop Keep_Trying;

      if Swap then
         Result := One - Result;
      end if;

   end Get_Beta;

   -----------------------------
   -- Chi_Squared_Probability --
   -----------------------------

   --  f(x) = (1 / 2) * (x / 2)**(s-1) * Exp (-x / 2)  / Gamma (s)
   --
   --   (where s = 0.5 * Degrees_of_Freedom).

   function Chi_Squared_Probability
     (Degrees_of_Freedom : in Real;
      x : in Real)
      return Real
   is
      s : constant Real := Half * Degrees_of_Freedom;
      Arg, Result : Real := Zero;
   begin

      if s <= Zero then   raise Constraint_Error;   end if;
      if x <= Zero then
         return Zero;
      end if;

      Arg :=  Log (x * Half) * (s - One)
             -x * Half
             -Log_Gamma (s);

      if Arg < -Two_to_the_Ninth then
         Result := Zero;
      else
         Result := Half * Exp (Arg);
      end if;

      return Result;

   end Chi_Squared_Probability;

   -----------------------
   -- Gamma_Probability --
   -----------------------

   --  f(x) = x**(s-1) * Exp (-x)  / Gamma (s)

   function Gamma_Probability
     (s : in Real;
      x : in Real)
      return Real
   is
      Arg, Result : Real := Zero;
   begin

      if s <= Zero then   raise Constraint_Error;   end if;
      if x <= Zero then
         return Zero;
      end if;

      Arg :=  Log (x) * (s - One)
             -x
             -Log_Gamma (s);

      if Arg < -Two_to_the_Ninth then
         Result := Zero;
      else
         Result := Exp (Arg);
      end if;

      return Result;

   end Gamma_Probability;

   -----------------
   -- Get_Gamma_1 --
   -----------------

   -- Uses the algorithm in
   -- Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
   -- gamma variables', Trans. on Math. Software (TOMS), vol.26(3), pp.363-372.
   -- translated from Alan Miller's Fortran 90.
   --
   --  Generates a random gamma deviate for shape parameter s >= 1.
   --  s := Shape parameter of Gamma distribution
   --  Must have s > 0.0
   --  Must have s < 1.0

   procedure Get_Gamma_1
     (s           : in     Real;
      Stream      : in out State;
      Result      :    out Real)
   is
      c, d : Real;
      u, v, x, y : Real;
   begin

      d := s - One / Three;
      c := One / (Three * Sqrt (d));

      loop
         loop
            Get_Normal_2 (Zero, One, Stream, x, y); -- Mean=0, Sigma=1
            v := (One + c*x)**3;
            exit when v > Zero;
            v := (One + c*y)**3; -- 1st try failed
            exit when v > Zero;
         end loop;

         Get_Random_Real (u, Stream);
         if u < One - 0.0331 * x**4 then
            Result := d*v;
            exit;
         elsif Log (u) < Half * x**2 + d * (One - v + Log (v)) then
            Result := d*v;
            exit;
         end if;
      end loop;

   end Get_Gamma_1;

   -----------------
   -- Get_Gamma_2 --
   -----------------

   -- Adapted from Fortran 77 code from the book:
   --     Dagpunar, J. 'Principles of random variate generation'
   --     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
   --
   --  Generates a random deviate in [0,infinity) from a gamma distribution.
   --  s = Shape parameter of distribution Gamma.
   --
   --  s := Shape parameter of Gamma distribution
   --  Must have s > 0.0
   --  Must have s < 1.0


   procedure Get_Gamma_2
     (s           : in     Real;
      Gamma_Init  : in out Gamma_Initializer;
      Stream      : in out State;
      Result      :    out Real)
   is
      a, r, x, w : Real;
      p, c, d, uf, vr : Real;
   begin

      if (s <= Zero or s >= One) then
         raise Constraint_Error;
      end if;

      a := One - s;

      if Gamma_Init.Uninitialized or Gamma_Init.s /= s then
         p := a / (a + s*Exp (-a));
         if s < Min_Allowed_Real then
            raise Constraint_Error;
         end if;
         c  := One / s;
         uf := p * (Min_Allowed_Real / a)**s;
         vr := One - Min_Allowed_Real;
         d  := a * Log (a);
         Gamma_Init.s  := s;
         Gamma_Init.p  := p;
         Gamma_Init.c  := c;
         Gamma_Init.d  := d;
         Gamma_Init.uf := uf;
         Gamma_Init.vr := vr;
         Gamma_Init.Uninitialized := False;
      else
         p  := Gamma_Init.p;
         c  := Gamma_Init.c;
         d  := Gamma_Init.d;
         uf := Gamma_Init.uf;
         vr := Gamma_Init.vr;
      end if;

      Keep_Trying: loop
         This_Attempt: loop
            Get_Random_Real (r, Stream);
            if (r >= vr) then
               exit This_Attempt;
            elsif (r > p) then
               x := a - Log ((One - r) / (One - p));
               w := a * Log (x) - d;
            elsif (r > uf) then
               x := a * (r / p)**c;
               w := x;
            else
               Result := Zero;
               return;
            end if;

            Get_Random_Real (r, Stream);
            if (One-r <= w and r > Zero) then
               exit This_Attempt when   r*(w + One) >= One;
               exit This_Attempt when  -Log(r) <= w;
            end if;
            exit Keep_Trying;
         end loop This_Attempt;
      end loop Keep_Trying;

      Result := x;

   end Get_Gamma_2;

   ---------------
   -- Get_Gamma --
   ---------------

   -- Adapted from Fortran 77 code from the book:
   -- Dagpunar, J. 'Principles of random variate generation'
   -- Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
   --
   --  s := Shape parameter of Gamma distribution
   --  Must have s > 0.

   procedure Get_Gamma
     (s          : in     Real;
      Gamma_Init : in out Gamma_Initializer;
      Stream     : in out State;
      Result     :    out Real)
   is
   begin
      if s <= Zero then
         raise Constraint_Error;
      end if;

      if s > One then
         Get_Gamma_1 (s, Stream, Result);
      elsif s < One then
         Get_Gamma_2 (s, Gamma_Init, Stream, Result);
      else
         Get_Exponential (One, Stream, Result); -- Mean = 1
      end if;

   end Get_Gamma;

   ---------------------
   -- Get_Chi_Squared --
   ---------------------

   --  Generates a random variate from the chi-squared distribution

   procedure Get_Chi_Squared
     (Degrees_of_Freedom : in     Real;
      Chi_Init           : in out Chi_Initializer; -- subtype of Gamma_Initializer
      Stream             : in out State;
      Result             :    out Real)
   is
      X: Real;
   begin
      if Degrees_of_Freedom <= Zero then raise Constraint_Error; end if;

      Get_Gamma (Half*Degrees_of_Freedom, Chi_Init, Stream, X);
      Result := Two * X;

   end Get_Chi_Squared;

   --------------
   -- LU_Solve --
   --------------

   --  Solve for X in:  L*X = b

   procedure LU_Solve
     (L : in     Matrix; -- L is lower triangle produced by Choleski decomp. 
      b : in     Vector;
      X :    out Vector) -- solve for this
   is
      Index_First : constant Natural := b'First; -- Index of Matrix is Pos.
      Index_Last  : constant Natural := b'Last;
      Z : Vector (Index_First..Index_Last);
      Sum : Real;
   begin
      if L'Length(1) <= 1          then  raise Constraint_Error;  end if;
      if L'First(2) /= Index_First then  raise Constraint_Error;  end if;
      if L'Last(2)  /= Index_Last  then  raise Constraint_Error;  end if;

      -- U' = L
      -- Solve for Z in  L*Z = b.
 
      Z(Index_First) := b(Index_First) / L(Index_First,Index_First);
      for Row in Index_First+1 .. Index_Last loop
         Sum := Zero;
         for Col in Index_First .. Row-1 loop
            Sum := Sum + L(Row,Col) * Z(Col);
         end loop;
         Z(Row) := (b(Row) - Sum) / L(Row,Row);
      end loop; --Row

      -- Solve for X in  U*X = Z.
 
      X(Index_Last) := Z(Index_Last) / L(Index_Last,Index_Last);
      for Row in reverse Index_First .. Index_Last-1 loop
         Sum := Zero;
         for Col in Row+1 .. Index_Last loop
            Sum := Sum + L(Col,Row) * X(Col); -- U = Tr(L) = Tr(L)
         end loop;
         X(Row) := (Z(Row) - Sum) / L(Row,Row);
      end loop; --Row

   end LU_Solve;

   ------------------------
   -- Choleski_Decompose --
   ------------------------

   --  At present adds a small number to diagonal to reduce unecessary failures.
   --  Result is it will not calculate inverse to full precision (15 digits).
   --  More like 13 digits w/ present setting.

   procedure Choleski_Decompose
     (Covariance       : in     Matrix;
      LU_of_Covariance :    out Matrix) -- Choleski Decomp of Covariance matrix.
   is
      A : Matrix renames LU_of_Covariance;
      Index_First : constant Natural := Covariance'First(1); -- Index of Matrix is Pos.
      Index_Last  : constant Natural := Covariance'Last(1);
      Reciprocal_Sqrt_Pivot, Pivot, Sum : Real;
      Max_Diag_Element : Real := Zero;
      Safety_Shift : Real;
      Safety_Fraction : constant Real := Two**(-(Real'Machine_Mantissa * 5) / 6);
      --  here Safety_Shift is about 10^(-13)
   begin
      if Covariance'Length(1) <= 1          then  raise Constraint_Error;  end if;
      if Covariance'First(2) /= Index_First then  raise Constraint_Error;  end if;
      if Covariance'Last(2)  /= Index_Last  then  raise Constraint_Error;  end if;

      A := Covariance;

      --  Positive semi-definite matrices have non-neg diag elements.

      --  If they're < 0, something went wrong only the client knows how to fix:

      for i in Index_First .. Index_Last loop
         if A(i,i) < Zero then  raise Constraint_Error;  end if;
      end loop;

      --  Otherwise, add a constant. First find max element of diagonal:

      Max_Diag_Element := Zero;
      for i in Index_First .. Index_Last loop
         if A(i,i) > Max_Diag_Element then Max_Diag_Element := A(i,i); end if;
      end loop;

      Safety_Shift := Max_Diag_Element * Safety_Fraction;
      for i in Index_First .. Index_Last loop
         A(i,i) := A(i,i) + Safety_Shift;
      end loop;

      for i in Index_First .. Index_Last loop
         Sum := Zero;
         for Col in Index_First .. i-1 loop
            Sum := Sum + A(i,Col) * A(i,Col);
         end loop;
         Pivot := A(i,i) - Sum;

         if (Pivot <= Zero) then  raise Constraint_Error;  end if;
         --  NOT Positive definite.

         A(i,i)                := Sqrt (Pivot);
         Reciprocal_Sqrt_Pivot := One / A(i,i);
         for Row in i+1 .. Index_Last loop
            Sum := Zero;
            for Col in Index_First .. i-1 loop
               Sum := Sum + A(i,Col) * A(Row,Col);
            end loop; -- Col
            A(Row,i) := (A(Row,i) - Sum) * Reciprocal_Sqrt_Pivot;
         end loop; -- Row
      end loop; -- i

   end Choleski_Decompose;

   -----------------------------
   -- Get_Multivariate_Normal --
   -----------------------------

   -- Adapted from Fortran 77 code from the book:
   --     Dagpunar, J. 'Principles of random variate generation'
   --     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
   --

   procedure Get_Multivariate_Normal
     (Mean             : in     Vector;
      LU_of_Covariance : in     Matrix;
      MV_Init          : in out MV_Normal_Initializer;
      Stream           : in out State;
      Result           :    out Vector)
   is
      Z : Real;
      Index_First : constant Natural := Mean'First; -- Index of Mat must be Positive.
      Index_Last  : constant Natural := Mean'Last;
   begin

      if Mean'Length < 1 then  raise Constraint_Error;  end if;
      if Index_First /= LU_of_Covariance'First(1) then  raise Constraint_Error;  end if;
      if Index_Last  /= LU_of_Covariance'Last(1)  then  raise Constraint_Error;  end if;

      Result := Mean;

      for Col in LU_of_Covariance'Range(2) loop
         Get_Normal (Zero, One, MV_Init, Stream, Z); --mean=0,sigma=1
         for Row in Col .. Index_Last loop
            Result(Row) := Result(Row) + LU_of_Covariance(Row,Col) * Z;
         end loop;
      end loop;

   end Get_Multivariate_Normal;

   -------------------------------------
   -- Multivariate_Normal_Probability --
   -------------------------------------

   --   f(X) = A * Exp (-(X - Mean)*Covariance_Mat_Inverse*(X - Mean) / 2)
   --
   --     A = 1.0 / Sqrt (Determinant(Covariance_Mat_Inverse) * (2*Pi)^N)
   --
   --   where N = number of variables (length of vector X).

   function Multivariate_Normal_Probability
     (Mean             : in Vector;   -- Means of random variables X
      LU_of_Covariance : in Matrix;   -- L of LU decomp of Covariance matrix
      X                : in Vector)
      return Real
   is
      Index_First : constant Natural := Mean'First;
      Index_Last  : constant Natural := Mean'Last;--Index of Matrix is type Positive.
      Y, Solution : Vector (Index_First..Index_Last);
      A, B, Arg, Sum, Prod, Log_Sqrt_Determinant : Real;
   begin
      Prod := One;
      for i in Index_First .. Index_Last loop
         Prod := Prod * LU_of_Covariance(i,i);
      end loop;
      Log_Sqrt_Determinant := Log (Prod + Min_Allowed_Real);

      for i in Index_First .. Index_Last loop
         Y(i) := X(i) - Mean (i);
      end loop;

      LU_Solve (LU_of_Covariance, Y, Solution);

      Sum := Zero;
      for i in Index_First .. Index_Last loop
         Sum := Sum + Y(i) * Solution (i);
      end loop;

      Arg :=  Half * Sum
             +Log_Sqrt_Determinant;

      if Arg > Two_to_the_Ninth then 
        return Zero;
      end if;

      B := Exp (-Arg);

      A := Sqrt_of_Two_Pi ** Integer (X'Length);

      return  B / A;

   end Multivariate_Normal_Probability;

   ------------
   -- Invert --
   ------------
  
   --  Get Inverse of Matrix M.  M must be positive definite.
   
   procedure Invert 
     (M              : in     Matrix;
      M_Inv          :    out Matrix;
      Max_Error      :    out Real)
   is 
      Index_Last      : constant Positive := M'Last(1);
      Index_First     : constant Positive := M'First(1);
      Solution_Vector : Vector(Index_First..Index_Last);
      Unit_Vector : Vector(Index_First..Index_Last) := (others => 0.0);
      Product     : Vector(Index_First..Index_Last) := (others => 0.0);
      Error       : Vector(Index_First..Index_Last) := (others => 0.0);
      M_LU        : Matrix(Index_First..Index_Last,Index_First..Index_Last);
      Sum : Real;
   begin
      Max_Error := 0.0;

      Choleski_Decompose 
        (Covariance       => M,
         LU_of_Covariance => M_LU);
   
      for I in Index_First .. Index_Last loop

         if I > Index_First then 
            Unit_Vector(I-1) := 0.0;
         end if;
         Unit_Vector(I) := 1.0;

         --  Solve  A*X = B:

         LU_Solve
	   (L => M_LU,
            B => Unit_Vector,
            X => Solution_Vector);

         -- get Error =  Unit_Vector - M*Solution_Vector

         -- Matrix vector product:
         for Row in Index_First .. Index_Last  loop
            Sum := Zero;
            for k in Index_First .. Index_Last  loop
               Sum := Sum + M(Row, k) * Solution_Vector(k);
            end loop;
	    Product(Row) := Sum;
         end loop;

         for Row in Index_First .. Index_Last  loop
            Error(Row) := Unit_Vector(Row) - Product(Row);
         end loop;

         for I in Index_First..Index_Last loop
            if Abs(Error(I)) > Max_Error then
               Max_Error := Abs(Error(I));
            end if;
         end loop;

         -- Solution vector is the I-th column of M_Inverse:

         for J in Index_First .. Index_Last loop
            M_Inv (J, I) := Solution_Vector(J);
         end loop;

      end loop;
   end Invert;

   -------------------
   -- Test_Choleski --
   -------------------
  
   procedure Test_Choleski 
   is 
      Index_First : constant := 1;
      Index_Last  : constant := 5;
      Test  : Matrix(Index_First .. Index_Last, Index_First .. Index_Last);
      T_Inv : Matrix(Index_First .. Index_Last, Index_First .. Index_Last);
      Max_Error, X : Real;
   begin
      X := 0.1234;
      for i in Index_First .. Index_Last loop
      for j in Index_First .. Index_Last loop
         X := X + 1.2345;
         Test(i,j) := X;
      end loop;
      end loop;

      for i in Index_First .. Index_Last loop
      for j in i .. Index_Last loop
         Test(i,j) := Test(j,i);
      end loop;
      end loop;

      for i in Index_First .. Index_Last loop
         Test(i,i) := Test(i,i) + 20.0 * Real (i);
      end loop;

      Invert 
        (M         => Test,
         M_Inv     => T_Inv,
         Max_Error => Max_Error);

      if Max_Error > 1.0e-8 then
        raise Constraint_Error;
      end if;

   end Test_Choleski;

end Disorderly.Random.Deviates;

