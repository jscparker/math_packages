
-------------------------------------------------------------------------------
-- package body Chi_Gaussian_CDF, CDF of Normal and Chi-square distributions
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

with Ada.Numerics.Generic_Elementary_Functions;
with Gamma;

package body Chi_Gaussian_CDF is

  package Math is new Ada.Numerics.Generic_Elementary_Functions(Real); use Math;
  package Gam is new Gamma(Real); use Gam; --provides log_gamma; for reals to 32 digits

  Real_Epsilon : constant Real := Real'Epsilon * 0.5;    -- 4.4*e-16 for ieee
  Max_Log      : constant Real := 666.0; 

  ---------------------
  -- Chi_Squared_CDF --
  ---------------------

  -- Cumulative distribution function for Chi^2 distribution. 

  function Chi_Squared_CDF 
   (True_Degrees_Freedom : Real; Chi_Squared : Real) 
    return Real is

    a : constant Real := True_Degrees_Freedom * 0.5;
    x : constant Real := Chi_Squared * 0.5;

  begin

    if not x'Valid then  raise Constraint_Error;  end if;
    --  In some cases, input of (say) inf or NaN will make numeric programs
    --  hang rather than crash .. very difficult to diagnose, so this seems
    --  best policy for function calls that are rarely found in time-sensitive
    --  inner loops

    return Incomplete_Gamma (a, x);

  end Chi_Squared_CDF;

  function Normal_CDF_x (x : in Real) return Real;

  -----------------
  -- Normal_CDF --
  -----------------

  -- Cumulative distribution function for Standard Normal Distribution. 
  --
  --   Normal (v) = exp(-v*v/2)  / Sqrt(2 pi)
  --
  -- Want to integrate this from -inf to x using the Incomplete Gamma function:
  --
  --   Inc_Gamma(a, w) = Int{from 0 to w} [exp(-t) t**(a-1) dt]  / Gamma(a)
  --
  --   Set a = 0.5.  Use Gamma(0.5) = Sqrt(pi) and let t = u*u/2:
  --
  --   Inc_Gamma(0.5, w) = Int{from 0 to sqrt(2w)} [2 exp(-v*v/2) dv] / Sqrt(2 pi)
  --
  --   Inc_Gamma(0.5, x*x/2) / 2 = Int{from 0 to x)} [exp(-v*v/2) dv] / Sqrt(2 pi)
  --
  --

  function Normal_CDF 
   (x : Real) 
    return Real is

    a : constant Real := 0.5;

  begin

    if not x'Valid then  raise Constraint_Error;  end if;
    --  In some cases, input of (say) inf or NaN will make numeric programs
    --  hang rather than crash .. very difficult to diagnose, so this seems
    --  best policy for function calls that are rarely found in time-sensitive
    --  inner loops.

    if x > 0.0 then
      return 0.5 * (1.0 + Incomplete_Gamma (a, x*x*0.5));
    elsif x <= -4.5 then
      return Normal_CDF_x (x); -- get better accuracy here out on tail.
    else
      return 0.5 * (1.0 - Incomplete_Gamma (a, x*x*0.5));
    end if;

  end Normal_CDF;


  ------------------------
  -- Incomplete_Gamma_C --
  ------------------------

  -- This is the complemented Incomplete_Gamma function:
  --
  -- Incomplete_Gamma_C (a,x) 
  --      = Integral {t=x to inf} [exp(-t) * t**(a-1) * dt]  / Gamma(a)
  --
  -- Notice that because the integrals are nomalized by 1 / Gamma(a):
  --
  --      Incomplete_Gamma (a,x) + Incomplete_Gamma_C (a,x) = 1.0
  --
  -- Uses Gauss' (c. 1813) continued fraction (see wikipedia for inc. gamma)
  --

  function Incomplete_Gamma_C (a : Real; x : Real) return Real is

    Exagam, Arg : Real;
    C_Fraction, N, g_2, g_1 : Real;
    r, t : Real;
    P, P_1, P_2 : Real;
    Q, Q_1, Q_2 : Real;

    Max_Allowed_Val     : constant Real := 2.0**48;
    Reciprocal_of_Max_Allowed_Val : constant Real := 2.0**(-48);
    Min_Allowed_Real  : constant Real := 2.0**(Real'Machine_Emin / 2);
  begin

    if not x'Valid then  raise Constraint_Error;  end if;

    if (x <= 0.0) or (a <= 0.0) then
      return  1.0;
    end if;
 
    -- Use series solution for small x:

    if (x < 1.0) or (x < a) then
      return  1.0 - Incomplete_Gamma(a,x);
    end if;
  
    -- Exagam :=  x**a * Exp(-x) / Gamma(a):
  
    Arg := a * Log (x) - x - Log_Gamma (a); -- notice never gets big > 0 if x>0
    if (Arg < -Max_Log) then
      return 0.0;
    else
      Exagam := Exp (Arg);
    end if;


    N   := 0.0;
    g_1 := x - a + 2.0;

    P_2 := 1.0;        --P(-2)
    Q_2 := x;          --Q(-2)

    P_1 := x + 1.0;    --P(-1)
    Q_1 := Q_2 * g_1;  --Q(-1)

    C_Fraction := P_1 / Q_1;
    
    Continued_Fraction_Evaluation: loop
      N   := N + 1.0;
      g_1 := g_1 + 2.0;
      g_2 := (N + 1.0 - a) * N;
  
      P := P_1 * g_1  -  P_2 * g_2;
      Q := Q_1 * g_1  -  Q_2 * g_2;

      if (Abs Q > Min_Allowed_Real) then
        r := P / Q;
        t := Abs ((C_Fraction - r) / r);
        C_Fraction := r;
      else
        t := 1.0;
      end if;
 
      -- scale P's and Q's identically with 2.0**n if P gets too large.  
      -- Final answer is P/Q.

      P_2 := P_1;
      P_1 := P;
      Q_2 := Q_1;
      Q_1 := Q;
      if (Abs (P) > Max_Allowed_Val) then
        P_2 := P_2 * Reciprocal_of_Max_Allowed_Val;
        P_1 := P_1 * Reciprocal_of_Max_Allowed_Val;
        Q_2 := Q_2 * Reciprocal_of_Max_Allowed_Val;
        Q_1 := Q_1 * Reciprocal_of_Max_Allowed_Val;
      end if;
  
      if (t < Real_Epsilon) then   exit Continued_Fraction_Evaluation;   end if;
  
    end loop Continued_Fraction_Evaluation;

    return C_Fraction * Exagam;

  end Incomplete_Gamma_C;

  ----------------------
  -- Incomplete_Gamma --
  ----------------------

  -- Incomplete_Gamma (a,x) = Integral {t=0 to x} [exp(-t) * t**(a-1) * dt]  / Gamma(a)
  --
  -- Notice as x -> inf, the Integral -> Gamma(a), and Incomplete_Gamma (a,x) -> 1.0
  --
  -- Abramowitz, Stegun 6.5.29:
  --
  --       = Exp(-x) * x**a * Sum {k=0 to inf} [x**k / Gamma(a+k+1)] 
  --
  --    use Gamma(a+1) = a * Gamma(a):
  --
  --       = Exp(-x) * x**a * Sum {k=0 to inf} [a! * (-x)**k / (a+k)!] / Gamma(a)
  --
  --       = Exp(-x) * x**a * [1/a + x/(a(a+1)) + x^2/(a(a+1)(a+2)) + ...] / Gamma(a)
  --

  function Incomplete_Gamma (a : Real;  x : Real)  return Real is

    Sum, Exagam, Arg, Next_Term, a_plus_n : Real;

  begin

    if not x'Valid then  raise Constraint_Error;  end if;

    if (x <= 0.0) or (a <= 0.0) then
      return 0.0;
    end if;
  
    if not ((x < 1.0) or (x < a)) then
      return (1.0 - Incomplete_Gamma_C (a,x));
    end if;

    -- Exagam :=  x**a * Exp(-x) / Gamma(a):

    Arg := a * Log (x) - x - Log_Gamma (a);
    if Arg < -Max_Log then
      return 0.0;
    else
      Exagam := Exp (Arg);
    end if;

    --
    --  (Exp(-x) * x**a / Gamma(a))* [1/a + x/(a(a+1)) + x^2/(a(a+1)(a+2)) + ...] 
    --
    -- factor out the 1/a:
    --
    --  = (Exagam / a) * [1 +  x/(a+1) + x^2/((a+1)(a+2)) + ...]
    --

    a_plus_n  := a;     -- n = 0
    Next_Term := 1.0;
    Sum       := 1.0;

    --  max number of allowed iterations arbitrarily set at 2**12:

    Series_Summation:
    for Iteration_id in 1 .. 2**12 loop
      a_plus_n  := a_plus_n + 1.0;
      Next_Term := Next_Term * x / a_plus_n;
      Sum       := Sum + Next_Term;
      if (Next_Term/Sum < Real_Epsilon) then   exit Series_Summation;  end if;
    end loop Series_Summation;

    return Sum * Exagam / a;

  end Incomplete_Gamma;


 ------------------
 -- Normal_CDF_x --
 ------------------

 --  Here just used for x out on tail, where it seems more accurate
 --  than the other method used above.
 --
 --  Evaluates the cumulative distribution function of a gaussian.
 --  Finds area of the standardized normal distribution
 --  (normalized gaussian with unit width) from -inf to x. 
 --  based on  Algorithm AS66 Applied Statistics (1973) vol.22, no.3.
 --    Usually  7 digits accuracy; better on x<<0 tail.

 function Normal_CDF_x (x : in Real) return Real is

   z, y, Result : Real;

   con : constant Real := 1.28;
   ltone  : constant Real := 7.0;
   utzero : constant Real := 40.0;
   p  : constant Real := 0.398942280444;
   q  : constant Real := 0.39990348504;
   r  : constant Real := 0.398942280385;
   a1 : constant Real := 5.75885480458;
   a2 : constant Real := 2.62433121679;
   a3 : constant Real := 5.92885724438;
   b1 : constant Real :=-29.8213557807;
   b2 : constant Real := 48.6959930692;
   c1 : constant Real :=-3.8052E-8;
   c2 : constant Real := 3.98064794E-4;
   c3 : constant Real :=-0.151679116635;
   c4 : constant Real := 4.8385912808;
   c5 : constant Real := 0.742380924027;
   c6 : constant Real := 3.99019417011;
   d1 : constant Real := 1.00000615302;
   d2 : constant Real := 1.98615381364;
   d3 : constant Real := 5.29330324926;
   d4 : constant Real :=-15.1508972451;
   d5 : constant Real := 30.789933034;

begin

   if not x'Valid then  raise Constraint_Error;  end if;

   z := Abs (x);

   if (z <= ltone  OR  (x < 0.0  AND  (z <= utzero)) ) then
     y := 0.5*z*z;
     if (z > con) then
       Result := r*Exp(-y) / (z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))));
     else
       Result := 0.5 - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))));
     end if;
   else
     Result := 0.0;
   end if;

   if (x >= 0.0) then   Result := 1.0 - Result;   end if;

   return Result;

 end Normal_CDF_x;


  -- rough test; more of working order than actual accuracy

  procedure Test_Normal_CDF (x_0 : in Real; Error : out Real) is

    Delta_x : constant Real := 2.0**(-6); -- 6 ok.

    type First_Deriv_Range_17 is range -8 .. 8;

    d_17 : constant array(First_Deriv_Range_17) of Real := 
   (9.71250971250971250971250971250971250E-6 ,
   -1.77600177600177600177600177600177600E-4 ,
    1.55400155400155400155400155400155400E-3 ,
   -8.70240870240870240870240870240870240E-3 ,
    3.53535353535353535353535353535353535E-2 ,
   -1.13131313131313131313131313131313131E-1 ,
    3.11111111111111111111111111111111111E-1 ,
   -8.88888888888888888888888888888888888E-1 ,
    0.0,
    8.88888888888888888888888888888888888E-1 ,
   -3.11111111111111111111111111111111111E-1 ,
    1.13131313131313131313131313131313131E-1 ,
   -3.53535353535353535353535353535353535E-2 ,
    8.70240870240870240870240870240870240E-3 ,
   -1.55400155400155400155400155400155400E-3 ,
    1.77600177600177600177600177600177600E-4 ,
   -9.71250971250971250971250971250971250E-6);

    One_over_sqrt_2_pi : constant := 0.398942280401432677939946;

    sum, x, Integral_of_Gaussian, Deriv : Real;

  begin

    sum := 0.0;
    for i in First_Deriv_Range_17 loop
      x := x_0 + Real(i)*Delta_x;
      Integral_of_Gaussian := Normal_CDF (x);
      sum := sum + d_17 (i) * Integral_of_Gaussian;
    end loop;
    Deriv := sum / Delta_x;

    -- relative
  --Error := (Deriv - Exp (-x_0**2 * 0.5) * One_over_sqrt_2_pi) / (Abs(Deriv)+1.0e-307);

    -- absolute
    Error := (Deriv - Exp (-x_0**2 * 0.5) * One_over_sqrt_2_pi);

  end Test_Normal_CDF;


end Chi_Gaussian_CDF;
