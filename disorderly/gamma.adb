
-------------------------------------------------------------------------------
-- package body Gamma, Log of Gamma Function
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

-- Thanks are due to Peter Luschny for the web page and discussion of
-- the Stieltjes continued fraction gamma function.
--
-- Parameter Error correction is based on visual inspection!
-- Parameter Error detection  is based on procedure Test_Stieltjes_Coefficients.

with ada.numerics.generic_elementary_functions;
with Gamma_1_to_3;
package body Gamma is

  package mathz is new ada.numerics.generic_elementary_functions(Real); use Mathz;
  package Remez_Gamma is new Gamma_1_to_3(Real); use Remez_Gamma;

  Zero  : constant Real := +0.0;
  One   : constant Real := +1.0;
  Three : constant Real := +3.0;

  Half_Log_Two_Pi : constant Real :=
          +0.91893853320467274178032973640561763986139747363778341;
        --+0.91893853320467274178032973640561763986139747363778341;
        --+0.91893853320467274178032973640561763986139747363778341;


  -- Coefficients for the Stieltjes continued fraction gamma function.

  b00 : constant Real :=  Zero;
  b01 : constant Real :=  One / (+12.0);
  b02 : constant Real :=  One / (+30.0);
  b03 : constant Real := (+53.0) / (+210.0);
  b04 : constant Real := (+195.0) / (+371.0);
  b05 : constant Real := (+22999.0) / (+22737.0);
  b06 : constant Real := (+29944523.0) / (+19733142.0);
  b07 : constant Real := (+109535241009.0) / (+48264275462.0);
  b08 : constant Real :=
          +3.00991738325939817007314073420771572737347407676490429;
        --+3.00991738325939817007314073420771572737347407676490429;
        --+3.00991738325939817007314073420771572737347407676490429;
  b09 : constant Real :=
          +4.02688719234390122616887595318144283632660221257447127;
        --+4.02688719234390122616887595318144283632660221257447127;
        --+4.02688719234390122616887595318144283632660221257447127;
  b10 : constant Real :=
          +5.00276808075403005168850241227618857481218316757980572;
        --+5.00276808075403005168850241227618857481218316757980572;
        --+5.00276808075403005168850241227618857481218316757980572;
  b11 : constant Real :=
          +6.28391137081578218007266315495164726427843982033318771;
        --+6.28391137081578218007266315495164726427843982033318771;
        --+6.28391137081578218007266315495164726427843982033318771;
  b12 : constant Real :=
          +7.49591912238403392975235470826750546574579869781542533;
        --+7.49591912238403392975235470826750546574579869781542533;
        --+7.49591912238403392975235470826750546574579869781542533;


  type CF_Coeff_range is range 0 .. 12;
  type CF_Coeff  is array(CF_Coeff_range) of Real;

  ---------------------------------
  -- Evaluate_Continued_Fraction --
  ---------------------------------

  -- no attempt here to prevent overflow

  -- CF_value := a0 + b1/(a1 + b2/(a2 + b3/( ...)))

  procedure Evaluate_Continued_Fraction
    (CF_value               :    out Real;
     Truncation_Error       :    out Real;
     a, b                   : in     CF_Coeff;
     Max_Term_in_Series     : in     CF_Coeff_Range := CF_Coeff_Range'Last;
     Error_Estimate_Desired : in     Boolean        := False)
  is
     P, Q : array (-1 .. Max_Term_in_Series) of Real;
  begin
     Q(-1) := Zero;
     Q(0)  := One;

     P(-1) := One;
     P(0)  := a(0);

     for j in 1 .. Max_Term_in_Series loop
       P(j) := a(j) * P(j-1)  +  b(j) * P(j-2);
       Q(j) := a(j) * Q(j-1)  +  b(j) * Q(j-2);
     end loop;

     CF_value := P(Max_Term_in_Series) / Q(Max_Term_in_Series);

     Truncation_Error := Zero;
     if Error_Estimate_Desired then
       Truncation_Error :=
            CF_value - P(Max_Term_in_Series-1) / Q(Max_Term_in_Series-1);
     end if;

  end Evaluate_Continued_Fraction;

  ---------------
  -- Log_Gamma --
  ---------------

  -- For x < 14 uses rational polynomial approximations.
  -- For x > 14, uses Stieltjes' continued fraction gamma:
  --
  -- gamma(N) := (N-1)!
  -- gamma(m+1) := m!
  --
  -- Stieltjes' continued fraction gamma:
  --
  -- good to 16 sig figures for x > 4.5
  -- good to 19 sig figures for x > 7
  -- good to 21 sig figures for x > 10
  --
  -- based on trunc_error, and on comparison with rational poly-gamma function,
  -- but very approximate!
  -- Expect 32 sig figures at not much higher than x=12, but not yet measured.
  --

  function Log_Gamma (x : in Real) return Real
  is
    a : constant CF_Coeff := (Zero, others => x);
    b : constant CF_Coeff :=
          (b00, b01, b02, b03, b04, b05, b06, b07, b08, b09, b10, b11, b12);
    CF, Trunc_Error, Log_Kernal_x, Log_Gamma_x : real;
  begin

    if not x'Valid then  raise Constraint_Error;  end if;
    --  In some cases, input of (say) inf or NaN will make numeric programs
    --  hang rather than crash .. very difficult to diagnose, so this seems
    --  best policy for function calls that are rarely found in time-sensitive
    --  inner loops. Very nice feature, 'Valid!

    if x <= Zero then
      raise Constraint_Error;  -- or is math arg error.
    end if;

    if x < 14.0 then
      return Log_Gamma_0_to_16 (x);
    end if;

    -- For testing. these 3 should give identical answers:
    -- CF :=   (1.0/12.0)/(x  + (1.0/30.0)/(x + (53.0/210.0)/(x  + (195.0/371.0)/x)));
    -- CF := a(0) +   b(1)/(a(1) +   b(2)/(a(2) +     b(3)/(a(3) +     b(4)/(a(4)))));
    -- Evaluate_Continued_Fraction (CF, Trunc_Error, a, b, 4, True);

    Evaluate_Continued_Fraction (CF, Trunc_Error, a, b, CF_Coeff_range'Last, True);
    --  text_io.put_line(Real'image(Trunc_Error));

    Log_Kernal_x    := Half_Log_Two_Pi + (x - 0.5)*Log (x) - x;
    Log_Gamma_x     := Log_Kernal_x + CF;

    return Log_Gamma_x;

  end Log_Gamma;


  -- reduce argument from x > 3 to x in [2,3]
  --
  --     gamma (x + 1) = x * gamma(x)      (gamma(i) = (i-1)!)
  --
  -- x in [3,4]:
  -- log_gamma (x + 1) = log(x) + log_gamma(x)
  -- log_gamma (x) = log((x-1)) + log_gamma(x-1)
  --
  -- x in [4,5]:
  -- log_gamma (x + 2) = log(x+1) + log_gamma(x+1) = log((x)*(x+1)) + log_gamma(x)
  -- log_gamma (x) = log((x-2)*(x-1)) + log_gamma(x-2)
  --
  -- x in [5,6]:
  -- log_gamma (x + 3) = log(x+2) + log_gamma(x+2) = log((x)*(x+1)*(x+2)) + log_gamma(x)
  -- log_gamma (x) = log((x-3)*(x-2)*(x-1)) + log_gamma(x-3)
  --

  function Log_Gamma_0_to_16 (x : in Real) return Real
  is
    Val, Arg, Prod : Real := Zero;
    Steps : Integer;
  begin

    if x <= Zero then
      raise Constraint_Error;
    end if;

    if x > +16.0 then
      raise Constraint_Error;
    end if;

    if x > Three then

      -- log_gamma (x) = log((x-1)) + log_gamma(x-1)              x in [3,4]
      -- log_gamma (x) = log((x-2)*(x-1)) + log_gamma(x-2)        x in [4,5]
      -- log_gamma (x) = log((x-3)*(x-2)*(x-1)) + log_gamma(x-3)  x in [5,6]

      Steps := 1 + Integer (x - Three);

      Arg  := x - Real(Steps);
      Prod := Arg;
      for i in 1..Steps-1 loop
        Prod := Prod * (x - Real(i));
      end loop;

      Val := Log (Prod) + Log_Gamma_1_to_3 (Arg);

    elsif x >= One then

      Val := Log_Gamma_1_to_3 (x);

    elsif x > Real_Epsilon then

      -- x in [eps/8, 1]:
      --
      --     gamma (x + 1) = x * gamma(x)
      -- log_gamma (x + 1) = Log(x) + log_gamma(x)
      -- log_gamma (x) = -Log(x) + log_gamma(x+1)

      Val := -Log (x) + Log_Gamma_1_to_3 (x + One);

    else

      Val := -Log (x);

    end if;

    return Val;

  end Log_Gamma_0_to_16;

  -- Just make sure that CF_Coeff's have not mutated:

  procedure Test_Stieltjes_Coefficients
  is
     Difference : Real;
 
     Numerator : constant CF_Coeff :=
       (0.0, 1.0, 1.0, 53.0, 195.0, 22999.0, 29944523.0, 109535241009.0,
        29404527905795295658.0,
        455377030420113432210116914702.0,
        26370812569397719001931992945645578779849.0,
        152537496709054809881638897472985990866753853122697839.0,
        100043420063777451042472529806266909090824649341814868347109676190691.0);
 
     Denominator : constant CF_Coeff :=
       (1.0, 12.0, 30.0, 210.0, 371.0, 22737.0, 19733142.0, 48264275462.0,
        9769214287853155785.0,
        113084128923675014537885725485.0,
        5271244267917980801966553649147604697542.0,
        24274291553105128438297398108902195365373879212227726.0,
        13346384670164266280033479022693768890138348905413621178450736182873.0);

     B_coeff : constant CF_Coeff :=
       (b00, b01, b02, b03, b04, b05, b06, b07, b08, b09, b10, b11, b12);

  begin
     
     for i in CF_Coeff_range loop
        Difference := B_coeff(i) - Numerator(i) / Denominator(i);
        if Abs Difference > 16.0 * Real_Epsilon then
          raise Program_Error;
        end if;
     end loop;

  end Test_Stieltjes_Coefficients;

end Gamma;

