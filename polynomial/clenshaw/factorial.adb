
-----------------------------------------------------------------------
-- package body Factorial. Use Stieltjes' continued fraction to get Log of N!
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
---------------------------------------------------------------------------

with Ada.Numerics.Generic_Elementary_Functions;
package body Factorial is

  package Maths is new Ada.Numerics.Generic_Elementary_Functions (Real); 
  use Maths;

  Zero  : constant Real := +0.0;
  One   : constant Real := +1.0;

  -- Coefficients for the Stieltjes continued fraction gamma function.
  -- The "+" functions can be used to translate the nums to an extended
  -- precision floating point.
  --
  -- A call to the test routine verifies that the b's have not mutated.

  b00 : constant Real := Zero;
  b01 : constant Real := One / (+12.0);
  b02 : constant Real := One / (+30.0);
  b03 : constant Real := (+53.0) / (+210.0);
  b04 : constant Real := (+195.0) / (+371.0);
  b05 : constant Real := (+22999.0) / (+22737.0);
  b06 : constant Real := (+29944523.0) / (+19733142.0);
  b07 : constant Real := (+109535241009.0) / (+48264275462.0);
  b08 : constant Real :=
    +3.009917383259398170073140734207715727373474076764904290594203015613;
  b09 : constant Real :=
    +4.026887192343901226168875953181442836326602212574471272142220773421;
  b10 : constant Real :=
    +5.002768080754030051688502412276188574812183167579805723221777492737;
  b11 : constant Real :=
    +6.283911370815782180072663154951647264278439820333187714994251481107;
  b12 : constant Real :=
    +7.495919122384033929752354708267505465745798697815425327954440696079;

  Half_Log_Two_Pi : constant Real :=
    +0.918938533204672741780329736405617639861397473637783412817151540483;
  --+0.918938533204672741780329736405617639861397473637783412817151540483;

  ---------------------------------
  -- Log_Factorial_Table_0_to_32 --
  ---------------------------------

  subtype Factorial_Table_Range is Natural range 0 .. 32;

  type Table is array (Factorial_Table_Range) of Real;

  -- A call to the test routine verifies that the table has not mutated.

  Log_Factorial_Table_0_to_32 : constant Table :=
    (0.0,
     0.0,
     0.69314718055994530941723212145817656807550013436025525412068000949,
     1.79175946922805500081247735838070227272299069218300470585537434313,
     3.17805383034794561964694160129705540887399096090351521409673436211,
     4.78749174278204599424770093452324304839959231517203293600938225359,
     6.57925121201010099506017829290394532112258300735503764186475659672,
     8.52516136106541430016553103634712505075966773693689883032414674666,
     10.60460290274525022841722740072165475498616814001766459268618677514,
     12.80182748008146961120771787456670616428114925566316349615557544241,
     15.10441257307551529522570932925107037188225074429193647218890334338,
     17.50230784587388583928765290721619967170395759822935364740747105251,
     19.98721449566188614951736238705507851250244842477261360738352540513,
     22.55216385312342288557084982862039711730771636953282072380257091580,
     25.19122118273868150009343469352175341502030123347493716638264107523,
     27.89927138384089156608943926367046675919339314556620434002998330034,
     30.67186010608067280375836774950317303149539368300722535651270333831,
     33.50507345013688888400790236737629956708359669559297014380994107619,
     36.39544520803305357621562496267952754445407794559872430140000975296,
     39.33988418719949403622465239456738108169145720689785311993796998937,
     42.33561646075348502965987597070992185736805882988688135009197789983,
     45.38013889847690802616047395107562729165263411729149199028606238341,
     48.47118135183522387963964965049893315954984110558916441962531010203,
     51.60667556776437357044640248230912927799222142042960016162394547952,
     54.78472939811231919009334408360618468686621238133311537572067984163,
     58.00360522298051993929486275005855996591741508987015081954597562458,
     61.26170176100200198476558231308205513879818316899061319008570114474,
     64.55753862700633105895131802384963225274065484245886154528978414565,
     67.88974313718153498289113501020916511852873984076123324199053431458,
     71.25703896716800901007440704257107672402325275468397732091220466622,
     74.65823634883016438548764373417796663627184480113549974868022690082,
     78.09222355331531063141680805872032384672178373161609172043694497330,
     81.55795945611503717850296866601120668709928440341736799104034502077);

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

  -------------------
  -- Log_Factorial --
  -------------------

  -- For x < 33,  uses table for log (factorial).
  -- For x >= 33, uses Stieltjes' continued fraction gamma:
  --
  -- gamma(N) := (N-1)!
  -- gamma(m+1) := m!
  --
  -- Stieltjes' continued fraction gamma:
  --
  --   good to 21 sig figures for x > 10
  --
  -- based on trunc_error, and on comparison with rational poly-gamma function,
  -- but very approximate!
  --

  function Log_Factorial (N : in Natural) return Real is

    x : constant Real := Real (N + 1); -- so gamma(x) = N!

    a : constant CF_Coeff := (Zero, others => x);
    b : constant CF_Coeff :=
       (b00, b01, b02, b03, b04, b05, b06, b07, b08, b09, b10, b11, b12);
    CF, Trunc_Error, Log_Kernal_x, Log_Gamma_x : Real;

  begin

    if N < 0  then  raise Constraint_Error; end if;

    if N in Factorial_Table_Range then
       return Log_Factorial_Table_0_to_32 (N);  
    end if;

    -- For testing. these 3 should give identical answers:
    -- CF :=   (1.0/12.0)/(x  + (1.0/30.0)/(x + (53.0/210.0)/(x  + (195.0/371.0)/x)));
    -- CF := a(0) +   b(1)/(a(1) +   b(2)/(a(2) +     b(3)/(a(3) +     b(4)/(a(4)))));
    -- Evaluate_Continued_Fraction (CF, Trunc_Error, a, b, 4, True);

    Evaluate_Continued_Fraction (CF, Trunc_Error, a, b, CF_Coeff_range'Last, False);

  --Evaluate_Continued_Fraction (CF, Trunc_Error, a, b, CF_Coeff_range'Last, True);
  --text_io.put_line(Real'image(Trunc_Error));

    Log_Kernal_x := Half_Log_Two_Pi + (x - 0.5)*Log (x) - x;
    Log_Gamma_x  := Log_Kernal_x + CF;

    return Log_Gamma_x;

  end Log_Factorial;

  -- Make sure that CF_Coeff's have not mutated:

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
        if Abs Difference > 4.0 * Real'Epsilon then
          raise Program_Error;
        end if;
     end loop;

  end Test_Stieltjes_Coefficients;

  procedure Test_Log_Factorial_Table is
     Factorial : Real := 1.0;
     Max_Err, Err : Real := 0.0;
  begin
     for i in Factorial_Table_Range loop
        Err := Abs (Log_Factorial_Table_0_to_32(i) - Log (Factorial));
        if Err > Max_Err then
           Max_Err := Err;
        end if;
        Factorial := Factorial * Real (i+1);
     end loop;
  
     if Max_Err > 16.0 * Real'Epsilon then
        raise Program_Error;
     end if;

  end Test_Log_Factorial_Table;

end Factorial;

