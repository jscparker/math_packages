
-----------------------------------------------------------------------
-- package Runge_Coeffs_PD_8, coefficients for Prince-Dormand Runge Kutta
-- Copyright (C) 2008-2018 Jonathan S. Parker.
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

-- package Runge_Coeffs_PD_8
--
-- Package contains coefficients for the Prince-Dormand
-- 8th order, 13 stage Runge Kutta method with error control.
--
-- You have a choice between the old 17 significant figure coefficients
-- and newer ones, given to about 28 significant figures.
-- The older ones are entered in rational form.
--
-- A test routine is provided to make make sure that they
-- have not been corrupted.  If the coefficients do not pass 
-- the tests, then they are given in redundant and alternate 
-- forms at the end of this package (commented out).
--
-- Reference: "High Order Embedded Runge-Kutta Formulae" by P.J.
-- Prince and J.R. Dormand, J. Comp. Appl. Math.,7, pp. 67-75, 1981.
--
--
-- NOTES ON THE ALGORITHM
--
-- Given a differential equ. dY/dt = F(t,Y) and an initial condition
-- Y(t), the Runge-Kutta formulas predict a new value of Y at
-- Y(t + h) with the formula
--
--                          max
--        Y(t + h) = Y(t) + SUM (K(j) * B(j))
--                          j=1
--
-- where max = Stages'Last, (which equals 13 here, since this is a 13
-- stage Runge-Kutta formula) where the quantities K are calculated from
--
--        K(1) = h * F(t, Y(t))
--
--                                         j-1
--        K(j) = h * F (t + C(j)*h, Y(t) + SUM (K(i)*A(j,i)) )
--                                         i=1
--
--  The Fehberg method, used here, provides two versions of the array
--  B, so that Y(t + h) may be calculated to both 7th order and to
--  8th order using the same K(i)'s.  This way, the local truncation
--  error can be estimated without calculating the K(i)'s twice.
--
--  If we apply this formula to a time-independent linear F, then we can
--  derive a condition that can be used to test the correctness of the
--  initialization values of the arrays A, B and C.  To derive such a
--  formula we use the fact that the RK prediction of Y(t + h) must equal
--  the Taylor series prediction up to the required order.  So,
--
--     K(1) = h*F*Y    (where F now is a matrix, and Y a vector)
--
--     K(2) = h*F*(Y + K(1)*A21)
--
--     K(3) = h*F*(Y + K(1)*A31 + K(2)*A32)
--
--     K(4) = h*F*(Y + K(1)*A41 + K(2)*A42 + K(3)*A43)
--
--  The linearity of F implies that F(a*Y + Z) = a*F*Y + F*Z so:
--
--     K(1) = h*F*Y
--
--     K(2) = h*F*Y +  A21*h^2*F^2*Y
--
--     K(3) = h*F*Y + (A31 + A32)*h^2*F^2*Y   +          A32*A21*h^3*F^3*Y
--
--     K(4) = h*F*Y + (A41+A42+A43)*h^2*F^2*Y + (A42*A21+A43*A31)h^3*F^3*Y
--                                              +       A43*A32*A21*h^4F^4*Y
--
--  Now we use the fact that we must have the RK prediction equal that
--  of the Taylor's series up to a certain order:
--
--     max
--     SUM (K(j) * B(j)) = h*F*Y + ... + (1/n!)*h^n*F^n*Y +  O(h^n+1)
--     j=1
--
-- Here n=8 for the coefficients B = B8 given below.  Its n=7 for B7.
-- The above formula gives us a relation between 1/n! and B and A.
-- This formula is used in the procedure TestRKPD given at the end
-- of this package.  We see immediately that we must have:
--
--  max
--  SUM (B(i))                  = 1/1!
--  i=1
--
--  max         i-1
--  SUM (B(i) * SUM(Aij)) = 1/2!
--  i=2         j=1
--
--***************************************************************
generic

   type Real is digits <>;

package Runge_Coeffs_PD_8 is

   subtype RK_Range  is Integer  range 0..12; 
   subtype Stages    is RK_Range range 0..12; -- always 0 .. 12
   type Coefficient  is array(RK_Range) of Real;
   type Coefficient_Array  is array(RK_Range) of Coefficient;

   --  At the bottom a renaming is used to choose which
   --  version of the Coefficients to use:
   --
   --   A_rational are good to about 17 decimal digits.
   --   A_extended are good to about 28 decimal digits.
   --
   --C  : Coefficient renames  C_extended;
   --B7 : Coefficient renames B7_extended;
   --B8 : Coefficient renames B8_extended;
   --A  : Coefficient_Array renames A_extended;
         

   procedure Test_Runge_Coeffs;

    
   --  Below the Coefficients A, B, and C are given to
   --  18 significant figures, in rational form.

   --  coefficients C(0..12) for getting Dt

   C_rational : constant Coefficient :=
     (0.0,
      1.0 / 18.0,
      1.0 / 12.0,
      1.0 / 8.0,
      5.0 / 16.0,
      3.0 / 8.0,
      59.0 / 400.0,
      93.0 / 200.0,
      5490023248.0 / 9719169821.0,
      13.0 / 20.0,
      1201146811.0 / 1299019798.0,
      1.0,
      1.0,
      others => 0.0);

   --  eighth order RKPD coefficients B8(0..12)
   --  so that Y(t + DeltaT) = Y(t) + SUM(B8*K)

   B8_rational : constant Coefficient :=
     (14005451.0 / 335480064.0,
      0.0,
      0.0,
      0.0,
      0.0,
     -59238493.0 / 1068277825.0,
      181606767.0 / 758867731.0,
      561292985.0 / 797845732.0,
     -1041891430.0 / 1371343529.0,
      760417239.0 / 1151165299.0,
      118820643.0 / 751138087.0,
     -528747749.0 / 2220607170.0,
      0.25,
      others => 0.0);


   --  Seventh order RKPD coefficients B7(0..12)
   --  so that Y(t + DeltaT) = Y(t) + SUM(B7*K)

   B7_rational : constant Coefficient :=
     (13451932.0 / 455176623.0,
      0.0,
      0.0,
      0.0,
      0.0,
     -808719846.0 / 976000145.0,
      1757004468.0 / 5645159321.0,
      656045339.0 / 265891186.0,
     -3867574721.0 / 1518517206.0,
      465885868.0 / 322736535.0,
      53011238.0 / 667516719.0,
      2.0 / 45.0,
      others => 0.0);



     --  The Runge-Kutta matrix for calculating K:

   A_rational : constant Coefficient_Array :=

     ((others => 0.0),

      --  coefficients A(1)(0..12) for getting K(1),
      --  so K(1) = DeltaT*F (Time + Dt(1), Y + SUM (A(1), K))
      --  A(1):

      (1.0  / 18.0,
       0.0,
       0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, others => 0.0),


      --  coefficients A(2)(0..12) for getting K(2),
      --  so K(2) = DeltaT*F (Time + Dt(2),Y+SUM(A2,K))
      --  A(2):

      (1.0  / 48.0,
       1.0  / 16.0,
       0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, others => 0.0),


      --  coefficients A(3)(0..12) for getting K(3)
      --  so K4 = DeltaT*F (Time + Dt(3), Y + SUM (A(3), K))

      (1.0 / 32.0,
       0.0,
       3.0 / 32.0,
       0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, others => 0.0),


      --  coefficients A(4)(0..12) for getting K(4)

      (5.0 / 16.0,
       0.0,
      -75.0 / 64.0,
       75.0 / 64.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, others => 0.0),


      -- coefficients A(5)(0..12) for getting K(5)

      (3.0 / 80.0,
       0.0,
       0.0,
       3.0 / 16.0,
       3.0 / 20.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, others => 0.0),


      -- coefficients A(6)(0..12) for getting K(6)

      (29443841.0 / 614563906.0,
       0.0,
       0.0,
       77736538.0 / 692538347.0,
      -28693883.0 / 1125000000.0,
       23124283.0 / 1800000000.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, others => 0.0),


      -- coefficients A(7)(0..12) for getting K(7)

      (16016141.0 / 946692911.0,
       0.0,
       0.0,
       61564180.0 / 158732637.0,
       22789713.0 / 633445777.0,
       545815736.0 / 2771057229.0,
      -180193667.0 / 1043307555.0,
       0.0,
       0.0, 0.0, 0.0, 0.0, others => 0.0),

       -- coefficients A(8)(0..12) for getting K(9)

      (39632708.0 /  573591083.0,
       0.0,
       0.0,
      -433636366.0 / 683701615.0,
      -421739975.0 / 2616292301.0,
       100302831.0 / 723423059.0,
       790204164.0 / 839813087.0,
       800635310.0 / 3783071287.0,
       0.0,
       0.0, 0.0, 0.0, others => 0.0),


      -- coefficients A(9)(0..12) for getting K(9)

      (246121993.0 / 1340847787.0,
       0.0,
       0.0,
      -37695042795.0 / 15268766246.0,
      -309121744.0 / 1061227803.0,
      -12992083.0 / 490766935.0,
       6005943493.0 / 2108947869.0,
       393006217.0 / 1396673457.0,
       123872331.0 / 1001029789.0,
       0.0, 0.0, 0.0, others => 0.0),


     -- coefficients A(10)(0..12) for getting K(10)

     (-1028468189.0 / 846180014.0,
       0.0,
       0.0,
       8478235783.0 / 508512852.0,
       1311729495.0 / 1432422823.0,
      -10304129995.0 / 1701304382.0,
      -48777925059.0 / 3047939560.0,
       15336726248.0 / 1032824649.0,
      -45442868181.0 / 3398467696.0,
       3065993473.0 / 597172653.0,
       0.0, 0.0, others => 0.0),


       -- coefficients A(11)(0..12) for getting K(11)

      (185892177.0 / 718116043.0,
       0.0,
       0.0,
      -3185094517.0 / 667107341.0,
      -477755414.0 / 1098053517.0,
      -703635378.0 / 230739211.0,
       5731566787.0 / 1027545527.0,
       5232866602.0 / 850066563.0,
      -4093664535.0 / 808688257.0,
       3962137247.0 / 1805957418.0,
       65686358.0 / 487910083.0,
       0.0, 
       others => 0.0),


       --  coefficients A(12)(0..12) for getting K(12)


      (403863854.0 / 491063109.0,
       0.0,
       0.0,
      -5068492393.0 / 434740067.0,
      -411421997.0 / 543043805.0,
       652783627.0 / 914296604.0,
       11173962825.0 / 925320556.0,
      -13158990841.0 / 6184727034.0,
       3936647629.0 / 1978049680.0,
      -160528059.0 / 685178525.0,
       248638103.0 / 1413531060.0,
       0.0, 
       others => 0.0),

       others => (others => 0.0));


   -- Below the Coefficients A, B, and C are given to (usually)
   -- 28 significant figures.

   --  coefficients C(0..12) for getting Dt

   C_extended : constant Coefficient :=
      (0.0,
       5.5555555555555555555555555555556E-02,
       8.3333333333333333333333333333333E-02,
       1.25E-01,
       3.125E-01,
       3.75E-01,
       1.475E-01,
       4.65E-01,
       5.64865451382259575398358501426E-01,
       6.5E-01,
       9.24656277640504446745013574318E-01,
       1.0,
       1.0,
       others => 0.0);

   --  eighth order RKPD coefficients B8(0..12)
   --  so that Y(t + DeltaT) = Y(t) + SUM(B8*K)
   B8_extended : constant Coefficient :=
      (4.17474911415302462220859284685E-02,
       0.0,
       0.0,
       0.0,
       0.0,
      -5.54523286112393089615218946547E-02,
       2.39312807201180097046747354249E-01,
       7.0351066940344302305804641089E-01,
      -7.59759613814460929884487677085E-01,
       6.60563030922286341461378594838E-01,
       1.58187482510123335529614838601E-01,
      -2.38109538752862804471863555306E-01,
       0.25,
       others => 0.0);

   --  Seventh order RKPD coefficients B7(0..12)
   --  so that Y(t + DeltaT) = Y(t) + SUM(B7*K)
   B7_extended : constant Coefficient :=
      (2.9553213676353496981964883112E-02,
       0.0,
       0.0,
       0.0,
       0.0,
      -8.28606276487797039766805612689E-01,
       3.11240900051118327929913751627E-01,
       2.46734519059988698196468570407E+00,
      -2.54694165184190873912738007542E+00,
       1.44354858367677524030187495069E+00,
       7.94155958811272872713019541622E-02,
       4.4444444444444444444444444444444E-02,
       0.0,
       others => 0.0);

   -- Here is the Runge-Kutta matrix for calculating K

   A_extended : constant Coefficient_Array :=

      -- coefficients A(1)(0..12) for getting K(1),
      -- so K(2) = DeltaT*F (Time + Dt(1), Y + SUM (A(1), K))
      -- A(2):
     ((others => 0.0),

      (1.0  / 18.0,
       0.0,
       0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       others => 0.0),


       -- coefficients A(2)(0..12) for getting K(2),
       -- so K(3) = DeltaT*F (Time + Dt(2),Y+SUM(A(2),K))
       -- A(3):
      (1.0  / 48.0,
       1.0  / 16.0,
       0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       others => 0.0),


       -- coefficients A(3)(0..12) for getting K(3)
       -- so K4 = DeltaT*F (Time + Dt(3), Y + SUM (A(3), K))
      (1.0 / 32.0,
       0.0,
       3.0 / 32.0,
       0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0,
       others => 0.0),


       -- coefficients A(4)(0..12) for getting K(4)
       -- so K(5) = DeltaT*F (Time + Dt(4), Y + SUM (A(4), K))
      (5.0 / 16.0,
       0.0,
      -75.0 / 64.0,
       75.0 / 64.0,
       0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0,
       others => 0.0),


       -- coefficients A(5)(0..12) for getting K(5)
       -- so K6 = DeltaT*F (Time + Dt(5),Y+SUM(A(5),K))
      (3.0 / 80.0,
       0.0,
       0.0,
       3.0 / 16.0,
       3.0 / 20.0,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0,
       others => 0.0),


       -- coefficients A(6)(0..12) for getting K(6)
       -- so K7 = DeltaT*F (Time + Dt(6),Y+SUM(A(6),K))
      (4.7910137111111111111111111111111E-02,
       0.0,
       0.0,
       1.1224871277777777777777777777778E-01,
      -2.5505673777777777777777777777778E-02,
       1.2846823888888888888888888888889E-02,
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0,
       others => 0.0),


       -- coefficients A(7)(0..12) for getting K(7)
       -- so K(8) = DeltaT*F (Time + Dt(7), Y + SUM (A(7), K))
      (1.6917989787292281181431107136E-02,
       0.0,
       0.0,
       3.87848278486043169526545744159E-01,
       3.59773698515003278967008896348E-02,
       1.96970214215666060156715256072E-01,
      -1.72713852340501838761392997002E-01,
       0.0,
       0.0, 0.0, 0.0, 0.0, 0.0,
       others => 0.0),


       -- coefficients A(8)(0..12) for getting K(8)
       -- so K(9) = DeltaT*F (Time + Dt(8),Y + SUM (A(8), K))
      (6.90957533591923006485645489846E-02,
       0.0,
       0.0,
      -6.34247976728854151882807874972E-01,
      -1.61197575224604080366876923982E-01,
       1.38650309458825255419866950133E-01,
       9.4092861403575626972423968413E-01,
       2.11636326481943981855372117132E-01,
       0.0,
       0.0, 0.0, 0.0, 0.0,
       others => 0.0),


       -- coefficients A(9)(0..12) for getting K(9)
       -- so K(10) = DeltaT*F (Time + Dt(9),Y + SUM (A(9), K))
      (1.83556996839045385489806023537E-01,
       0.0,
       0.0,
      -2.46876808431559245274431575997E+00,
      -2.91286887816300456388002572804E-01,
      -2.6473020233117375688439799466E-02,
       2.84783876419280044916451825422E+00,
       2.81387331469849792539403641827E-01,
       1.23744899863314657627030212664E-01,
       0.0, 0.0, 0.0, 0.0,
       others => 0.0),


       -- coefficients A(10)(0..12) for getting K(10)
       -- so K(11) = DeltaT*F (Time + Dt(10),Y + SUM (A(10), K))
     (-1.21542481739588805916051052503E+00,
       0.0,
       0.0,
       1.66726086659457724322804132886E+01,
       9.15741828416817960595718650451E-01,
      -6.05660580435747094755450554309E+00,
      -1.60035735941561781118417064101E+01,
       1.4849303086297662557545391898E+01,
      -1.33715757352898493182930413962E+01,
       5.13418264817963793317325361166E+00,
       0.0, 0.0, 0.0,
       others => 0.0),


       -- coefficients A(11)(0..12) for getting K(11)
       -- so K(12) = DeltaT*F (Time + Dt(11),Y + SUM (A(11), K))
      (2.58860916438264283815730932232E-01,
       0.0,
       0.0,
      -4.77448578548920511231011750971E+00,
      -4.3509301377703250944070041181E-01,
      -3.04948333207224150956051286631E+00,
       5.57792003993609911742367663447E+00,
       6.15583158986104009733868912669E+00,
      -5.06210458673693837007740643391E+00,
       2.19392617318067906127491429047E+00,
       1.34627998659334941535726237887E-01,
       0.0, 0.0,
       others => 0.0),


       -- coefficients A(12)(0..12) for getting K(12)
       -- so K(13) = DeltaT*F (Time + Dt(12),Y + SUM (A(12), K))
      (8.22427599626507477963168204773E-01,
       0.0,
       0.0,
      -1.16586732572776642839765530355E+01,
      -7.57622116690936195881116154088E-01,
       7.13973588159581527978269282765E-01,
       1.20757749868900567395661704486E+01,
      -2.12765911392040265639082085897E+00,
       1.99016620704895541832807169835E+00,
      -2.34286471544040292660294691857E-01,
       1.7589857770794226507310510589E-01,
       0.0, 0.0,
       others => 0.0),

       others => (others => 0.0));


 --C  : Coefficient renames  C_rational;
 --B7 : Coefficient renames B7_rational;
 --B8 : Coefficient renames B8_rational;
 --A  : Coefficient_Array renames A_rational;

   -- These are preferred for any Real'Digits > 15:

   C  : Coefficient renames  C_extended;
   B7 : Coefficient renames B7_extended;
   B8 : Coefficient renames B8_extended;
   A  : Coefficient_Array renames A_extended;
 

  ----------------------------------------------------------------------- 
  --              -- coefficients A(2)(0..12) for getting K(2),
  --              -- so K(2) = DeltaT*F (Time + Dt(2), Y + SUM (A(2), K))
  --              -- A(2):
  --             ((1.0  / 18.0,
  --               0.0,
  --               0.0, 0.0, 0.0, 0.0, 0.0,
  --               0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  --
  --   
  --              -- coefficients A(3)(0..12) for getting K(3),
  --              -- so K(3) = DeltaT*F (Time + Dt(3),Y+SUM(A3,K))
  --              -- A(3):
  --              (1.0  / 48.0,
  --               1.0  / 16.0,
  --               0.0, 0.0, 0.0, 0.0, 0.0,
  --               0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  --
  --   
  --              -- coefficients A(4)(0..12) for getting K(4)
  --              -- so K4 = DeltaT*F (Time + Dt(4), Y + SUM (A(4), K))
  --               (1.0 / 32.0,
  --                0.0,
  --                3.0 / 32.0,
  --                0.0, 0.0, 0.0, 0.0, 0.0,
  --                0.0, 0.0, 0.0, 0.0, 0.0),
  --
  --   
  --              -- coefficients A(5)(0..12) for getting K(5)
  --              -- so K(5) = DeltaT*F (Time + Dt(5), Y + SUM (A(5), K))
  --               (5.0 / 16.0,
  --               0.0,
  --              -75.0 / 64.0,
  --               75.0 / 64.0,
  --               0.0, 0.0, 0.0, 0.0,
  --               0.0, 0.0, 0.0, 0.0, 0.0),
  --
  --   
  --              -- coefficients A(6)(0..12) for getting K(6)
  --              -- so K6 = DeltaT*F (Time + Dt(6),Y+SUM(A6,K))
  --              (3.0 / 80.0,
  --               0.0,
  --               0.0,
  --               3.0 / 16.0,
  --               3.0 / 20.0,
  --               0.0, 0.0, 0.0,
  --               0.0, 0.0, 0.0, 0.0, 0.0),
  --
  --   
  --              -- coefficients A(7)(0..12) for getting K(7)
  --              -- so K7 = DeltaT*F (Time + Dt(7),Y+SUM(A7,K))
  --               (29443841.0 / 614563906.0,
  --               0.0,
  --               0.0,
  --                77736538.0 / 692538347.0,
  --               -28693883.0 / 1125000000.0,
  --                23124283.0 / 1800000000.0,
  --               0.0, 0.0, 0.0,
  --               0.0, 0.0, 0.0, 0.0),
  --
  --   
  --              -- coefficients A(8)(0..12) for getting K(8)
  --              -- so K(8) = DeltaT*F (Time + Dt(8), Y + SUM (A(8), K))
  --               (16016141.0 / 946692911.0,
  --               0.0,
  --               0.0,
  --               61564180.0 / 158732637.0,
  --               22789713.0 / 633445777.0,
  --              545815736.0 / 2771057229.0,
  --             -180193667.0 / 1043307555.0,
  --               0.0,
  --               0.0, 0.0, 0.0, 0.0, 0.0),
  --
  --   
  --              -- coefficients A(9)(0..12) for getting K(9)
  --              -- so K(9) = DeltaT*F (Time + Dt(9),Y + SUM (A(9), K))
  --               (39632708.0 /  573591083.0,
  --               0.0,
  --               0.0,
  --              -433636366.0 / 683701615.0,
  --              -421739975.0 / 2616292301.0,
  --               100302831.0 / 723423059.0,
  --               790204164.0 / 839813087.0,
  --               800635310.0 / 3783071287.0,
  --               0.0,
  --               0.0, 0.0, 0.0, 0.0),
  --
  --   
  --              -- coefficients A(10)(0..12) for getting K(10)
  --              -- so K(10) = DeltaT*F (Time + Dt(10),Y + SUM (A(10), K))
  --               (246121993.0 / 1340847787.0,
  --                0.0,
  --                0.0,
  --             -37695042795.0 / 15268766246.0,
  --               -309121744.0 / 1061227803.0,
  --                -12992083.0 / 490766935.0,
  --               6005943493.0 / 2108947869.0,
  --                393006217.0 / 1396673457.0,
  --                123872331.0 / 1001029789.0,
  --                 0.0, 0.0, 0.0, 0.0),
  --
  --   
  --              -- coefficients A(11)(0..12) for getting K(11)
  --              -- so K(11) = DeltaT*F (Time + Dt(11),Y + SUM (A(11), K))
  --               (-1028468189.0 / 846180014.0,
  --                0.0,
  --                0.0,
  --                 8478235783.0 / 508512852.0,
  --                 1311729495.0 / 1432422823.0,
  --               -10304129995.0 / 1701304382.0,
  --               -48777925059.0 / 3047939560.0,
  --                15336726248.0 / 1032824649.0,
  --               -45442868181.0 / 3398467696.0,
  --                 3065993473.0 / 597172653.0,
  --                 0.0, 0.0, 0.0),
  --
  --   
  --              -- coefficients A(12)(0..12) for getting K(12)
  --              -- so K(12) = DeltaT*F (Time + Dt(12),Y + SUM (A(12), K))
  --              (185892177.0 / 718116043.0,
  --               0.0,
  --               0.0,
  --             -3185094517.0 / 667107341.0,
  --              -477755414.0 / 1098053517.0,
  --              -703635378.0 / 230739211.0,
  --              5731566787.0 / 1027545527.0,
  --              5232866602.0 / 850066563.0,
  --             -4093664535.0 / 808688257.0,
  --              3962137247.0 / 1805957418.0,
  --                65686358.0 / 487910083.0,
  --                0.0, 0.0),
  --
  --   
  --              -- coefficients A(13)(0..12) for getting K(13)
  --              -- so K(13) = DeltaT*F (Time + Dt(13),Y + SUM (A(13), K))
  --               (403863854.0 / 491063109.0,
  --               0.0,
  --               0.0,
  --              -5068492393.0 / 434740067.0,
  --               -411421997.0 / 543043805.0,
  --                652783627.0 / 914296604.0,
  --              11173962825.0 / 925320556.0,
  --             -13158990841.0 / 6184727034.0,
  --               3936647629.0 / 1978049680.0,
  --               -160528059.0 / 685178525.0,
  --               248638103.0 / 1413531060.0,
  --               0.0, 0.0));
  --
  --  This pair is from "High Order Embedded Runge-Kutta Formulae" by P.J.
  --  Prince and J.R. Dormand, J. Comp. Appl. Math.,7, pp. 67-75, 1981.
  --  These coefficients were taken from a public domain RK Suite written
  --  R. W. Brankin,I. Gladwell, and L. F. Shampine.  They were calculated by
  --  Prince and Dormand.
  --        A(2,1) := 5.55555555555555555555555555556E-2
  --        A(3,1) := 2.08333333333333333333333333333E-2
  --        A(3,2) := 6.25E-2
  --        A(4,1) := 3.125E-2
  --        A(4,2) := 0.0
  --        A(4,3) := 9.375E-2
  --        A(5,1) := 3.125E-1
  --        A(5,2) := 0.0
  --        A(5,3) := -1.171875E0
  --        A(5,4) := 1.171875E0
  --        A(6,1) := 3.75E-2
  --        A(6,2) := 0.0
  --        A(6,3) := 0.0
  --        A(6,4) := 1.875E-1
  --        A(6,5) := 1.5E-1
  --        A(7,1) := 4.79101371111111111111111111111E-2
  --        A(7,2) := 0.0
  --        A(7,3) := 0.0
  --        A(7,4) := 1.12248712777777777777777777778E-1
  --        A(7,5) := -2.55056737777777777777777777778E-2
  --        A(7,6) := 1.28468238888888888888888888889E-2
  --        A(8,1) := 1.6917989787292281181431107136E-2
  --        A(8,2) := 0.0
  --        A(8,3) := 0.0
  --        A(8,4) := 3.87848278486043169526545744159E-1
  --        A(8,5) := 3.59773698515003278967008896348E-2
  --        A(8,6) := 1.96970214215666060156715256072E-1
  --        A(8,7) := -1.72713852340501838761392997002E-1
  --        A(9,1) := 6.90957533591923006485645489846E-2
  --        A(9,2) := 0.0
  --        A(9,3) := 0.0
  --        A(9,4) := -6.34247976728854151882807874972E-1
  --        A(9,5) := -1.61197575224604080366876923982E-1
  --        A(9,6) := 1.38650309458825255419866950133E-1
  --        A(9,7) := 9.4092861403575626972423968413E-1
  --        A(9,8) := 2.11636326481943981855372117132E-1
  --        A(10,1) := 1.83556996839045385489806023537E-1
  --        A(10,2) := 0.0
  --        A(10,3) := 0.0
  --        A(10,4) := -2.46876808431559245274431575997E0
  --        A(10,5) := -2.91286887816300456388002572804E-1
  --        A(10,6) := -2.6473020233117375688439799466E-2
  --        A(10,7) := 2.84783876419280044916451825422E0
  --        A(10,8) := 2.81387331469849792539403641827E-1
  --        A(10,9) := 1.23744899863314657627030212664E-1
  --        A(11,1) := -1.21542481739588805916051052503E0
  --        A(11,2) := 0.0
  --        A(11,3) := 0.0
  --        A(11,4) := 1.66726086659457724322804132886E1
  --        A(11,5) := 9.15741828416817960595718650451E-1
  --        A(11,6) := -6.05660580435747094755450554309E0
  --        A(11,7) := -1.60035735941561781118417064101E1
  --        A(11,8) := 1.4849303086297662557545391898E1
  --        A(11,9) := -1.33715757352898493182930413962E1
  --        A(11,10) := 5.13418264817963793317325361166E0
  --        A(12,1) := 2.58860916438264283815730932232E-1
  --        A(12,2) := 0.0
  --        A(12,3) := 0.0
  --        A(12,4) := -4.77448578548920511231011750971E0
  --        A(12,5) := -4.3509301377703250944070041181E-1
  --        A(12,6) := -3.04948333207224150956051286631E0
  --        A(12,7) := 5.57792003993609911742367663447E0
  --        A(12,8) := 6.15583158986104009733868912669E0
  --        A(12,9) := -5.06210458673693837007740643391E0
  --        A(12,10) := 2.19392617318067906127491429047E0
  --        A(12,11) := 1.34627998659334941535726237887E-1
  --        A(13,1) := 8.22427599626507477963168204773E-1
  --        A(13,2) := 0.0
  --        A(13,3) := 0.0
  --        A(13,4) := -1.16586732572776642839765530355E1
  --        A(13,5) := -7.57622116690936195881116154088E-1
  --        A(13,6) := 7.13973588159581527978269282765E-1
  --        A(13,7) := 1.20757749868900567395661704486E1
  --        A(13,8) := -2.12765911392040265639082085897E0
  --        A(13,9) := 1.99016620704895541832807169835E0
  --        A(13,10) := -2.34286471544040292660294691857E-1
  --        A(13,11) := 1.7589857770794226507310510589E-1
  --        A(13,12) := 0.0
  --
  --        B8 : constant Coefficient :=
  --          (4.17474911415302462220859284685E-2,
  --          0.0,
  --          0.0,
  --          0.0,
  --          0.0,
  --         -5.54523286112393089615218946547E-2,
  --          2.39312807201180097046747354249E-1,
  --          7.0351066940344302305804641089E-1,
  --         -7.59759613814460929884487677085E-1,
  --          6.60563030922286341461378594838E-1,
  --          1.58187482510123335529614838601E-1,
  --         -2.38109538752862804471863555306E-1,
  --          2.5E-1);
  --
  --        B7 : constant Coefficient :=
  --        (2.9553213676353496981964883112E-2,
  --         0.0,
  --         0.0,
  --         0.0,
  --         0.0,
  --         -8.28606276487797039766805612689E-1,
  --         3.11240900051118327929913751627E-1,
  --         2.46734519059988698196468570407E0,
  --         -2.54694165184190873912738007542E0,
  --         1.44354858367677524030187495069E0,
  --         7.94155958811272872713019541622E-2,
  --         4.44444444444444444444444444445E-2,
  --         0.0);
  --
  --        C : constant Coefficient :=
  --         (0.0,
  --         5.55555555555555555555555555556E-2,
  --         8.33333333333333333333333333334E-2,
  --         1.25E-1,
  --         3.125E-1,
  --         3.75E-1,
  --         1.475E-1,
  --         4.65E-1,
  --         5.64865451382259575398358501426E-1,
  --         6.5E-1,
  --         9.24656277640504446745013574318E-1,
  --         1.0,
  --         1.0);

--**************coefficients C(0..12) for getting Dt*******
--    C : constant Coefficient :=
--            (0.0,
--             5.5555555555555555555555555555556E-02,
--             8.3333333333333333333333333333333E-02,
--             1.25E-01,
--             3.125E-01,
--             3.75E-01,
--             1.475E-01,
--             4.65E-01,
--             5.64865451382259575398358501426E-01,
--             6.5E-01,
--             9.24656277640504446745013574318E-01,
--             1.0,
--             1.0,
--           others => 0.0);
--
--        --*******eighth order RKPD coefficients B8(0..12)*******
--        --******so that Y(t + DeltaT) = Y(t) + SUM(B8*K)*****
--    B8 : constant Coefficient :=
--           (4.17474911415302462220859284685E-02,
--           0.0,
--           0.0,
--           0.0,
--           0.0,
--           -5.54523286112393089615218946547E-02,
--            2.39312807201180097046747354249E-01,
--            7.0351066940344302305804641089E-01,
--           -7.59759613814460929884487677085E-01,
--            6.60563030922286341461378594838E-01,
--            1.58187482510123335529614838601E-01,
--           -2.38109538752862804471863555306E-01,
--            0.25, others => 0.0);
--
--        --*******Seventh order RKPD coefficients B7(0..12)*******
--        --******so that Y(t + DeltaT) = Y(t) + SUM(B7*K)*****
--    B7 : constant Coefficient :=
--            (2.9553213676353496981964883112E-02,
--             0.0,
--             0.0,
--             0.0,
--             0.0,
--            -8.28606276487797039766805612689E-01,
--             3.11240900051118327929913751627E-01,
--             2.46734519059988698196468570407E+00,
--            -2.54694165184190873912738007542E+00,
--             1.44354858367677524030187495069E+00,
--             7.94155958811272872713019541622E-02,
--             4.4444444444444444444444444444444E-02,
--             others => 0.0);
--
--  
--            -- Here is the Runge-Kutta matrix for calculating K
--
--        A : constant Coefficient_Array :=
--
--  
--              -- coefficients A(2)(0..12) for getting K(2),
--              -- so K(2) = DeltaT*F (Time + Dt(2), Y + SUM (A(2), K))
--              -- A(2):
--             ((1.0  / 18.0,
--               0.0,
--               0.0, 0.0, 0.0, 0.0, 0.0,
--               0.0, 0.0, 0.0, 0.0, 0.0, others => 0.0),
--
--  
--              -- coefficients A(3)(0..12) for getting K(3),
--              -- so K(3) = DeltaT*F (Time + Dt(3),Y+SUM(A3,K))
--              -- A(3):
--              (1.0  / 48.0,
--               1.0  / 16.0,
--               0.0, 0.0, 0.0, 0.0, 0.0,
--               0.0, 0.0, 0.0, 0.0, 0.0, others => 0.0),
--
--  
--              -- coefficients A(4)(0..12) for getting K(4)
--              -- so K4 = DeltaT*F (Time + Dt(4), Y + SUM (A(4), K))
--               (1.0 / 32.0,
--                0.0,
--                3.0 / 32.0,
--                0.0, 0.0, 0.0, 0.0, 0.0,
--                0.0, 0.0, 0.0, 0.0, others => 0.0),
--
--  
--              -- coefficients A(5)(0..12) for getting K(5)
--              -- so K(5) = DeltaT*F (Time + Dt(5), Y + SUM (A(5), K))
--               (5.0 / 16.0,
--               0.0,
--              -75.0 / 64.0,
--               75.0 / 64.0,
--               0.0, 0.0, 0.0, 0.0,
--               0.0, 0.0, 0.0, 0.0, others => 0.0),
--
--  
--              -- coefficients A(6)(0..12) for getting K(6)
--              -- so K6 = DeltaT*F (Time + Dt(6),Y+SUM(A6,K))
--              (3.0 / 80.0,
--               0.0,
--               0.0,
--               3.0 / 16.0,
--               3.0 / 20.0,
--               0.0, 0.0, 0.0,
--               0.0, 0.0, 0.0, 0.0, others => 0.0),
--
--  
--              -- coefficients A(7)(0..12) for getting K(7)
--              -- so K7 = DeltaT*F (Time + Dt(7),Y+SUM(A7,K))
--               (4.7910137111111111111111111111111E-02,
--               0.0,
--               0.0,
--               1.1224871277777777777777777777778E-01,
--              -2.5505673777777777777777777777778E-02,
--               1.2846823888888888888888888888889E-02,
--               0.0, 0.0, 0.0,
--               0.0, 0.0, 0.0, others => 0.0),
--
--  
--              -- coefficients A(8)(0..12) for getting K(8)
--              -- so K(8) = DeltaT*F (Time + Dt(8), Y + SUM (A(8), K))
--               (1.6917989787292281181431107136E-02,
--                0.0,
--                0.0,
--                3.87848278486043169526545744159E-01,
--                3.59773698515003278967008896348E-02,
--                1.96970214215666060156715256072E-01,
--               -1.72713852340501838761392997002E-01,
--                0.0,
--                0.0, 0.0, 0.0, 0.0, others => 0.0),
--
--  
--              -- coefficients A(9)(0..12) for getting K(9)
--              -- so K(9) = DeltaT*F (Time + Dt(9),Y + SUM (A(9), K))
--               (6.90957533591923006485645489846E-02,
--                0.0,
--                0.0,
--                -6.34247976728854151882807874972E-01,
--                -1.61197575224604080366876923982E-01,
--                1.38650309458825255419866950133E-01,
--                9.4092861403575626972423968413E-01,
--                2.11636326481943981855372117132E-01,
--                0.0,
--                0.0, 0.0, 0.0, others => 0.0),
--
--  
--              -- coefficients A(10)(0..12) for getting K(10)
--              -- so K(10) = DeltaT*F (Time + Dt(10),Y + SUM (A(10), K))
--                (1.83556996839045385489806023537E-01,
--                 0.0,
--                 0.0,
--                -2.46876808431559245274431575997E+00,
--                -2.91286887816300456388002572804E-01,
--                -2.6473020233117375688439799466E-02,
--                 2.84783876419280044916451825422E+00,
--                 2.81387331469849792539403641827E-01,
--                 1.23744899863314657627030212664E-01,
--                 0.0, 0.0, 0.0, others => 0.0),
--
--  
--              -- coefficients A(11)(0..12) for getting K(11)
--              -- so K(11) = DeltaT*F (Time + Dt(11),Y + SUM (A(11), K))
--               (-1.21542481739588805916051052503E+00,
--                 0.0,
--                 0.0,
--                 1.66726086659457724322804132886E+01,
--                 9.15741828416817960595718650451E-01,
--                -6.05660580435747094755450554309E+00,
--                -1.60035735941561781118417064101E+01,
--                1.4849303086297662557545391898E+01,
--                -1.33715757352898493182930413962E+01,
--                 5.13418264817963793317325361166E+00,
--                 0.0, 0.0, others => 0.0),
--
--  
--              -- coefficients A(12)(0..12) for getting K(12)
--              -- so K(12) = DeltaT*F (Time + Dt(12),Y + SUM (A(12), K))
--              (2.58860916438264283815730932232E-01,
--                0.0,
--                0.0,
--               -4.77448578548920511231011750971E+00,
--               -4.3509301377703250944070041181E-01,
--               -3.04948333207224150956051286631E+00,
--                5.57792003993609911742367663447E+00,
--                6.15583158986104009733868912669E+00,
--               -5.06210458673693837007740643391E+00,
--                2.19392617318067906127491429047E+00,
--                1.34627998659334941535726237887E-01,
--                0.0, others => 0.0),
--
--  
--              -- coefficients A(13)(0..12) for getting K(13)
--              -- so K(13) = DeltaT*F (Time + Dt(13),Y + SUM (A(13), K))
--               (8.22427599626507477963168204773E-01,
--                0.0,
--                0.0,
--               -1.16586732572776642839765530355E+01,
--               -7.57622116690936195881116154088E-01,
--                7.13973588159581527978269282765E-01,
--                1.20757749868900567395661704486E+01,
--               -2.12765911392040265639082085897E+00,
--                1.99016620704895541832807169835E+00,
--               -2.34286471544040292660294691857E-01,
--                1.7589857770794226507310510589E-01,
--                0.0, others => 0.0));
end Runge_Coeffs_PD_8;
