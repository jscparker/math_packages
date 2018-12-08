

with e_Derivs;
with Ada.Numerics.Generic_Elementary_Functions;
with Extended_Real;
with Extended_Real.Elementary_Functions;
with Extended_Real.IO;
with Text_IO; use Text_IO;

procedure e_deriv_tst_1 is

   type Real_8 is digits 15;

   package mth is new Ada.Numerics.Generic_Elementary_Functions (Real_8);
   use mth;
   package ext is new Extended_Real (Real_8);
   use ext;
   package fnc is new ext.Elementary_Functions (Sqrt, Log, Exp, Arcsin);
   use fnc;
   package eio is new ext.IO; -- extented IO
   use eio;


   subtype Real is e_Real;

   Max_Order_Of_Deriv : constant Positive := 40;

   package dif is new e_Derivs (Max_Order_Of_Deriv, Real, Real_8);
   use dif;

   w, Phase     : Real    := +1.0;
   Time, t      : Real;

   FullDuration : constant Real    := +1.0;
   No_Of_Steps  : constant Integer := 100;

   DeltaT       : Real    := FullDuration / (+Real_8(No_Of_Steps));

   G, H : Derivatives := (others => +0.0);

   Deriv0, Deriv1, Deriv2, Deriv3, Deriv4, Deriv5, Deriv6 : Real;
   Error : Array(Deriv_Index) of Real := (others => +0.0);
   Max_Error : Real;

   Three  : constant Real := +3.0;

begin

-- REMEMBER, Functions return REDUCED derivatives.

--*********************************************************************
-- Test 1.  H =  Exp_d(Sin_d(t))
--*********************************************************************
--  H = f(g(t)) =         Exp(Sin(t))
-- d^1 H        =  Cos(t)*Exp(Sin(t))
-- d^2 H        = (-Sin(t) + Cos(t)**2) * Exp(Sin(t))
-- d^3 H        = (-Cos(t) -3*Sin(t)*Cos(t) + Cos(t)**3) * Exp(Sin(t))
-- d^4 H        = ((-Cos(t) -3*Sin(t)*Cos(t) + Cos(t)**3) * Cos(t) +
--           (Sin(t) - 3*Cos(t)**2 + 3*Sin(t)**2 - 3*Sin(t)*Cos(t)**2))*Exp(Sin(t))

    Max_Error := Zero;

    for I in 0..No_Of_Steps loop

        Time := (+Real_8 (I)) * DeltaT;
        t    := Time;

        G := Sin_d (t);
        H := Compose (Exp_d (G(0)), G);  -- Exp (Sin (t))

        Deriv0 := Exp (Sin (t));

        Deriv1 := Cos(t) * Exp(Sin(t));

        Deriv2 :=  (-Sin(t) + Cos(t)**2) * Exp(Sin(t));

        Deriv3   := (-Cos(t) - Three*Sin(t)*Cos(t) + Cos(t)**3) * Exp(Sin(t));

        Deriv4 := ((-Cos(t) -Three*Sin(t)*Cos(t) + Cos(t)**3) * Cos(t) +
         (Sin(t) - Three*Cos(t)**2 + Three*Sin(t)**2
	         - Three*Sin(t)*Cos(t)**2))*Exp(Sin(t));

        Un_Reduce (H);

        Error(0) := Abs (H(0)  - Deriv0);
        Error(1) := Abs (H(1)  - Deriv1);
        Error(2) := Abs (H(2)  - Deriv2);
        Error(3) := Abs (H(3)  - Deriv3);
        Error(4) := Abs (H(4)  - Deriv4);

        for I in 0..4 loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
        end loop;

     end loop;
     new_line; put("Max error, test 1:  "); 
     put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 2.  H = Sin_d(t**3)
--*********************************************************************
--  H = f(g(t)) =  Sin(t**3)
-- d^1 H        = 3*t**2  * Cos(t**3)
-- d^2 H        = 6*t     * Cos(t**3) - 9*t**4  *  Sin(t**3)

    Max_Error := Zero;

    for I in 0..No_Of_Steps loop

       Time := (+Real_8 (I)) * DeltaT;

       G := Time ** 3;
       H := Compose (Sin_d (G(0)), G);   -- Sin (Time**3)

       Deriv0 := Sin(Time**3);
       Deriv1 := (+3.0)*Time**2 * COS(Time**3);
       Deriv2 := (+6.0)*Time*Cos(Time**3) - (+9.0)*Time**4*Sin(Time**3);

       Un_Reduce (H);

       Error(0) := Abs (H(0) -  Deriv0);
       Error(1) := Abs (H(1) -  Deriv1);
       Error(2) := Abs (H(2) -  Deriv2);

       for I in 0..2 loop
        if Error(I) > Max_Error then
           Max_Error := Error(I);
        end if;
       end loop;

    end loop;
    new_line; put("Max error, test 2:  "); 
    put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 3.  H = Sin_d(t)**5
--*********************************************************************
--  H = g(t)**5 =  Sin(t)**5
-- d^1 H        =  5 * Cos(t) * Sin(t)**4
-- d^2 H        =  20 * Cos(t)**2 * Sin(t)**3 - 5 * Sin(t)**5
-- d^3 H        = -40 * Cos(t) * Sin(t)**4  + 60 * Cos(t)**3 * Sin(t)**2
--                                               - 25 * Cos(t) * Sin(t)**4

    Max_Error := Zero;

    for I in 0..No_Of_Steps loop

      Time := (+Real_8 (I)) * DeltaT; t := Time;

      H := Sin_d(Time) ** 5;

      Deriv0 := Sin(t)**5;
      Deriv1 :=  (+5.0) * Cos(t) * Sin(t)**4 ;
      Deriv2 := (+20.0) * Cos(t)**2 * Sin(t)**3 - (+5.0) * Sin(t)**5;
      Deriv3 := (-40.0) * Cos(t) * Sin(t)**4 + (+60.0) * Cos(t)**3 * Sin(t)**2
                                             - (+25.0) * Cos(t) * Sin(t)**4;

      Un_Reduce (H);

      Error(0) := Abs (H(0) -  Deriv0);
      Error(1) := Abs (H(1) -  Deriv1);
      Error(2) := Abs (H(2) -  Deriv2);
      Error(3) := Abs (H(3) -  Deriv3);

      for I in 0..3 loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
      end loop;

   end loop;
   new_line; put("Max error, test 3:  ");
   put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 4.  H = t**4
--*********************************************************************
--  H =         =  t**4
-- d^1 H        =  4 * t**3
-- d^2 H        =  12 * t**2
-- d^3 H        =  24 * t
-- d^4 H        =  24
-- d^5 H        =  0
-- d^6 H        =  0

    Max_Error := Zero;

    for I in 0..No_Of_Steps loop

      Time := (+Real_8 (I)) * DeltaT; t := Time;

      H := t ** 4;

      Deriv0 := t ** 4;
      Deriv1 := (+4.0) * t**3;
      Deriv2 := (+12.0) * t**2;
      Deriv3 := (+24.0) * t;
      Deriv4 := +24.0;
      Deriv5 := Zero;
      Deriv6 := Zero;

      Un_Reduce (H);

      Error(0) := Abs (H(0) -  Deriv0);
      Error(1) := Abs (H(1) -  Deriv1);
      Error(2) := Abs (H(2) -  Deriv2);
      Error(3) := Abs (H(3) -  Deriv3);
      Error(4) := Abs (H(4) -  Deriv4);
      Error(5) := Abs (H(5) -  Deriv5);
      Error(6) := Abs (H(6) -  Deriv6);

      for I in 0..6 loop
       if Error(I) > Max_Error then
          Max_Error := Error(I);
       end if;
      end loop;

   end loop;
   new_line; put("Max error, test 4:  "); 
   put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 5.   H = Sin_d (w * t + phase)
--*********************************************************************
--     H        =       Sin (w * t + phase)
-- d^1 H        =   w**1 * Cos (w * t + phase)
-- d^2 H        =  -w**2 * Sin (w * t + phase)
-- d^3 H        =  -w**3 * Cos (w * t + phase)
-- d^4 H        =   w**4 * Sin (w * t + phase)
-- d^5 H        =   w**5 * Cos (w * t + phase)
-- d^6 H        =  -w**6 * Sin (w * t + phase)

    Max_Error := Zero;
    w         := +0.92345;
    phase     := +0.34567;

    for I in 0..No_Of_Steps loop

        Time := (+Real_8 (I)) * DeltaT; t := Time;

        H := Sin_d (t, w, phase);

        Deriv0 := Sin (w * t + phase);
        Deriv1 :=  w**1 * Cos (w * t + phase);
        Deriv2 := -w**2 * Sin (w * t + phase);
        Deriv3 := -w**3 * Cos (w * t + phase);
        Deriv4 :=  w**4 * Sin (w * t + phase);
        Deriv5 :=  w**5 * Cos (w * t + phase);
        Deriv6 := -w**6 * Sin (w * t + phase);

        Un_Reduce (H);

        Error(0) := Abs (H(0) -  Deriv0);
        Error(1) := Abs (H(1) -  Deriv1);
        Error(2) := Abs (H(2) -  Deriv2);
        Error(3) := Abs (H(3) -  Deriv3);
        Error(4) := Abs (H(4) -  Deriv4);
        Error(5) := Abs (H(5) -  Deriv5);
        Error(6) := Abs (H(6) -  Deriv6);

        for I in 0..6 loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
        end loop;

     end loop;
     new_line; put("Max error, test 5:  "); 
     put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 6.   H = Log_d (w * t + phase).
-- This also test Reciprocal : Real -> Derivative.
--*********************************************************************

    Max_Error := Zero;
    w         := +0.22345;
    phase     := +0.34567;

    for I in 0..No_Of_Steps loop

        Time := (+Real_8 (I)) * DeltaT; t := Time;

        H := Log_d (t, w, phase);

        Deriv0 := Log (w * t + phase);
        Deriv1 :=            w**1  / (w * t + phase);
        Deriv2 :=           -w**2  / (w * t + phase)**2;
        Deriv3 :=   (+2.0) * w**3  / (w * t + phase)**3;
        Deriv4 :=   (-6.0) * w**4  / (w * t + phase)**4;
        Deriv5 :=  (+24.0) * w**5  / (w * t + phase)**5;
        Deriv6 := (-120.0) * w**6  / (w * t + phase)**6;

        Un_Reduce (H);

        Error(0) := Abs (H(0) -  Deriv0);
        Error(1) := Abs (H(1) -  Deriv1);
        Error(2) := Abs (H(2) -  Deriv2);
        Error(3) := Abs (H(3) -  Deriv3);
        Error(4) := Abs (H(4) -  Deriv4);
        Error(5) := Abs (H(5) -  Deriv5);
        Error(6) := Abs (H(6) -  Deriv6);

        for I in 0..6 loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
        end loop;

     end loop;
     new_line; put("Max error, test 6:  "); 
     put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 7.   H = Cos_d (w * t + phase)
--*********************************************************************
    Max_Error := Zero;
    w         := +0.92345;
    phase     := +0.34567;

    for I in 0..No_Of_Steps loop

        Time := (+Real_8 (I)) * DeltaT; t := Time;

        H := Cos_d (t, w, phase);

        Deriv0 := Cos (w * t + phase);
        Deriv1 := -w**1 * Sin (w * t + phase);
        Deriv2 := -w**2 * Cos (w * t + phase);
        Deriv3 :=  w**3 * Sin (w * t + phase);
        Deriv4 :=  w**4 * Cos (w * t + phase);
        Deriv5 := -w**5 * Sin (w * t + phase);
        Deriv6 := -w**6 * Cos (w * t + phase);

        Un_Reduce (H);

        Error(0) := Abs (H(0) -  Deriv0);
        Error(1) := Abs (H(1) -  Deriv1);
        Error(2) := Abs (H(2) -  Deriv2);
        Error(3) := Abs (H(3) -  Deriv3);
        Error(4) := Abs (H(4) -  Deriv4);
        Error(5) := Abs (H(5) -  Deriv5);
        Error(6) := Abs (H(6) -  Deriv6);

        for I in 0..6 loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
        end loop;

     end loop;
     new_line; put("Max error, test 7:  "); 
     put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 8.   H = Sin_d / Cos_d;
--*********************************************************************
    Max_Error := Zero;

    for I in 0..No_Of_Steps loop

        Time := (+Real_8 (I)) * DeltaT; t := Time;

        H := Sin_d (t) / Cos_d (t);

        Deriv0 := Sin(t) / Cos(t);
        Deriv1 := (+1.0) / Cos(t)**2;
        Deriv2 := (+2.0) * Sin(t) / Cos(t)**3;
        Deriv3 := (+2.0) / Cos(t)**2 + (+6.0) * Sin(t)**2 / Cos(t)**4;
        Deriv4 := (+4.0) * Sin(t) / Cos(t)**3 + (+12.0) * Sin(t) / Cos (t)**3
                    + (+24.0) * Sin(t)**3 / Cos(t)**5;

        Un_Reduce (H);

        Error(0) := Abs (H(0) -  Deriv0);
        Error(1) := Abs (H(1) -  Deriv1);
        Error(2) := Abs (H(2) -  Deriv2);
        Error(3) := Abs (H(3) -  Deriv3);
        Error(4) := Abs (H(4) -  Deriv4);

        for I in 0..4 loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
        end loop;

     end loop;
     new_line; put("Max error, test 8:  "); 
     put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 9.   H = Sqrt_d;
--*********************************************************************
    Max_Error := Zero;
    DeltaT    := One / (+(No_Of_Steps + 1));

    for I in 2..No_Of_Steps loop

        Time := (+Real_8 (I)) * DeltaT; t := Time;

        H := Sqrt_d (t);

        Deriv0 := Sqrt(t);
        Deriv1 := (+0.5) / Sqrt (t);
        Deriv2 := -(+0.5)*(+0.5)/ ((t) * Sqrt (t));
        Deriv3 :=  (+0.5)*(+0.5)*(+0.5)*(+3.0)/ (t * t * Sqrt (t));
      --Deriv4 := -(+0.5)*(+0.5)*(+0.5)*(+0.5)*(+3.0)*(+5.0)/ (t * t * t * Sqrt (t));

        Un_Reduce (H);

        Error(0) := Abs (H(0) -  Deriv0);
        Error(1) := Abs (H(1) -  Deriv1);
        Error(2) := Abs (H(2) -  Deriv2);
        Error(3) := Abs (H(3) -  Deriv3);
      --Error(4) := Abs (H(4) -  Deriv4);

        for I in 0..3 loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
        end loop;

     end loop;
     new_line; put("Max error, test 9:  "); 
     put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 10.   H = Arcsin_d; tests Compose and "/".
--*********************************************************************
    Max_Error := Zero;
    DeltaT    := One / (+Real_8(No_Of_Steps + 7));

    for I in 1..No_Of_Steps loop

        Time := (+Real_8 (I)) * DeltaT;
        t    := Time;

        H := Arcsin_d (t);

        Deriv0 := Arcsin(t);
        Deriv1 := One / Sqrt (One - t*t);
        Deriv2 :=  t  / ((One - t*t) * Sqrt (One - t*t));
        Deriv3 := One / ((One - t*t) * Sqrt (One - t*t))
                  + (+3.0)*t*t / ((One - t*t)*(One - t*t) * Sqrt (One - t*t));

        Deriv4 := (+9.0)* t / ((One - t*t)*(One - t*t)*Sqrt (One - t*t))
       + (+3.0)*(+5.0)*t*t*t / ((One - t*t)*(One - t*t)*(One - t*t)*Sqrt (One - t*t));

        Un_Reduce (H);

        Error(0) := Abs (H(0) -  Deriv0);
        Error(1) := Abs (H(1) -  Deriv1);
        Error(2) := Abs (H(2) -  Deriv2);
        Error(3) := Abs (H(3) -  Deriv3);
      --Error(4) := Abs (H(4) -  Deriv4);

        for I in 0..3 loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
        end loop;

     end loop;
     new_line; put("Max error, test 10: "); 
     put (e_Real_Image (Max_Error, Aft => 20));


--*********************************************************************
-- Test 11.   H = Arctan_d; tests Compose and "/".
--*********************************************************************
    Max_Error := Zero;
    DeltaT    := One / (+(No_Of_Steps + 7));

    for I in 1..No_Of_Steps loop

        Time := (+Real_8(I)) * DeltaT;
        t    := Time;

        H := Arctan_d (t);

        Deriv0 := Arcsin (t / Sqrt (One + t*t));
        Deriv1 := One / (One + t*t);
        Deriv2 := -(+2.0) * t  / ((One + t*t)**2);
        Deriv3 := -(+2.0)  / ((One + t*t)**2)
                  + (+8.0)*t*t / ((One + t*t)**3);

        Un_Reduce (H);

        Error(0) := Abs (H(0) -  Deriv0);
        Error(1) := Abs (H(1) -  Deriv1);
        Error(2) := Abs (H(2) -  Deriv2);
        Error(3) := Abs (H(3) -  Deriv3);

        for I in 0..3 loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
        end loop;

     end loop;
     new_line; put("Max error, test 11: "); 
     put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 12.  H =  Exp_d (Log_d (t))
--
-- Test the reduced derivs.  Un-reduced get too large.
--*********************************************************************
--  H = f(g(t)) =  Exp_d (Log_d (t))
-- d^1 H        =  One
-- d^2 H        =  Zero

    Max_Error := Zero;

    for I in 0..No_Of_Steps loop

        Time := (+Real_8 (I)) * (+0.5) + (+5.0);
        t    := Time;

        G := Log_d (t);
        H := Compose (Exp_d (G(0)), G);  -- Exp (Log (t))

        Deriv0 := t;

        Deriv1 := One;

        Deriv2 := Zero;

        --Un_Reduce (H);

        Error(0) := Abs (H(0)  - Deriv0);
        Error(1) := Abs (H(1)  - Deriv1);
        Error(2) := Abs ((+2.0)*H(2)  - Deriv2);

        for I in 3..Max_Order_Of_Deriv loop
           Error(I) := Abs (H(I));
        end loop;

        for I in 0..Max_Order_Of_Deriv loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
        end loop;

     end loop;
     new_line; put("Max error, test 12: "); 
     put (e_Real_Image (Max_Error, Aft => 20));

--*********************************************************************
-- Test 13.  H =  Sin_d (Arcsin_d (t))
--
-- Test the reduced derivs.  Un-reduced get too large.
--*********************************************************************
--  H = f(g(t)) =  Sin_d (Arcsin_d (t))
--  d^1 H       =  One
--  d^2 H       =  Zero

    Max_Error := Zero;

    for I in 0..No_Of_Steps loop

        Time := (+Real_8 (I + 1)) * DeltaT / (+Real_8 (I + 2));
        t    := Time;

        G := Arcsin_d (t);
        H := Compose (Sin_d (G(0)), G);  -- Sin_d (Arcsin_d (t))

        Deriv0 := t;

        Deriv1 := One;

        Deriv2 := Zero;

        --Un_Reduce (H);

        Error(0) := Abs (H(0)  - Deriv0);
        Error(1) := Abs (H(1)  - Deriv1);
        Error(2) := Abs ((+2.0)*H(2)  - Deriv2);

        for I in 3..Max_Order_Of_Deriv loop
           Error(I) := Abs (H(I));
        end loop;

        for I in 0..Max_Order_Of_Deriv loop
         if Error(I) > Max_Error then
            Max_Error := Error(I);
         end if;
        end loop;

     end loop;
     new_line; put("Max error, test 13: "); 
     put (e_Real_Image (Max_Error, Aft => 20));

end;
