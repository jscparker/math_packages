
-----------------------------------------------------------------------
-- package body Runge_8th, 8th order Prince and Dormand Runge-Kutta
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
-----------------------------------------------------------------------

with Ada.Numerics.Generic_Elementary_Functions;
with Runge_Coeffs_PD_8;

package body Runge_8th is

 Real_Small : constant Real  := Two * Real'Small;

 package mth is new Ada.Numerics.Generic_Elementary_Functions(Real);
 use mth;

 package Prince_Dormand_Coeffs is new Runge_Coeffs_PD_8 (Real);
 use Prince_Dormand_Coeffs;
 --  Package contains coefficients for the Prince-Dormand 8th
 --  order, 13 stage Runge Kutta method with error control.

 type K_type1 is array (RK_Range) of Real;

 type K_type is array (Dyn_Index) of K_type1;

 ---------------
 -- Runge_Sum --
 ---------------

 --  Function to multiply vector G times matrix K.
 --
 --  This works only for the RKPD coefficients given above.
 --  It's optimized to take into account all of
 --  Zero's in array A; review A in the package above.
 --  The general formula is:
 --
 --       Result := A(Stage)(0) * K(0);
 --       for j in 1..Stage-1 loop
 --          Result := Result + A(Stage)(j) * K(j);
 --       end loop;
 --
 --  The above loop can be used for arbitrary RK coefficients.
 --  See Runge_Sum0.
 --
 --  What follows is optimized for the Prince-Dormond Coefficients.

 procedure Runge_Sum
   (Y      : in  Dynamical_Variable;
    Next_Y : out  Dynamical_Variable;
    Stage  : in Stages;
    K      : in K_type)
 is
    Stage2 : Stages;
    Sum  : Real;
 begin
    case Stage is
    when 1 =>

       Stage2 := 1;
       for l in Dyn_Index loop
          Sum :=  A(Stage2)(0) * K(l)(0);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 2 =>

       Stage2 := 2;
       for l in Dyn_Index loop
          Sum :=  A(Stage2)(0) * K(l)(0) + A(Stage2)(1) * K(l)(1);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 3  =>

       Stage2 := 3;
       for l in Dyn_Index loop
          Sum :=  A(Stage2)(0) * K(l)(0) +  A(Stage2)(2) * K(l)(2);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 4 =>

       Stage2 := 4;
       for l in Dyn_Index loop
         Sum := A(Stage2)(0)  * K(l)(0)  +
       -- A(Stage2)(1)  * K(l)(1) +
          A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 5 =>

       Stage2 := 5;
       for l in Dyn_Index loop
         Sum :=                     A(Stage2)(0)  * K(l)(0)  +
       -- A(Stage2)(1)  * K(l)(1) + A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3) + A(Stage2)(4)  * K(l)(4);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 6 =>

       Stage2 := 6;
       for l in Dyn_Index loop
         Sum :=                     A(Stage2)(0)  * K(l)(0)  +
       -- A(Stage2)(1)  * K(l)(1) + A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3) + A(Stage2)(4)  * K(l)(4)  +
          A(Stage2)(5)  * K(l)(5);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 7 =>

       Stage2 := 7;
       for l in Dyn_Index loop
         Sum :=                     A(Stage2)(0)  * K(l)(0)  +
       -- A(Stage2)(1)  * K(l)(1) + A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3) + A(Stage2)(4)  * K(l)(4)  +
          A(Stage2)(5)  * K(l)(5) + A(Stage2)(6)  * K(l)(6);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 8 =>

       Stage2 := 8;
       for l in Dyn_Index loop
         Sum :=                     A(Stage2)(0)  * K(l)(0)  +
       -- A(Stage2)(1)  * K(l)(1) + A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3) + A(Stage2)(4)  * K(l)(4)  +
          A(Stage2)(5)  * K(l)(5) + A(Stage2)(6)  * K(l)(6)  +
          A(Stage2)(7)  * K(l)(7);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 9 =>

       Stage2 := 9;
       for l in Dyn_Index loop
         Sum :=                     A(Stage2)(0)  * K(l)(0)  +
       -- A(Stage2)(1)  * K(l)(1) + A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3) + A(Stage2)(4)  * K(l)(4)  +
          A(Stage2)(5)  * K(l)(5) + A(Stage2)(6)  * K(l)(6)  +
          A(Stage2)(7)  * K(l)(7) + A(Stage2)(8)  * K(l)(8);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 10 =>

       Stage2 := 10;
       for l in Dyn_Index loop
         Sum :=                     A(Stage2)(0)  * K(l)(0)  +
       -- A(Stage2)(1)  * K(l)(1) + A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3) + A(Stage2)(4)  * K(l)(4)  +
          A(Stage2)(5)  * K(l)(5) + A(Stage2)(6)  * K(l)(6)  +
          A(Stage2)(7)  * K(l)(7) + A(Stage2)(8)  * K(l)(8)  +
          A(Stage2)(9)  * K(l)(9);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 11 =>

       Stage2 := 11;
       for l in Dyn_Index loop
         Sum :=                     A(Stage2)(0)  * K(l)(0)  +
       -- A(Stage2)(1)  * K(l)(1) + A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3) + A(Stage2)(4)  * K(l)(4)  +
          A(Stage2)(5)  * K(l)(5) + A(Stage2)(6)  * K(l)(6)  +
          A(Stage2)(7)  * K(l)(7) + A(Stage2)(8)  * K(l)(8)  +
          A(Stage2)(9)  * K(l)(9) + A(Stage2)(10) * K(l)(10);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 12 =>

       Stage2 := 12;
       for l in Dyn_Index loop
         Sum :=                     A(Stage2)(0)  * K(l)(0)  +
       -- A(Stage2)(1)  * K(l)(1) + A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3) + A(Stage2)(4)  * K(l)(4)  +
          A(Stage2)(5)  * K(l)(5) + A(Stage2)(6)  * K(l)(6)  +
          A(Stage2)(7)  * K(l)(7) + A(Stage2)(8)  * K(l)(8)  +
          A(Stage2)(9)  * K(l)(9) + A(Stage2)(10) * K(l)(10);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when others =>

       null;

    end case;

 end Runge_Sum;

 pragma Inline (Runge_Sum);

 ----------------
 -- Runge_Sum1 --
 ----------------

 procedure Runge_Sum1
   (Y      : in   Dynamical_Variable;
    Next_Y : out  Dynamical_Variable;
    Stage  : in   Stages;
    K      : in   K_type)
 is
    Stage2 : Stages;
    Sum  : Real;
 begin

    case Stage is
    when 1 =>

       Stage2 := 1;
       for l in Dyn_Index loop
          Sum :=  A(Stage2)(0) * K(l)(0);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 2 =>

       Stage2 := 2;
       for l in Dyn_Index loop
          Sum :=  A(Stage2)(0) * K(l)(0) + A(Stage2)(1) * K(l)(1);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 3  =>

       Stage2 := 3;
       for l in Dyn_Index loop
          Sum :=  A(Stage2)(0) * K(l)(0) +  A(Stage2)(2) * K(l)(2);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 4 =>

       Stage2 := 4;
       for l in Dyn_Index loop
          Sum := Zero;
          for n in reverse RK_Range range 0..3 loop
             Sum := Sum + A(Stage2)(n) * K(l)(n);
          end loop;
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 5 =>

       Stage2 := 5;
       for l in Dyn_Index loop
          Sum := Zero;
          for n in reverse RK_Range range 0..4 loop
             Sum := Sum + A(Stage2)(n) * K(l)(n);
          end loop;
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 6 =>

       Stage2 := 6;
       for l in Dyn_Index loop
          Sum := Zero;
          for n in reverse RK_Range range 0..5 loop
             Sum := Sum + A(Stage2)(n) * K(l)(n);
          end loop;
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 7 =>

       Stage2 := 7;
       for l in Dyn_Index loop
          Sum := Zero;
          for n in reverse RK_Range range 0..6 loop
             Sum := Sum + A(Stage2)(n) * K(l)(n);
          end loop;
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 8 =>

       Stage2 := 8;
       for l in Dyn_Index loop
          Sum := Zero;
          for n in reverse RK_Range range 0..7 loop
             Sum := Sum + A(Stage2)(n) * K(l)(n);
          end loop;
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 9 =>

       Stage2 := 9;
       for l in Dyn_Index loop
          Sum := Zero;
          for n in reverse RK_Range range 0..8 loop
             Sum := Sum + A(Stage2)(n) * K(l)(n);
          end loop;
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 10 =>

       Stage2 := 10;
       for l in Dyn_Index loop
          Sum := Zero;
          for n in reverse RK_Range range 0..9 loop
             Sum := Sum + A(Stage2)(n) * K(l)(n);
          end loop;
          Next_Y(l) := Y(l) + Sum;
       end loop;
    when 11 =>

       Stage2 := 11;
       for l in Dyn_Index loop
          Sum := Zero;
          for n in reverse RK_Range range 0..10 loop
             Sum := Sum + A(Stage2)(n) * K(l)(n);
          end loop;
          Next_Y(l) := Y(l) + Sum;
          --Sum := Zero; Sum2  := 0.0;   n := RK_Range'first;
          --for p in RK_Range range 0..4 loop
           --  Sum := Sum + A(Stage2)(n) * K(l)(n);
           --  Sum2  := Sum2 + A(Stage2)(n+1)  * K(l)(n+1);
           --  n := n+2;
          --end loop;
          --Next_Y(l) := Y(l) + Sum +  Sum2 + A(Stage2)(10)  * K(l)(10);

         -- Sum := Zero; Sum2  := 0.0;
         -- for n in reverse RK_Range range 0..10 loop
         --    Sum := Sum + A(Stage2)(n) * K(l)(n);
         --    Sum2  := Sum2 + A(Stage2)(n) * K(l)(m+1)(n);
         -- end loop;
         -- Next_Y(l) := Y(l) + Sum;
         -- Result(l)(m+1) :=  Sum2;
         -- m := m+2;
       end loop;

    when 12 =>

       Stage2 := 12;
       for l in Dyn_Index loop
          Sum := Zero;
          for n in reverse RK_Range range 0..10 loop
             Sum := Sum + A(Stage2)(n) * K(l)(n);
          end loop;
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when others =>

       null;

    end case;

 end Runge_Sum1;

 pragma Inline (Runge_Sum1);

 ----------------
 -- Runge_Sum0 --
 ----------------

 --  General formula for getting corrections K.
 --  The Next_Y is the Y at which F(t, Y) is next evaluated.

 procedure Runge_Sum0
   (Y    : in  Dynamical_Variable;
    Next_Y : out  Dynamical_Variable;
    Stage : in Stages;
    K     : in K_type)
 is
    Stage2 : Stages;
    Sum  : Real;
 begin
    Stage2 := Stage;
    for l in Dyn_Index loop
       Sum := 0.0;
       for n in reverse RK_Range range 0..Stage2-1 loop -- Sum small 1st
          Sum := Sum + A(Stage2)(n) * K(l)(n);
       end loop;
       Next_Y(l) := Y(l) + Sum;
    end loop;
 end Runge_Sum0;

 ---------------
 -- Integrate --
 ---------------

 --  Integrate to eighth order using the Prince-Dormond Coefficients.

 procedure Integrate
   (Final_State           :    out Dynamical_Variable;
    Final_Time            : in     Real;
    Initial_State         : in     Dynamical_Variable;
    Initial_Time          : in     Real;
    No_Of_Steps           : in     Step_Integer;
    Error_Control_Desired : in     Boolean       := False;
    Error_Tolerance       : in     Real          := +1.0E-10)
 is
    N : constant Real := Real (No_Of_Steps);
    Static_Delta_t : constant Real := (Final_Time - Initial_Time) / N;

    Delta_t, Error : Real;
    Present_t, Time1 : Real;

    Y        : Dynamical_Variable;
    Error_Y  : Dynamical_Variable;
    Delta_Y  : Dynamical_Variable;
    This_Is_The_Final_Time_Step : Boolean := False;

    --  Increments of Independent variable
    --  so K(13) = Delta_t*F (Time + Dt(13), Y + SUM (A(13), K))

    K : K_type;

    Dt : Coefficient;

    ---------------------------
    -- Seventh_Order_Delta_Y --
    ---------------------------

    --  function to Sum Series For Delta Y efficiently
    --  Force it sum small terms 1st.  It hardly matters.

    function Seventh_Order_Delta_Y
       return Dynamical_Variable
    is
       Result : Dynamical_Variable;
       Sum : Real;
    begin
       for l in Dyn_Index loop
          Sum :=       B7(12) * K(l)(12);
          Sum := Sum + B7(11) * K(l)(11);
          Sum := Sum + B7(10) * K(l)(10);
          Sum := Sum + B7(9) * K(l)(9);
          Sum := Sum + B7(8) * K(l)(8);
          Sum := Sum + B7(5) * K(l)(5);
          Sum := Sum + B7(7) * K(l)(7);
          Sum := Sum + B7(6) * K(l)(6);
          Sum := Sum + B7(0) * K(l)(0);
          Result(l) := Sum;
       end loop;
       return Result;
    end Seventh_Order_Delta_Y;

    --------------------------
    -- Eighth_Order_Delta_Y --
    --------------------------

    --  function to Sum Series For Delta Y efficiently
    --  Force it sum small terms 1st.  It hardly matters.

    function Eighth_Order_Delta_Y
       return Dynamical_Variable
    is
       Result : Dynamical_Variable;
       Sum : Real;
    begin
       for l in Dyn_Index loop
          Sum :=       B8(12) * K(l)(12);
          Sum := Sum + B8(11) * K(l)(11);
          Sum := Sum + B8(10) * K(l)(10);
          Sum := Sum + B8(9) * K(l)(9);
          Sum := Sum + B8(8) * K(l)(8);
          Sum := Sum + B8(5) * K(l)(5);
          Sum := Sum + B8(7) * K(l)(7);
          Sum := Sum + B8(6) * K(l)(6);
          Sum := Sum + B8(0) * K(l)(0);
          Result(l) := Sum;
       end loop;
       return Result;
    end Eighth_Order_Delta_Y;

    --------------------------------
    --  Get_New_Y_to_Eighth_Order --
    --------------------------------

    --  Force it sum small terms 1st.  It hardly matters.

    procedure Get_New_Y_to_Eighth_Order
      (Y : in out Dynamical_Variable)
    is
       Sum : Real;
    begin
       for l in Dyn_Index loop
          Sum :=       B8(12) * K(l)(12);
          Sum := Sum + B8(11) * K(l)(11);
          Sum := Sum + B8(10) * K(l)(10);
          Sum := Sum + B8(9) * K(l)(9);
          Sum := Sum + B8(8) * K(l)(8);
          Sum := Sum + B8(5) * K(l)(5);
          Sum := Sum + B8(7) * K(l)(7);
          Sum := Sum + B8(6) * K(l)(6);
          Sum := Sum + B8(0) * K(l)(0);
          Y(l) := Y(l) + Sum;
       end loop;
    end Get_New_Y_to_Eighth_Order;

    -----------------------------
    --  Get_New_Delta_t_and_Dt --
    -----------------------------

    --  Modifies Delta_t and Dt(i) as global variables.

    Inverse_Error_Tolerance : constant Real := One / Error_Tolerance;

    procedure Get_New_Delta_t_and_Dt
      (Error   : in     Real;
       Delta_t : in out Real;
       Dt      :    out Coefficient)
    is
       Error_Epsilon : constant Real := Error_Tolerance * 1.0E-6 + Real_Small;
       New_Delta_t_Factor : Real := One;
       Delta_t_Fractional_Change : Real := One;
    begin
       Delta_t_Fractional_Change := Nine_Tenths * Exp (-One_Eighth * Log (
            Inverse_Error_Tolerance * (Error + Error_Epsilon)));

       if Delta_t_Fractional_Change < One_Eighth then
          New_Delta_t_Factor := One_Eighth;
       elsif  Delta_t_Fractional_Change > Four then
          New_Delta_t_Factor := Four;
       else
          New_Delta_t_Factor := Delta_t_Fractional_Change;
       end if;

       Delta_t := New_Delta_t_Factor * Delta_t;
       for i in Stages Loop
          Dt(i) := Delta_t * C(i);
       end loop;
    end Get_New_Delta_t_and_Dt;

 begin

    Y         := Initial_State;
    Present_t := Initial_Time;
    Delta_t   := Static_Delta_t;

    for i in Stages loop
       Dt(i) := Delta_t * C(i);
    end loop;

    Time_Steps: loop

       --  First get DeltaY to 8th Order by calculating the
       --  Runge-Kutta corrections K.
       --
       --  K(Stages'First)  := Delta_t * F (Time, Y);
       --  for Stage in Stages'First+1 .. Stages'Last loop
       --    K(Stage) := Delta_t * F (Time + Dt(Stage),  Y + Sum (Stage));
       --  end loop;
 
       Make_New_Corrections_K:
       declare
          Next_t : Real := Present_t;
          Next_Deriv, Next_Y : Dynamical_Variable;
       begin
          Next_Deriv := F (Next_t, Y);
          for l in Dyn_Index loop
             K(l)(Stages'First) := Delta_t * Next_Deriv(l);
          end loop;
    
          for Stage in Stages'First+1 .. Stages'Last loop
             Runge_Sum (Y, Next_Y, Stage, K);
             Next_t     := Present_t + Dt(Stage);
             Next_Deriv := F (Next_t, Next_Y);
             for l in Dyn_Index loop
                K(l)(Stage) := Delta_t * Next_Deriv(l);
             end loop;
          end loop;
       end Make_New_Corrections_K;
    
       if This_Is_The_Final_Time_Step then --used only if error control is enabled.
          Get_New_Y_to_Eighth_Order (Y);
          exit Time_Steps;
       end if;
 
       --  Now increment Y and Time, and if desired, Delta_t.
       --  There are two algorithms below: with and
       --  without error control.
 
       if not Error_Control_Desired then
 
          --  Increment time and sum the Runge-Kutta series to increment Y.
          --  With the new Y, we can exit if Time is very near Final_Time.
          --  (Notice that Delta_t is negative if we integrate back in time.)
 
           -- use globally updated K to get new Y:

          Get_New_Y_to_Eighth_Order (Y);
 
          Present_t := Present_t + Static_Delta_t;
          exit Time_Steps when Abs (Final_Time-Present_t) < Abs (0.125*Static_Delta_t);
 
       else
 
          --  Error control desired, so first calculate error.
 
          Delta_Y := Eighth_Order_Delta_Y;
          Error_Y := Delta_Y - Seventh_Order_Delta_Y;
          Error   := Norm (Error_Y);
          --  Error in 7th order Y really; scales as dt**8
 
          --  Next increment Y and Time if error is OK.
 
          if Error <= Error_Tolerance then  --  error is OK.
 
             Time1 := Present_t + Delta_t;
 
             if Abs (Time1-Initial_Time) < Abs (Final_Time-Initial_Time) then
 
                --  Increment both Time and Y. Afterwards
                --  get the new step size Delta_t.
 
                Present_t := Time1;
                Y         := Y + Delta_Y;
 
             elsif Abs (Time1-Initial_Time) > Abs (Final_Time-Initial_Time) then
 
                --  Have to go through the loop again even
                --  though the error was small, because overshot
                --  the end.  Decrease Delta_t here, so there's
                --  no need to check error again after final Loop.
 
                This_Is_The_Final_Time_Step := True;
 
                Delta_t  := Final_Time - Present_t;
                for i in Stages loop
                   Dt(i) := Delta_t * C(i);
                end loop;
 
             else
 
                 --  Time1 = Final_Time, to just about the maximum accuracy
                 --  of the floating point, so get the final Y and exit.
 
                 Y := Y + Delta_Y;
                 exit Time_Steps;
 
             end if;
          end if;
 
          --  If this isn't the final loop, then want to adjust Delta_t.
          --  We want to make Delta_t smaller if necessary for accuracy, and
          --  larger if possible for speed.  If function ** is too slow then
          --  modify code so that this stuff is done only when we flunk
          --  the error test above.  But as is, it is done each time step.
 
          if not This_Is_The_Final_Time_Step then
             Get_New_Delta_t_and_Dt (Error, Delta_t, Dt);
          end if;
 
       end if; -- not error_control_desired
 
    end loop Time_Steps;

    Final_State := Y;

 end Integrate;

begin
   if Test_of_Runge_Coeffs_Desired then
      Prince_Dormand_Coeffs.Test_Runge_Coeffs;
   end if;
end Runge_8th;
