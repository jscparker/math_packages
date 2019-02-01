
-----------------------------------------------------------------------
-- package body Runge_5th, 5th order Runge-Kutta
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
with Runge_Coeffs_pd_5;

package body Runge_5th is

 Real_Small : constant Real  := Two * Real'Small;

 package mth is new Ada.Numerics.Generic_Elementary_Functions(Real);
 use mth;

 package Prince_Dormand_Coeffs is new Runge_Coeffs_pd_5 (Real);
 use Prince_Dormand_Coeffs;

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
 --  Notice uses only the lower triangle of A; no diagonal elements of A.
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
          Sum :=  A(Stage2)(0) * K(l)(0)
                + A(Stage2)(1) * K(l)(1);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 3  =>

       Stage2 := 3;
       for l in Dyn_Index loop
          Sum :=  A(Stage2)(0) * K(l)(0)
               +  A(Stage2)(1) * K(l)(1)
               +  A(Stage2)(2) * K(l)(2);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 4 =>

       Stage2 := 4;
       for l in Dyn_Index loop
         Sum := A(Stage2)(0)  * K(l)(0) +
                A(Stage2)(1)  * K(l)(1) +
                A(Stage2)(2)  * K(l)(2) +
                A(Stage2)(3)  * K(l)(3);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 5 =>

       Stage2 := 5;
       for l in Dyn_Index loop
         Sum :=                     A(Stage2)(0)  * K(l)(0)  +
          A(Stage2)(1)  * K(l)(1) + A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3) + A(Stage2)(4)  * K(l)(4);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when 6 =>

       Stage2 := 6;
       for l in Dyn_Index loop
         Sum :=                     A(Stage2)(0)  * K(l)(0)  +
          A(Stage2)(1)  * K(l)(1) + A(Stage2)(2)  * K(l)(2)  +
          A(Stage2)(3)  * K(l)(3) + A(Stage2)(4)  * K(l)(4)  +
          A(Stage2)(5)  * K(l)(5);
          Next_Y(l) := Y(l) + Sum;
       end loop;

    when others =>

       null;

    end case;

 end Runge_Sum;

 pragma Inline (Runge_Sum);

 ----------------
 -- Runge_Sum0 --
 ----------------

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
       Sum := Zero;
       for n in RK_Range range 0..Stage2-1 loop --check: 0..Stages2-2?
          Sum := Sum + A(Stage2)(n) * K(l)(n);
       end loop;
       Next_Y(l) := Y(l) + Sum;
    end loop;
 end Runge_Sum0;

 ---------------
 -- Integrate --
 ---------------

 --  Integrate to 5th order.

 procedure Integrate
   (Final_State           :    out Dynamical_Variable;
    Final_Time            : in     Real;
    Initial_State         : in     Dynamical_Variable;
    Initial_Time          : in     Real;
    No_Of_Steps           : in     Step_Integer;
    Error_Control_Desired : in     Boolean       := False;
    Error_Tolerance       : in     Real          := 1.0E-6)
 is
    N : constant Real := Real (No_Of_Steps);
    Static_Delta_t    : constant Real := (Final_Time - Initial_Time) / N;
    Final_Time_Test     : constant Real := Final_Time - Static_Delta_t*Half;
    Error_Tolerance_Sqr : constant Real := Error_Tolerance**2;

    type Step_Type is new Natural;
    Step_id : Step_Type;

    Delta_t, Error      : Real;
    Present_t, Time1    : Real;
    Next_Deriv, Next_Y  : Dynamical_Variable;
    Y, Error_Y, Delta_Y : Dynamical_Variable;
    This_Is_The_Final_Time_Step : Boolean := False;

    --  Increments of Independent variable
    --  so K(13) = Delta_t*F (Time + Dt(13), Y + SUM (A(13), K))

    K : K_type;

    Dt  : Coefficient;

    --------------------------
    -- Fourth_Order_Delta_Y --
    --------------------------

    --  function to Sum Series For Delta Y efficiently

    function Fourth_Order_Delta_Y
       return Dynamical_Variable
    is
       Result : Dynamical_Variable;
    begin
       for l in Dyn_Index loop
       Result(l) :=
          B4(0)  * K(l)(0)  +
        --B4(1)  * K(l)(1)  +  -- B4(1) = 0
          B4(2)  * K(l)(2)  +
	  B4(3)  * K(l)(3)  +
          B4(4)  * K(l)(4)  +
	  B4(5)  * K(l)(5)  +
          B4(6)  * K(l)(6);
       end loop;
       return Result;
    end Fourth_Order_Delta_Y;

    -------------------------
    -- Fifth_Order_Delta_Y --
    -------------------------

    --  function to Sum Series For Delta Y efficiently

    function Fifth_Order_Delta_Y
       return Dynamical_Variable
    is
       Result : Dynamical_Variable;
    begin
       for l in Dyn_Index loop
       Result(l) :=
          B5(0)  * K(l)(0)  + 
	--B5(1)  * K(l)(1)  + -- 0 for the pd coeffs
          B5(2)  * K(l)(2)  + 
	  B5(3)  * K(l)(3)  +
          B5(4)  * K(l)(4)  + 
	  B5(5)  * K(l)(5);
        --B5(6)  * K(l)(6);   -- 0 for the pd coeffs
       end loop;
       return Result;
    end Fifth_Order_Delta_Y;

    -------------------------------
    --  Get_New_Y_to_Fifth_Order --
    -------------------------------

    procedure Get_New_Y_to_Fifth_Order
      (Y : in out Dynamical_Variable)
    is
    begin
       for l in Dyn_Index loop
       Y(l) := Y(l) +
          B5(0)  * K(l)(0)  + 
	--B5(1)  * K(l)(1)  + -- 0 for the pd coeffs
          B5(2)  * K(l)(2)  + 
	  B5(3)  * K(l)(3)  +
          B5(4)  * K(l)(4)  + 
	  B5(5)  * K(l)(5);
        --B5(6)  * K(l)(6);   -- 0 for the pd coeffs
       end loop;
    end Get_New_Y_to_Fifth_Order;

    -----------------------------
    --  Get_New_Delta_t_and_Dt --
    -----------------------------

    --  Modifies Delta_t and Dt(i).
    --  Make Error_Tolerance 4 times smaller than requested for safety:

    Inverse_Error_Tolerance : constant Real := Four / Error_Tolerance;

    procedure Get_New_Delta_t_and_Dt
      (Error     : in     Real;
       Delta_t   : in out Real;
       Dt        :    out Coefficient)
    is
       Error_Epsilon : constant Real := Error_Tolerance * 1.0E-6 + Real_Small;
       New_Delta_t_Factor : Real := One;
       Delta_t_Fractional_Change : Real := One;
    begin
       Delta_t_Fractional_Change := Four_Fifths * Exp (-One_Fifth * Log (
            Inverse_Error_Tolerance * (Error + Error_Epsilon)));

       if Delta_t_Fractional_Change < One_Sixteenth then
          New_Delta_t_Factor := One_Sixteenth;
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
    Step_id   := Step_Type'First;

    for i in Stages loop
       Dt(i) := Delta_t * C(i);
    end loop;

    Time_Steps: loop

      --  First get Delta_Y to 5th Order by calculating the
      --  Runge-Kutta corrections K.
      --
      --  K(Stages'First)  := Delta_t * F (Time, Y);
      --  for Stage in Stages'First+1 .. Stages'Last loop
      --    K(Stage) := Delta_t * F (Time + Dt(Stage),  Y + Sum (Stage));
      --  end loop;

      -- Special Case:    Stage = Stages'First;
      -- In the initial stage the prince-dormand method allows us to
      -- reuse the most recently calculated Next_Deriv (from the previous
      -- step). That's why it's effectively a 6 stage method, rather than 7.

       Make_New_Corrections_K:
       declare
          Next_t : Real := Present_t; --local to K calculation.
       begin
          --  optimization only for prince-dormand 5th order coeffs:
          if Error_Control_Desired or else Step_id = Step_Type'First then 
             Next_Deriv := F (Next_t, Y); 
          end if;
          for l in Dyn_Index loop
             K(l)(Stages'First) := Delta_t * Next_Deriv(l);
          end loop;
          Step_id := Step_id + 1;
    
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
          Get_New_Y_to_Fifth_Order (Y);
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

          Get_New_Y_to_Fifth_Order (Y);

          Present_t := Present_t + Static_Delta_t;
          exit Time_Steps when Abs (Final_Time-Present_t) < Abs (0.125*Static_Delta_t);
 
       else
 
          --  Error control desired, so first calculate error.
 
          Delta_Y := Fifth_Order_Delta_Y;               -- used below
          Error_Y := Delta_Y - Fourth_Order_Delta_Y;
          Error   := Norm (Error_Y);
          --  Error in 4th order Y really; scales as dt**5 (per step)

          --  Next increment Y and Time if error is OK.
 
          if Error <= Error_Tolerance then  --  error is OK.
 
             Time1 := Present_t + Delta_t;
 
             if Abs (Time1-Initial_Time) < Abs (Final_Time-Initial_Time) then
 
                --  Increment both Time and Y. Afterwards
                --  get the new step size Delta_t.
 
                Present_t := Time1;
                Y         := Y + Delta_Y; -- 5th order
 
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
end Runge_5th;
