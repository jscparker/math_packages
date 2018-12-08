
-----------------------------------------------------------------------
-- package body Predictor_1, 17th order Predictor-Corrector integrator
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
with pc_1_interface;
with Runge_8th;

package body Predictor_1 is

 package Coeffs is new pc_1_interface.Predictor_Corrector_Rules (Real);
 use Coeffs;

 package Maths is new Ada.Numerics.Generic_Elementary_Functions (Real);
 use Maths;

 package Runge is new
    Runge_8th (Real, Dyn_Index, Dynamical_Variable, F, "*", "+", "-", Norm);

 type Derivative_History is array(PC_Rule_Range) of Dynamical_Variable;

 ---------------
 -- Integrate --
 ---------------

 procedure Integrate
   (Final_State   :    out Dynamical_Variable;
    Final_Time    : in     Real;
    Initial_State : in     Dynamical_Variable;
    Initial_Time  : in     Real;
    No_Of_Steps   : in     Step_Integer)
 is
    Static_dt       : constant Real := (Final_Time-Initial_Time) / Real (No_Of_Steps);
    Final_Time_Test : constant Real := Abs (Final_Time-Initial_Time + 0.5*Static_dt);

    Delta_t, New_t, Previous_t : Real;
    New_Y_Corrector, Predict_Last, Correct_Last : Real;
    Predictor, Corrector : Integration_Rule;

    Previous_Y : Dynamical_Variable;
    New_Y      : Dynamical_Variable;

    Index_Of_Oldest_Deriv : PC_Rule_Range := PC_Rule_Range'First;
    Derivative_Storage : Derivative_History;

    --------------------------
    -- Vectors_times_Matrix --
    --------------------------

    subtype Matrix     is Derivative_History;
    subtype Matrix_Row is Dynamical_Variable;
    subtype Vector     is Integration_Rule;

    procedure Vectors_times_Matrix
      (Mat       : in Matrix;
       Vector1   : in Vector;
       Vector2   : in Vector;
       Product1 : out Matrix_Row;
       Product2 : out Matrix_Row) is

       Row_ID  : PC_Rule_Range'Base;
       Sum1, Sum2 : Real;

    begin

       for Col in Dyn_Index loop

          Row_ID  :=  PC_Rule_Range'First;

          Sum1 := 0.0; Sum2 := 0.0;
          for i in PC_Rule_Range loop
             Sum1 := Sum1 + Mat(Row_ID)(Col) * Vector1(Row_ID);
             Sum2 := Sum2 + Mat(Row_ID)(Col) * Vector2(Row_ID);
             Row_ID := Row_ID + 1;
          end loop;
          Product1(Col) := Sum1; Product2(Col) := Sum2;

       end loop;

    end Vectors_times_Matrix;

 begin

    Delta_t := Static_dt;

    for i in PC_Rule_Range loop
       Predictor(i) := Delta_t * Predictor_Rule(i);
       Corrector(i) := Delta_t * Corrector_Rule(i);
    end loop;
    New_Y_Corrector := Delta_t * Final_Step_Corrector;


    --  use Runge-Kutta to Initialize Derivative_Storage:

    Previous_Y := Initial_State;
    Previous_t := Initial_Time;
    Derivative_Storage(PC_Rule_Range'First) := F (Previous_t, Previous_Y);

    for I in PC_Rule_Range'First+1 .. PC_Rule_Range'Last loop

       New_t := Previous_t + Delta_t;

       if Abs (New_t - Initial_Time) > Final_Time_Test then -- tried to go too far
          Final_State := Previous_Y;
          return;
       end if;

       Runge.Integrate
         (Final_State        => New_Y,  -- output latest Y
          Final_Time         => New_t,
          Initial_State      => Previous_Y,
          Initial_Time       => Previous_t,
          No_Of_Steps        => 2);

       Previous_Y := New_Y;
       Previous_t := New_t;
       Derivative_Storage(i) := F (Previous_t, Previous_Y);

    end loop;
    Index_Of_Oldest_Deriv := PC_Rule_Range'First;

    --  Starting pt. is Previous_Y, and Previous_t:

    Time_Steps: loop

       New_t := Previous_t + Delta_t;

       exit Time_Steps when Abs (New_t - Initial_Time) > Final_Time_Test;

       Predict_and_Correct:
       declare
          Predictor_Delta_Y     : Dynamical_Variable;
          Corrector_Delta_Y     : Dynamical_Variable;
          Predicted_New_Y       : Dynamical_Variable;
          Almost_New_Y          : Dynamical_Variable;
          Error_in_New_Y        : Dynamical_Variable;
          Latest_Derivative_Val : Dynamical_Variable;
       begin
   
          Vectors_times_Matrix
            (Mat      => Derivative_Storage,
             Vector1  => Predictor,
             Vector2  => Corrector,
             Product1 => Predictor_Delta_Y,
             Product2 => Corrector_Delta_Y);
   
          --  Use Predicted Y to evaluate the derivative of Y:
   
          Predicted_New_Y       := Previous_Y + Predictor_Delta_Y;
          Latest_Derivative_Val := F (New_t, Predicted_New_Y);
   
          --  Complete calculation of the Corrected Y: New_Y. It can be done
          --  iteratively, but usually just one Iterations is best:
   
          Almost_New_Y := Previous_Y + Corrector_Delta_Y;
   
          for I in 1 .. No_Of_Corrector_Evaluate_Iterations loop
             New_Y := Almost_New_Y + New_Y_Corrector * Latest_Derivative_Val;
             Latest_Derivative_Val := F (New_t, New_Y);
          end loop;
   
          if Final_Corrector_Desired then
             New_Y := Almost_New_Y + New_Y_Corrector * Latest_Derivative_Val;
          end if;
   
          if Extrapolation_Desired then  -- Subtract away estimated error:
             Error_in_New_Y := Extrap_Factor*(New_Y - Predicted_New_Y);
             New_Y := New_Y - Error_in_New_Y;
          end if;
   
          --  Save F(t,Y) in Derivative_Storage (write over oldest F),
          --  and adjust the predictor-corrector rules to the right
          --  starting points:
   
          Derivative_Storage(Index_Of_Oldest_Deriv) := Latest_Derivative_Val;
   
       end Predict_and_Correct;


       Previous_Y := New_Y;
       Previous_t := New_t;
   
       if Index_Of_Oldest_Deriv = PC_Rule_Range'Last then
          Index_Of_Oldest_Deriv := PC_Rule_Range'First;
       else
          Index_Of_Oldest_Deriv := Index_Of_Oldest_Deriv +  1;
       end if;

       Predict_Last := Predictor(PC_Rule_Range'Last);
       for I in reverse PC_Rule_Range'First+1 .. PC_Rule_Range'Last loop
          Predictor(I) := Predictor(I-1);
       end loop;
       Predictor(PC_Rule_Range'First) := Predict_Last;

       Correct_Last := Corrector(PC_Rule_Range'Last);
       for I in reverse PC_Rule_Range'First+1 .. PC_Rule_Range'Last loop
          Corrector(I) := Corrector(I-1);
       end loop;
       Corrector(PC_Rule_Range'First) := Correct_Last;

    end loop Time_Steps;

    -- All done, go home:

    Final_State := Previous_Y;

 end Integrate;

end Predictor_1;
