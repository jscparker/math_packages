
-----------------------------------------------------------------------
-- package body Predictor_2, 20th order Predictor-Corrector
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

with pc_2_interface;
with Runge_pc_2;

package body Predictor_2 is

 package Runge_pc is new Runge_pc_2 (Real, Dyn_Index, Dynamical_Variable, F);

 package Coeffs is new pc_2_interface.Predictor_Corrector_Rules (Real);
 use Coeffs;

 type Second_Deriv_History is array(PC_Rule_Range) of Dynamical_Variable;

 ---------------
 -- Integrate --
 ---------------

 procedure Integrate
   (Final_Y              :    out Dynamical_Variable;
    Final_deriv_Of_Y     :    out Dynamical_Variable;
    Final_Time           : in     Real;
    Initial_Y            : in     Dynamical_Variable;
    Initial_deriv_Of_Y   : in     Dynamical_Variable;
    Initial_Time         : in     Real;
    No_Of_Steps          : in     Real)
 is
    Delta_t         : constant Real := (Final_Time - Initial_Time) / No_Of_Steps;
    Final_Time_Test : constant Real := Abs (Final_Time-Initial_Time + 0.9*Delta_t);
    --  Notice Delta_t can be negative, (when itegrating backward in time).

    New_t, Previous_t   : Real;

    New_Y_Corrector, Predict_Last, Correct_Last, Integral_Last : Real;
    Predictor, Corrector             : Integration_Rule;
    Integrate_to_End, First_Integral : Integration_Rule;

    Previous_Deriv_of_Y_x_dt : Dynamical_Variable;
    Previous_Deriv_of_Y      : Dynamical_Variable;
    New_Deriv_of_Y_x_dt      : Dynamical_Variable;
    New_Deriv_of_Y           : Dynamical_Variable;
    New_Y, Previous_Y        : Dynamical_Variable;

    Index_Of_Oldest_Deriv : PC_Rule_Range := PC_Rule_Range'First;
    Second_Deriv_Storage : Second_Deriv_History;

    --------------------------
    -- Vectors_x_Matrix --
    --------------------------

    --  in applications where Predictor_Corrector's are appropriate,
    --  Second_Deriv_History is usually a giant array that spills out
    --  cache.  Fetching it from mem is slow, so below we fetch it
    --  once, rather than 3 times.

    subtype Matrix     is Second_Deriv_History;
    subtype Matrix_Row is Dynamical_Variable;
    subtype Vector     is Integration_Rule;

    procedure Vectors_x_Matrix
      (Mat       : in Matrix;
       Vector1   : in Vector;
       Vector2   : in Vector;
       Vector3   : in Vector;
       Product1 : out Matrix_Row;
       Product2 : out Matrix_Row;
       Product3 : out Matrix_Row) is

       Row_ID  : PC_Rule_Range'Base;
       Sum1, Sum2, Sum3 : Real;

    begin

       for Col in Dyn_Index loop

          Row_ID  :=  PC_Rule_Range'First;

          Sum1 := 0.0; Sum2 := 0.0; Sum3 := 0.0;
          for i in PC_Rule_Range loop
             Sum1 := Sum1 + Mat(Row_ID)(Col) * Vector1(Row_ID);
             Sum2 := Sum2 + Mat(Row_ID)(Col) * Vector2(Row_ID);
             Sum3 := Sum3 + Mat(Row_ID)(Col) * Vector3(Row_ID);
             Row_ID := Row_ID + 1;
          end loop;

          Product1(Col) := Sum1; Product2(Col) := Sum2; Product3(Col) := Sum3;

       end loop;

    end Vectors_x_Matrix;

    ---------------------------------
    -- First_Integral_of_2nd_deriv --
    ---------------------------------

    function First_Integral_of_2nd_deriv
      (New_Deriv_of_Y_x_dt : in Dynamical_Variable;
       Index_Of_Oldest_Deriv    : in PC_Rule_Range;
       Delta_t                  : in Real)
       return  Dynamical_Variable
    is
       Product1 : Dynamical_Variable;
       Shift, i  : PC_Rule_Range'Base;
       Sum1 : Real;
       Inverse_Delta_t : Real := 1.0 / Delta_t;
    begin
       Shift := PC_Rule_Range'First;
       for k in Index_Of_Oldest_Deriv .. PC_Rule_Range'Last loop
          Integrate_to_End(k) := Delta_t * Center_to_End_Integration(Shift);
          Shift := Shift + 1;
       end loop;

       if Index_Of_Oldest_Deriv > PC_Rule_Range'First then
       for k in PC_Rule_Range'First .. Index_Of_Oldest_Deriv-1 loop
          Integrate_to_End(k) := Delta_t * Center_to_End_Integration(Shift);
          Shift := Shift + 1;
       end loop;
       end if;

       for Col in Dyn_Index loop

          i  :=  PC_Rule_Range'First;

          Sum1 := 0.0;
          for k in PC_Rule_Range loop
             Sum1 := Sum1 + Second_Deriv_Storage(i)(Col) * Integrate_to_End(i);
             i := i + 1;
          end loop;

          Product1(Col) := Sum1;

       end loop;

       return Product1 + Inverse_Delta_t * New_Deriv_of_Y_x_dt;

    end First_Integral_of_2nd_deriv;

 begin -- Integrate

    -- always init out params:
    Final_Y          := Initial_Y;
    Final_deriv_Of_Y := Initial_deriv_Of_Y;

    for i in PC_Rule_Range loop
       Predictor(i) := Delta_t**2 * Predictor_Rule(i);
    end loop;
    for i in PC_Rule_Range loop
       Corrector(i) := Delta_t**2 * Corrector_Rule(i);
    end loop;
    New_Y_Corrector := Delta_t**2 * Final_Step_Corrector;

    for i in PC_Rule_Range loop
       First_Integral(i) := Delta_t**2 * Center_Integration(i);
    end loop;


    Previous_Y          := Initial_Y;
    Previous_deriv_Of_Y := Initial_deriv_Of_Y;
    Previous_t          := Initial_Time;

    Index_Of_Oldest_Deriv                     := PC_Rule_Range'First;
    Second_Deriv_Storage(PC_Rule_Range'First) := F (Previous_t, Previous_Y);

    for I in PC_Rule_Range'First+1 .. PC_Rule_Range'Last loop

       New_t := Previous_t + Delta_t;

       if Abs (New_t - Initial_Time) > Final_Time_Test then
          Final_Y          := Previous_Y;
          Final_deriv_Of_Y := Previous_deriv_Of_Y;
          return;
       end if;

       Runge_pc.Integrate
         (Final_Y            => New_Y,
          Final_deriv_Of_Y   => New_Deriv_of_Y,
          Final_Time         => New_t,
          Initial_Y          => Previous_Y,
          Initial_deriv_Of_Y => Previous_Deriv_of_Y,
          Initial_Time       => Previous_t,
          No_Of_Steps        => 4.0);
       --  Use 7th order method and small Delta_t for better accuracy?

       if I = Starting_Id_of_First_Deriv_of_Y then
          Previous_Deriv_of_Y_x_dt :=  Delta_t * New_Deriv_of_Y;
       end if;

       Second_Deriv_Storage(i) := F (New_t, New_Y);

       Previous_Y          := New_Y;
       Previous_deriv_Of_Y := New_deriv_Of_Y;
       Previous_t          := New_t;

    end loop;


    Time_Steps: Loop

       New_t :=  Previous_t + Delta_t;

       if Abs (New_t - Initial_Time) > Final_Time_Test then -- tried to go too far
          Final_Y          := Previous_Y;
          Final_deriv_Of_Y :=
            First_Integral_of_2nd_deriv (Previous_Deriv_of_Y_x_dt,
                                    Index_Of_Oldest_Deriv, Delta_t);
          exit Time_Steps;
       end if;
         
       Predict_and_Correct:
       declare
          Predictor_Delta_Y     : Dynamical_Variable;
          Corrector_Delta_Y     : Dynamical_Variable;
          Delta_Deriv_of_Y_x_dt : Dynamical_Variable;
          Predicted_New_Y       : Dynamical_Variable;
          Almost_New_Y          : Dynamical_Variable;
          Error_in_New_Y        : Dynamical_Variable;
          New_2nd_Deriv_Val     : Dynamical_Variable;
       begin

          Vectors_x_Matrix
            (Mat      => Second_Deriv_Storage,
             Vector1  => Predictor,
             Vector2  => Corrector,
             Vector3  => First_Integral,
             Product1 => Predictor_Delta_Y,
             Product2 => Corrector_Delta_Y,
             Product3 => Delta_Deriv_of_Y_x_dt);
   
          New_Deriv_of_Y_x_dt := Previous_Deriv_of_Y_x_dt + Delta_Deriv_of_Y_x_dt;
   
          Predicted_New_Y := Predictor_Delta_Y + Previous_Y + New_Deriv_of_Y_x_dt;
          New_2nd_Deriv_Val := F (New_t, Predicted_New_Y);
   
          Almost_New_Y := Corrector_Delta_Y + Previous_Y + New_Deriv_of_Y_x_dt;
   
          for I in 1 .. No_Of_Corrector_Evaluate_Iterations loop
             New_Y := Almost_New_Y + New_Y_Corrector * New_2nd_Deriv_Val;
             New_2nd_Deriv_Val := F (New_t, New_Y);
          end loop;
   
          if Final_Corrector_Desired then
             New_Y := Almost_New_Y + New_Y_Corrector * New_2nd_Deriv_Val;
          end if;
   
          if Extrapolation_Desired then
             Error_in_New_Y := Extrap_Factor * (New_Y - Predicted_New_Y);
             New_Y := New_Y - Error_in_New_Y;
          end if;
   
          Second_Deriv_Storage(Index_Of_Oldest_Deriv) := New_2nd_Deriv_Val;
   
       end Predict_and_Correct;


       --  rotate the rule arrays.  seems as fast as shifting the index
       --  during the dot product with the derivatives.

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

       Integral_Last := First_Integral(PC_Rule_Range'Last);
       for I in reverse PC_Rule_Range'First+1 .. PC_Rule_Range'Last loop
          First_Integral(I) := First_Integral(I-1);
       end loop;
       First_Integral(PC_Rule_Range'First) := Integral_Last;


       Previous_Deriv_of_Y_x_dt := New_Deriv_of_Y_x_dt;
       Previous_Y := New_Y;
       Previous_t := New_t;

    end loop Time_Steps;

 end Integrate;

end Predictor_2;
