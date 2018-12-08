
-----------------------------------------------------------------------
-- package body Runge_pc_2, Runge-Kutta integrator.
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

with runge_coeffs_pd_8;

package body Runge_pc_2 is

  package Coeffs is new Runge_Coeffs_pd_8 (Real);
  use Coeffs;

  type K_type1 is array (RK_Range) of Real;

  type K_type2 is array (Dyn_Index)  of K_type1;

  type K_type is array (0..1) of K_type2;

  type RK_Dynamical_Variable is array (0..1) of Dynamical_Variable;

  procedure Runge_Sum
    (Y      : in   RK_Dynamical_Variable;
     Next_Y : out  RK_Dynamical_Variable;
     Stage  : in Stages;
     K      : in K_type)
  is
     Stage2 : Stages;
     Sum  : Real;
  begin
     Stage2 := Stage;

     for j in 0..1 loop
     for l in Dyn_Index loop
        Sum  := 0.0;
        for n in RK_Range range 0..Stage2-1 loop
           Sum  := Sum + A(Stage2)(n)  * K(j)(l)(n);
        end loop;
        Next_Y(j)(l) :=  Y(j)(l) + Sum;
     end loop;
     end loop;
  end Runge_Sum;

  pragma Inline (Runge_Sum);

  function F1
    (Time : Real;
     Y1   : RK_Dynamical_Variable)
     return RK_Dynamical_Variable
  is
     Result : RK_Dynamical_Variable;
  begin
     Result(0) := Y1(1);
     Result(1) := F (Time, Y1(0));
     return Result;
  end F1;

  --  Integrate to eighth order using RKPD

  procedure Integrate
    (Final_Y              :    out Dynamical_Variable;
     Final_deriv_Of_Y     :    out Dynamical_Variable;
     Final_Time           : in     Real;
     Initial_Y            : in     Dynamical_Variable;
     Initial_deriv_Of_Y   : in     Dynamical_Variable;
     Initial_Time         : in     Real;
     No_of_Steps          : in     Real)
  is
     Delta_t : constant Real := (Final_Time - Initial_Time) / No_of_Steps;

     Present_t : Real;
     Y : RK_Dynamical_Variable;

     K : K_type;

     --  Increments of Independent variable
     --  so K(13) = Delta_t*F (Time + Dt(13), Y + SUM (A(13), K))

     Dt  : Coefficient;

     --  function to Sum Series For Delta Y efficiently


     function Seventh_Order_Delta_Y return RK_Dynamical_Variable is

       Result : RK_Dynamical_Variable;

     begin

        for j in 0..1 loop
        for l in Dyn_Index loop

        Result(j)(l) :=
          B7(0)  * K(j)(l)(0)   + B7(5)  * K(j)(l)(5)  +
          B7(6)  * K(j)(l)(6)   + B7(7)  * K(j)(l)(7)  +
          B7(8)  * K(j)(l)(8)   + B7(9)  * K(j)(l)(9)  +
          B7(10) * K(j)(l)(10)  + B7(11) * K(j)(l)(11);

        end loop;
        end loop;

        return Result;
     end Seventh_Order_Delta_Y;

      --  function to Sum Series For Delta Y efficiently

     function Eighth_Order_Delta_Y return RK_Dynamical_Variable is

       Result : RK_Dynamical_Variable;

     begin

        for j in 0..1 loop
        for l in Dyn_Index loop

        Result(j)(l) :=
          B8(0)  * K(j)(l)(0)   + B8(5)  * K(j)(l)(5)  +
          B8(6)  * K(j)(l)(6)   + B8(7)  * K(j)(l)(7)  +
          B8(8)  * K(j)(l)(8)   + B8(9)  * K(j)(l)(9)  +
          B8(10) * K(j)(l)(10)  + B8(11) * K(j)(l)(11) +
          B8(12) * K(j)(l)(12);

        end loop;
        end loop;

        return Result;
     end Eighth_Order_Delta_Y;

     --  function to Sum Series For Delta Y efficiently

     procedure Get_New_Y_to_Seventh_Order
       (Y : in out RK_Dynamical_Variable)
     is
     begin

        for j in 0..1 loop
        for l in Dyn_Index loop

        Y(j)(l) := Y(j)(l) +
          B7(0)  * K(j)(l)(0)   + B7(5)  * K(j)(l)(5)  +
          B7(6)  * K(j)(l)(6)   + B7(7)  * K(j)(l)(7)  +
          B7(8)  * K(j)(l)(8)   + B7(9)  * K(j)(l)(9)  +
          B7(10) * K(j)(l)(10)  + B7(11) * K(j)(l)(11);

        end loop;
        end loop;

     end Get_New_Y_to_Seventh_Order;

     procedure Get_New_Y_to_Eighth_Order
       (Y : in out RK_Dynamical_Variable) is
     begin

        for j in 0..1 loop
        for l in Dyn_Index loop

        Y(j)(l) := Y(j)(l) +
          B8(0)  * K(j)(l)(0)   + B8(5)  * K(j)(l)(5)  +
          B8(6)  * K(j)(l)(6)   + B8(7)  * K(j)(l)(7)  +
          B8(8)  * K(j)(l)(8)   + B8(9)  * K(j)(l)(9)  +
          B8(10) * K(j)(l)(10)  + B8(11) * K(j)(l)(11) +
          B8(12) * K(j)(l)(12);

        end loop;
        end loop;

     end Get_New_Y_to_Eighth_Order;

  begin

     Y(0)      := Initial_Y;
     Y(1)      := Initial_deriv_Of_Y;
     Present_t := Initial_Time;

     for i in Stages Loop
        Dt(i) := Delta_t * C(i);
     end loop;

     Time_Steps:
     for step in 1 .. Integer(No_Of_Steps) loop

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
           Next_Deriv, Next_Y : RK_Dynamical_Variable;
        begin
           Next_Deriv := F1 (Next_t, Y);
           for j in 0..1 loop
           for l in Dyn_Index loop
              K(j)(l)(Stages'First) := Delta_t * Next_Deriv(j)(l);
           end loop;
           end loop;
   
           for Stage in Stages'First+1 .. Stages'Last loop
   
              Runge_Sum (Y, Next_Y, Stage, K);
              Next_t     := Present_t + Dt(Stage);
              Next_Deriv := F1 (Next_t, Next_Y);
              for j in 0..1 loop
              for l in Dyn_Index loop
                 K(j)(l)(Stage)  := Delta_t * Next_Deriv(j)(l);
              end loop;
              end loop;
   
           end loop;
        end Make_New_Corrections_K;

        -- use globally updated K to get new Y:
        Get_New_Y_to_Eighth_Order (Y);  
      --Get_New_Y_to_Seventh_Order (Y);
        Present_t := Present_t + Delta_t; -- new time

     end loop Time_Steps;

     Final_Y          := Y(0);
     Final_deriv_Of_Y := Y(1);

  end Integrate;

end Runge_pc_2;
