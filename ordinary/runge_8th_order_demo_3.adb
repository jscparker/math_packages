
with Runge_8th;
with Text_IO; use Text_IO;
with Quadrature;
with Ada.Numerics.Generic_Elementary_Functions;

procedure runge_8th_order_demo_3 is

   type Real is digits 15;

   package mth is new Ada.Numerics.Generic_Elementary_Functions(Real);
   use mth; -- for Sin

   package quadr is new Quadrature (Real, Sin);
   use quadr;

   package Area_Under_the_Curve is new 
      Runge_8th (Real, Dyn_Index, Dynamical_Variable, F, "*", "+", "-", Norm);
   use Area_Under_the_Curve;

   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Step_Integer);
   use iio;

   Initial_Condition         : Dynamical_Variable;
   Previous_Y, Final_Y       : Dynamical_Variable;
   Final_Time, Starting_Time : Real;
   Previous_t, Final_t       : Real;
   Delta_t                   : Real;

   Steps	         : Step_Integer;
   Error_Tolerance       : Real          := 1.0E-10;
   Error_Control_Desired : Boolean       := False;
   Answer                : Character     := 'n';

begin

   --  choose initial conditions

   new_line;
   put ("The test integrates  (d/dt) Y(t) = Sin(t)  for Y(t).");
   new_line(2);
   put ("Ordinary quadrature is just an area-under-the-curve problem, ");
   new_line;
   put ("where you solve: (d/dt) Y(t) = Sin(t) for Y.");
   new_line(2);
   put ("Input number of steps (try 16 with and without error control)");
   new_line;
   get (Steps);
   new_line;
   put ("Every time the integration advances this number of steps, ERROR is printed.");
   new_line;
   new_line;
   put ("Use error control? Enter y or n:"); new_line; get (Answer);
   if (Answer = 'Y') or (Answer = 'y') then
      Error_Control_Desired := True;
      put ("Error control it is."); new_line;
   else
      Error_Control_Desired := False;
      put ("OK, no error control."); new_line;
   end if;

   Initial_Condition(0) :=  -1.0;

   Starting_Time   :=  0.0;
   Final_Time      :=  32.0;

   Previous_Y := Initial_Condition;
   Previous_t := Starting_Time;
   Delta_t    := Final_Time - Starting_Time;

   for i in 1..20 loop

      Final_t := Previous_t + Delta_t;

      Integrate
        (Final_State           => Final_Y,    -- the result (output).
         Final_Time            => Final_t,    -- end integration here.
         Initial_State         => Previous_Y, -- the initial condition (input).
         Initial_Time          => Previous_t, -- start integrating here.
         No_Of_Steps           => Steps,      -- if no err control, uses this.
         Error_Control_Desired => Error_Control_Desired,
         Error_Tolerance       => Error_Tolerance);

      Previous_Y := Final_Y;
      Previous_t := Final_t;

      new_line;
      put ("Time = t =");
      put (Final_t, Aft => 7);
      put (",     Error = (Cos (t) - Integrated Cos) = ");
      put (Abs (-Cos(Final_t) - Final_Y(0)), Aft => 7);

   end loop;

   if (Answer = 'Y') or (Answer = 'y') then
      new_line(2);
      put ("With error control enabled, program attempted to reduce");
      new_line;
      put ("error *per step* to (well) under:  "); 
      put (Error_Tolerance, Aft => 6);
      new_line;
      put ("Over thousands of steps, accumulated error will be much larger than that.");
      new_line(2);
   end if;

end;
