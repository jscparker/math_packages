with Runge_8th;
with Text_IO; use Text_IO;
with Sinu;
with Ada.Numerics.Generic_Elementary_Functions;

procedure runge_8th_order_demo_1 is

   type Real is digits 15;

   package Sinusoid is new Sinu (Real);
   use Sinusoid;

   package Sin_Integrate is new 
      Runge_8th (Real, Dyn_Index, Dynamical_Variable, F, "*", "+", "-", Norm);
   use Sin_Integrate;

   package math is new Ada.Numerics.Generic_Elementary_Functions(Real);
   use math;
   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Step_Integer);
   use iio;

   Initial_Condition   : Dynamical_Variable;
   Previous_Y, Final_Y : Dynamical_Variable;
   Previous_t, Final_t : Real;
   Final_Time          : Real;
   Starting_Time       : Real;
   Delta_t             : Real;
   Steps	       : Step_Integer;
begin
   new_line;
   put ("Test integrates an equation whose solution is Sin(t).");
   new_line;
   put ("The equation is (d/dt)**2 (Y) = -Y. (Y = Sin(t) has a period of 2*Pi.)");
   new_line(2);
   put ("Enter number of time steps (try 1024): ");
   get (Steps);
   new_line;
   put ("Every time the integration advances this number of steps, ERROR is printed.");
   new_line;

   --  choose initial conditions  

   Initial_Condition(0) :=  0.0;
   Initial_Condition(1) :=  1.0;

   Starting_Time        :=  0.0;
   Final_Time           :=  64.0;

   Previous_Y := Initial_Condition;
   Previous_t := Starting_Time;
   Delta_t    := Final_Time - Starting_Time;

   for i in 1..28 loop

      Final_t := Previous_t + Delta_t;

      Integrate
         (Final_State        => Final_Y,   -- the result (output).
          Final_Time         => Final_t,   -- integrate out to here (input).
          Initial_State      => Previous_Y, -- input an initial condition (input).
          Initial_Time       => Previous_t, -- start integrating here (input).
          No_Of_Steps        => Steps);     -- if no err control, use this no_of_steps

      Previous_Y := Final_Y;
      Previous_t := Final_t;
      --  Start over.


    --put ("Final value of numerically integrated sin:  ");
    --put (Final_Y(0), Aft => 7);
      new_line;
      put ("Time = t =");
      put (Final_t, Aft => 7);
      put (",     Error = Sin(t)-Y = ");
      put (Sin (Final_t) - Final_Y(0), Aft => 7);

   end loop;

   new_line(2);
   put ("Total elapsed time is:");
   put (Final_t / (2.0*3.14159265), Aft => 7);
   put (" in units of 2 pi.");

end;
