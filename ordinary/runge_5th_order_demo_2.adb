
with Runge_5th;
with Text_IO; use Text_IO;
with Orbit;

procedure runge_5th_order_demo_2 is

   type Real is digits 15;

   package One_Electron_Atom is new Orbit (Real);
   use One_Electron_Atom;

   package Orb_Integrate is new 
      Runge_5th (Real, Dyn_Index, Dynamical_Variable, F, "*", "+", "-", Norm);
   use Orb_Integrate;

   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Step_Integer);
   use iio;

   Initial_Condition         : Dynamical_Variable;
   Previous_Y, Final_Y       : Dynamical_Variable;
   Final_Time, Starting_Time : Real;
   Previous_t, Final_t       : Real;
   Delta_t                   : Real;
   Initial_Energy            : Real;

   Steps	         : Step_Integer  := 200;
   Error_Tolerance       : Real          := 1.0E-7;
   Error_Control_Desired : Boolean       := False;
   Answer                : Character     := 'n';

begin

   --  choose initial conditions

   new_line;
   put ("The test calculates the trajectory of a body in a highly elliptical orbit.");
   new_line;
   put ("During most of the orbit a large step-size is fine.  During the near-");
   new_line;
   put ("collision of the 2 bodies, a tiny step-size is necessary. The test");
   new_line;
   put ("demonstrates that the error control option (which uses variable step-size)");
   new_line;
   put ("is more efficient.");
   new_line(2);
   put ("Input number of steps (try 400_000 with and without error control): "); 
   new_line;
   get (Steps);
   new_line;
   put ("Every time the integration advances Delta_t = 4, ERROR is printed.");
   new_line;
   put ("Use error control? Enter y or n:"); new_line; get (Answer);
   if (Answer = 'Y') or (Answer = 'y') then
      Error_Control_Desired := True;
      put ("Error control it is. Program attempts to reduce error *per step* to: "); 
      put (Error_Tolerance, Aft => 6);
      new_line;
   else
      Error_Control_Desired := False;
      put ("OK, no error control."); new_line;
   end if;
 --if Error_Control_Desired then
    --put ("Input error tolerance:"); 
    --new_line; get (Error_Tolerance);
 --end if;

   Initial_Condition(0) :=  0.0;
   Initial_Condition(1) :=  1.0;
 --Initial_Condition(2) :=  2.0;  -- circular orbit
   Initial_Condition(2) :=  0.20; -- highly ellitical orbit
   Initial_Condition(3) :=  0.0;

   Starting_Time   :=  0.0;
   Final_Time      :=  4.0;
   Initial_Energy  :=  Energy (Initial_Condition);

   Previous_Y := Initial_Condition;
   Previous_t := Starting_Time;
   Delta_t    := Final_Time - Starting_Time;

   for i in 1..30 loop

      Final_t := Previous_t + Delta_t;

      Integrate
        (Final_State           => Final_Y,   -- the result (output).
         Final_Time            => Final_t,   -- end integration here.
         Initial_State         => Previous_Y, -- the initial condition (input).
         Initial_Time          => Previous_t, -- start integrating here.
         No_Of_Steps           => Steps,      -- if no err control, this is no_of_steps
         Error_Control_Desired => Error_Control_Desired,
         Error_Tolerance       => Error_Tolerance);

      Previous_Y := Final_Y;
      Previous_t := Final_t;

      new_line;
      put ("Time = t =");
      put (Final_t, Aft => 7);
      put (",     Error = (True Energy - Integrated) = ");
      put (Abs (Initial_Energy - Energy (Final_Y)), Aft => 7);

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
