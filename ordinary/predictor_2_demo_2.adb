
with Predictor_2;
with Text_IO; use Text_IO;
with Orbit_2;

procedure predictor_2_demo_2 is

   type Real is digits 15;

   package One_Electron_Atom is new Orbit_2 (Real);
   use One_Electron_Atom;

   package Orb_Integrate is new 
      Predictor_2 (Real, Dyn_Index, Dynamical_Variable, F, "*", "+", "-");
   use Orb_Integrate;

   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Integer);
   use iio;

   Initial_Y, Initial_Y_dot    : Dynamical_Variable;
   Previous_Y, Final_Y         : Dynamical_Variable;
   Previous_Y_dot, Final_Y_dot : Dynamical_Variable;
   Final_Time, Starting_Time : Real;
   Previous_t, Final_t       : Real;
   Delta_t                   : Real;
   Initial_Energy            : Real;

   Steps : Real  := 128.0;

begin

   --  choose initial conditions

   new_line;
   put ("The test calculates the trajectory of a body in a circular orbit.");
   new_line(2);
   put ("Enter number of time steps (try 512): ");
   get (Steps);
   new_line;
   put ("Every time the integration advances this no of steps, ERROR is printed.");
   new_line;

   Initial_Y(0)     :=  0.0; -- x
   Initial_Y(1)     :=  1.0; -- z
 --Initial_Y_dot(0) :=  2.0;  -- circular orbit, x_dot
   Initial_Y_dot(0) :=  1.8;  -- near circular orbit, x_dot
   Initial_Y_dot(1) :=  0.0;  -- z_dot

   Starting_Time   :=  0.0;
 --Final_Time      :=  3.14159_26535_89793_23846 * 10.0; --10 orbits, if orbit circular.
   Final_Time      :=  32.0;
   Initial_Energy  :=  Energy (Initial_Y, Initial_Y_dot);

   Previous_Y     := Initial_Y;
   Previous_Y_dot := Initial_Y_dot;
   Previous_t     := Starting_Time;

   Delta_t        := Final_Time - Starting_Time;

   for i in 1..30 loop

      Final_t := Previous_t + Delta_t;

      Integrate
        (Final_Y            => Final_Y,       -- the result (output).
         Final_deriv_Of_Y   => Final_Y_dot,
         Final_Time         => Final_t,       -- integrate out to here (input).
         Initial_Y          => Previous_Y,    -- input an initial condition (input).
         Initial_deriv_Of_Y => Previous_Y_dot,
         Initial_Time       => Previous_t,    -- start integrating here (input).
         No_Of_Steps        => Steps);        -- use this no_of_steps

      Previous_Y     := Final_Y;
      Previous_Y_dot := Final_Y_dot;
      Previous_t     := Final_t;

      put ("Time = t =");
      put (Final_t, Aft => 7);
      new_line;
      put ("(True Energy - Integrated Energy) = ");
      put (Abs (Initial_Energy - Energy (Final_Y, Final_Y_dot)), Aft => 7);
    --put (", (x - Integrated x) = "); -- need to integrate integer no of orbits.
    --put (Abs (0.0 - Final_Y(0)), Aft => 7);
      new_line;

   end loop;

end;

