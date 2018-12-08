
-- F defines a differential equation whose solution is Exp (i*t).
-- dY/dt = F(Y).
-- For testing.
generic

   type Real is digits <>;

package Sinu_2 is

   type Dyn_Index is range 0..1;  -- the 2nd component is just ignored in the tests.

   type Dynamical_Variable is array(Dyn_Index) of Real;

   DynZero : constant Dynamical_Variable := (others => 0.0);

   function F 
     (Time : Real;
      Y    : Dynamical_Variable)
      return Dynamical_Variable;

   -- Defines the equation to be integrated,
   -- dY/dt = F (t, Y).  Even if the equation is t or Y
   -- independent, it must be entered in this form.

   function "*" 
     (Left  : Real;
      Right : Dynamical_Variable) 
      return Dynamical_Variable;

   function "+" 
     (Left  : Dynamical_Variable;
      Right : Dynamical_Variable) 
      return Dynamical_Variable;

   function "-" 
     (Left  : Dynamical_Variable;
      Right : Dynamical_Variable) 
      return Dynamical_Variable;

   function Norm (Y : Dynamical_Variable) return Real;

   pragma Inline (F, "*", "+", "-", Norm);

end Sinu_2;
