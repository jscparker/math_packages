
-- package Quadrature
--
-- The equation is dY/dt = F (t).
-- The solution is just the area under the curve - quadrature.
-- F(t) is input as a generic formal function called Integrand.
--
-- If you want (d/dt)**N Y = F(t), then N is written as
--
--    N = Order_of_Equation
--
-- Then the integrated Y(0) is solution of (d/dt)**N   Y = F(t),
-- Then the integrated Y(1) is solution of (d/dt)**N-1 Y = F(t),
-- and so on.

generic

   type Real is digits <>;

   with function Integrand (X : Real) return Real is <>;

package Quadrature is

   Order_of_Equation : constant Positive := 1;

   type Dyn_Index is range 0 .. Order_of_Equation-1;

   type Dynamical_Variable is array(Dyn_Index) of Real;

   DynZero : constant Dynamical_Variable := (others => (+0.0));

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

end Quadrature;
