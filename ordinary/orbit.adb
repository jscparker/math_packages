
with Ada.Numerics.Generic_Elementary_Functions;

package body Orbit is

  package mth is new Ada.Numerics.Generic_Elementary_Functions(Real);
  use mth;

  function Norm (Y : Dynamical_Variable) 
     return Real
  is
     Sum : Real := +0.0;
  begin
     for i in Dyn_Index loop
        Sum := Sum + Abs Y(i);
      --Sum := Sum + Y(i) * Y(i);
     end loop;
     return Sum;
  end Norm;

  function "-" 
    (Left  : Dynamical_Variable;
     Right : Dynamical_Variable)
     return Dynamical_Variable 
  is
     Result : Dynamical_Variable;
  begin
     for I in Dyn_Index loop
        Result(I) := Left(I) - Right(I);
     end loop;
     return Result;
  end "-";

  function "*" 
    (Left  : Real;
     Right : Dynamical_Variable)
     return Dynamical_Variable 
  is
     Result : Dynamical_Variable;
  begin
     for I in Dyn_Index loop
        Result(I) := Left * Right(I);
     end loop;
     return Result;
  end "*";

  function "+" 
    (Left  : Dynamical_Variable;
     Right : Dynamical_Variable)
     return Dynamical_Variable 
  is
     Result : Dynamical_Variable;
  begin
     for I in Dyn_Index loop
        Result(I) := Left(I) + Right(I);
     end loop;
     return Result;
  end "+";

  --  The differential equation is dY/dt = F (t, Y)

  function F
    (Time : Real;
     Y : Dynamical_Variable)
     return Dynamical_Variable
  is
     Deriv : Dynamical_Variable;
     R2, R3, R6, x, z  : Real;
     One  : constant Real := +1.0;
     Four : constant Real := +4.0;
     Min_Allowed_Real : constant Real := 2.0**(Real'Machine_Emin / 2);
  begin
     x := Y(Dyn_Index'First + 0);
     z := Y(Dyn_Index'First + 1);
     R2 := x*x + z*z;
     R6 := R2*R2*R2;
     R3 := Sqrt (One / (R6 + Min_Allowed_Real));
 
     Deriv(Dyn_Index'First + 0) := Y(Dyn_Index'First + 2);
     Deriv(Dyn_Index'First + 1) := Y(Dyn_Index'First + 3);
     Deriv(Dyn_Index'First + 2) := -Four * x * R3;
     Deriv(Dyn_Index'First + 3) := -Four * z * R3;
 
     return Deriv;
  end F;
 
  function Energy 
    (Y : in Dynamical_Variable) 
     return Real 
  is
     x  : constant Real := Y(Dyn_Index'First + 0);
     z  : constant Real := Y(Dyn_Index'First + 1);
     u  : constant Real := Y(Dyn_Index'First + 2);
     v  : constant Real := Y(Dyn_Index'First + 3);
     R : Real;
     Result : Real;
  begin
     R      := Sqrt (X*X + Z*Z);
     Result := 0.25 * (U*U + V*V) - 2.0 / R;
     return Result;
  end Energy;

end Orbit;
