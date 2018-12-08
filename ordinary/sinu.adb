
package body Sinu is

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
     Y    : Dynamical_Variable)
     return Dynamical_Variable 
  is
     Deriv            : Dynamical_Variable;
     SpringConstant   : Constant Real := 1.0;
   --MaxForce         : Constant Real := 5.0;
   --WellWidthInverse : Constant Real := 10.0;
   --WellCenter       : Constant Real := 0.0;
  begin
     Deriv(0) := Y(1);
     Deriv(1) := -SpringConstant*Y(0);
     return Deriv;
  end F;

end Sinu;
