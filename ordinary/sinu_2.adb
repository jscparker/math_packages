
package body Sinu_2 is

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

  --  The differential equation is (d/dt)^2 Y = F (t, Y)

  function F 
    (Time : Real;
     Y    : Dynamical_Variable)
     return Dynamical_Variable 
  is
     Second_Deriv : Dynamical_Variable;
  begin
     Second_Deriv(0) := -Y(0);
     Second_Deriv(1) := 0.0; -- just ignore 2nd component.
     return Second_Deriv;
  end F;

end Sinu_2;
