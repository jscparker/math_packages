
package body Quadrature is

  function Norm (Y : Dynamical_Variable) 
     return Real
  is
     Sum : Real := +0.0;
  begin
     for i in Dyn_Index loop
        Sum := Sum + Abs Y(i);
      --Sum := Sum + Y(i) * Y(i); -- return Sqrt (Sum);
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
  begin
     for I in Dyn_Index range Dyn_Index'First .. Dyn_Index'Last-1 loop
       Deriv(I) := Y(I+1);
     end loop;
     Deriv(Dyn_Index'Last) := Integrand (Time);
     return Deriv;
  end F;
 
end Quadrature;
