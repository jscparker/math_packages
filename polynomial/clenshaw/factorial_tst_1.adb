
-- gamma(1+x) = x! for x=0,1,2 ...

-- Calls 2 important Test procedures from Factorial, then quick test.

with Factorial;
with text_io; use text_io;
with ada.numerics.generic_elementary_functions;

procedure Factorial_tst_1 is

  type Real is digits 15;

  package F is new Factorial (Real); use F;
  package Math is new ada.numerics.generic_elementary_functions(Real); use Math;

  function Factorial (n : in Natural) return Real is
    Prod : Real := 1.0;
   begin
    if n < 2 then
      return 1.0;
    end if;
    for j in 2 .. n loop
      Prod := Prod * Real(j);
    end loop;
    return Prod;
  end;

begin

   Test_Stieltjes_Coefficients;
   Test_Log_Factorial_Table;

   new_line;
   new_line; put(real'image( (log_factorial(1) - Log (Factorial(1)) )));
   new_line; put(real'image( (log_factorial(2) - Log (Factorial(2)) )));
   new_line; put(real'image( (log_factorial(6) - Log (Factorial(6)) )));
   new_line; put(real'image( (log_factorial(10) - Log (Factorial(10)) )));
   new_line; put(real'image( (log_factorial(14) - Log (Factorial(14)) )));
   new_line; put(real'image( (log_factorial(20) - Log (Factorial(20)) )));
   new_line; put(real'image( (log_factorial(21) - Log (Factorial(21)) )));
   new_line; put(real'image( (log_factorial(27) - Log (Factorial(27)) )));
   new_line; put(real'image( (log_factorial(37) - Log (Factorial(37)) )));
   new_line; put(real'image( (log_factorial(47) - Log (Factorial(47)) )));
   new_line; put(real'image( (log_factorial(52) - Log (Factorial(52)) )));
   new_line; put(real'image( (log_factorial(57) - Log (Factorial(57)) )));
   new_line; put(real'image( (log_factorial(87) - Log (Factorial(87)) )));

end;
