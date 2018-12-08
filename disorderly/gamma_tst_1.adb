
-- gamma(1+x) = x! for x=0,1,2 ...
-- gamma(0.5) = sqrt(pi)
-- ln (sqrt(pi)) = 0.57236494292470008707171367567652935582354993281988
--     sqrt(pi)  = 1.77245385090551602729816748334114518279737668859774

with Gamma;
with text_io; use text_io;
with ada.numerics.generic_elementary_functions;

procedure gamma_tst_1 is

  type Real is digits 15;

  x : real;

  package G is new Gamma (Real); use G;
  package math is new ada.numerics.generic_elementary_functions(Real); use math;
  --for exp

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

   for j in 0 .. 159 loop
     x := 0.1 + Real(j) * 0.1;
     new_line;
     put(Real'image(x)); 
     put(real'image(log_gamma(x)-log_gamma_0_to_16(x)));
   end loop;

   new_line;
   put(real'image(log_gamma(0.5)-0.57236494292470008707171367567652935582355));
   new_line;
   put(real'image(Exp (log_gamma(0.5))-1.7724538509055160272981674833411451828));
   new_line;
   new_line; put(real'image( (log_gamma(7.0) - Log (Factorial(6)) )));
   new_line; put(real'image( (log_gamma(11.0) - Log (Factorial(10)) )));
   new_line; put(real'image( (log_gamma(15.0) - Log (Factorial(14)) )));
   new_line; put(real'image( (log_gamma(28.0) - Log (Factorial(27)) )));

end;
