
with text_io;
procedure rank_distr is

  type Real is digits 15;

  L : constant Integer := 61; -- length of vector in bits
  k : constant Integer := L;  -- no of vectors in matrix whose rank is to be calculated.
  Sum, Product : Real;
  Prob : array(1..k) of Real := (others => 0.0);

begin

 for x in k-7 .. k loop

  product := 1.0;
  for i in 0 .. x-1 loop
    product := product*(1.0 - 2.0**(i-L))*(1.0 - 2.0**(i-k)) / (1.0 - 2.0**(i-x));
  end loop;

  prob(x) := 2.0**(x*(L + k - x) - L*k) * product;
  
 end loop;

 for x in k-7 .. k loop
   Sum := 0.0;
   for j in x .. k loop
     Sum := Sum + Prob(j);
   end loop;
 end loop;

 Sum := 0.0;
 for j in k-7 .. k-4 loop
   Sum := Sum + Prob(j);
 end loop;

 text_io.new_line;
 text_io.put (Integer'Image(k-4)); text_io.put (Real'Image(Sum));
 text_io.new_line;

for x in k-3 .. k loop

  text_io.new_line;
  text_io.put (Integer'Image(x)); text_io.put (Real'Image (Prob(x)));
  --text_io.put (Real'Image(-Log (Prob (x))));

 end loop;

 text_io.new_line(2);
 Sum := 0.0;
 for j in k-7 .. k loop
   Sum := Sum + Prob(j);
 end loop;
 text_io.put ("Total probability"); text_io.put (Real'Image(Sum));

end;
