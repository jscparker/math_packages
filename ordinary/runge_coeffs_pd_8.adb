
-----------------------------------------------------------------------
-- package body Runge_Coeffs_PD_8, coefficients for Prince-Dormand Runge Kutta
-- Copyright (C) 2008-2018 Jonathan S. Parker.
--
-- Permission to use, copy, modify, and/or distribute this software for any
-- purpose with or without fee is hereby granted, provided that the above
-- copyright notice and this permission notice appear in all copies.
-- THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
-- WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
-- MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
-- ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
-- WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
-- ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
-- OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
-------------------------------------------------------------------------------
--
-----------------------------------------------------------------------
with Text_IO; use Text_IO;

package body Runge_Coeffs_PD_8 is

   Print_to_Screen_Desired : constant Boolean := False;

   -----------------------
   -- Test_Runge_Coeffs --
   -----------------------

   procedure Test_Runge_Coeffs
   is
      Factors : array(Stages, Stages) of Real := (others => (others => 0.0));
      Factorial, Sum, Error, Max_Error : Real := +0.0;
   begin
      for i in RK_Range loop
      for j in RK_Range loop
         Error := Abs Real((A_rational(i)(j)) - (A_extended(i)(j)));
         if Error > Max_Error then Max_Error := Error; end if;
      end loop;
      end loop;

      if Print_to_Screen_Desired then
         new_line(1);
         put ("Error in A coeffients: "); put (Real'Image (Max_Error));
      end if;

      if Max_Error > 1.0e-14 then
         put ("Problem with the Runge Kutta Coeffs in Runge_Coeffs_PD_8 (0).");
         raise Program_Error;
      end if;
      

      for I in Stages loop
         Factors(I,Stages'First) := 1.0;
      end loop;
      --    K(I = 1) = h*F*Y                         J = 1..0
      --
      --    K(I = 2) = h*F*(Y + K(1)*A21)               J = 1..1
      --
      --    K(I = 3) = h*F*(Y + K(1)*A31 + K(2)*A32)          J = 1..2
      --
      --    K(I = 4) = h*F*(Y + K(1)*A41 + K(2)*A42 + K(3)*A43)     J = 1..3

      for I in stages range 1..Stages'Last loop

         for J in Stages range Stages'First..I-1 loop
         for Order in Stages range Stages'First+1..J+1 loop
            Factors(I,Order) := Factors(I,Order) + A(I)(J)*Factors(J,Order-1);
         end loop;
         end loop;

      end loop;

      -- Now should have SUM(B(i) * Factors(i,n)) = 1/n! up to the order of the
      -- the Taylors series.

      for Order in Stages range 1..8 loop
         Sum := 0.0;
         for I in Stages loop
            Sum := Sum + B8(I) * Factors (I, Order-1);
         end loop;

         --  Calculate 1.0 / Order! to get
         Factorial := 1.0;
         for N in Stages range 2..Order loop
            Factorial := Factorial * Real(N);
         end loop;
         Factorial := 1.0 / Factorial;

         if Print_to_Screen_Desired then
            new_line(2);
            put ("           "); put (Real'Image (Sum));
            new_line;
            put ("Should be: "); put (Real'Image (Factorial));
            new_line;
         end if;

         Error := Abs (Sum - Factorial);
         if Error > 1.0e-14 then
            put ("Problem with the Runge Kutta Coeffs in Runge_Coeffs_PD_8 (1).");
            raise Program_Error;
         end if;
      
      end loop;

      for Order in Stages range 1..7 loop
         Sum := 0.0;
         for I in Stages loop
            Sum := Sum + B7(I) * Factors (I, Order-1);
         end loop;

         --  Calculate 1.0 / Order! to get
         Factorial := 1.0;
         for N in Stages range 2..Order loop
            Factorial := Factorial * Real(N);
         end loop;
         Factorial := 1.0 / Factorial;

         if Print_to_Screen_Desired then
            new_line(2);
            put ("           "); put (Real'Image (Sum));
            new_line;
            put ("Should be: "); put (Real'Image (Factorial));
            new_line;
         end if;

         Error := Abs (Sum - Factorial);
         if Error > 1.0e-14 then
            put ("Problem with the Runge Kutta Coeffs in Runge_Coeffs_PD_8 (2).");
            raise Program_Error;
         end if;
      
      end loop;

   end Test_Runge_Coeffs;

end Runge_Coeffs_PD_8;

