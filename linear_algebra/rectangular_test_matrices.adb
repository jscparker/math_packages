
---------------------------------------------------------------------------
-- package body rectangular_test_matrices, test matrix generator
-- Copyright (C) 2018 Jonathan S. Parker.
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
---------------------------------------------------------------------------

with Ada.Numerics.Discrete_Random;

package body Rectangular_Test_Matrices is

 type Unsigned32 is mod 2**32;

 package Discrete_32_bit is new Ada.Numerics.Discrete_Random (Unsigned32);

 Random_Stream_id : Discrete_32_bit.Generator;

 function Real_Random return Real is
   X : Real := 0.0;
 begin
  
   -- up to 64 random bits in X:

   --X := 2.0**(-32) * Real (Discrete_32_bit.Random (Random_Stream_id)) 
      --+ 2.0**(-64) * Real (Discrete_32_bit.Random (Random_Stream_id));

   -- 16 random bits in X:

   --X := 2.0**(32-16) * Real (Discrete_32_bit.Random (Random_Stream_id) / 2**(32-16));
     X := 2.0**(-16) * Real (Discrete_32_bit.Random (Random_Stream_id) / 2**(32-16));

   -- 1 random bit in X:

   --X := 2.0**(-1) * Real (Discrete_32_bit.Random (Random_Stream_id) / 2**(32-1));

   return X;

 end Real_Random;


   
 procedure Init_Matrix 
   (M               : out Matrix;
    Desired_Matrix  : in  Matrix_Id := Random;
    Final_Row       : in  R_Index   := R_Index'Last;
    Final_Col       : in  C_Index   := C_Index'Last;
    Starting_Row    : in  R_Index   := R_Index'First;
    Starting_Col    : in  C_Index   := C_Index'First) is
 
    Denominator : Real;
   
 begin

   M := (others => (others => 0.0));     --essential init.
   --the 0.0 is essential.

   case Desired_Matrix is

   when Zero_Diagonal  =>

      declare
         Row : R_Index := Starting_Row;
      begin
         for Col in C_Index range Starting_Col+1 .. C_Index'Last loop
            M(Row, Col) := 1.0;
	    if Row < R_Index'Last then  Row := Row + 1; else exit; end if;
         end loop;

         Row := Starting_Row+1;
         for Col in C_Index range Starting_Col .. C_Index'Last-1 loop
            M(Row, Col) := 1.0;
	    if Row < R_Index'Last then  Row := Row + 1; else exit; end if;
         end loop;
      end;

   when Upper_Ones  =>

      M := (others => (others => 1.0));     --essential init.

      declare
         Row_Starting_Index : R_Index := Starting_Row;
      begin
      for Col in C_Index range Starting_Col .. Final_Col loop
         for Row in Row_Starting_Index .. R_Index'Last loop
            M(Row, Col) := 0.0; 
         end loop;
	 if Row_Starting_Index < R_Index'Last then   
	    Row_Starting_Index := Row_Starting_Index + 1;
	 else
	    exit;
	 end if;
      end loop;
      end;

      declare
         Row : R_Index := Starting_Row;
      begin
         for Col in C_Index range Starting_Col .. Final_Col loop
            M(Row, Col) := 1.0;
	    if Row < R_Index'Last then 
	       Row := Row + 1;
	    else
	       exit;
	    end if;
         end loop;
      end;

   when Lower_Ones  =>

      M := (others => (others => 0.0));     --essential init.

      declare
         Row_Starting_Index : R_Index := Starting_Row;
      begin
      for Col in C_Index range Starting_Col .. C_Index'Last loop
         for Row in Row_Starting_Index .. R_Index'Last loop
            M(Row, Col) := 1.0; 
         end loop;
	 if Row_Starting_Index < R_Index'Last then   
	    Row_Starting_Index := Row_Starting_Index + 1;
	 else
	    exit;
	 end if;
      end loop;
      end;

   when Random  =>

      Discrete_32_bit.Reset (Random_Stream_id);

      for Row in R_Index loop
      for Col in C_Index loop
         M(Row, Col) := Real_Random; 
      end loop;
      end loop;
   
   when Zero_Cols_and_Rows  =>

      Discrete_32_bit.Reset (Random_Stream_id);

      for Row in R_Index loop
      for Col in C_Index loop
         M(Row, Col) := Real_Random; 
      end loop;
      end loop;
   
      for Row in Starting_Row .. Starting_Row+2 loop
      for Col in Starting_Col .. C_Index'Last loop
         M(Row, Col) := 0.0; 
      end loop;
      end loop;
   
      for Col in Starting_Col .. Starting_Col+2 loop
      for Row in Starting_Row .. R_Index'Last loop
         M(Row, Col) := 0.0; 
      end loop;
      end loop;
   
   when Frank  =>

      declare
        c, r : Real;
      begin
         for Row in R_Index loop
         for Col in C_Index loop
            c := Real(Col) - Real(C_Index'First);
            r := Real(Row) - Real(R_Index'First);
            M(Row, Col) := Real'Min (c,r) + 1.0; 
         end loop;
         end loop;
      end;

   when Ding_Dong  =>
  
      for Row in Starting_Row .. R_Index'Last loop
      for Col in Starting_Col .. C_Index'Last loop
          Denominator :=  Real(R_Index'Last) 
	             + Real(Starting_Row) - Real(Row) - Real(Col) + 0.5;
	  M(Row, Col) :=  1.0 / Denominator;
      end loop;
      end loop;
       
   when Vandermonde  =>

   --A := 1.0 / 16.0; -- up to 256 x 256 matrix without over/under flow (if digits 15).
   --A := 1.0; -- up to 128 x 128 matrix ? without over/under flow (if digits 15).
 
     Vandermondes_Matrix:
     declare
       Exp, X_Count : Integer;
       X : Real;
       Half_No_Of_Rows : constant Integer 
                   := (Integer(Final_Row) - Integer(Starting_Row) + 1)/2;
       B : constant Real := 1.0 / Real (Half_No_Of_Rows);
       A : constant Real := 2.0 ** (Real'Exponent(B) - 1);
     begin
       for Row in Starting_Row .. R_Index'Last loop
       for Col in Starting_Col .. C_Index'Last loop
         Exp         := Integer(Col) - Integer(Starting_Col);
	 X_Count     := Integer(Row) - Integer(Starting_Row); 
         X           := A * (Real (X_Count - Half_No_Of_Rows));
         M(Row, Col) := X ** Exp;
       end loop;
       end loop;
     end Vandermondes_Matrix;

   when All_Ones  =>
 
      M := (others => (others => 1.0));
 
   when All_Zeros  =>
 
      M := (others => (others => 0.0));
 
   end case;

 end Init_Matrix;

end Rectangular_Test_Matrices;
