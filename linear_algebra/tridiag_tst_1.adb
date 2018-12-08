
-- Test Eigendecomposition of real valued square matrices.

with Ada.Numerics.Generic_elementary_functions;
with Tridiagonal;
with Test_Matrices;
with Text_IO; use Text_IO;

procedure tridiag_tst_1 is

   type Real is digits 15;
 --type Real is digits 18;

 --Desired_Size_Of_Matrix : constant := 2291;
 --Desired_Size_Of_Matrix : constant := 1024;
 --Desired_Size_Of_Matrix : constant := 478;
 --Desired_Size_Of_Matrix : constant := 254;
   Desired_Size_Of_Matrix : constant := 136;
 --Desired_Size_Of_Matrix : constant := 87;
 --Desired_Size_Of_Matrix : constant := 95;

   pragma Assert (Desired_Size_Of_Matrix < 2**24 and Desired_Size_Of_Matrix > 1);

   -- Next: make sure Matrix_Storage_Size is never odd, and never a power of 2.
   --
   -- This doesn't help for small matrices (Size < 128 on some machines).
   --
   -- But odd sizes and power of 2 sizes are a performance disaster on certain
   -- processors (when the matrix is too large for certain cache sizes).
   -- Lapack, compiled by gfortran, seems to have similar problems.

   type Unsigned_32 is mod 2**32;
   Padding_0      : constant Unsigned_32 := 2 + Desired_Size_Of_Matrix mod 2;
   Storage_0      : constant Unsigned_32 := Padding_0 + Desired_Size_Of_Matrix;
   Power_of_2_Tst : constant Unsigned_32 := Storage_0 and (Storage_0 - 1);
   More_Padding   : constant Unsigned_32 := 2 * (1 / (1 + Power_of_2_Tst));

   -- Power_of_2_Tst = 0 iff Storage = 2**n.   More_Padding = 2 iff Storage = 2**n

   Matrix_Storage_Size : constant Integer := Integer (Storage_0 + More_Padding);

   -- the test matrix is square-shaped matrix on:  Index x Index.
   -- eg Hilbert's matrix is a square matrix with unique elements on the range
   -- Index'First .. Index'Last.  However, you  have the option or using any 
   -- diagonal sub-block of the matrix defined by Index x Index

   subtype Index is Integer range 0 .. Matrix_Storage_Size-1;

   Starting_Col : constant Index := Index'First;
   Final_Col    : constant Index := Starting_Col + Desired_Size_Of_Matrix - 1;

   -- You have the option or using any diagonal sub-block
   -- of the matrix defined by Index x Index.  Since it's square, the
   -- corners of this diagonal sub-block are defined by the 2 numbers
   -- defined above, Starting_Col and Final_Col.

   type Matrix is array(Index, Index) of Real;

   pragma Convention (Fortran, Matrix);
   -- use Convention Fortran.  Over twice as fast as Ada convention.

   package Tridiag is new Tridiagonal
     (Real     => Real, 
      Index    => Index, 
      A_Matrix => Matrix);

   -- Create a package of test matrices:

   package Square_Matrices is new Test_Matrices (Real, Index, Matrix);
   use Square_Matrices;

   package rio is new Float_IO(Real);
   use rio;

 --subtype Real_e is Real;   -- general case, works fine
   type Real_e is digits 18;   -- 18 ok on intel

   package Math_e is new Ada.Numerics.Generic_Elementary_Functions (Real_e);
   use Math_e;

   type Matrix_e is array(Index, Index) of Real_e;

   Min_Exp          : constant Integer := Real'Machine_Emin - Real'Machine_Emin/32;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (A            : in     Matrix_e;
      Starting_Col : in     Index;
      Final_Col    : in     Index)
      return Real_e
    is
      Min_Allowed_Real : constant Real_e := 2.0**Min_Exp;  
      Max_A_Val : Real_e := +0.0;
      Sum, Scaling, tmp : Real_e := +0.0;
    begin
 
      Max_A_Val := +0.0;
      for Row in Starting_Col .. Final_Col  loop
      for Col in Starting_Col .. Final_Col  loop
         if Max_A_Val < Abs A(Row, Col) then Max_A_Val := Abs A(Row, Col); end if;
      end loop;
      end loop;

      Max_A_Val := Max_A_Val + Min_Allowed_Real;
      Scaling   := +1.0 / Max_A_Val;

      Sum := +0.0;
      for Row in Starting_Col .. Final_Col  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * A(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;

 
   -- Q * T * Q'
 
   function Product
     (Q, T         : in Matrix;
      Starting_Col : in Index;
      Final_Col    : in Index)
      return Matrix_e 
   is
      E, Result : Matrix_e := (others => (others => 0.0));
      Sum : Real_e;
   begin

      for Row in Starting_Col .. Final_Col loop
      for Col in Starting_Col .. Final_Col loop
         Sum := +0.0;
         for k in Starting_Col .. Final_Col loop
            Sum := Sum + Real_e (Q(Row, k)) * Real_e (T(k, Col));
         end loop;
         E (Row, Col) := Sum;
      end loop;
      end loop;

      for Row in Starting_Col .. Final_Col loop
      for Col in Starting_Col .. Final_Col loop
         Sum := +0.0;
         for k in Starting_Col .. Final_Col loop
            Sum := Sum + E(Row, k) * Real_e (Q(Col, k));  -- Q_transpose
         end loop;
         Result (Row, Col) := Sum;

         if not Sum'valid then  put("!!!"); end if;

      end loop;
      end loop;
      return Result;

   end Product;

   ----------------------------
   -- Error_in_Decomposition --
   ----------------------------

   -- rough estimate, relative err.

   procedure Error_in_Decomposition
     (A            : in     Matrix;  -- A_true
      U, T         : in     Matrix;
      Starting_Col : in     Index;
      Final_Col    : in     Index;
      Max_Error    :    out Real)
   is
      Err, Arc : Real_e := 0.0;
      D : Matrix_e;
      A_e : Matrix_e;
      Min_Allowed_Real : constant Real_e := 2.0**Min_Exp;  
   begin

      -- Want Err = A - U * T * U'

      Max_Error := 0.0;

      -- D = U * T * U'

      D := Product (U, T, Starting_Col, Final_Col);

      -- resuse D as D = A - U * T * U'

      for r in Starting_Col .. Final_Col loop
      for c in Starting_Col .. Final_Col loop
         Arc       := Real_e(A(r, c));
         A_e(r, c) := Arc;
         D(r, c)   := D(r, c) - Arc;
      end loop;
      end loop;
       
      Err :=  Frobenius_Norm (D, Starting_Col, Final_Col) / 
         (Frobenius_Norm (A_e, Starting_Col, Final_Col) + Min_Allowed_Real);

      Max_Error := Real (Err);

   end Error_in_Decomposition;

   -----------
   -- Pause --
   -----------

   procedure Pause (s0,s1,s2,s3,s4,s5,s6,s7,s8,s9 : string := "") is
     Continue : Character := ' ';
   begin
     new_line;
     if S0 /= "" then put_line (S0); end if;
     if S1 /= "" then put_line (S1); end if;
     if S2 /= "" then put_line (S2); end if;
     if S3 /= "" then put_line (S3); end if;
     if S4 /= "" then put_line (S4); end if;
     if S5 /= "" then put_line (S5); end if;
     if S6 /= "" then put_line (S6); end if;
     if S7 /= "" then put_line (S7); end if;
     if S8 /= "" then put_line (S8); end if;
     if S9 /= "" then put_line (S9); end if;
     new_line;
     begin
        put ("Type a character to continue: ");
        get_immediate (Continue);
     exception
        when others => null;
     end;
   end pause;

   x, Max_Error : Real := 0.0;

   A_true, A : Matrix := (others => (others => 0.0));
   U : Matrix := (others => (others => 0.0));
   N : constant Integer := Integer(Final_Col) -  Integer(Starting_Col) + 1;
begin

   Pause( 
     "Test1: Tri-diagonalization of matrix A.",
     "Error is estimated by measuring Max value of",
     " ",
     "         ||Q * A_tridiag * Q' - A|| / ||A||",
     " ",
     "where ||M|| denotes the Frobenius norm of matrix M. If 15 digit Reals are",
     "used and if all goes well, then expect error to be a few parts per 10**15."
     );

   new_line(2);
   put ("More padding = "); put (Unsigned_32'Image(more_padding));
   new_line(2);
   put ("Testing Tridiagonal on the following"); 
   put (Integer'Image (N)); put (" x"); put (Integer'Image (N)); 
   put (" symmetrized matrices:"); 
   new_line(1);

   for Chosen_Matrix in Matrix_Id loop
 --for Chosen_Matrix in pas_fib .. pas_fib loop
 --for Chosen_Matrix in vandermonde .. companion_0 loop
 --for Chosen_Matrix in vandermonde .. vandermonde loop
 --for Chosen_Matrix in companion_1 .. companion_1 loop
 --for Chosen_Matrix in clustered .. clustered loop
 --for Chosen_Matrix in fiedler_1 .. fiedler_1 loop
 --for Chosen_Matrix in laguerre .. laguerre loop
  
      Square_Matrices.Init_Matrix (A, Chosen_Matrix, Starting_Col, Final_Col);
 
      transpose (A);

      -- Tridi only works on symmetric matrices:
 
      for c in Starting_Col .. Final_Col loop
      for r in Starting_Col .. Final_Col loop
         A_true(r, c) :=  A(r, c) + A(c, r);
      end loop;
      end loop;
  
      A := A_true; -- A will be destroyed, so save it in A_true.
  
      --  A_tridiag = Q_tr * A_true * Q.
      --
      --  A_true    = Q  * A_tridiag * Q_tr.

      Tridiag.Tridiagonalize 
        (A                    => A,   -- A is replaced with:  A_tridiagonal
         Q                    => U,   -- A_true  =  U * A_tridiagonal * U_transpose
         Starting_Col         => Starting_Col,
         Final_Col            => Final_Col,
         Initial_Q            => Tridiag.Identity,
         Q_Matrix_Desired     => True);
 
      x := x + A(2,2) + U(2,2);

      Error_in_Decomposition 
        (A_true, 
         U, A,                         -- A should be tridiagonal now
         Starting_Col, Final_Col,
         Max_Error);

      for c in Starting_Col .. Final_Col-1 loop
      for r in c .. Final_Col  loop
         if A(r, c) /= A(c, r) then
            new_line; 
            put("Symmetry error for matrix"); 
            put (Matrix_id'Image(Chosen_Matrix));
         end if;
      end loop;
      end loop;
  
      new_line; 
      put ("Error in decomposition ~"); 
      put (Max_Error);
      put (" for matrix "); 
      put (Matrix_id'Image(Chosen_Matrix));
  
   end loop;
 --put (x); 
  
end tridiag_tst_1;
