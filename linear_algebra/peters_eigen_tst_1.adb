
-- Test Eigendecomposition of real valued square matrices.
--
-- Good optimization switches:
--
--   gnatmake -gnatnp -O3 -march=native -funroll-loops peters_eigen_tst_1
--
-- For initial test:
--
--   gnatmake -Wall -gnat12 -gnatwa -gnatVa -gnata -gnato -fstack-check -gnateE 
--
-- Large matrices may run out of stack space.  On linux Bash
-- shell the command for more stack is:  ulimit -s unlimited.

with Ada.Numerics.Generic_elementary_functions;
with Peters_Eigen;
with Test_Matrices;
with Text_IO; use Text_IO;

procedure peters_eigen_tst_1 is

   type Real is digits 15;
 --type Real is digits 18; -- useful for estimating err in the 'digits 15' results.

 --Size_Of_Test_Matrix : constant := 57;
   Size_Of_Test_Matrix : constant := 137;
 --Size_Of_Test_Matrix : constant := 7010;

   Index_First : constant := 1; -- matrix corner is (Index_First, Index_First)

   -- Sometimes it's faster if you use a storage matrix that's a little
   -- larger than the original matrix. (It depends on the hardware, the problem,
   -- and the size of the matrix.) A rule of thumb: Width of storage matrix
   -- should be an even number, greater than matrix size, and never 2**n.
   --
   -- Most important thing is to avoid 2**n matrix width.
   --
   -- Padding usually doesn't help much, but just in case let's add 1 or 2:

   Padding : constant := 2 - (Size_Of_Test_Matrix - Index_First + 1) mod 2;

   Matrix_Storage_Size : constant := Size_Of_Test_Matrix + Padding - Index_First + 1;

   pragma Assert (Matrix_Storage_Size < 2**24 and Matrix_Storage_Size > 1);

   -- You can perform the calculation on any square diagonal sub-block of
   -- the matrix defined above by Index x Index.  Since it's square, the
   -- corners of this diagonal sub-block are defined by the 2 numbers
   -- defined below: Starting_Col and Final_Col.

   subtype Index is Integer range Index_First .. Matrix_Storage_Size;

   Starting_Col : constant Index := Index'First;
   Final_Col    : constant Index := Starting_Col + Size_Of_Test_Matrix - 1;

   type Matrix is array(Index, Index) of Real;
   pragma Convention (Fortran, Matrix);

   package Math is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use Math;

   -- Eigen exports Col_Vector:

   package Eigen is new Peters_Eigen
     (Real   => Real, 
      Index  => Index, 
      Matrix => Matrix);
   use Eigen;

   -- Create a package of test matrices:

   package Square_Matrices is new Test_Matrices (Real, Index, Matrix);
   use Square_Matrices;

   package rio is new Float_IO(Real);
   use rio;

   package iio is new Integer_IO(Integer);
   use iio;

 --subtype Real_e is Real;   -- general case, works fine
   type Real_e is digits 18; -- ok on intel, a bit more accuracy in err estimation.

   type Col_Vector_e is array (Index) of Real_e;

   Min_Exp          : constant Integer := Real'Machine_Emin/2 + Real'Machine_Emin/4;
   Min_Allowed_Real : constant Real    := 2.0**Min_Exp;  

   ---------
   -- "*" --
   ---------
  
   function Product
     (A            : in Matrix;
      Vector       : in Col_Vector_e;
      Starting_Col : in Index;
      Final_Col    : in Index)
      return Col_Vector_e is
      Result : Col_Vector_e := (others => 0.0);
      Sum : Real_e;
   begin

      for Row in Starting_Col .. Final_Col loop
         Sum := +0.0;
         for Col in Starting_Col .. Final_Col loop
            Sum := Sum + Real_e (A(Row, Col)) * Vector(Col);
         end loop;
         Result (Row) := Sum;
      end loop;

      return Result;

   end Product;

   ----------------
   -- Hypotenuse --
   ----------------

   -- Hypotenuse = Sqrt (a*a + b*b)

   function Hypotenuse (a, b : Real) return Real is
      Result, r : Real := 0.0;
   begin
      if Abs (a) = 0.0 and then Abs (b) = 0.0 then
         Result := 0.0;
      elsif Abs (b) < Abs (a) then
         r      := Abs (b) / Abs (a);
         Result := Abs (a) * Sqrt (1.0 + r*r);
      else
         r      := Abs (a) / Abs (b);
         Result := Abs (b) * Sqrt (1.0 + r*r);
      end if;
      return Result;
   end Hypotenuse;

   ----------------------------
   -- Error_in_Decomposition --
   ----------------------------

   -- rough estimate, relative err.

   procedure Error_in_Decomposition
     (A            : in     Matrix;
      W_r, W_i     : in     Col_Vector;
      Z_r, Z_i     : in     Matrix;
      Starting_Col : in     Index;
      Final_Col    : in     Index;
      Max_Err      :    out Real)
   is
      Err, Max_Eig, tst : Real;
      Err_2 : Real_e := 0.0;
      A_times_Eig_r, A_times_Eig_i : Col_Vector_e;
      Eigvec_r, Eigvec_i : Col_Vector_e;
      Err_Vec_r, Err_Vec_i : Col_Vector;
      --ZZ_r, ZZ_i : Col_Vector;
      W_times_EigVec_r, W_times_EigVec_i : Col_Vector_e;
   begin

      Max_Err  := 0.0;

      -- Make sure all uncalculated Eigs are set to 0 in W?

      Max_Eig  := 0.0;
      for Col_id in Starting_Col .. Final_Col loop
         tst := Hypotenuse (W_r(Col_id), W_i(Col_id));
         if tst > Max_Eig then  Max_Eig := tst;  end if;
      end loop;

      for Col_id in Starting_Col .. Final_Col loop

         for j in Starting_Col .. Final_Col loop
            Eigvec_r(j) := Real_e (Z_r(j, Col_id));
            Eigvec_i(j) := Real_e (Z_i(j, Col_id));
            --ZZ_r(j) := (Z_r(j, Col_id));
            --ZZ_i(j) := (Z_i(j, Col_id));
         end loop;

         A_times_Eig_r := Product (A, Eigvec_r, Starting_Col, Final_Col);
         A_times_Eig_i := Product (A, Eigvec_i, Starting_Col, Final_Col);

         Err_2 := 0.0;
         for j in Starting_Col .. Final_Col loop
            W_times_EigVec_r(j) := 
               Real_e(W_r(Col_id)) * Eigvec_r(j) - Real_e(W_i(Col_id)) * Eigvec_i(j);
            W_times_EigVec_i(j) := 
               Real_e(W_i(Col_id)) * Eigvec_r(j) + Real_e(W_r(Col_id)) * Eigvec_i(j);

            Err_2 := Err_2 + (A_times_Eig_r(j) - W_times_EigVec_r(j))**2
                           + (A_times_Eig_i(j) - W_times_EigVec_i(j))**2;

            Err_Vec_r(j) := Real (A_times_Eig_r(j) - W_times_EigVec_r(j));
            Err_Vec_i(j) := Real (A_times_Eig_i(j) - W_times_EigVec_i(j));
         end loop;
      
         -- The vectors Eigvec are already normalized, even if A was balanced.
 
         Err := Norm (Err_Vec_r, Err_Vec_i, Starting_Col, Final_Col) / 
            (Max_Eig + Min_Allowed_Real);

         -- ||A*v - lambda*v|| / ||lambda_max||  over all normalized v

         --Err := Sqrt (Real(Err_2)) / (Max_Eig + Min_Allowed_Real);
 
         if Err > Max_Err then  Max_Err := Err;  end if;

       end loop;

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

   failed_eig : Integer;
   Max_Error : Real;

   A_true, A : Matrix;
   Z_r, Z_i  : Matrix;
   Eigvals_r : Col_Vector;
   Eigvals_i : Col_Vector;
   N : constant Integer := Integer(Final_Col) -  Integer(Starting_Col) + 1;

   No_Of_Repetitions : constant Integer := 1;
   Err_Sum : Real := 0.0;
begin

   Pause( 
     "Test1: Eigen-decomposition of matrix A.",
     "Error is estimated by measuring Max value of",
     " ",
     "         |A*v - lambda*v| / |max_eigenvalue_of_A|",
     " ",
     "over all normalized eigenvectors v. If 15 digit Reals are used and if",
     "all goes well, then expect the error to be a few parts per 10**15.",
     " ",
     "Reminder: not all matrices can be diagonalized. Some of the matrices in",
     "the test suite will demonstrate this."
     );

   new_line(2);
   put ("Testing Peters_Eigen on the following"); 
   put (Integer'Image (N)); put (" x"); put (Integer'Image (N)); 
   put (" matrices:"); 
   new_line(1);

   for repetitions in 1 .. No_Of_Repetitions loop

   for Chosen_Matrix in Matrix_Id loop
 --for Chosen_Matrix in Hilbert .. Gear_1 loop
 --for Chosen_Matrix in Hilbert .. Hilbert loop
 --for Chosen_Matrix in Lesp .. Lesp loop
 --for Chosen_Matrix in Chow1 .. Chow1 loop
 --for Chosen_Matrix in COMPANION_1 .. COMPANION_1 loop
 --for Chosen_Matrix in Vandermonde .. Vandermonde loop
 --for Chosen_Matrix in Random_32_Bit .. Random_32_Bit loop
 --for Chosen_Matrix in Random_1_Bit .. Random_1_Bit loop
 --for Chosen_Matrix in Pascal .. Pascal loop
 --for Chosen_Matrix in Pas_Fib .. Pascal loop
 --for Chosen_Matrix in kahan .. kahan loop
  
      Square_Matrices.Init_Matrix (A, Chosen_Matrix, Starting_Col, Final_Col);

      --  For symmetric matrices:
      --  
      --      for c in Starting_Col .. Final_Col loop
      --      for r in Starting_Col .. Final_Col loop
      --         A_true(r, c) := 0.5 * (A(c, r) + A(r, c));
      --      end loop;
      --      end loop;
      --
      --      A := A_true;

      A_true := A; -- Save original A

      Eigen.Decompose 
        (A                    => A,   -- A is destroyed
         Z_r                  => Z_r, -- Columns of output Z_r = Re part of Eigvecs
         Z_i                  => Z_i, -- Columns of output Z_i = Im part of Eigvecs
         W_r                  => Eigvals_r,
         W_i                  => Eigvals_i,
         Id_of_Failed_Eig     => Failed_eig,
         Starting_Col         => Starting_Col,
         Final_Col            => Final_Col,
         Eigenvectors_Desired => True,
         Balance_Policy       => Disabled);
         --Balance_Policy       => Partial);
         --Balance_Policy       => Full);
  
      new_line(1);
      if Failed_Eig >= Integer(Starting_Col) then
         put ("Failed to converge! Converged Eigs start at index:");
         put (Integer'Image (Failed_Eig+1));
      end if;
  
      if false then
         Sort_Eigs_And_Vecs
           (Z_r           => Z_r,
            Z_i           => Z_i,
            W_r           => Eigvals_r,
            W_i           => Eigvals_i,
            Starting_Col  => Starting_Col,
            Final_Col     => Final_Col);
  
         -- the  **2  operation can overflow, but OK here for a quick test.
         for j in Starting_Col .. Final_Col-1 loop
            if Sqrt (Eigvals_r(j+0)**2 + Eigvals_i(j+0)**2)   < 
               Sqrt (Eigvals_r(j+1)**2 + Eigvals_i(j+1)**2) 
            then
               put_line ("Detected Error in Sort_Eigs_And_Vecs"); 
            end if;
         end loop;
      end if;
  
      Error_in_Decomposition 
        (A_true, 
         Eigvals_r, Eigvals_i, 
         Z_r, Z_i, 
         Starting_Col, Final_Col,
         Max_Error);
  
      put ("Error in decomposition ~"); 
      put (Max_Error);
      put (" for matrix "); 
      put (Matrix_id'Image(Chosen_Matrix));

      Err_Sum := Err_Sum + Max_Error;
  
      if false then
         New_Line;
         for r in Starting_Col .. Final_Col loop
         for c in Starting_Col .. Final_Col loop
            put (A(r, c), Aft => 4);
         end loop;
         new_line(1);
         end loop;
      end if;
  
      if false then
         New_Line;
         for i in Starting_Col .. Final_Col loop
            put (Eigvals_r(i));  put (Eigvals_i(i)); new_line(1);
          --put (Hypotenuse (Eigvals_r(i), Eigvals_i(i)));
         end loop;
      end if;
  
      if true then
         declare 
            Sum : Real := +0.0; 
          --Product, R : Real := +1.0; 
         begin
            for I in Starting_Col .. Final_Col loop
             --R := Hypotenuse (Eigvals_r(I), Eigvals_i(I));
             --Sum := Sum + R;
             --Product := Product * R;
             --Sum := Hypotenuse (Sum, Hypotenuse (Eigvals_r(I), Eigvals_i(I)));
             --Sum := Sum + Abs Eigvals_r(I) + Abs Eigvals_i(I);
               Sum := Sum + Eigvals_r(I);
            end loop;
            New_Line; put ("Eigval sum: "); put (Sum); 
         end;
      end if;
  
   end loop;
   end loop;

   new_line; put("Average = "); put (real'Image(Err_Sum / Real(No_Of_Repetitions)));
  
end peters_eigen_tst_1;
