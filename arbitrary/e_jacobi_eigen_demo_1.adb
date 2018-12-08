
-----------------------------------------------------------------------
-- procedure e_jacobi_eigen_tst_1, test of extended precision Jacobi
-- Copyright (C) 2008-2009 Jonathan S. Parker.
--
-- This program is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
-- You should have received a copy of the GNU General Public License
-- along with this program.  If not, see  http://www.gnu.org/licenses/ 
--
-- As a special exception, if other files instantiate generics from
-- this unit, or you link this unit with other files to produce an
-- executable, this  unit  does not  by itself cause  the resulting
-- executable to be covered by the GNU General Public License. This
-- exception does not however invalidate any other reasons why the
-- executable file  might be covered by the  GNU Public License.
-----------------------------------------------------------------------

-- Test Jacobi Eigendecomposition of real valued square matrices.

with Ada.Numerics.Generic_Elementary_Functions;
with Extended_Real;
with Extended_Real.Elementary_Functions;
with Extended_Real.IO;
with E_Jacobi_Eigen;
With Text_IO; use Text_IO;

procedure e_jacobi_eigen_demo_1 is

   type Real_8 is digits 15;

   package mth is new Ada.Numerics.Generic_Elementary_Functions (Real_8);
   use mth;
   package ext is new Extended_Real (Real_8);
   use ext; 
   package fnc is new ext.Elementary_Functions (Sqrt, Log, Exp, Arcsin);
   use fnc;
   package eio is new ext.IO; -- extented IO
   use eio;
   package rio is new Text_IO.Float_IO (Real_8);
   use rio;
   package iio is new Integer_IO (Integer);
   use iio;

   subtype Real is e_Real;

   subtype Index is Integer range 1..36;

   -- the test matrix is square-shaped matrix on:  Index x Index.
   -- eg Hilbert's matrix is a square matrix with unique elements on the range
   -- Index'First .. Index'Last.  However, you  have the option or using any 
   -- diagonal sub-block of the matrix defined by Index x Index

   subtype Row_Index is Index;
   subtype Col_Index is Index;

   Starting_Col : constant Index := Index'First + 0;
   Final_Col    : constant Index := Index'Last  - 0;

   -- Can't change:

   Starting_Row : constant Index := Starting_Col;
   Final_Row    : constant Index := Final_Col;

   type Matrix is array(Index, Index) of Real;

   --pragma Convention (Fortran, Matrix); --No! prefers Ada convention.

   package eig is new e_Jacobi_Eigen (e_Real, Index, Matrix, Real_8, e_Integer);
   use eig;

   -- Eig exports Col_Vector

   subtype Real_Extended is Real;     -- general case, works fine

   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;

   A, A_true, Q_tr : Matrix;
   Eigenvals : Col_Vector;
   Frobenius_QtrQ_Err, Frobenius_QQtr_Err, Frobenius_QEQ_Err : Real;

   No_of_Sweeps_Done, No_of_Rotations : Natural;
   
   N : constant Real_8 := Real_8 (Starting_Col) - Real_8 (Final_Col) + 1.0;

   -----------------
   -- Rayleigh_Eig --
   -----------------
  
   function Rayleigh_Eig 
     (A            : in     Matrix;
      Q_tr         : in     Matrix;
      j            : in     Index)
      return Real 
   is
      V : Col_Vector := (others => Zero);
      Sum : Real_Extended;
   begin

      for Row in Row_Index range Starting_Row .. Final_Row loop
         Sum := Zero;
         for Col in Col_Index range Starting_Col .. Final_Col loop
	   Sum := Sum + Real_Extended (A(Row, Col)) * Real_Extended (Q_tr(j, Col));
         end loop;
	 V (Row) := Real (Sum);
      end loop;

      Sum := Zero;
      for Col in Col_Index range Starting_Col .. Final_Col loop
         Sum := Sum + Real_Extended (V(Col)) * Real_Extended (Q_tr(j, Col));
      end loop;

      return Real (Sum);

   end Rayleigh_Eig;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (A : in Matrix)
      return Real
    is
      Max_A_Val : Real := Zero;
      Sum, Scaling, tmp : Real := Zero;
    begin
 
      Max_A_Val := Zero;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         if Max_A_Val < Abs A(Row, Col) then Max_A_Val := Abs A(Row, Col); end if;
      end loop;
      end loop;

      Max_A_Val := Max_A_Val + Two ** Integer (e_Real_Machine_Emin + 8);
      Scaling := One / Max_A_Val;

      Sum := Zero;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * A(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;

   ------------------------------------
   -- Get_Err_in_Reassembled_Q_and_A --
   ------------------------------------
  
   --  check that A = V*E*Q_tr
   --  E is diagonal with the eigs along the diag.
   --  V is orthogonal with the eig vecs as columns.

   
   procedure Get_Err_in_Reassembled_Q_and_A 
     (A                  : in     Matrix;  -- true original A
      Q_tr               : in     Matrix;
      E                  : in     Col_Vector;
      Final_Col          : in     Col_Index;
      Starting_Col       : in     Col_Index;
      Frobenius_QtrQ_Err :    out Real;
      Frobenius_QQtr_Err :    out Real;
      Frobenius_QEQ_Err  :    out Real)
   is
      Err, S          : Real;
      Min_Real : constant Real := Two ** Integer (e_Real_Machine_Emin + 8);
    
      Sum : Real_Extended;

      Identity    : Matrix := (others => (others => Zero));
      Product_QQ  : Matrix := (others => (others => Zero));
      Product_QEQ : Matrix := (others => (others => Zero));

      subtype Index_Subrange is Index range Starting_Col .. Final_Col;

   begin
  
      for r in Index_Subrange loop
        Identity(r, r) := One;
      end loop;

      -- Find error in I - Q*Q' etc.
      -- Notation: Q' == Q_tr == transpose of Q.

      for Col in Index_Subrange loop
      for Row in Index_Subrange loop
        Sum := Zero;
	for j in Index_Subrange loop
           Sum := Sum + Real_Extended (Q_tr(j, Row)) * Real_Extended (Q_tr(j, Col)); 
        end loop;
	Product_QQ(Row, Col) := Real (Sum);
      end loop;
      end loop;

      -- Get Frobenius norm of:  Product_QQ - I:

      S := Zero;
      for Col in Index_Subrange loop
      for Row in Index_Subrange loop
         Err := Identity(Row, Col) - Product_QQ(Row, Col);
         S := S + Err * Err;
      end loop;
      end loop;

      -- Get fractional Frobenius :  Frobenius (Product_QQ - I) / Frobenius (I):

      Frobenius_QQtr_Err := Sqrt(S) / 
              Sqrt (Make_Extended (Real_8(Starting_Col)+Real_8(Final_Col)+1.0));
  

      -- Find error in I - Q'*Q.
      -- reuse array Product_QQ:

      for Col in Index_Subrange loop
      for Row in Index_Subrange loop
        Sum := Zero;
	for j in Index_Subrange loop
           Sum := Sum + Real_Extended (Q_tr(Row, j)) * Real_Extended (Q_tr(Col, j)); 
        end loop;
	Product_QQ(Row, Col) := Real (Sum);
      end loop;
      end loop;

      -- Get Frobenius norm of:  Product_QQ - I:

      S := Zero;
      for Col in Index_Subrange loop
      for Row in Index_Subrange loop
         Err := Identity(Row, Col) - Product_QQ(Row, Col);
         S := S + Err * Err;
      end loop;
      end loop;

      -- Get fractional Frobenius :  Frobenius (Product_QQ - I) / Frobenius (I):

      Frobenius_QtrQ_Err := Sqrt(S) / 
              Sqrt (Make_Extended (Real_8(Starting_Col)+Real_8(Final_Col)+1.0));
   
      --  check that A = Q*E*Q_tr
      --  E is diagonal with the eigs along the diag.
      --  Q is orthogonal with the eig vecs as columns.

      -- explicitly calculate Q*E*Q_tr:

      for Col in Index_Subrange loop
      for Row in Index_Subrange loop
        Sum := Zero;
	for j in Index_Subrange loop
           Sum := Sum +
              Real_Extended (Q_tr(j, Row)) *           --  Q(Row, j)
	      Real_Extended (E(j))         *           --  j-th eig is const along Q col
	      Real_Extended (Q_tr(j, Col)); 
        end loop;
	Product_QEQ(Row, Col) := Real (Sum);
      end loop;
      end loop;

      -- resuse array Product_QEQ to get Error Matrix := Product_QEQ - A:

      for Col in Starting_Col .. Final_Col loop
      for Row in Starting_Row .. Final_Row loop
         Product_QEQ(Row, Col) := A(Row, Col) - Product_QEQ(Row, Col);
      end loop;
      end loop;

      Frobenius_QEQ_Err := Frobenius_Norm (Product_QEQ) / 
                                     (Frobenius_Norm (A) + Min_Real);
   
   end Get_Err_in_Reassembled_Q_and_A;

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
     new_line;
   end pause;

  ------------------------
  -- Get_Hilbert_Matrix --
  ------------------------
   
   procedure Get_Hilbert_Matrix 
     (A : out Matrix)
   is
      Prime_Factors, Denominator : Real;
   begin

    --Prime_Factors := 3.0*3.0*3.0*5.0*5.0*7.0*11.0*13.0*17.0*19.0*23.0*29.0*31.0*37.0;
      -- so  Prime_Factors / D  is exactly represented in 15 digit floating point
      -- up to D = 39 (allowing an 20x20 matrix).  Prime_Factors = 166966608033225.0

    --Prime_Factors := (+166966608033225.0) * (+580027.0) * Two**(-68);
    --Prime_Factors := One;
      Prime_Factors := (+0.25);
  
      for Row in Starting_Col .. Final_Col loop
      for Col in Starting_Col .. Final_Col loop
         Denominator := +(Real_8(Row) + Real_8(Col) - 2.0*Real_8(Starting_Col) + 1.0);
         A(Row, Col) := Prime_Factors / Denominator;
      end loop;
      end loop;
   end Get_Hilbert_Matrix;

   pragma Inline (Get_Hilbert_Matrix);

   ----------------------------------
   -- Print_Extended_Real_Settings --
   ----------------------------------

   procedure Print_Extended_Real_Settings
   is
      Bits_In_Radix : constant := Desired_No_Of_Bits_In_Radix;
   begin
      new_line(1); 
      put ("    Desired_Decimal_Digit_Precision =");
      put (Integer'Image(Desired_Decimal_Digit_Precision)); 
      new_line(1); 
      new_line(1); 
      put ("Number of decimal digits of precision requested:    ");
      put (Integer'Image(Desired_Decimal_Digit_Precision));
      new_line(1); 
      put ("Number of digits in use (including 2 guard digits): "); 
      put (e_Integer'Image(e_Real_Machine_Mantissa));
      new_line(1); 
      put ("These digits are not decimal; they have Radix:       2**(");
      put (e_Integer'Image(Bits_In_Radix)); put(")");
      new_line(1); 
      put ("In other words, each of these digits is in range:    0 .. 2**(");
      put (e_Integer'Image(Bits_In_Radix)); put(")"); put (" - 1.");
      new_line(1); 
      put ("Number of decimal digits per actual digit is approx: 9"); 
      new_line(2);
      put("Guard digits (digits of extra precision) are appended to the end of");
      new_line(1); 
      put("each number.  There are always 2 guard digits. This adds up to 18"); 
      new_line(1); 
      put("decimal digits of extra precision. The arithmetic operators, (""*"",");
      new_line(1); 
      put("""/"", ""+"" etc) usually produce results that are correct to all");
      new_line(1); 
      put("digits except the final (guard) digit.");
      new_line(2);
      put("If a number is correct to all digits except the final (guard) digit,");
      new_line(1); 
      put("expect errors of the order:");
      new_line(2);
      put(e_Real_Image (e_Real_Model_Epsilon / (One+One)**Bits_In_Radix, Aft => 20));
      new_line(2); 
      put("If you lose 2 digits of accuracy (i.e. both guard digits) instead");
      new_line(1); 
      put("of 1 (as in the above case) then you lose another 9 decimal digits");
      new_line(1); 
      put("of accuracy.  In this case expect errors of the order:");
      new_line(2);
      put(e_Real_Image (e_Real_Model_Epsilon, Aft => 20));
      new_line(1);

      Pause;

   end Print_Extended_Real_Settings;

begin

   Print_Extended_Real_Settings;

   Pause( 
     "Test 1: Jacobi Eigendecomposition of matrix A.",
     "   ",
     "The Jacobi Eigendecomposition of A is successful if the identities Q*Q' = I,",
     "Q'*Q = I and Q*E*Q' = A are satisfied. Here Q' denotes the transpose of Q, and",
     "E is any diagonal matrix. (E will hold the Eigvals.) If 15 digit Reals are used,",
     "then we expect the error in the calculation of A = Q*E*Q' to be (hopefully)",
     "a few parts per 10**15. In other words ||Q*E*Q' - A|| / ||A|| should be a few",
     "multiples of 10**(-15). Here ||*|| denotes the Frobenius Norm. Other matrix",
     "norms give slightly different answers, so its an order of magnitude estimate."
     );

   -- 
   --  Get A = Q*E*Q_tr, or E = Q_tr*A*Q.  
   --  E is diagonal with the eigs along the diag.
   --  V is orthogonal with the eig vecs as columns.

   Get_Hilbert_Matrix (A);

   --  Usually A is not symmetric.  Eigen_Decompose doesn't care about
   --  that. It uses the upper triangle of A, and pretends that A is
   --  symmetric. But for subsequent analysis,  we symmetrize:

   if false then -- use lower triangle of A to make a fully symmetric A:

     for Col in Starting_Col .. Final_Col  loop
     for Row in Col .. Final_Row  loop
        A(Col, Row) := A(Row, Col);  -- write lower triangle of A to upper triangle
     end loop;
     end loop;

   else                -- use upper triangle of A to make a fully symmetric A:

     for Col in Starting_Col .. Final_Col  loop
     for Row in Starting_Col .. Col  loop
        A(Col, Row) := A(Row, Col);  -- write lower triangle of A to upper triangle
     end loop;
     end loop;

   end if;

   A_true := A; -- Save original A

   Eigen_Decompose 
     (A                      => A,   -- A is destroyed
      Q_tr                   => Q_tr,
      Eigenvals              => Eigenvals,
      No_of_Sweeps_Performed => No_of_Sweeps_Done,
      Total_No_of_Rotations  => No_of_Rotations,
      Final_Col              => Final_Col,
      Start_Col              => Starting_Col,
      Eigenvectors_Desired   => True);

   Sort_Eigs
     (Eigenvals          => Eigenvals,
      Q_tr               => Q_tr,
      Start_Col          => Starting_Col,
      Final_Col          => Final_Col,
      Sort_Eigvecs_Also  => True);


   Get_Err_in_Reassembled_Q_and_A 
     (A                 => A_True, 
      Q_tr              => Q_tr, 
      E                 => Eigenvals,
      Final_Col         => Final_Col, 
      Starting_Col      => Starting_Col, 
      Frobenius_QtrQ_Err => Frobenius_QtrQ_Err, 
      Frobenius_QQtr_Err => Frobenius_QQtr_Err, 
      Frobenius_QEQ_Err  => Frobenius_QEQ_Err);

   --  Froebenius norm fractional error:
   --                Max_Error_F  = ||Err_Matrix|| / ||A||

   new_line;
   new_line;
   put(" No of sweeps performed, and Total_No_of_Rotations / (N*(N-1)/2) ="); 
   new_line;
   put(Real_8 (No_of_Sweeps_Done)); 
   put(Real_8 (No_of_Rotations) / (N*(N-1.0)/2.0));
   new_line;
   put(" Err in I-Q*Q'   (Q = orthogonal) is  ||I-Q*Q'|| / ||I|| ="); 
   put(e_Real_Image (Frobenius_QQtr_Err));
   new_line;
   put(" Err in I-Q'*Q   (Q = orthogonal) is  ||I-Q'*Q|| / ||I|| ="); 
   put(e_Real_Image (Frobenius_QtrQ_Err));
   new_line;
   put(" Err in A-Q*E*Q' (E = eigenvals) is ||A-Q*E*Q'|| / ||A|| ="); 
   put(e_Real_Image (Frobenius_QEQ_Err));
   new_line;

   Pause( 
    "Test 2: Estimate of error in eigenvalues.",
    "   ",
    "Actual error can be calculated accurately by rerunning the calculation",
    "with higher precision, but the following method works very well. Simply ",
    "subtract the calculated eigenvalues from the Rayleigh eigs. In other",
    "words, calculate Err = E - Q'*A*Q, where the column vectors of Q are the",
    "calculated Eigenvectors of A, and the diagonal of E contains the", 
    "calculated Eigenvalues of A. The diagonal elements of Err follow:"
    );

   new_line (1);
   for I in Starting_Col .. Final_Col loop
      new_line;
      put("Eigenvalue number"); put (Index'Image (I));
      new_line;
      put("Eigenvalue      = "); 
      put(e_Real_Image (Eigenvals(I))); 
      new_line (1);
      put("Err in Eigenval = "); 
      put(e_Real_Image ((Rayleigh_Eig(A_True, Q_tr, I)-Eigenvals(I)))); 
   end loop; 
   new_line (2);

   Pause( 
    "Final comment.",
    "Errors were unusually small here.  (Got the 1st guard digit right.)",
    "The Jacobi Eigendecomposition is good at minimizing accumulation of error.",
    "In large flt. pt. calculations, you should expect to get both guard digits",
    "wrong (i.e. lose an additional 9 decimal digits of accuracy).",
    "That's why 2 guard digits are used."
    );
   Pause( 
    "Also worth noting that no rounding was done in the Eigendecomposition",
    "arithmetic.  Rounding would lower accuracy and slow down the execution.",
    "In other algorithms rounding may be beneficial or necessary. It can be",
    "done intermittantly, not at all, or frequently. Calls to routines that",
    "perform rounding are inserted by the user, entirely at his discretion."
    );

end;
