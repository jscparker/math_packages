
-- Demonstrate Jacobi Eigendecomposition of real valued square matrices.

with Ada.Numerics.Generic_Elementary_Functions;
with extended_real;
with extended_real.elementary_functions;
with extended_real.io;
with e_jacobi_eigen;
with text_io; use text_io;

procedure e_jacobi_eigen_tst_1 is

   type Real_8 is digits 15;

   package mth is new Ada.Numerics.Generic_Elementary_Functions (Real_8);
   use mth;
   package ext is new Extended_Real (Real_8);
   use ext; 
   package eio is new ext.IO;
   use eio;
   package fnc is new Ext.Elementary_Functions (Sqrt, Log, Exp, Arcsin);
   use fnc;

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
   
   N : constant Real_8 := Real_8 (Final_Col) - Real_8 (Starting_Col) + 1.0;

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

      Prime_Factors := (+166966608033225.0) * (+580027.0) * Two**(-68);
  
      for Row in Starting_Col .. Final_Col loop
      for Col in Starting_Col .. Final_Col loop
         Denominator := +(Real_8(Row) + Real_8(Col) - 2.0*Real_8(Starting_Col) + 1.0);
         A(Row, Col) := Prime_Factors / Denominator;
      end loop;
      end loop;
   end Get_Hilbert_Matrix;

   pragma Inline (Get_Hilbert_Matrix);

begin

   -- 
   --  Get A = Q*E*Q_tr, or E = Q_tr*A*Q.  
   --  E is diagonal with the eigs along the diag.
   --  V is orthogonal with the eig vecs as columns.

   Get_Hilbert_Matrix (A);

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

   new_line;
   put ("Matrix size (N):              ");
   put (Real_8'Image(N));
   new_line;
   put ("First Eigenvalue:             ");
   put (e_Real_Image(Eigenvals(Starting_Col)));
   new_line;
   put ("Smallest Eigenvalue:          ");
   put (e_Real_Image(Eigenvals(Final_Col)));
   new_line;
   put ("No_of_Rotations / (N(N-1)/2): ");
   put (Real_8'Image(Real_8(No_of_Rotations) / (0.5*N*(N-1.0))));

end;
