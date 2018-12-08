
-- PACKAGE e_Jacobi_Eigen
--
-- Extended precision version of Jacobi's iterative algorithm for 
-- eigen-decomposition of square real-valued symmetric matrices.
-- 
-- PROCEDURE Eigen_Decompose 
--
-- Works on arbitrary diagonal blocks of input matrix. For other blocks just
-- copy the matrix to desired position; copy overhead is negligable compared
-- to the O(N^3) running time of the decomposition.
--

generic

   type e_Real is private;

   type Index is range <>;
   type Matrix is array (Index, Index) of e_Real;

   type Real_8 is digits <>; 
   --  Real_8 is for easy communication with e_Real. 
   --  Must be digits 15 or more. 
   --  The function "+" below translates  Real_8  to  e_Real.
   
   --  Exported by   Extended_Real:
   
   type e_Integer is range <>;

   with function e_Real_Model_Epsilon return e_Real is <>;
   with function e_Real_Machine_Emin  return e_Integer is <>;
   with function e_Real_Machine_Emax return e_Integer is <>;
   with function e_Real_Machine_Radix return Real_8 is <>;
 --with function Exponent (X : e_Real) return e_Integer is <>;

   with function "*" (X : e_Real; Y : e_Real) return e_Real is <>;
   with function "/" (X : e_Real; Y : e_Real) return e_Real is <>;
   with function "-" (X : e_Real; Y : e_Real) return e_Real is <>;
   with function "+" (X : e_Real; Y : e_Real) return e_Real is <>;
   with function "**" (X : e_Real; I : Integer) return e_Real is <>;
   with function "<=" (X : e_Real; Y : e_Real) return Boolean is <>;
   with function ">=" (X : e_Real; Y : e_Real) return Boolean is <>;
   with function  "<" (X : e_Real; Y : e_Real) return Boolean is <>;
   with function  ">" (X : e_Real; Y : e_Real) return Boolean is <>;
   with function  "-" (X : e_Real) return e_Real is <>;
   with function  "+" (X : Real_8) return e_Real is <>;
   with function "Abs" (X : e_Real) return e_Real is <>;
 --with function  "=" (X : e_Real; Y : e_Real) return Boolean is <>;
   
   --  Exported by   Extended_Real.Elementary_Functions:
   
   with function Sqrt (X : e_Real) return e_Real is <>;
   with function Reciprocal_Sqrt (X : e_Real) return e_Real is <>; --1/Sqrt(X)
 --with function Divide (X, Y : e_Real) return e_Real is <>; -- X/Y
 --with function "**" (X : e_Real; Y : e_Real) return e_Real is <>;
 --with function  Cos (X : e_Real) return e_Real is <>;
 --with function  Sin (X : e_Real) return e_Real is <>;

package e_Jacobi_Eigen is

   subtype Real is e_Real;

   type Col_Vector is array(Index) of Real;


   -- PROCEDURE Eigen_Decompose 
   --
   -- Standard Jacobi iterative eigendecomposition. The routine returns 
   -- eigenvectors and eigenvalues of any real-valued square symmetric matrix.
   --
   -- The orthonormal (unordered) eigenvectors are the Columns of Q.
   -- The orthonormal (unordered) eigenvectors are returned as the Rows of Q'=Q_tr.
   -- Eigenvals (returned in array Eigenvals) are ordered the same as Eigvecs in Q.
   -- So A = QEQ'. The diagonal elements of diagonal matrix E are the eigvals.
   -- The routine performs the eigen-decomposition on arbitrary square
   -- diagonal blocks of matrix A. 
   -- It is assumed the blocks are symmetric.
   -- The upper left corner of the square matrix is (Start_Col, Start_Col).
   -- The lower rgt  corner of the square matrix is (Final_Col, Final_Col).
   -- Matrix A doesn't need to be positive definite, or semi-definite.
   -- If Eigenvectors_Desired = False, then Q_tr is not calculated.
   --
   -- Routine only sees and operates on the upper triangle of matrix.
   --
   -- Input matrix A is destroyed. Save a copy of A if you need it.
   --
   -- Eigenvectors of A are returned as the ROWS of matrix:  Q_tr
   --
   -- so   Q_tr * A * Q = Diagonal_Eigs
   --
   procedure Eigen_Decompose
     (A                      : in out Matrix;     -- destroyed
      Q_tr                   : out Matrix;        -- rows of Q_tr are the eigvectors
      Eigenvals              : out Col_Vector;
      No_of_Sweeps_Performed : out Natural;
      Total_No_of_Rotations  : out Natural;
      Start_Col              : in Index        := Index'First;
      Final_Col              : in Index        := Index'Last;
      Eigenvectors_Desired   : in Boolean      := False);


   procedure Sort_Eigs
     (Eigenvals         : in out Col_Vector;
      Q_tr              : in out Matrix;          -- rows of Q_tr are the eigvectors
      Start_Col         : in     Index   := Index'First;
      Final_Col         : in     Index   := Index'Last;
      Sort_Eigvecs_Also : in     Boolean := False);

private

   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;
   Half : constant Real := +0.5;

   e_Radix : constant Real := +e_Real_Machine_Radix;
   e_Real_Safe_Min : constant Real := e_Radix**Integer(e_Real_Machine_Emin / 4);
   e_Real_Safe_Max : constant Real := e_Radix**Integer(e_Real_Machine_Emax / 4);

   Min_Exp : constant Integer := Integer (e_Real_Machine_Emin);
   Min_Allowed_Real : constant Real := e_Radix**(Min_Exp/2 + Min_Exp/4);

end e_Jacobi_Eigen;

