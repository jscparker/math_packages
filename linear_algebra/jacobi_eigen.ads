
--------------------------------------------------------------------------
-- package Jacobi_Eigen, Jacobi iterative eigen-decomposition
-- Copyright (C) 2008-2018 Jonathan S. Parker
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

-- package Jacobi_Eigen
--
-- Jacobi's iterative algorithm for eigen-decomposition of
-- square real-valued symmetric matrices.
-- 
-- The Jacobi method converges quadratically and with high
-- reliability, but it is usually several times slower than the 
-- standard Golub-Reinsch algorithm (unless the matrix is small).
-- Jacobi is usually best if accuracy and reliability are more
-- important than speed.
--
-- procedure Eigen_Decompose 
--
-- Works on arbitrary diagonal blocks of input matrix. For other blocks just
-- copy the matrix to desired position; copy overhead is negligable compared
-- to the O(N^3) running time of the decomposition.
--
-- Procedure Eigen_Decompose is based on the Heinz Rutishauser ALGOL routine.
-- If you want to see the original program in its full glory, you should 
-- be able to find it by googling "Rutishauser, Jacobi, Algol, ethistory".
-- It was written over half a century ago! Changes are few:
-- Slightly different calculation of the rotation angles.
-- Also a bit more care in avoiding overflows; the test suite
-- caught a few of them in the original routine.
--

generic

   type Real is digits <>;
   type Index is range <>;
   type Matrix is array (Index, Index) of Real;

package Jacobi_Eigen is

   type Col_Vector is array(Index) of Real;


   -- procedure Eigen_Decompose 
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
      Q_tr                   :    out Matrix;     -- rows of Q_tr are the eigvecs
      Eigenvals              :    out Col_Vector;
      No_of_Sweeps_Performed :    out Natural;
      Total_No_of_Rotations  :    out Natural;
      Start_Col              : in     Index        := Index'First;
      Final_Col              : in     Index        := Index'Last;
      Eigenvectors_Desired   : in     Boolean      := True);


   procedure Sort_Eigs
     (Eigenvals         : in out Col_Vector;
      Q_tr              : in out Matrix;          -- rows of Q_tr are the eigvecs
      Start_Col         : in     Index   := Index'First;
      Final_Col         : in     Index   := Index'Last;
      Sort_Eigvecs_Also : in     Boolean := True);


   Standard_Threshold_Policy : constant Boolean := True;
   --  True is faster.
   --  False sometimes improves accuracy if the matrix is badly scaled.

end Jacobi_Eigen;

