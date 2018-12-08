
----------------------------------------------------------------------
-- package Peters_Eigen, eigendecomposition
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
---------------------------------------------------------------------------

-- PACKAGE Peters_Eigen
--
-- Performs an eigen-decomposition on general square real matrices. All 
-- arithmetic is real. The set of (potentially) complex eigenvalues is 
-- returned as a pair of vectors W_r, W_i. (If the Matrix is symmetric, 
-- then the eigenvalues should be real valued.)
--
-- Peters_Eigen is based on the original Eispack hqr2.f [1] routine with 
-- several changes. The Hessenberg reduction was replaced with one based 
-- on Givens rotations. The Eispack QR algorithm is retained with small 
-- changes (guards to prevent floating point problems that occur with near-
-- zero numbers, and slightly modified inner loops to improve accuracy).
--
-- Eigenvalue error can be estimated if you have a higher precision floating 
-- point available. Most Intel CPUs provide 18 digit floating point, so
-- instantiate Peters_Eigen with  
--
--                 "type Real is digits 18".
--
-- 1. Peters and Wilkinson, Num. Math. 16, 181-204 (1970)
--
-- 2. J G F Francis, "The QR Transformation, I", The Computer Journal,
--    vol. 4, no. 3, pages 265-271 (1961)
--
-- 3. Vera N Kublanovskaya, "On some algorithms for the solution of the
--    complete eigenvalue problem," USSR Computational Mathematics and
--    Mathematical Physics, vol. 1, no. 3, pages 637657

generic

   type Real is digits <>;
   type Index is range <>;
   type Matrix is array (Index, Index) of Real;

package Peters_Eigen is

   type Col_Vector is array(Index) of Real;

   type Balance_Code is (Disabled, Partial, Full);

   -- The default disables balancing. Balancing isn't recommended for general use.
   -- If balancing is enabled, then eigenvectors are not orthogonal in general.

   procedure Decompose
     (A                    : in out Matrix;
      Z_r, Z_i             :    out Matrix;
      W_r, W_i             :    out Col_Vector;
      Id_of_Failed_Eig     :    out Integer;
      Starting_Col         : in     Index        := Index'First;
      Final_Col            : in     Index        := Index'Last;
      Eigenvectors_Desired : in     Boolean      := True;
      Balance_Policy       : in     Balance_Code := Disabled);
   --
   --  If you set "Eigenvectors_Desired := False" then only Eigenvalues are
   --  calculated. 
   --
   --  Decomposition is performed on arbitrary diagonal blocks of matrix A.
   --
   --  A is destroyed (overwritten) by procedure Decompose.
   --
   --  Upper left  corner of the diagonal block is (r,c) = (Starting_Col,Starting_Col)
   --  Lower right corner of the diagonal block is (r,c) = (Final_Col, Final_Col)
   --
   --  Convergence is judged to be OK in range: Id_of_Failed_Eig+1 .. Final_Col.
   --
   --  If Eigenvectors_Desired = True, then
   --
   --     Real      parts of eigvecs are returned as columns of Z_r.
   --     Imaginary parts of eigvecs are returned as columns of Z_i.
   --
   --     Eigenvectors are normalized.
   --     Computes right eigenvectors z: A*z = lambda*z
   --
   --  If Eigenvectors_Desired = False, then
   --
   --     Z_i is *not* initialized, (which frees up space if matrices are large).
   --
   --  Balancing is disabled by default.
   --
   --     The Q matrices won't be orthogonal if the matrix is balanced.
   --     Sometimes balancing improves eigenvalue accuracy, sometimes it
   --     degrades accuracy.


   procedure Sort_Eigs_And_Vecs
     (Z_r, Z_i     : in out Matrix;   -- eigvecs are columns of Z
      W_r, W_i     : in out Col_Vector;
      Starting_Col : in     Index   := Index'First;
      Final_Col    : in     Index   := Index'Last);
   --
   -- Sorted according to size of: Sqrt (W_r**2 + W_i**2), largest first.
   -- 
   -- Notice that the eigvectors (Z_r, Z_i) are in out, so must be initialized.
   -- 
   -- Largest Eigs and their associated vectors are placed at: 
   -- Index'First, Index'First+1, ...
   -- Notice that Sqrt (W_r**2 + W_i**2) is used, so numerical noise may mean
   -- that the eigs are not sorted by magnitude of W_r**2 + W_i**2 alone.

   function Norm
     (V_r, V_i     : in Col_Vector;
      Starting_Col : in Index := Index'First;
      Final_Col    : in Index := Index'Last)
      return Real;
   --
   --  Norm of complex vector V with components V(j) =  (V_r(j), V_i(j))

end Peters_Eigen;
