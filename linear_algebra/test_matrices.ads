   
-- Collection of square matrices.  Most are ill conditioned, singular
-- or badly scaled.
-- Several matrices are from John Burkardt's fortran 90 Test_Mat.

generic
  type real is digits <>;
  type Index is range <>;
  type Matrix is array(Index, Index) of Real;
package Test_Matrices is

   -- In notes below, N is the Order of the Matrix.

   type Matrix_Id is 
     (Laguerre,        -- eig vec calc fails w/ balancing for lower tri version.
      Ding_Dong,       -- eigs clustered very near (18+ sig. figs) +/- pi, symmetric
      Lesp,             
      -- Sensitive eigs., tridiag, not symmetric, transpo hard, triggers underflows,
      -- Balancing v. bad., hessenberg pivoting solves difficult eigen calc.

      Vandermonde,   -- non-symmetric
      Combin,        -- solutions with maximally high error for given condition num?
      Kahan,
      Kahan_2,
      Kahan_Col_Scaled,
      Kahan_Row_Scaled,
      Sampling,    -- ill conditioned eigs, eigs = 0..N-1, ok up to ~ 21x21
      Sampling_1,  -- well behaved, eigs = 0..N-1
      Companion_2,      -- Eigs are roots of companion polynomial.
      Companion_1,      -- Badly scaled, difficult eigs, balancing essential for eigs.
      Companion_0,      -- Balancing important.
      Companion,        -- Eigs are roots of companion polynomial
      Lower_Integers,   -- almost triangular, eigs = 1, 2, 3 ...
      Pas_Fib,          -- ill-cond N > 9; rescaled if big
      Pascal_Symmetric,
      Pascal,           
      -- ill-cond N > 9; sing. vals ~ x, 1/x, all>0, rescaled if big, balancing bad
      -- if lower triangular

      Pascal_Row_Scaled,
      Pascal_Col_Scaled,
      Frank_0,     -- upper Hess, ill cond eigs, eig prod = 1, balance bad, 17x17max
      Frank_1,         -- upper Hess, ill cond eigs, eig prod = N, balance bad
      Frank_2,         -- symmetric, easy, product of eigs = 1
      Fiedler_0,       -- symm, non singular
      Fiedler_1,       -- non symm, poor-cond eigs
      Fibonacci,       -- non symm, non singular, tri diag, max defective
      Wilkinson_Minus, -- symm. W- only if odd N: then just one 0.0 eigval
      Hilbert,     -- exactly represented in 15 digit Real only up to ~ 20x20
      Lotkin,      -- like Hilbert, but not symmetric
      Clustered,   -- eigs clustered, poorly scaled
      Zielke_0,    -- symmetric, clustered eigs, ill conditioned
      Zielke_1,    -- symmetric, clustered eigs, singular
      Zielke_2,    -- symmetric, clustered eigs, singular
      Gear_0,      -- non symm, non singular, one small sing val
      Gear_1,      -- non symm, singular
      Diag_Test,   -- one small singular val.
      Moler,       -- 1 small eig, prod. of eigs = 1, eig sum = N(N+1)/2
      Peters,      -- 1 small eig, like Givens_Moler, lower triangular
      Peters_0,    -- 1 small eig, like Givens_Moler, upper triangular
      Peters_1,    -- ok w. row pivoting, LU
      Peters_2,    -- bad w/ row pivoting, LU
      Gregory,     -- symmetric, all but 2 of the eigs = 1.
      Anti_Hadamard_Upper_Tri, -- 1 small sing.val, non singular
      Anti_Hadamard_Lower_Tri, -- 1 small sing.val, non singular
      Wilkinson_Plus,     -- strictly its W+ only for odd N.
      Wilkinson_Plus_2I,
      U_Hard,             -- slow convergence on Golub SVD, all 1's but 1st col=0
      Zero_Cols_and_Rows, -- 1st 3 rows and cols all 0; else its 1.
      Easy_Matrix, 
      Symmetric_Banded,   -- not diagonally dominant, moderately ill-conditioned
      Small_Diagonal,     -- but not symmetric
      Trench,
      Trench_1,
      Forsythe_0,         -- eigs on circle of radius 1, center=A=0.
      Forsythe_1,         -- eigs on circle of radius 1, center=A=1.
      Forsythe_Symmetric,
      Zero_Diagonal,      -- if odd then N singular, eg N X N = 13x13
      QR_Test,            -- 4x4 hard on QR
      Ring_Adjacency_0,   -- Symm.
      Ring_Adjacency_1,   -- Non-Symm.
      Upper_Tri_K,
      Lower_Tri_K,        -- lower tri, non-singular, v high condition num
      Upper_Ones,         -- triangular
      Lower_Ones,         -- triangular
      All_Ones,
      Redheff, -- sensitive eigs, hard tst - slow convergence on Peters_Eigen.
      Chow,    -- ill-conditioned eigs, N / 2 eigs = 0, balancing bad, lower hess.
      Chow1,   -- chow w/ alpha=-1.05
      Chow2,   -- chow w/ alpha=gamma
      Chow3,   -- non-sym, full
      Lehmer,  -- easy
      Random_1_bit_anti, -- anti-sym, singular, on [0,1] new seed each run, hard tst
      Random_1_bit,      -- new seed each run.
      Random_32_bit,     -- on [0,1), new seed each run. 
      All_Zeros);

   procedure Init_Matrix 
     (M              : out Matrix;
      Desired_Matrix : in  Matrix_Id := Easy_Matrix;
      Starting_Index : in  Index := Index'First;
      Max_Index      : in  Index := Index'Last);
     
   -- Can optionally add constant Matrix_Addend to some of the Matrices.
   -- Increases the variety of tests, and usually makes them harder to 
   -- decompose.  The matrices are:
   --
   --   Anti_Hadamard_Lower_Tri, Anti_Hadamard_Upper_Tri, Lower_Tri_K, Upper_Tri_K
   --   Pascal_Col_Scaled, Pascal_Row_Scaled, Pascal
   --   Lower_Ones, Upper_Ones, Lower_Integers
   --   Frank_0, Frank_1, Fibonacci, Peters, Peters_0
   --   Kahan_Row_Scaled, Kahan_Col_Scaled, Kahan_Col_Scaled_2, Kahan

 --Matrix_Addend : constant Real := +1.0e-5;
   Matrix_Addend : constant Real := +0.0;

   procedure Transpose 
     (A              : in out Matrix;
      Starting_Index : in     Index     := Index'First;
      Max_Index      : in     Index     := Index'Last);

  procedure Symmetrize
    (A : in out Matrix);

end Test_Matrices;
