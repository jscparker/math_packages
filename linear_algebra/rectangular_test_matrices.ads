
generic
  type Real is digits <>;
  type R_Index is range <>;
  type C_Index is range <>;
  type Matrix  is array(R_Index, C_Index) of Real;
package Rectangular_Test_Matrices is

   type Matrix_Id is 
     (Zero_Diagonal,    -- all 1's with 0 down the diagonal
      Upper_Ones,       -- 
      Lower_Ones,       --
      Zero_Cols_and_Rows,  -- 1st 4 rows and cols all 0; else its Upper_Ones
      Random,       -- new seed each run. 16 bit elements (1 bit => 4x4 often singular)
      Frank, 
      Ding_Dong,   -- eigs clustered ~ pi; 99x99 matrix has lots of pi's to 18 sig. figs.
      Vandermonde, -- max allowed matrix size is 256 x 256 without over/under flow
      All_Zeros,
      All_Ones);

   procedure Init_Matrix 
     (M               : out Matrix;
      Desired_Matrix  : in  Matrix_Id := Random;
      Final_Row       : in  R_Index   := R_Index'Last;
      Final_Col       : in  C_Index   := C_Index'Last;
      Starting_Row    : in  R_Index   := R_Index'First;
      Starting_Col    : in  C_Index   := C_Index'First);
     
end Rectangular_Test_Matrices;
