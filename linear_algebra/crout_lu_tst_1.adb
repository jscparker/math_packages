
-- Test LU decomposition on a real valued square matrix.

with Ada.Numerics.Generic_elementary_functions;
with Crout_LU;
with Text_IO; use Text_IO;
with Test_Matrices;

procedure crout_lu_tst_1 is

   type Real is digits 15;

   type Index is range 0..69;

   Starting_Index : constant Index := Index'First + 0;
   Final_Index    :          Index := Index'Last  - 0;

   type Matrix is array(Index, Index) of Real;

   type Real_Extended is digits 18;
   -- use 18 on Intel for better error estimates.

   package math is new Ada.Numerics.Generic_elementary_functions(Real); --for Sqrt
   use math;
   package lu is new crout_lu (Real, Index, Matrix);
   use lu;
   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Integer);
   use iio;
   package Make_Square_Matrix is new Test_Matrices (Real, Index, Matrix);
   use Make_Square_Matrix;


   e_Sum : Real_Extended;
   Sum : Real;

   Permute : Rearrangement;
   Scale   : Scale_Vectors;
     
   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;
   Min_Allowed_Real : constant Real := Two**(Real'Machine_Emin - Real'Machine_Emin / 8);

   Zero_Vector : constant Row_Vector := (others => Zero);

   A, A_LU, Err, L, U : Matrix := (others => (others => Zero));
   Relative_Err   : Real;
   IO_Final_Index : Integer := 4;

   Scale_the_Matrix : constant Boolean := True; 

   -----------
   -- Pause --
   -----------

   procedure Pause (s1,s2,s3,s4,s5,s6,s7,s8 : string := "") is
     Continue : Character := ' ';
   begin
     New_Line;
     if S1 /= "" then put_line (S1); end if;
     if S2 /= "" then put_line (S2); end if;
     if S3 /= "" then put_line (S3); end if;
     if S4 /= "" then put_line (S4); end if;
     if S5 /= "" then put_line (S5); end if;
     if S6 /= "" then put_line (S6); end if;
     if S7 /= "" then put_line (S7); end if;
     if S8 /= "" then put_line (S8); end if;
     new_line;
     begin
        put ("Enter a character to continue: ");
        get_immediate (Continue);
        new_line;
     exception
	when others => null;
     end;
   end pause;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (A            : in     Matrix;
      Final_Row    : in     Index := Final_Index; 
      Final_Col    : in     Index := Final_Index;
      Starting_Row : in     Index := Starting_Index; 
      Starting_Col : in     Index := Starting_Index)
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

      Scaling := One / (Max_A_Val + Min_Allowed_Real);

      Sum := Zero;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * A(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;

begin

   put("Maximum matrix size is "& 
      Integer'Image (Zero_Vector'length-(Integer(Starting_Index)-Integer(Index'First))));
      new_Line;
   put("Input Size Of Matrix To Invert (enter an Integer)"); new_Line;
   get(IO_Final_Index);
   Final_Index := Index (Integer(Starting_Index) + IO_Final_Index - 1);
   
   Pause(
     "Test1: LU Decomposition of matrix A = P*L*U, where P is the permutation",
     "matrix that comes from Row pivoting. If 15 digit Reals are used, then we",
     "expect the error in the calculation of A = P*L*U to be (hopefully) a few",
     "parts per 10**15. In other words ||P*L*U - A|| / ||A|| should be a few",
     "multiples of 10**(-15). Here |*| denotes the Frobenius Norm. Other matrix",
     "norms give slightly different answers, so its an order of magnitude estimate.",
     " ",
     "Note: LU Decomp. with row pivoting always fails on (large) Peters_2 matrices."
     );

  new_line;
  for Chosen_Matrix in Matrix_id loop
  
     put("For matrix of type "); put(Matrix_id'Image(Chosen_Matrix)); Put(":"); 
     new_line; 
 
     -- Get A:

     Init_Matrix (A, Chosen_Matrix, Starting_Index, Final_Index);

     -- Get A = P * L * U:

     A_LU := A;
 
     LU_decompose 
       (A               => A_LU, 
        Scalings        => Scale, 
        Row_Permutation => Permute,
        Final_Index     => Final_Index,
        Starting_Index  => Starting_Index,
        Scaling_Desired => Scale_the_Matrix);

     -- if matrix A was scaled, then the L*U equals the scaled A.

     if Scale_the_Matrix then 
        Scale_Cols_Then_Rows
          (A              => A,
           Scalings       => Scale,
           Final_Index    => Final_Index,
           Starting_Index => Starting_Index);
     end if;
  
     for Row in Index loop
     for Col in Row .. Index'Last loop
        U(Row, Col) := A_LU(Row, Col);
     end loop;
     end loop;
  
     for Col in Index loop
     for Row in Col .. Index'Last loop
        L(Row, Col) := A_LU(Row, Col);
     end loop;
     end loop;

     for Col in Index loop
        L(Col, Col) := One;
     end loop;
  
     --  Multiply Original L and U as test.  Get Max error:
    
     for I in Starting_Index .. Final_Index loop
     for J in Starting_Index .. Final_Index loop
    
        e_Sum := +0.0;
        for K in Starting_Index .. Final_Index loop
           e_Sum := e_Sum + Real_Extended (L(I, K)) * Real_Extended (U(K, J));
        end loop;
        Sum:= Real(e_Sum);
    
        --  Product(I,J) := Sum;
        --  Calculate the error:
    
        Err(i, j) :=  Sum - A(Permute(i), j);
  
     end loop;
     end loop;
     
     Relative_Err := Frobenius_Norm (Err) / (Frobenius_Norm (A) + Min_Allowed_Real);
  
     put(" Err in A - P*L*U is  ~  ||A - P*L*U|| / ||A|| ="); 
     put(Relative_Err);
     new_line; 
  
  end loop; -- Matrix_id
     
end;
