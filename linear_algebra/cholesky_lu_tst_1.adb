
-- Test LU decomposition on a real valued square matrix.

with Ada.Numerics.Generic_elementary_functions;
with Cholesky_LU;
with Text_IO; use Text_IO;
with Test_Matrices;

procedure cholesky_lu_tst_1 is

   type Real is digits 15;

   subtype Index is Integer range 0..191;

   Starting_Index : constant Index := Index'First + 0;
   Final_Index    :          Index := Index'Last  - 0;


   package Math is new Ada.Numerics.Generic_elementary_functions(Real);
   use Math; --for sqrt
   package lu is new Cholesky_LU (Real, Index, Matrix);
   use lu;
   package rio is new Float_IO(Real);
   use rio;
   package iio is new Integer_IO(Integer);
   use iio;
   package Make_Square_Matrix is new Test_Matrices (Real, Index, Matrix);
   use Make_Square_Matrix;


   type Matrix is array(Index, Index) of Real;

   type Real_Extended  is digits 15;   -- or 18 on intel

   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;
   Min_Allowed_Real : constant Real := Two**(Real'Machine_Emin / 4);

   e_Sum : Real_Extended;
   Sum : Real;

   Diag_Inverse : Row_Vector := (others => Zero);
   Zero_Vector : constant Row_Vector := (others => Zero);

   A, A_LU, Err, L, U : Matrix := (others => (others => Zero));
   Trace, Relative_Err   : Real;
   IO_Final_Index : Integer := 4;

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

   -----------------------------------
   -- Transpose_of_Left_Times_Right --
   -----------------------------------
  
   function Transpose_of_Left_Times_Right
     (A, B         : in     Matrix;
      Final_Row    : in     Index := Final_Index; 
      Final_Col    : in     Index := Final_Index;
      Starting_Row : in     Index := Starting_Index; 
      Starting_Col : in     Index := Starting_Index)
      return Matrix
   is
      Sum : Real := Zero;
      Result : Matrix := (others => (others => Zero));
   begin
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         Sum := Zero;
         for k in Starting_Col .. Final_Col  loop
            Sum := Sum + A(k, Row) * B(k, Col);
         end loop;
	 Result(Row, Col) := Sum;
      end loop;
      end loop;
      return Result;
   end Transpose_of_Left_Times_Right;

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
     Integer'Image(Zero_Vector'length-(Integer(Starting_Index)-Integer(Index'First))));
     new_Line;
   put("Input Size Of Matrix To Invert (enter an Integer)"); new_Line;
   get(IO_Final_Index);
   Final_Index := Starting_Index + Index (IO_Final_Index-1);
   
   Pause( 
     "Test 1: Cholesky's LU Decomposition of matrix A = L*U. A positive definite",
     "A is obtained from A = B'*B + eps*I, where B' = Transpose (B). If 15 digit",
     "Reals are used, then expect error in calculation of A = L*U to be a few",
     "parts per 10**15. In other words ||L*U - A|| / ||A|| should be a few",
     "multiples of 10**(-15). Here |*| denotes the Frobenius Norm. Other matrix",
     "norms give slightly different answers, so its an order of magnitude estimate."
     );

  new_line;
  for Chosen_Matrix in Matrix_id loop
  
     put("For matrix A = B'*B + eps*I, where B is "); 
     put(Matrix_id'Image(Chosen_Matrix)); Put(":"); 
     new_line; 
 
     -- Get A:

     Init_Matrix (A, Chosen_Matrix);

     A := Transpose_of_Left_Times_Right (A, A);

     -- A = A'A is now positive semi-definite. Shift all eigs +Epsilon, where
     -- Epsilon = 2**(-Real'Machine_Mantissa / 2) * Upper_Bound_of_Largest_Eig:

     Trace := Min_Allowed_Real;
     for Col in Starting_Index .. Final_Index loop
        Trace := Trace + A(Col, Col);
     end loop;

     for Col in Starting_Index .. Final_Index loop
        A(Col, Col) := Two ** (-Real'Machine_Mantissa / 2) * Trace + A(Col, Col);
     end loop;

     -- Get A = L * U:

     A_LU := A;
 
     LU_decompose 
       (A               => A_LU, 
        Diag_Inverse    => Diag_Inverse, 
        Final_Index     => Final_Index,
        Starting_Index  => Starting_Index);
 
     --  L, U initialized to 0:

     for Row in Starting_Index .. Index'Last loop
     for Col in Row .. Index'Last loop
        U(Row, Col) := A_LU(Row, Col);
     end loop;
     end loop;
  
     for Col in Starting_Index .. Index'Last loop
     for Row in Col .. Index'Last loop
        L(Row, Col) := A_LU(Row, Col);
     end loop;
     end loop;
  
     --  Multiply Original L and U as test.  Get Max error:
    
     for I in Starting_Index..Final_Index loop
     for J in Starting_Index..Final_Index loop
    
        e_Sum := +0.0;
        for K in Starting_Index .. Final_Index loop
           e_Sum := e_Sum + Real_Extended (L(I, K)) * Real_Extended (U(K, J));
        end loop;
        Sum:= Real(e_Sum);
    
        --  Product(I,J) := Sum;
        --  Calculate the error:
    
        Err(i, j) :=  Sum - A(i, j);
  
     end loop;
     end loop;
     
     Relative_Err := Frobenius_Norm (Err) / (Frobenius_Norm (A) + Min_Allowed_Real);
  
     put("        Err in A - L*U is  ~  ||A - L*U|| / ||A|| ="); 
     put(Relative_Err);
     new_line; 
  
  end loop; -- Matrix_id
     
end;
