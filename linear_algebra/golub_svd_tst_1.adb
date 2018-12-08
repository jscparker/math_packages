
with Ada.Numerics.Generic_Elementary_Functions;
with Golub_SVD;
with Test_Matrices;
with Ada.Text_IO; 
use  Ada.Text_IO;

--  Square matrix test. Demonstrates use of SVD when:  Final_Col = Final_Row.
--
--  It depends on the hardware and compiler, but in some tests procedure
--  SVD_Decompose runs about twice as fast if you use gcc's -ffast-math switch.
--  The switch -Ofast has the same affect. Best optimization on GNAT seems to be
--
--     gnatmake -O3 -funroll-loops -ffast-math  golub_svd_tst_1.adb

procedure golub_svd_tst_1 is

   type Real is digits 15;
 --type Real is digits 18;
 
 --Size_Of_Test_Matrix : constant := 19;
   Size_Of_Test_Matrix : constant := 137;
 --Size_Of_Test_Matrix : constant := 7010;

   Index_First : constant := 1;

   -- Sometimes it's faster if you use a storage matrix that's a little
   -- larger than the original matrix. (Depends on the hardware, the problem,
   -- and the size of the matrix.) A rule of thumb: Width of storage matrix
   -- should be an even number, greater than matrix size, and never 2**n.
   --
   -- Most important thing is to avoid 2**n matrix width.
   --
   -- Padding usually doesn't help a lot, but just in case let's add 1 or 2:

   Padding : constant := 2 - (Size_Of_Test_Matrix - Index_First + 1) mod 2;

   Matrix_Storage_Size : constant := Size_Of_Test_Matrix + Padding - Index_First + 1;

   pragma Assert (Matrix_Storage_Size < 2**24 and Matrix_Storage_Size > 1);

   subtype Col_Index is Integer range Index_First .. Matrix_Storage_Size;
   subtype Row_Index is Integer range Index_First .. Matrix_Storage_Size;
   
   type A_Matrix is array (Row_Index, Col_Index) of Real;
   pragma Convention (Fortran, A_Matrix); -- convention doesn't matter here.
 
   package lin_svd is new golub_svd (Real, Row_Index, Col_Index, A_Matrix);
   use lin_svd;
 
   package Make_Square_Matrices is new Test_Matrices (Real, Col_Index, A_Matrix);
   use Make_Square_Matrices;
 
   package Math is new Ada.Numerics.Generic_Elementary_Functions(Real);
   use Math;
 
   package rio is new Ada.Text_IO.Float_IO(Real); use rio;
 
   Final_Row    : constant Row_Index := Size_Of_Test_Matrix - 0;
   Final_Col    : constant Col_Index := Col_Index (Final_Row) - 0;
   -- in this test do only squares
 
   Starting_Col : constant Col_Index := Col_Index'First + 0;
   Starting_Row : constant Row_Index := Row_Index (Starting_Col); -- for test routines
   -- Starting_Row = Starting_Col is a requirement of SVD routine.
 
 
   A_true, A, SVD_Product : A_Matrix  := (others => (others => 0.0));
   V, VV_Product : V_Matrix;
   U, UU_Product : U_matrix;
 
   Singular_Vals : Singular_Vector := (others => 0.0);
 
   Col_id_of_1st_Good_S_Val : Col_Index;
 
   Del : Real;
   Max_Error_F, Frobenius_Norm_of_Identity : Real; 
 
   type Real_Extended is digits 18;
   -- subtype Real_Extended is Real;
   -- digits 15 is fine here, and the subtype Real_Extended is Real is also fine,
   -- but if you have an extended available, (eg 18 on intel) then it slightly 
   -- improves accuracy of the test measurements below.
 
   Sum : Real_Extended;
 
   Min_Allowed_Real : constant Real := 2.0 ** (Real'Machine_Emin / 8);

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (M            : in   A_Matrix;
      Final_Row    : in  Row_Index; 
      Final_Col    : in  Col_Index;
      Starting_Row : in  Row_Index; 
      Starting_Col : in  Col_Index)
      return Real
    is
      Max_A_Val : Real := 0.0;
      Sum, Scaling, tmp : Real := 0.0;
    begin
 
      Max_A_Val := 0.0;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         if Max_A_Val < Abs M(Row, Col) then Max_A_Val := Abs M(Row, Col); end if;
      end loop;
      end loop;

      Max_A_Val := Max_A_Val + 2.0 ** (Real'Machine_Emin + 4);
      Scaling := 1.0 / Max_A_Val;

      Sum := 0.0;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * M(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (M            : in   U_Matrix;
      Final_Row    : in  Row_Index; 
      Starting_Row : in  Row_Index)
      return Real
    is
      Max_A_Val : Real := 0.0;
      Sum, Scaling, tmp : Real := 0.0;
    begin
 
      Max_A_Val := 0.0;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Row .. Final_Row  loop
         if Max_A_Val < Abs M(Row, Col) then Max_A_Val := Abs M(Row, Col); end if;
      end loop;
      end loop;

      Max_A_Val := Max_A_Val + 2.0 ** (Real'Machine_Emin + 4);
      Scaling := 1.0 / Max_A_Val;

      Sum := 0.0;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Row .. Final_Row  loop
         tmp := Scaling * M(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (M            : in   V_Matrix;
      Final_Col    : in  Col_Index;
      Starting_Col : in  Col_Index)
      return Real
    is
      Max_A_Val : Real := 0.0;
      Sum, Scaling, tmp : Real := 0.0;
    begin
 
      Max_A_Val := 0.0;
      for Row in Starting_Col .. Final_Col  loop
      for Col in Starting_Col .. Final_Col  loop
         if Max_A_Val < Abs M(Row, Col) then Max_A_Val := Abs M(Row, Col); end if;
      end loop;
      end loop;

      Max_A_Val := Max_A_Val + 2.0 ** (Real'Machine_Emin + 4);
      Scaling := 1.0 / Max_A_Val;

      Sum := 0.0;
      for Row in Starting_Col .. Final_Col  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * M(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;

   -----------
   -- Pause --
   -----------

   procedure Pause (s0,s1,s2,s3,s4,s5,s6,s7 : string := "") is
     Continue : Character := ' ';
   begin
     New_Line;
     if S0 /= "" then put_line (S0); end if;
     if S1 /= "" then put_line (S1); end if;
     if S2 /= "" then put_line (S2); end if;
     if S3 /= "" then put_line (S3); end if;
     if S4 /= "" then put_line (S4); end if;
     if S5 /= "" then put_line (S5); end if;
     if S6 /= "" then put_line (S6); end if;
     if S7 /= "" then put_line (S7); end if;
     new_line;
     begin
	put ("Type a character to continue: ");
	get_immediate (Continue);
     exception
	when others => null;
     end;
   end pause;

   x : Real := 0.0;
  
begin

   Pause( 
     "Test1: Singular Value Decomposition of matrix A.",
     "The Singular Value Decomposition of A is successful if the identities V'*V = I,",
     "U'*U = I and U*S*V' = A are satisfied. Here V' denotes the transpose of V, and",
     "S is any diagonal matrix. (S will hold the sing. vals.) If 15 digit Reals are ",
     "used, then we expect the error in the calculation of A = U*S*V' to be (hopefully)",
     "a few parts per 10**15. In other words ||U*S*V' - A|| / ||A|| should be a few",
     "multiples of 10**(-15). Here ||*|| denotes the Frobenius Norm. Other matrix",
     "norms give slightly different answers, so it's an order of magnitude estimate."
     );

   New_Line(2);

   for Reps in 1 .. 1 loop

   for Chosen_Matrix in Matrix_id loop
 --for Chosen_Matrix in Hilbert .. Gear_1 loop
 --for Chosen_Matrix in Clustered .. Clustered loop
 --for Chosen_Matrix in Sampling_1 .. Sampling_1 loop
 --for Chosen_Matrix in Pascal_row_scaled .. Pascal_row_scaled loop
 --for Chosen_Matrix in Moler .. Moler loop
 --for Chosen_Matrix in U_hard .. U_hard loop
 --for Chosen_Matrix in pas_fib .. pas_fib loop
 --for Chosen_Matrix in kahan .. kahan_row_scaled loop
 --for Chosen_Matrix in hilbert .. hilbert loop
 --for Chosen_Matrix in vandermonde .. vandermonde loop
 --for Chosen_Matrix in lotkin .. lotkin loop
 --for Chosen_Matrix in wilkinson_plus .. wilkinson_plus loop
 --for Chosen_Matrix in wilkinson_minus .. wilkinson_minus loop
  
     put("For matrix of type "); put(Matrix_id'Image(Chosen_Matrix)); Put(":"); 
 
     Init_Matrix (A, Chosen_Matrix, Starting_Col, Final_Col);
   
--     for i in Col_Index range Starting_Col .. Final_Col loop
--     for j in Col_Index range Starting_Col .. Final_Col loop
--        A(i,j) := A(i,j) + 1.0e-7;
--     end loop;
--     end loop;

     -- Get A = U * W * V'  where W is diagonal containing the singular vals
   
     A_true := A;

     SVD_Decompose
       (A                => A, -- A is destroyed, so we saved a copy of A in A_true 
        U                => U,
        V                => V,
        S                => Singular_Vals,
        Id_of_1st_S_Val  => Col_id_of_1st_Good_S_Val,
        Starting_Col     => Starting_Col,
        Final_Col        => Final_Col,
        Final_Row        => Final_Row, 
        Matrix_U_Desired => True,
        Matrix_V_Desired => True);

     if Col_id_of_1st_good_S_Val /= Starting_Col then
        new_line; 
        put ("QR failed to fully converge at column = ");
        put(Col_Index'Image(Col_id_of_1st_good_S_Val));
        new_line(2); 
     else
        declare Max_Sing_Val, Min_Sing_Val : Real; begin
           new_line; 
           put ("Max Sing. Val =");
           Max_Sing_Val := Singular_Vals(Col_id_of_1st_good_S_Val);
           put(Real'Image(Max_Sing_Val)); put(", ");
           put ("Min Sing. Val =");
           Min_Sing_Val := Singular_Vals(Final_Col);
           put(Real'Image(Min_Sing_Val));
         --put ("Ratio of Max to Min (Condition Number) =");
         --put(Real'Image(Max_Sing_Val / (Min_Sing_Val + Min_Allowed_Real)));
        end;
     end if;
   
     if true then
        new_line;
        if Integer(Final_Col - Starting_Col + 1) > 2 then
      --for i in Starting_Col .. Final_Col loop
      --for i in Starting_Col .. Starting_Col+2 loop
        for i in Final_Col-2 .. Final_Col loop
           put (Singular_Vals(i));
        end loop;
        new_line(1);
        end if;
      --new_line; put ("Singular_Vals printed above.");
     end if;

     x := x + U(1,1) + V(1,1) + A(1,1);

     if false then
        new_line (1);
        put ("Condition number: ");
        put (Singular_Vals(Starting_Col) / (Singular_Vals(Final_Col)+2.0**(Real'Machine_Emin+16)));
        new_line (2);
     end if;

     -- Bypass the following tests for benchmarking.

     -- goto Endies;
 
     -- Get Product = V'* V = Transpose(V)*V
   
     for i in Col_Index range Starting_Col .. Final_Col loop
     for j in Col_Index range Starting_Col .. Final_Col loop
        Sum := 0.0;
        for k in Col_Index range Starting_Col .. Final_Col loop
          Sum := Sum + Real_Extended (V(k, i)) * Real_Extended (V(k, j)); 
 	       -- V' * V for convention Fortran is fastest.
        end loop;
        VV_Product(i, j) := Real (Sum);
     end loop;
     end loop;
   
     -- reuse VV_Product to store error = I - V'*V
   
     for i in Col_Index range Starting_Col .. Final_Col loop
     for j in Col_Index range Starting_Col .. Final_Col loop
        if i = j then
           Del := Abs (1.0 - VV_Product(i, j));
        else
           Del := Abs (0.0 - VV_Product(i, j));
        end if;
        VV_Product(i, j) := Del;
     end loop;
     end loop;
   
     Frobenius_Norm_of_Identity := Sqrt (Real (Final_Col) - Real (Starting_Col) + 1.0);
 
     Max_Error_F := 
       Frobenius_Norm (VV_Product, Final_Col, Starting_Col) /
          Frobenius_Norm_of_Identity; 
 
     new_line; 
     put(" Err in I -  V'*V  is ~ ||I -  V'*V || / ||I|| ="); 
     put(Max_Error_F);
 
     -- Get Product = U' * U =  transpose(U) * U
   
     for i in Row_Index range Starting_Row .. Final_Row loop
     for j in Row_Index range Starting_Row .. Final_Row loop
        Sum := 0.0;
        for k in Starting_Row .. Final_Row loop
          Sum := Sum + Real_Extended (U(k, i)) * Real_Extended (U(k, j));
        end loop;
        UU_Product(i, j) := Real (Sum);
     end loop;
     end loop;
   
     --  reuse UU_Product to store error I - U'*U
   
     for i in Row_Index range Starting_Row .. Final_Row loop
     for j in Row_Index range Starting_Row .. Final_Row loop
        if i = j then
           Del := Abs (1.0 - UU_Product(i, j));
        else
           Del := Abs (0.0 - UU_Product(i, j));
        end if;
        UU_Product(i, j) := Del;
     end loop;
     end loop;
   
     Frobenius_Norm_of_Identity := Sqrt (Real (Final_Col) - Real (Starting_Col) + 1.0);
 
     Max_Error_F := 
       Frobenius_Norm (UU_Product, Final_Row, Starting_Row) /
          Frobenius_Norm_of_Identity; 
 
     new_line; 
     put(" Err in I -  U'*U  is ~ ||I -  U'*U || / ||I|| ="); 
     put(Max_Error_F);
     new_line; 
 
     --  Get Product = U * W * V' = U * W * transpose(V)
   
     --  U is  usually Row_Index x Row_Index but can maybe use Col_Index for its 2nd
     --  V is  Col_Index x Col_Index
     --  S is  conceptually Row_Index x Col_Index, with Sing_Vals on its diagonal, 
     --  but its really in a vector on Col_Index.
     --  So the No_of_Singular vals is MIN of No_of_Rows, No_Of_Cols
     --  Singular vals are sorted so: largest are at beginning of Singular_Vals(k).
     --
   
     --  Get matrix product S * V'; place in matrix V':
   
       for k   in Col_Index range Starting_Col .. Final_Col loop
       for Col in Col_Index range Starting_Col .. Final_Col loop
          V(Col, k) := V(Col, k) * Singular_Vals(k);
       end loop;
       end loop;
   
       --  V' = S*V'  is now conceptually Row_Index x Col_Index:
   
       for i in Starting_Row .. Final_Row loop  -- U  is     Row_Index x Row_Index
       for j in Starting_Col .. Final_Col loop  -- V' is now Row_Index x Col_Index
   
          Sum := 0.0;
          for k in Col_Index range Starting_Col .. Final_Col loop
             Sum := Sum + Real_Extended (U(i, k)) * Real_Extended (V(j, k));
          end loop;
          SVD_Product(i, j) := Real (Sum);

       end loop;
       end loop;
     
       --  Product - A = 0?
       --  Reuse SVD_Product as the error matrix.
     
       for i in Starting_Row .. Final_Row loop
       for j in Starting_Col .. Final_Col loop
          SVD_Product(i, j) := Abs (A_true(i, j) - SVD_Product(i, j));
       end loop;
       end loop;
     
       --  Frobenius norm fractional error = ||Err_Matrix|| / ||A||
       Max_Error_F := 
         Frobenius_Norm (SVD_Product, Final_Row, Final_Col, Starting_Row, Starting_Col) /
         (Frobenius_Norm (A_true, Final_Row, Final_Col, Starting_Row, Starting_Col) 
          + Min_Allowed_Real);
   
       put(" Err in A - U*W*V' is ~ ||A - U*W*V'|| / ||A|| ="); 
       put(Max_Error_F);
       new_line(2); 
   
       <<Endies>> null;
   
   end loop; -- over Chosen_Matrix
   end loop; -- over No_Of_Repetitions

   put (x);
  
end golub_svd_tst_1;
