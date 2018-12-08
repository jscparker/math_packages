
with Golub_SVD;
with Rectangular_Test_Matrices;
with Ada.Text_IO; 
use  Ada.Text_IO;

-- Demonstrates use of SVD when  No_of_Rows > No_of_Cols.
-- So you have more equations than unknowns in equation solving.
-- (Can't have No_of_Rows < No_of_Cols.)

procedure golub_svd_tst_2 is

   type Real is digits 15;
 
   Start   : constant := 1;
   Limit_r : constant := 222;
   Limit_c : constant := 100;
 
   --  You can make them different types if you want:
 
   type Row_Index is new Integer range Start .. Limit_r;
   type Col_Index is new Integer range Start .. Limit_c;
   
   --subtype Row_Index is Integer range Start .. Limit_r;
   --subtype Col_Index is Integer range Start .. Limit_c;
   
   type A_Matrix is array (Row_Index, Col_Index) of Real;
   pragma Convention (Fortran, A_Matrix);
 
   package lin_svd is new golub_svd (Real, Row_Index, Col_Index, A_Matrix);
   use lin_svd;
 
   package Rect_Matrix is new Rectangular_Test_Matrices(Real, Row_Index, Col_Index, A_Matrix);
   use Rect_Matrix;
 
   package rio is new Ada.Text_IO.Float_IO(Real); use rio;
 
 
   Final_Row    : constant Row_Index := Row_Index'Last  - 3;
   Final_Col    : constant Col_Index := Col_Index'Last  - 3;
   Starting_Col : constant Col_Index := Col_Index'First + 2;
   Starting_Row : constant Row_Index := Row_Index (Starting_Col);-- requirement of SVD
 
   A0, A, SVD_Product : A_Matrix  := (others => (others => 0.0));
   V : V_Matrix;
   VV_Product : V_Matrix;
   U : U_Matrix;
   UU_Product : U_matrix;
 
   Singular_Vals : Singular_Vector := (others => 0.0);
 
   Col_id_of_1st_Good_S_Val : Col_Index;
 
   Max_Singular_Val, Sum, Error_Off_Diag, Error_Diag, Del : Real;
 
begin

   new_line(1);
   new_line; put ("Error in the singular value decomposition is estimated by finding");
   new_line; put ("The max element of  A - U*W*V', and dividing it  by the largest");
   new_line; put ("singular value of A. Should be somewhere around Real'Epsilon"); 
   new_line; put ("in magnitude if all goes well.");
   new_line(1);
   new_line; put ("Error in the calculation of V is estimated by finding the max element");
   new_line; put ("of I - V'*V = I - Transpose(V)*V");
   new_line; put ("I - V'*V should be near Real'Epsilon.");
   new_line(1);
   new_line; put ("Error in the calculation of U is estimated by finding the max element");
   new_line; put ("of I - U'*U = I - Transpose(U)*U");
   new_line; put ("I - U'*U should be near Real'Epsilon.");
   new_line(1);
   new_line; put ("Notice that U'U = I implies A'A = VW'U'UWV' = VEV' where E is a");
   new_line; put ("diagonal matrix (E = W'W) containing the squares of the singular");
   new_line; put ("values. So U'U = I and V'V = I imply that the cols of V are the "); 
   new_line; put ("eigvectors of A'A. (Multiply the above equality by V to get A'AV = VE.)");
   new_line(1);

   for Chosen_Matrix in Matrix_id loop
  
      Init_Matrix (A, Chosen_Matrix);
    
      -- Get A = U * W * V'  where W is diagonal containing the singular vals
    
      for i in 1..1 loop
         A0 := A;
         SVD_Decompose 
           (A                => A0, 
            U                => U,
            V                => V,
            S                => Singular_Vals,
            Id_of_1st_S_Val  => Col_id_of_1st_Good_S_Val,
            Starting_Col     => Starting_Col,
            Final_Col        => Final_Col,
            Final_Row        => Final_Row, 
            Matrix_U_Desired => True,
            Matrix_V_Desired => True);
      end loop;
    
      new_line(2); 
      put ("Testing SVD on matrix A of type: "); 
      put (Matrix_id'Image(Chosen_Matrix));
  
      if Col_id_of_1st_good_S_Val /= Starting_Col then
         new_line; 
         put ("QR failed to fully converge at column = ");
         put(Col_Index'Image(Col_id_of_1st_good_S_Val));
         new_line; 
      else
         new_line; 
         put ("Max Singular Val = ");
         put(Real'Image(Singular_Vals(Col_id_of_1st_good_S_Val)));
      end if;
  
      Max_Singular_Val := Singular_Vals(Col_id_of_1st_good_S_Val);
      Max_Singular_Val := Abs Max_Singular_Val + 2.0**(Real'Machine_Emin + 4*Real'Digits);
  
      --new_line;
      --for i in Singular_Vals'Range loop
        --put (Singular_Vals(i));
      --end loop;
      --new_line(1);
      --new_line; put ("Singular_Vals printed above.");
    
    
      --  Get Product = V' * V = transpose(V) * V
    
      for i in Col_Index range Starting_Col .. Final_Col loop
      for j in Col_Index range Starting_Col .. Final_Col loop
         Sum := 0.0;
         for k in Col_Index range Starting_Col .. Final_Col loop
            Sum := Sum + V(k, i) * V(k, j);
         end loop;
         VV_Product(i,j) := Sum;
      end loop;
      end loop;
    
      --  Product - I = 0?
    
      Error_Diag     := 0.0;
      Error_Off_Diag := 0.0;
    
      for i in Col_Index range Starting_Col .. Final_Col loop
      for j in Col_Index range Starting_Col .. Final_Col loop
         if i = j then
            Del := Abs (1.0 - VV_Product(i,j));
            if Del > Error_Diag then
               Error_Diag := Del;
            end if;
         else
            Del := Abs (0.0 - VV_Product(i,j));
            if Del > Error_Off_Diag then
               Error_Off_Diag := Del;
            end if;
         end if;
      end loop;
      end loop;
    
      new_line; put ("Max err in I - V'*V,   On_Diagonal  ="); put (Error_Diag);
      new_line; put ("Max err in I - V'*V,   Off_Diagonal ="); put (Error_Off_Diag);
    
      --  Get Product = U * W * V' = U * W * transpose(V)
      --  U is  usually Row_Index x Row_Index but can probably use Col_Index for its 2nd
      --  V is  Col_Index x Col_Index
      --
      --  S is  conceptually Row_Index x Col_Index, with the Sing_Vals on its diagonal, 
      --  but its really in a vector on Col_Index.
      --
      --  So the No_of_Singular vals is MIN of No_of_Rows, No_Of_Cols
      --  Singular vals are sorted so that largest are at beginning of Singular_Vals(k).
      --
      --  Get matrix product S * V'; place in matrix V':
  
      for k in Col_Index range Starting_Col .. Final_Col loop
      for c in Col_Index range Starting_Col .. Final_Col loop
         V(c, k) := V(c, k) * Singular_Vals(k);
      end loop;
      end loop;
  
      --  V' = S*V'  is now conceptually Row_Index x Col_Index, but really it
      --  is all zeros for Rows > Final_Col:
  
      for r in Starting_Row .. Final_Row loop  -- U  is     Row_Index x Row_Index
      for c in Starting_Col .. Final_Col loop  -- V' is     Row_Index x Col_Index
  
         Sum := 0.0;
         for k in Col_Index range Starting_Col .. Final_Col loop
            Sum := Sum + U(r, Row_Index (k)) * V(c, k);
         end loop;
         SVD_Product(r, c) := Sum;
  
      end loop;
      end loop;
    
    
      Error_Diag     := 0.0;
      Error_Off_Diag := 0.0;
    
      for i in Starting_Row .. Final_Row loop
      for j in Starting_Col .. Final_Col loop
         if i = Row_Index'Base (j) then
            Del := Abs (A(i,j) - SVD_Product(i,j));
            if Del > Error_Diag then
               Error_Diag := Del;
            end if;
         else
            Del := Abs (A(i,j) - SVD_Product(i,j));
            if Del > Error_Off_Diag then
               Error_Off_Diag := Del;
            end if;
         end if;
      end loop;
      end loop;
    
      new_line(1);
      put ("Max err in A - U*W*V', On_Diagonal  ="); put (Error_Diag / Max_Singular_Val);
      new_line; 
      put ("Max err in A - U*W*V', Off_Diagonal ="); put (Error_Off_Diag / Max_Singular_Val);
    
      -- Get Product = U' * U =  transpose(U) * U
    
      for i in Row_Index range Starting_Row .. Final_Row loop
      for j in Row_Index range Starting_Row .. Final_Row loop
         Sum := 0.0;
         for k in Starting_Row .. Final_Row loop
           Sum := Sum + U(k, i) * U(k, j);
         end loop;
         UU_Product(i,j) := Sum;
      end loop;
      end loop;
    
      -- Product - I = 0?
    
      Error_Diag     := 0.0;
      Error_Off_Diag := 0.0;
    
      for i in Row_Index range Starting_Row .. Final_Row loop
      for j in Row_Index range Starting_Row .. Final_Row loop
         if i = j then
            Del := Abs (1.0 - UU_Product(i,j));
            if Del > Error_Diag then
               Error_Diag := Del;
            end if;
         else
            Del := Abs (0.0 - UU_Product(i,j));
            if Del > Error_Off_Diag then
               Error_Off_Diag := Del;
            end if;
         end if;
      end loop;
      end loop;
    
      new_line; put ("Max err in I - U'*U,   On_Diagonal  ="); put (Error_Diag);
      new_line; put ("Max err in I - U'*U,   Off_Diagonal ="); put (Error_Off_Diag);
    
   end loop; -- over Chosen_Matrix
  
end golub_svd_tst_2;
