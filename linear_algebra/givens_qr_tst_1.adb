
-- Test QR decomposition real valued square matrices.

with Ada.Numerics.Generic_elementary_functions;
with Givens_QR;
with Test_Matrices;
With Text_IO; use Text_IO;

procedure givens_qr_tst_1 is

   type Real is digits 15;

   subtype Index is Integer range 1..137;

   -- in this test, matrix is a square-shaped matrix on  Index x Index.
   -- eg Hilbert's matrix is a square matrix with unique elements on the range
   -- Index'First .. Index'Last.  However, you  have the option to QR any rectangular
   -- sub-block of the matrix that is defined on Index x Index (provided 
   -- number of rows is >= number of cols).
   -- To do that you choose new values for Starting_Row, Starting_Col, Final_Row
   -- Final_Col just below.

   subtype Row_Index is Index;
   subtype Col_Index is Index;

   Starting_Row : constant Row_Index := Index'First + 0;
   Starting_Col : constant Col_Index := Index'First + 0;
   Final_Row : constant Row_Index := Index'Last- 0;
   Final_Col : constant Col_Index := Index'Last- 0;

   type Matrix is array(Row_Index, Col_Index) of Real;

   --pragma Convention (Fortran, Matrix); --No! This QR prefers Ada convention.

   package math is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use math;

   package QR is new Givens_QR
     (Real     => Real, 
      R_Index  => Index, 
      C_Index  => Index, 
      A_Matrix => Matrix);
   use QR;

   -- QR exports Row_Vector and Col_Vector

   package Make_Square_Matrix is new test_matrices (Real, Index, Matrix);
   use Make_Square_Matrix;

   package rio is new Float_IO(Real);
   use rio;

 --subtype Real_Extended  is Real;     -- general case, and for best speed
   type Real_Extended  is digits 18;   -- 18 ok on intel

   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;


   A, R : Matrix;
   Q : Q_Matrix;
   Max_Error : Real;
   Frobenius_QR_Err_0 : Real;
   Max_Error_qq, Max_Error_qr : Real;
   Frobenius_QQ_Err, Frobenius_QR_Err : Real;
   Scale   :  Col_Vector;
   Permute :  Permutation;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (A            : in     Matrix)
    --Final_Row    : in     Index; 
    --Final_Col    : in     Index;
    --Starting_Row : in     Index; 
    --Starting_Col : in     Index)
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

      Max_A_Val := Max_A_Val + Two ** (Real'Machine_Emin + 4);
      Scaling := One / Max_A_Val;

      Sum := Zero;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * A(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;

   ------------------------------------
   -- Get_Error_in_QR_Decomposition --
   ------------------------------------
   
   procedure Get_Error_in_QR_Decomposition 
     (A              : in     Matrix;
      R              : in     Matrix;
      Q              : in     Q_Matrix;
      Scale          : in     Col_Vector;
      Permute        : in     Permutation;
      Max_Error      :    out Real;
      Max_Error_F    :    out Real) is
    
      Err             : Real;
      Min_Real        : constant Real := 2.0 ** (Real'Machine_Emin + 2);
      Product_Vector, Col_of_R : Col_Vector := (others => Zero);
      Err_Matrix : Matrix;
   begin
      -- find error in A - Q*R
      -- The Columns of R have been permuted; unpermute before comparison of A with Q*R
      -- The Columns of R have been scaled. Before comparison of A with Q*R
      -- Must unscale each col of R by multiplying them with 1/Scale(Permute(Col).

      Max_Error := Zero;
     
      for Col in Starting_Col .. Final_Col loop

         Col_of_R := (others => Zero);

         for Row in Starting_Row .. Final_Row loop
            Col_of_R(Row) := R(Row, Col);
         end loop;

         Product_Vector := Q_x_Col_Vector (Q, Col_of_R);

         for Row in Starting_Row .. Final_Row loop
            Err := Abs (A(Row, Permute(Col)) - 
                      --Product_Vector(Row) / (Scale(Permute(Col)) + Min_Real));
                        Product_Vector(Row) / (Scale(Row) + Min_Real));
            if Err > Max_Error then
               Max_Error := Err;
            end if;
            Err_Matrix(Row, Col) := Err;
         end loop;
      end loop;

      --  Froebenius norm fractional error = ||Err_Matrix|| / ||M||

      Max_Error_F := Frobenius_Norm (Err_Matrix) /
                       (Frobenius_Norm (A) + Min_Real); 

   end Get_Error_in_QR_Decomposition;

   ------------------------------------
   -- Get_Err_in_Reassembled_Q_and_A --
   ------------------------------------
  
   --  Get an explicit matrix version of the matrix Q.  Call it V.
   
   procedure Get_Err_in_Reassembled_Q_and_A 
     (A              : in     Matrix;
      R              : in     Matrix;
      Q              : in     Q_Matrix;
      Scale          : in     Col_Vector;
      Permute        : in     Permutation;
      Final_Row      : in     Row_Index;
      Final_Col      : in     Col_Index;
      Starting_Row   : in     Row_Index;
      Starting_Col   : in     Col_Index;
      Frobenius_QQ_Err :  out Real;
      Frobenius_QR_Err :  out Real;
      Max_Error_QQ     :  out Real;
      Max_Error_QR     :  out Real)
   is
      Err, S          : Real;
      Min_Real : constant Real := +2.0 **(Real'Machine_Emin + 4);
    
      V, V_tr, Identity, Product_QQ : V_Matrix := (others => (others => Zero));
      Sum : Real_Extended;

      Product_A : Matrix;

      subtype Row_Index_Subrange is Row_Index range Starting_Row .. Final_Row;

   begin
  
      for r in Row_Index_Subrange loop
        Identity(r, r) := 1.0;
      end loop;

      -- Find error in I - V*V' etc.
      -- Start V as identity matrix (V should be square. MxM == R_Index x R_Index)
      -- Turn V    into and explicit version of Q    by calling Q_x_V_Matrix
      -- (Q is an array of 2x2 rotation matrices.)

      V := Identity;
      Q_x_V_Matrix (Q, V);
 
      -- Turn V_tr into and explicit version of Q_tr by calling Q_trans..x_V_Matrix

      V_tr := Identity;
      Q_transpose_x_V_Matrix (Q, V_tr);

      -- Usually find that orthonormality of *Rows of V* and *Cols of V_tr* is best.
      -- Notation: V' == V_tr == transpose of V.

      for Col in Row_Index_Subrange loop
      for Row in Row_Index_Subrange loop
        Sum := 0.0;
        for j in Row_Index_Subrange loop
         --Sum := Sum + Real_Extended(V(j, Row)) * Real_Extended(V(j, Col));  --V'*V 
           Sum := Sum + Real_Extended(V(Row, j)) * Real_Extended(V(Col, j));  --V*V' 
         --Sum := Sum + V(Row, j) * V(Col, j);     --V*V' has least err; rows of Q ortho
         --Sum := Sum + V_tr(j, Row)*V_tr(j, Col); --V_tr'*V_tr also least err; cols V_tr
         --Sum := Sum + V(Row, j) *  V_tr(j, Col); --least err; 
         --Sum := Sum + V_tr(Row, j)  * V(j, Col);
        end loop;
        Product_QQ(Row, Col) := Real (Sum);
      end loop;
      end loop;

      -- get   Product_QQ - Identity

      Max_Error_QQ := Zero;
     
      for Col in Row_Index_Subrange loop
      for Row in Row_Index_Subrange loop
         Err := Abs (Identity(Row, Col) - Product_QQ(Row, Col));
         if Err > Max_Error_QQ then
            Max_Error_QQ := Err;
         end if;
      end loop;
      end loop;
  
      -- Get Frobenius norm of:  Product_QQ - I:

      S := Zero;
      for Col in Row_Index_Subrange loop
      for Row in Row_Index_Subrange loop
         Err := Identity(Row, Col) - Product_QQ(Row, Col);
         S := S + Err * Err;
      end loop;
      end loop;

      Frobenius_QQ_Err := Sqrt(S) / Sqrt (-Real(Starting_Row) + Real(Final_Row) + 1.0);
   
     -- explicitly calculate Q*R by getting  V*R == Q*R (V==Q).

      for Col in Starting_Col .. Final_Col loop
      for Row in Starting_Row .. Final_Row loop
        Sum := 0.0;
        for j in Starting_Row .. Final_Row loop
           Sum := Sum + Real_Extended(V(Row, j)) * Real_Extended(R(j, Col));  --V*R
        end loop;
        Product_A(Row, Col) := Real (Sum);
      end loop;
      end loop;

     -- recall that the actual decomposition is:  A*Scale*Permute = Q*R

      Max_Error_QR := Zero;
     
      for Col in Starting_Col .. Final_Col loop
      for Row in Starting_Row .. Final_Row loop
         Err := Abs (A(Row, Permute(Col))
           -- - Product_A(Row, Col) / (Scale(Permute(Col)) + Min_Real)); -- Scale > 0
              - Product_A(Row, Col) / (Scale(Row) + Min_Real)); -- Scale > 0
         if Err > Max_Error_QR then
            Max_Error_QR := Err;
         end if;
      end loop;
      end loop;

      -- resuse array Product_A to get error matrix Error := Product_A - A:

      for Col in Starting_Col .. Final_Col loop
      for Row in Starting_Row .. Final_Row loop
         Product_A(Row, Col) := A(Row, Permute(Col)) 
           -- - Product_A(Row, Col) / (Scale(Permute(Col)) + Min_Real); -- Scale > 0
              - Product_A(Row, Col) / (Scale(Row) + Min_Real); -- Scale > 0
      end loop;
      end loop;

      Frobenius_QR_Err := Frobenius_Norm (Product_A) / 
                                   (Frobenius_Norm (A) + Min_Real);
   
   end Get_Err_in_Reassembled_Q_and_A;

   -----------
   -- Pause --
   -----------

   procedure Pause (s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12 : string := "") is
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
      if S8 /= "" then put_line (S8); end if;
      if S9 /= "" then put_line (S9); end if;
      if S10 /= "" then put_line (S10); end if;
      if S11 /= "" then put_line (S11); end if;
      if S12 /= "" then put_line (S12); end if;
      new_line;
      begin
          Put ("Enter a character to continue: ");
          Get_Immediate (Continue);
      exception
          when others => null;
      end;
   end pause;

begin

   Pause( 
     "Test 1: QR decomposition of matrix A. The QR decomposition of A is",
     "successful if two identities are satisfied: Q'*Q = I and Q*R = A. If 15",
     "digit Reals are used, then we expect the error in the calculation of A = Q*R to",
     "be about 1 part in 10**15. In other words ||Q*R - A|| / ||A|| should be",
     "around 10**(-15). Here ||*|| denotes the Frobenius Norm. Other matrix norms",
     "give slightly different answers, so its an order of magnitude estimate. The",
     "Q matrix is a list of 2x2 Givens rotations. Most operations (eg. equation",
     "solving) use this representation of Q for best accuracy.",
     " ",
     "The tests are repeated using an explicit matrix for Q. The error in the",
     "QR decomposition of A, i.e. ||Q*R - A|| / ||A|| is recalculated and printed",
     "below so that you can see the increase in error."
     );


   for Chosen_Matrix in Matrix_Id loop
 --for Chosen_Matrix in kahan .. kahan loop
 --for Chosen_Matrix in kahan_col_scaled_2 .. kahan_col_scaled_2 loop
 --for Chosen_Matrix in kahan .. kahan_col_scaled_2 loop
 --for Chosen_Matrix in kahan .. kahan_row_scaled loop

      Init_Matrix (A, Chosen_Matrix, Index'First, Index'Last);

      R := A; -- A remains original A. Only R is input.
 
      QR_Decompose 
        (A                  => R,   -- A has now been tranformed into the R matrix
         Q                  => Q,
         Row_Scalings       => Scale,
         Col_Permutation    => Permute,
         Final_Row          => Final_Row,
         Final_Col          => Final_Col,
         Starting_Row       => Starting_Row,
         Starting_Col       => Starting_Col);

      --declare Row : Row_Index := Starting_Row; begin
      --for Col in Starting_Col .. Final_Col loop
        --put(R(Col, Col));  -- this is the R matrix but it has not been unscaled yet.
        --if Row < Row_Index'Last then Row := Row + 1; end if;
      --end loop;
      --end;

      if true then
         new_line;
         put("For matrix A of type  "); put(Matrix_id'Image(Chosen_Matrix));  put(":");
         new_line; put ("Min 3 diag elements:");
       --for i in Starting_Col .. Final_Col loop
         for i in Final_Col-2 .. Final_Col loop
            put (r(i,i));
         end loop;
         new_line(1);
      end if;

      --goto endies;

      Get_Err_in_Reassembled_Q_and_A 
        (A, R, Q, Scale, Permute,
         Final_Row, Final_Col, Starting_Row, Starting_Col, 
         Frobenius_QQ_Err, Frobenius_QR_Err,
         Max_Error_qq, Max_Error_qr);

      Get_Error_in_QR_Decomposition 
        (A, R, Q, Scale, Permute,
         Max_Error, Frobenius_QR_Err_0);

      --  Froebenius norm fractional error:
      --      Max_Error_F  = ||Err_Matrix|| / ||A||

      new_line;
      put("For matrix A of type  "); put(Matrix_id'Image(Chosen_Matrix));  put(":");
      new_line;
      put(" Err in I-Q'*Q (Q = explicit matrix) is ||I-Q'*Q|| / ||I|| ="); 
      put(Frobenius_QQ_Err);
      new_line;
      put(" Err in A-Q*R  (Q = explicit matrix) is ||A-Q*R || / ||A|| ="); 
      put(Frobenius_QR_Err);
      new_line;
      put(" Err in A-Q*R  (Q = Givens rotation) is ||A-Q*R || / ||A|| ="); 
      put(Frobenius_QR_Err_0);
      new_line;

      <<endies>> null;

   end loop;

end givens_qr_tst_1;
