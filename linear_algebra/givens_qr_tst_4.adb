
-- Test QR decomposition on a real valued rectangular matrix.
-- Here we calculate Q and Q' in explicit matrix form, and
-- verify that Q'*Q = Q*Q' = I, and Q*R = A
-- (Error greater this way; previously just left Q as a list of
-- 2 x 2 Givens rotation matrices.)

with Givens_QR;
with Ada.Numerics.Generic_Elementary_Functions;
with Rectangular_Test_Matrices;
with Text_IO; use Text_IO;

procedure givens_qr_tst_4 is

   type Real is digits 15;

   type Row_Index is new Integer range 0..2**8-9;
   type Col_Index is new Integer range 0..2**7-10;
   -- must have fewer Columns than Rows, (or same number).

   type Matrix is array(Row_Index, Col_Index) of Real;

   package Math is new Ada.Numerics.Generic_Elementary_Functions(Real);
   use Math;

   -- Can QR an arbitrary rectangle of the matrix, 
   -- but must have:  No_of_Cols <= No_of_Rows
   -- These 4 numbers define the sub-block to be transformed:

   Starting_Row : constant Row_Index := Row_Index'First + 2;
   Starting_Col : constant Col_Index := Col_Index'First + 3;
   Final_Row    : constant Row_Index := Row_Index'Last-7;
   Final_Col    : constant Col_Index := Col_Index'Last-17;

   package QR is new Givens_QR
     (Real     => Real, 
      R_Index  => Row_Index, 
      C_Index  => Col_Index, 
      A_Matrix => Matrix);
   use QR;

   package Make_Matrix is new Rectangular_Test_Matrices (Real, Row_Index, Col_Index, Matrix);
   use Make_Matrix;

   package rio is new Float_IO(Real);
   use rio;

   type Real_Extended  is digits 15;   -- general case, works fine
 --type Real_Extended  is digits 18;   -- 18 ok on intel

   Zero : constant Real := +0.0;

   C : Matrix := (others => (others => 0.0));
   Q : Q_Matrix;
   Max_Error_qq, Max_Error_qr : Real;
   Frobenius_QQ_Err, Frobenius_QR_Err : Real;

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (A            : in Matrix)
      --Final_Row    : in Row_Index;
      --Final_Col    : in Col_Index;
      --Starting_Row : in Row_Index;
      --Starting_Col : in Col_Index)
      return Real
    is
      Max_A_Val : Real := 0.0;
      Sum, Scaling, tmp : Real := 0.0;
    begin
 
      Max_A_Val := 0.0;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         if Max_A_Val < Abs A(Row, Col) then Max_A_Val := Abs A(Row, Col); end if;
      end loop;
      end loop;

      Max_A_Val := Max_A_Val + 2.0 ** (Real'Machine_Emin + 4);
      Scaling := 1.0 / Max_A_Val;

      Sum := 0.0;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * A(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;

   --------------------------------
   -- Get_Error_in_reassembled_Q --
   --------------------------------
  
   --  Get an explicit matrix version of the matrix Q.  Call it V.
   
   procedure Get_Error_in_reassembled_Q 
     (A              : in     Matrix;
      Final_Row      : in     Row_Index;
      Final_Col      : in     Col_Index;
      Starting_Row   : in     Row_Index;
      Starting_Col   : in     Col_Index;
      Frobenius_QQ_Err :  out Real;
      Frobenius_QR_Err :  out Real;
      Max_Error_QQ     :  out Real;
      Max_Error_QR     :  out Real)
   is
      A_qr            : Matrix;
      Scale           : Col_Vector;
      Permute         : Permutation;
      Err, S          : Real;
      Min_Real : constant Real := +2.0 **(Real'Machine_Emin + 4);
    
      V, V_tr, Identity, Product_QQ : V_Matrix := (others => (others => 0.0));
      Sum : Real_Extended;

      Product_A : Matrix;

      subtype Row_Index_Subrange is Row_Index range Starting_Row .. Final_Row;

   begin
  
      for r in Row_Index_Subrange loop
        Identity(r, r) := 1.0;
      end loop;

      -- find error in I - V*V' etc

      A_qr := A;
 
      QR_Decompose 
        (A                  => A_qr,   -- A has now been tranformed into R == A_qr
         Q                  => Q,
         Row_Scalings       => Scale,
         Col_Permutation    => Permute,
         Final_Row          => Final_Row,
         Final_Col          => Final_Col,
         Starting_Row       => Starting_Row,
         Starting_Col       => Starting_Col);

      -- start V as identity matrix (V should be square. MxM == R_Index x R_Index)
      -- Q is an array of 2x2 rotation matrices.
      -- Turn V    into and explicit version of Q    by calling Q_x_V_Matrix

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
	   Sum := Sum + Real_Extended(V(j, Row)) * Real_Extended(V(j, Col));  --V'*V 
	 --Sum := Sum + V(Row, j) * V(Col, j);     --V*V' has least err; rows of Q ortho
         --Sum := Sum + V_tr(j, Row)*V_tr(j, Col); --V_tr'*V_tr also least err; cols V_tr
	 --Sum := Sum + V(Row, j) *  V_tr(j, Col); --least err; 
	 --Sum := Sum + V_tr(Row, j)  * V(j, Col);
        end loop;
	Product_QQ(Row, Col) := Real (Sum);
      end loop;
      end loop;

      -- get   Product_QQ - Identity

      Max_Error_QQ := 0.0;
     
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
   
     -- explicitly calculate Q*R by getting V*A_qr  (V==Q and A_qr==R).

      for Col in Starting_Col .. Final_Col loop
      for Row in Starting_Row .. Final_Row loop
        Sum := 0.0;
	for j in Starting_Row .. Final_Row loop
	   Sum := Sum + Real_Extended(V(Row, j)) * Real_Extended(A_qr(j, Col));  --V*R
	end loop;
	Product_A(Row, Col) := Real (Sum);
      end loop;
      end loop;

     -- recall that the actual decomposition is:  A*Scale*Permute = Q*R

      Max_Error_QR := 0.0;
     
      for Col in Starting_Col .. Final_Col loop
      for Row in Starting_Row .. Final_Row loop
         Err := Abs (A(Row, Permute(Col))
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
             - Product_A(Row, Col) / (Scale(Row) + Min_Real); -- Scale > 0
      end loop;
      end loop;

      Frobenius_QR_Err := Frobenius_Norm (Product_A) / 
                                   (Frobenius_Norm (A) + Min_Real);
   
   end Get_Error_in_reassembled_Q;

   -----------
   -- Pause --
   -----------

   procedure Pause (s1,s2,s3,s4 : string := "") is
     Continue : Character := ' ';
   begin
     New_Line;
     if S1 /= "" then put_line (S1); end if;
     if S2 /= "" then put_line (S2); end if;
     if S3 /= "" then put_line (S3); end if;
     if S4 /= "" then put_line (S4); end if;
     new_line;
     begin
	Put ("Enter a character to continue: ");
	Get_Immediate (Continue);
     exception
	when others => null;
     end;
   end pause;

begin

   Pause ("Tests QR decomposition on a real valued rectangular matrix.",
          "Here we calculate Q and Q' in explicit matrix form, and then",
          "verify that Q'*Q = I, and Q*R = A.");

   for Chosen_Matrix in Matrix_id loop

      Init_Matrix 
        (C, Chosen_Matrix, Final_Row, Final_Col, Starting_Row, Starting_Col);

      Get_Error_in_reassembled_Q 
       (C, Final_Row, Final_Col, Starting_Row, Starting_Col, 
        Frobenius_QQ_Err, Frobenius_QR_Err,
        Max_Error_qq, Max_Error_qr);

      new_line;
      put("For matrix of type  "); put(Matrix_id'Image(Chosen_Matrix));  put(":");
      new_line;
      put(" Err in A - Q*R  is ~ ||A - Q*R || / ||A|| = "); 
      put(Frobenius_QR_Err);
      new_line;
      put(" Err in I - Q'*Q is ~ ||I - Q'*Q|| / ||I|| = "); 
      put(Frobenius_QQ_Err);
      new_line;

   end loop;

end;
