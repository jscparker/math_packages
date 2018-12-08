
--------------------------------------------------------------------------
-- package body Givens_QR, QR decomposition using Givens rotations
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

with Ada.Numerics;
with Ada.Numerics.Generic_Elementary_Functions;

package body Givens_QR is

  package math is new Ada.Numerics.Generic_Elementary_Functions (Real);  use math;

  Zero : constant Real := +0.0;
  One  : constant Real := +1.0;
  Two  : constant Real := +2.0;

  Exp_Min          : constant Integer := Real'Machine_Emin;
  Min_Allowed_Real : constant Real    := Two**(Exp_Min - Exp_Min / 8);

  type Int64 is range -2**63+1 .. 2**63-1;

  -- Make these True:

  Scaling_Desired  : constant Boolean := True;
  Pivoting_Desired : constant Boolean := True;
  --  Scaling Prevents overflow during pivoting. With pivoting and
  --  scaling enabled, the QR usually reveals the rank of matrix A.

  -------------------------------
  -- Inverse_Rotate_Col_Vector --
  -------------------------------
 
  -- for getting Q*B 
 
  procedure Inverse_Rotate_Col_Vector  
    (B              : in out Col_Vector;
     lower_right_id : in     R_Index;  -- low diagonal element (cos) of rot. matrix
     upper_left_id  : in     R_Index;  -- top diagonal element (cos) of rot. matrix
     Rot            : in     Rotation;
     P_bigger_than_L : in Boolean) 
  is
     B_pvt, B_Low : Real;
     c : real renames Rot.Cosine;
     s : real renames Rot.Sine;
  begin
 
     -- all you do is flip the sign of Sine.  s in part 1, 
     -- and in part 2, 1 + (s-1) becomes -1 - (s-1), where the s-1 is given as s.
 
     if P_bigger_than_L then
        if abs s > Zero then -- optimization uses fact that c is always positive.
           B_pvt := B(upper_left_id);
           B_Low := B(lower_right_id);
           B(upper_left_id)  := B_pvt - s*(-c*B_pvt +   B_Low);
           B(lower_right_id) := B_Low - s*(  -B_pvt - c*B_Low);
        end if;
     else
        B_pvt := B(upper_left_id);
        B_Low := B(lower_right_id);
        B(upper_left_id)  :=-B_Low + c*(   B_Pvt - s*B_Low);
        B(lower_right_id) := B_Pvt + c*(+s*B_Pvt +   B_Low);
     end if;
 
  end Inverse_Rotate_Col_Vector;
 
  pragma Inline (Inverse_Rotate_Col_Vector);
 
  -----------------------
  -- Rotate_Col_Vector --
  -----------------------
 
  -- for getting Q'*B . Uses actual rotation stored in Q, which is the one
  -- used to create R .. (so Q really contains Q_transpose.)
 
  procedure Rotate_Col_Vector  
    (B              : in out Col_Vector;
     lower_right_id : in     R_Index;  -- low diagonal element (cos) of rot. matrix
     upper_left_id  : in     R_Index;  -- top diagonal element (cos) of rot. matrix
     Rot            : in     Rotation;
     P_bigger_than_L : in Boolean) 
  is
     B_pvt, B_Low : Real;
     c : real renames Rot.Cosine;
     s : real renames Rot.Sine;
  begin
     if P_bigger_than_L then
        if abs s > Zero then -- rot not identity: use fact c is always positive.
           B_pvt := B(upper_left_id);
           B_Low := B(lower_right_id);
           B(upper_left_id)  := B_pvt + s*( c*B_pvt +   B_Low);
           B(lower_right_id) := B_Low + s*(  -B_pvt + c*B_Low);
        end if;
     else
        B_pvt := B(upper_left_id);
        B_Low := B(lower_right_id);
        B(upper_left_id)  := B_Low + c*(   B_Pvt + s*B_Low);
        B(lower_right_id) :=-B_Pvt + c*(-s*B_Pvt +   B_Low);
     end if;
 
  end Rotate_Col_Vector;
 
  pragma Inline (Rotate_Col_Vector);
 
  -----------------------------------
  -- Find_Id_Of_Max_Element_Of_Row --
  -----------------------------------

  procedure Get_Id_Of_Max_Element_Of_Row
    (The_Row            : in  Row_Vector;
     Starting_Col       : in  C_Index;
     Final_Col          : in  C_Index;
     Id_Of_Max_Element  : out C_Index;
     Val_of_Max_Element : out Real)
  is
     Val : Real;
  begin
     Id_Of_Max_Element  := Starting_Col;
     Val_of_Max_Element := Abs The_Row(Starting_Col);

     if Final_Col > Starting_Col then
     for k in Starting_Col+1 .. Final_Col loop
        Val := The_Row(k);
        if Abs Val > Val_of_Max_Element then
           Val_of_Max_Element     := Val;
           Id_Of_Max_Element      := k;
        end if;
     end loop;
     end if;
  end Get_Id_Of_Max_Element_Of_Row;

  ------------------
  -- QR_Decompose --
  ------------------
 
  procedure QR_Decompose 
    (A                  : in out A_Matrix;
     Q                  :    out Q_Matrix;
     Row_Scalings       :    out Col_Vector;
     Col_Permutation    :    out Permutation;
     Final_Row          : in     R_Index := R_Index'Last;
     Final_Col          : in     C_Index := C_Index'Last;
     Starting_Row       : in     R_Index := R_Index'First;
     Starting_Col       : in     C_Index := C_Index'First)
  is 
     No_of_Rows : constant Real 
                 := Real (Final_Row) - Real (Starting_Row) + One;
     No_of_Cols : constant Real 
                 := Real (Final_Col) - Real (Starting_Col) + One;
 
     Pivot_Row : R_Index := Starting_Row; -- essential init.
 
     Norm_Squ_of_Col_of_A : Row_Vector;
 
     ------------------------------------------------
     -- Rotate_Rows_to_Zero_Out_Element_Lo_of_pCol --
     ------------------------------------------------
 
     -- cos = P/r,   sin = L/r,   r = sqrt(P*P + L*L)
     --
     -- (P is for Pivot, L is for Low.
     --
     -- clockwise rotation: notice all rotations are centered on the diagonal
     --
     --   1  0  0  0  0       0       0
     --   0  c  s  0  0       P       r
     --   0 -s  c  0  0   x   L   =   0
     --   0  0  0  1  0       0       0
     --   0  0  0  0  1       0       0
     --
     --
     -- if |L| >= |P| then   tn = P/L  <= 1
     -- if |P| >  |L| then   tn = L/P  <  1
     --
     --
     --  let t = smaller / larger
     --
     --  let u = 1/sqrt(1+t*t) and w = sqrt(t*t/(1+t*t)) with
     --
     --  use (1 - 1/sqrt(1+t*t)) * (1 + 1/sqrt(1+t*t)) = t*t / (1 + t*t)
     --
     --  1/sqrt(1+t*t) - 1 = - t*t/(sqrt(1+t*t) + 1+t*t) 
     --                    = -(t/(sqrt(1+t*t))*(t/(1 + Sqrt(1+t*t)))
     --                    = - Abs (w) * Abs (t)/(1+ sqrt(1+t*t))
     --                    =   u_lo
     --
     --  u_hi = 1   =>   u_lo + u_hi  =  1/sqrt(1+t*t) = u = Abs (cs) if a<b
     --
     --
     --  hypot = |L| * sqrt(1+t*t) = |L| * (1  + t*t/(1+sqrt(1+t*t)) = hi + lo if t<1
     --
     --
 
     procedure Rotate_Rows_to_Zero_Out_Element_Lo_of_pCol
       (pCol   : in C_Index; 
        Hi_Row : in R_Index;
        Lo_Row : in R_Index)
     is
        Pivot_Col : C_Index renames pCol;
        Pivot_Row : R_Index renames Hi_Row;
        Low_Row   : R_Index renames Lo_Row;
 
        tn, sq : real;
        sn, cn : Real;
        A_pvt, A_low : real;
        cn_minus_1_over_sn, sn_minus_1_over_cn : real;
        P : Real := A(Pivot_Row, Pivot_Col); -- P is for Pivot
        L : Real := A(Low_Row, Pivot_Col);
        Abs_P : constant Real := Abs (P);
        Abs_L : constant Real := Abs (L);
        P_bigger_than_L : Boolean := True;
 
     begin
 
       -- use Identity if no rotation is performed:
 
       Q.Rot(Low_Row, Pivot_Col).Cosine := Zero;
       Q.Rot(Low_Row, Pivot_Col).Sine   := Zero;
       Q.P_bigger_than_L(Low_Row, Pivot_Col)  := True;
 
       if (not L'Valid) or (not P'Valid) then
         return; -- having updated Q w. identity
       end if;
 
       if Abs_P > Abs_L then   -- |s| < |c|
 
         if Abs_P > Zero then
 
            P_bigger_than_L := True;
 
            -- Make  P>0 always, so if P<0 then flip signs of both P and L
            -- to ensure zero'd out element.
            if P < Zero then  
               L := -L; P := Abs_P;
            end if;
            tn                 := L / P;                --  P>0 so t has sign of L
            sq                 := Sqrt (One + tn*tn);
            sn                 := tn / sq;              --  sn has sign of L
            cn_minus_1_over_sn :=-tn / (One + sq); 
            -- factored out the sn from cn_minus_1
 
            for Col in Pivot_Col .. Final_Col  loop 
               A_pvt := A(Pivot_Row, Col);
               A_low := A(Low_Row, Col);
               A(Pivot_Row, Col) := A_pvt + sn * (cn_minus_1_over_sn*A_pvt +  A_low);
               A(Low_Row, Col)   := A_low + sn * (-A_pvt + cn_minus_1_over_sn*A_low);
            end loop;
 
            Q.Rot(Low_Row, Pivot_Col).Cosine       := cn_minus_1_over_sn;
            Q.Rot(Low_Row, Pivot_Col).Sine         := sn;
            Q.P_bigger_than_L(Low_Row, Pivot_Col)  := P_bigger_than_L;
 
         end if;
 
       else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1
 
         if Abs_L > Zero then
 
            P_bigger_than_L := False;
 
            -- Make L>0 always, so if L<0 then flip signs of both P and L
            -- to ensure zero'd out element.
            if L < Zero then  
               L := Abs_L; P := -P;
            end if;
  
            tn                 := P / L;               -- L>0 now so, t has sign of P.
            sq                 := Sqrt (One + tn*tn);
            cn                 := tn / sq;             -- has sign of P
            sn_minus_1_over_cn :=-tn / (One + sq);
 
            -- factored out the cn from sn_minus_1
 
            for Col in Pivot_Col .. Final_Col  loop 
               A_pvt := A(Pivot_Row, Col);
               A_low := A(Low_Row, Col);
               A(Pivot_Row, Col) := A_low + cn * (A_pvt +  sn_minus_1_over_cn*A_low);
               A(Low_Row, Col)   :=-A_pvt + cn * (-sn_minus_1_over_cn*A_pvt + A_low);
            end loop;
 
            Q.Rot(Low_Row, Pivot_Col).Cosine       := cn;
            Q.Rot(Low_Row, Pivot_Col).Sine         := sn_minus_1_over_cn;
            Q.P_bigger_than_L(Low_Row, Pivot_Col)  := P_bigger_than_L;
 
          end if; -- Abs_L > Zero
 
       end if;
 
     end Rotate_Rows_to_Zero_Out_Element_Lo_of_pCol;
 
     procedure Scale_Rows
       (Row_Scalings : out Col_Vector;
        Starting_Col : in C_Index := C_Index'First;
        Final_Col    : in C_Index := C_Index'Last;
        Starting_Row : in R_Index := R_Index'First;
        Final_Row    : in R_Index := R_Index'Last)
     is
        Norm : Real;
        Scale_Factor, Abs_A : Real;
     begin

        Row_Scalings := (others => One); -- important init
 
        for Row in Starting_Row .. Final_Row loop
    
           --Norm := Zero;
           --Get_1_Norm:
           --for Col in Starting_Col .. Final_Col loop 
             --if not A(Row, Col)'valid then raise constraint_error with "7"; end if;
             --Norm := Norm + abs A(Row, Col);
           --end loop Get_1_Norm;
          
           Norm := Zero;
           Get_Infinity_Norm:
           for Col in Starting_Col .. Final_Col loop 
              Abs_A := Abs A(Row, Col);
              if  Abs_A > Norm  then   Norm := Abs_A;   end if;
           end loop Get_Infinity_Norm;
    
           if not Norm'Valid then
              raise Ada.Numerics.Argument_Error with "Elements of input matrix invalid.";
           end if;
     
           if Norm < Two ** (Real'Machine_Emin - Real'Machine_Emin/16) then
              Scale_Factor := Zero;
           else
              declare 
                 Norm_Exp : constant Integer := Real'Exponent (Norm);
              begin
                 Scale_Factor := Two ** (-Norm_Exp);
              end;
           end if;
 
           Row_Scalings(Row) := Scale_Factor;
     
        end loop;
     
        for Row in Starting_Row .. Final_Row loop 
        for Col in Starting_Col .. Final_Col loop 
           A(Row, Col) := Row_Scalings(Row) * A(Row, Col);
        end loop;
        end loop;
 
     end Scale_Rows;

     -- Find sub Col_Vector with max norm, swap Columns if necessary. 

     procedure Do_Column_Pivoting
       (Starting_Col     : in     C_Index;
        Final_Col        : in     C_Index;
        Starting_Row     : in     R_Index;
        Final_Row        : in     R_Index;
        Pivot_Col        : in     C_Index;
        Col_Permutation  : in out Permutation;
        Col_Norm_Squared : in out Row_Vector)
     is
        Pivot_Row : constant R_Index := R_Index (Pivot_Col);
        tmp_index, New_Pivot_Col : C_Index;
        Pivot_Val, tmp : Real;
     begin

        if Pivot_Col >= Final_Col then return; end if;

        if (Pivot_Col = Starting_Col) or else
           ((Int64(Pivot_Col) - Int64(Starting_Col)) mod 6 = 0) then

           --  Without row scaling, we get overflows if matrix is too large.
           --  Get the right answer accurately every N steps:
 
           Col_Norm_Squared := (others => Zero);
           for r in Pivot_Row .. Final_Row loop 
           for c in Pivot_Col .. Final_Col loop
              Col_Norm_Squared(c) := Col_Norm_Squared(c) + A(r, c)**2;
           end loop;
           end loop;  -- fastest with Ada array convention rather than fortran.
      
        else

           --  Optimization: length of Col is invariant under rotations, so just
           --  subtract away the unwanted element in Col from the old sum for Norm_Squ:
           --  Norm_Squared (Col) := Norm_Squared (Col) - A(Pivot_Row-1, Col)**2;
           --  Has bad accuracy.

           for c in Pivot_Col .. Final_Col loop 
              Col_Norm_Squared(c) := Col_Norm_Squared(c) - A(Pivot_Row-1, c)**2;
           end loop;

        end if;

        Get_Id_Of_Max_Element_Of_Row
          (The_Row            => Col_Norm_Squared,
           Starting_Col       => Pivot_Col,
           Final_Col          => Final_Col,
           Id_Of_Max_Element  => New_Pivot_Col,
           Val_of_Max_Element => Pivot_Val);
 
        if New_Pivot_Col > Pivot_Col then

           for r in Starting_Row .. Final_Row loop
              tmp                 := A(r, New_Pivot_Col);
              A(r, New_Pivot_Col) := A(r, Pivot_Col);
              A(r, Pivot_Col)     := tmp;
           end loop;
 
           tmp                             := Col_Norm_Squared (New_Pivot_Col);
           Col_Norm_Squared(New_Pivot_Col) := Col_Norm_Squared (Pivot_Col);
           Col_Norm_Squared(Pivot_Col)     := tmp;
 
           tmp_index                      := Col_Permutation(New_Pivot_Col);
           Col_Permutation(New_Pivot_Col) := Col_Permutation(Pivot_Col);
           Col_Permutation(Pivot_Col)     := tmp_index;

        end if;
 
     end Do_Column_Pivoting;

  begin
 
   -- Important Inits.
 
   Q.Rot             := (others =>  (others => Identity));
   Q.P_bigger_than_L := (others =>  (others => True));
   Q.Final_Row       := Final_Row;
   Q.Final_Col       := Final_Col;
   Q.Starting_Row    := Starting_Row;
   Q.Starting_Col    := Starting_Col;
 
   Row_Scalings := (others => One); -- important init
 
   for j in C_Index loop
      Col_Permutation(j) := j;
   end loop;
 
   if No_Of_Cols > No_of_Rows then
      raise Ada.Numerics.Argument_Error with "Can't have more columns than rows.";
   end if;
   if No_Of_Cols < 2.0 then
      raise Ada.Numerics.Argument_Error with "Need at least 2 columns in matrix.";
   end if;
 
   -- Step 1.  Column Scaling.  Scaling prevents overflow when 2-norms are calculated.
 
   if Scaling_Desired then
     Scale_Rows
       (Row_Scalings => Row_Scalings,
        Starting_Col => Starting_Col,
        Final_Col    => Final_Col,
        Starting_Row => Starting_Row,
        Final_Row    => Final_Row);
   end if;
 
   --  Start the QR decomposition
 
   Pivot_Row := Starting_Row; -- essential init.
 
   for Pivot_Col in Starting_Col .. Final_Col loop
 
      if Pivoting_Desired then
         Do_Column_Pivoting
           (Starting_Col     => Starting_Col,
            Final_Col        => Final_Col,
            Starting_Row     => Starting_Row,
            Final_Row        => Final_Row,
            Pivot_Col        => Pivot_Col,
            Col_Permutation  => Col_Permutation,
            Col_Norm_Squared => Norm_Squ_of_Col_of_A);
      end if; -- pivoting_desired
 
      if Pivot_Row < Final_Row then
 
      for Low_Row in Pivot_Row+1 .. Final_Row loop
 
         Rotate_Rows_to_Zero_Out_Element_Lo_of_pCol
           (pCol   => Pivot_Col, -- Element to zero out is A(Lo_Row, pCol)
            Hi_Row => Pivot_Row, -- High to the eye; its id is lower than Lo_Row's.
            Lo_Row => Low_Row); 
 
         A(Low_Row, Pivot_Col)  := Zero;
 
      end loop;  --  over Low_Row
 
      end if; -- Pivot_Row < Final_Row
 
 
      if Pivot_Row < Final_Row then  Pivot_Row := Pivot_Row + 1; end if; 
      --  We are always finished if Pivot_Row = Final_Row because have M rows >= N cols
 
    end loop; -- over Pivot_Col
 
    -- To unscale the R matrix (A) you scale the Cols of A=R w/ Inv_Row_Scalings(Col).
   -- Bad idea because the scaled R=A has diag elements in monotonic decreasing order.
 
  end QR_Decompose;
 
  --------------------
  -- Q_x_Col_Vector --
  --------------------
 
  -- Get Product = Q*B
  --
  -- Apply the inverses of all the 2x2 rotation matrices in the order opposite to the
  -- order in which they were created: Start w/ most recent and wrk backwards in time.
 
  function Q_x_Col_Vector
    (Q : in Q_Matrix;
     B : in Col_Vector)
     return Col_Vector
  is
     Product : Col_Vector;
 
     Final_Row    : R_Index renames Q.Final_Row;
     Final_Col    : C_Index renames Q.Final_Col;
     Starting_Row : R_Index renames Q.Starting_Row;
     Starting_Col : C_Index renames Q.Starting_Col;
 
     Pivot_Row : R_Index := Starting_Row;
     P_is_the_Bigger : Boolean;
  begin
 
     -- The loop bounds reflect the design of Q. Q is stored in
     -- matrix A, which was defined to be M x N via Starting_Row,
     -- Starting_Col, etc.
 
     Product := B;
     Pivot_Row := R_Index (Int64(Final_Col) + Int64(Starting_Row) - Int64(Starting_Col));
 
     for Pivot_Col in reverse Starting_Col .. Final_Col loop 
        if Pivot_Row < Final_Row then
           for Low_Row in reverse Pivot_Row+1 .. Final_Row loop
              P_is_the_bigger := Q.P_bigger_than_L (Low_Row, Pivot_Col); 
              Inverse_Rotate_Col_Vector
               (Product, Low_Row, Pivot_Row, Q.Rot(Low_Row, Pivot_Col), P_is_the_bigger);
           end loop;
        end if;
        if Pivot_Row > Starting_Row then   Pivot_Row := Pivot_Row - 1;   end if;
     end loop;
 
     return Product;
 
  end Q_x_Col_Vector;
 
  ------------------------------
  -- Q_transpose_x_Col_Vector --
  ------------------------------
 
  -- Get Product = Q'*B
  --
  -- Apply the the 2x2 rotation matrices in the order in which they were created.
  -- Start w/ oldest and wrk forwards in time.
 
  function Q_transpose_x_Col_Vector
    (Q : in Q_Matrix;
     B : in Col_Vector)
     return Col_Vector
  is
     Product : Col_Vector;
 
     Final_Row    : R_Index renames Q.Final_Row;
     Final_Col    : C_Index renames Q.Final_Col;
     Starting_Row : R_Index renames Q.Starting_Row;
     Starting_Col : C_Index renames Q.Starting_Col;
 
     Pivot_Row : R_Index := Starting_Row;
     P_is_the_bigger : Boolean;
 
  begin
 
     -- The loop bounds reflect the design of Q. Q is stored in
     -- matrix A, which was defined to be M x N via Starting_Row,
     -- Starting_Col, etc.
 
     Product := B;
     Pivot_Row := Starting_Row;
 
     for Pivot_Col in Starting_Col .. Final_Col loop 
        if Pivot_Row < Final_Row then
           for Low_Row in Pivot_Row+1 .. Final_Row loop
              P_is_the_bigger := Q.P_bigger_than_L (Low_Row, Pivot_Col); 
              Rotate_Col_Vector (Product, Low_Row, Pivot_Row,
                        Q.Rot(Low_Row, Pivot_Col), P_is_the_bigger);
           end loop;
           Pivot_Row := Pivot_Row + 1;
        end if;
     end loop;
 
     return Product;
 
  end Q_transpose_x_Col_Vector;
 
  --------------
  -- QR_Solve --
  --------------
 
  -- A*X = B.  
  -- S*A*P = Q*R where S is diagonal matrix of row scalings; P permutes cols.
  -- So X = P*R^(-1)*Q'*S*B where Q is orthogonal, so Q' = Q^(-1).
 
  procedure QR_Solve 
    (X                  :    out Row_Vector;
     B                  : in     Col_Vector;
     R                  : in     A_Matrix; -- We overwrote A with R, so input A here.
     Q                  : in     Q_Matrix; 
     Row_Scalings       : in     Col_Vector;
     Col_Permutation    : in     Permutation;
     Singularity_Cutoff : in     Real := Default_Singularity_Cutoff)
  is
     Final_Row    : R_Index renames Q.Final_Row;
     Final_Col    : C_Index renames Q.Final_Col;
     Starting_Row : R_Index renames Q.Starting_Row;
     Starting_Col : C_Index renames Q.Starting_Col;
 
     Max_R_Diagonal_Val : constant Real := Abs R(Starting_Row, Starting_Col);
     R_Diagonal_Inv     : Row_Vector;
 
     Sum : Real;
     B1, B2  : Col_Vector;
     B0, X0 : Row_Vector;
     Row : R_Index := Starting_Row;
  begin
 
     -- X = P * R_inverse * Q' * S * B  where  diagonal matrix S scaled the rows of A
 
     for i in R_Index loop
        B2(i) := Row_Scalings(i) * B(i); -- Row_Scalings was initialized to 1.
     end loop;
 
     -- X is out, so init:
 
     for j in C_Index loop
        X(j)  := Zero;
        X0(j) := Zero;
     end loop;
 
     -- Get B1 = Q'*B2
 
     B1 := Q_transpose_x_Col_Vector (Q, B2);
 
     -- put Column B1 into a row B0 of A (ie indexed by C_Index):
 
     for i in Starting_Col .. Final_Col loop
        Row   := R_Index(Int64 (i-Starting_Col) + Int64(Starting_Row));
        B0(i) := B1(Row);
     end loop;
 
     -- Now R*X0 = B0; solve for X0.
 
     -- No_of_Rows >= No_of_Cols.  We stop at Final_Col, so
     -- effectively we are assuming that R is square (upper triangular actually).
 
     Row := Starting_Row;
     for Col in Starting_Col .. Final_Col loop
        if Abs R(Row,Col) < Abs (Singularity_Cutoff * Max_R_Diagonal_Val) or
           Abs R(Row,Col) < Min_Allowed_Real then
           R_Diagonal_Inv(Col) :=  Zero;
        else
           R_Diagonal_Inv(Col) := One / R(Row,Col);
        end if;
        if Row < Final_Row then Row := Row + 1; end if;
     end loop;
 
     X0(Final_Col) := B0(Final_Col) * R_Diagonal_Inv(Final_Col);
 
     if Final_Col > Starting_Col then
     
        for Col in reverse Starting_Col .. Final_Col-1 loop
           Row := R_Index (Int64 (Col-Starting_Col) + Int64 (Starting_Row));
           Sum := Zero;
           for j in Col+1 .. Final_Col loop
              Sum := Sum + R(Row, j) * X0(j);
           end loop;
           X0(Col) := (B0(Col) - Sum) * R_Diagonal_Inv(Col);
        end loop;
        
     end if;
 
     for i in C_Index loop
        X(Col_Permutation(i)) := X0(i);
     end loop;
 
     -- x = S * P * R_inverse * Q' * B
 
  end QR_Solve;
 
  ------------------
  -- Q_x_V_Matrix --
  ------------------
 
  -- Get Product = Q*V
  -- The columns of V are vectors indexed by R_index.
  -- Not Q' times Matrix V but Q*V, so
  -- apply the inverses of all the 2x2 rotation matrices in the order opposite to the
  -- order in which they were created: Start w/ most recent and wrk backwards in time.
  -- The columns of V must be same length as Columns of original matrix A, (hence Q).
  -- The number of these Cols is unimportant (set statically: V is a generic formal).
  -- Q is MxM: No_of_Rows x No_of_Rows, since Q*R = A, and R has shape of A.
 
  procedure Q_x_V_Matrix
    (Q : in Q_Matrix;
     V : in out V_Matrix) 
  is
     Final_Row    : R_Index renames Q.Final_Row;
     Final_Col    : C_Index renames Q.Final_Col;
     Starting_Row : R_Index renames Q.Starting_Row;
     Starting_Col : C_Index renames Q.Starting_Col;
     s, c         : Real;
     V_pvt, V_Low : Real;
 
     Pivot_Row : R_Index := Starting_Row;
  begin
 
     -- The loop bounds reflect the design of Q. Q is stored in
     -- matrix A, which was defined to be M x N via Starting_Row,
     -- Starting_Col, etc.
 
     Pivot_Row := R_Index (Integer(Final_Col)
                            + Integer (Starting_Row) - Integer(Starting_Col));
 
     -- The following 2 loops range over all the 2x2 rotation matrices in Q.
     -- Each 2x2 rotates 2 Row_Vectors of V.
     -- The columns of V must be same length as Columns of original matrix A, (hence Q).
 
     for Pivot_Col in reverse Starting_Col .. Final_Col loop 
        if Pivot_Row < Final_Row then
        for Low_Row in reverse Pivot_Row+1 .. Final_Row loop
 
           s := Q.Rot(Low_Row, Pivot_Col).Sine;
           c := Q.Rot(Low_Row, Pivot_Col).Cosine;
 
           if Q.P_bigger_than_L(Low_Row, Pivot_Col) then --abs A(piv,piv )> abs A(low,piv)
 
              --if abs s > Zero then -- rot isn't identity: use fact c is always positive.
              for Col in V'Range(2) loop
                 V_pvt := V(Pivot_Row, Col);
                 V_Low := V(Low_Row, Col);
                 V(Pivot_Row, Col)  := V_pvt - s*(-c*V_pvt +   V_Low);
                 V(Low_Row, Col)    := V_Low - s*(  -V_pvt - c*V_Low);
              end loop;
              --end if;
 
           else
 
              for Col in V'Range(2) loop
                 V_pvt := V(Pivot_Row, Col);
                 V_Low := V(Low_Row, Col);
                 V(Pivot_Row, Col)  :=-V_Low + c*(   V_Pvt - s*V_Low);
                 V(Low_Row, Col)    := V_Pvt + c*(+s*V_Pvt +   V_Low);
              end loop;
 
           end if;
 
        end loop; -- Low_Row  in   Final_Row  to Pivot_Row+1 
        end if;
        if Pivot_Row > Starting_Row then   Pivot_Row := Pivot_Row - 1;   end if;
     end loop;
 
  end Q_x_V_Matrix;
 
  ----------------------------
  -- Q_transpose_x_V_Matrix --
  ----------------------------
 
  -- Get Product = Q'*V
  --
  -- Apply the the 2x2 rotation matrices in the order in which they were created.
  -- Start w/ oldest and wrk forwards in time.  Use the same rotation used to
  -- create R (ie the actual rotation stored on Q, not its inverse.
  -- Q really contains Q_transpose.)
  --
  -- The columns of V are vectors indexed by R_index.
  -- The columns of V must be same length as Columns of original matrix A, (hence Q).
  -- The number of these Cols is unimportant (set statically: V is a generic formal).
  -- Q is MxM: No_of_Rows x No_of_Rows, (since Q*R = A, and R has shape of A).
  --
  procedure Q_transpose_x_V_Matrix
    (Q : in Q_Matrix;
     V : in out V_Matrix) 
  is
     Final_Row    : R_Index renames Q.Final_Row;
     Final_Col    : C_Index renames Q.Final_Col;
     Starting_Row : R_Index renames Q.Starting_Row;
     Starting_Col : C_Index renames Q.Starting_Col;
     s, c         : Real;
     V_pvt, V_Low : Real;
 
     Pivot_Row : R_Index := Starting_Row;
  begin
 
     -- The loop bounds reflect the design of Q. Q is stored in
     -- matrix A, which was defined to be M x N via Starting_Row,
     -- Starting_Col, etc.
 
     -- The following 2 loops range over all the 2x2 rotation matrices in Q.
     -- Each 2x2 rotates 2 Row_Vectors of V.
     -- The columns of V must be same length as Columns of original matrix A, (hence Q).
 
     Pivot_Row := Starting_Row;
 
     for Pivot_Col in Starting_Col .. Final_Col loop 
        if Pivot_Row < Final_Row then
        for Low_Row in Pivot_Row+1 .. Final_Row loop
 
           s := Q.Rot(Low_Row, Pivot_Col).Sine;
           c := Q.Rot(Low_Row, Pivot_Col).Cosine;
 
           if Q.P_bigger_than_L(Low_Row, Pivot_Col) then --abs A(piv,piv )> abs A(low,piv)
 
              --if abs s > Zero then -- rot not identity: use fact c is always positive.
              for Col in V'Range(2) loop
                 V_pvt := V(Pivot_Row, Col);
                 V_Low := V(Low_Row, Col);
                 V(Pivot_Row, Col)  := V_pvt + s*(c*V_pvt +   V_Low);
                 V(Low_Row, Col)    := V_Low + s*( -V_pvt + c*V_Low);
              end loop;
              --end if;
 
           else
 
              for Col in V'Range(2) loop
                 V_pvt := V(Pivot_Row, Col);
                 V_Low := V(Low_Row, Col);
                 V(Pivot_Row, Col)  := V_Low + c*(   V_Pvt + s*V_Low);
                 V(Low_Row, Col)    :=-V_Pvt + c*(-s*V_Pvt +   V_Low);
              end loop;
 
           end if;
 
        end loop;
        Pivot_Row := Pivot_Row + 1;
        end if;
     end loop;
 
  end Q_transpose_x_V_Matrix;
 
  --------------------
  -- Q_x_Matrix --
  --------------------
 
  -- Get Product = Q*A
 
  procedure Q_x_Matrix
    (Q : in Q_Matrix;
     A : in out A_Matrix)
  is
     B : Col_Vector := (others => Zero);
     C : Col_Vector;
 
     Final_Row    : R_Index renames Q.Final_Row;
     Final_Col    : C_Index renames Q.Final_Col;
     Starting_Row : R_Index renames Q.Starting_Row;
     Starting_Col : C_Index renames Q.Starting_Col;
  begin
 
     for Col in Starting_Col .. Final_Col loop
 
       for Row in Starting_Row .. Final_Row loop
          B(Row) := A(Row, Col);
       end loop;
 
       C := Q_x_Col_Vector (Q, B);
 
       for Row in Starting_Row .. Final_Row loop
          A(Row, Col) := C(Row);
       end loop;
 
     end loop;
 
  end Q_x_Matrix;
 
end Givens_QR;
