
-----------------------------------------------------------------------
-- package body Bidiagonal; bidiagonalizes real valued matrices.
-- Copyright (C) 2015-2018 Jonathan S. Parker
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
-----------------------------------------------------------------------

with Ada.Numerics;
with Ada.Numerics.Generic_Elementary_Functions;
with Givens_Rotation;

package body Bidiagonal is

  package math is new Ada.Numerics.Generic_Elementary_Functions (Real);  use math;

  package Rotate is new Givens_Rotation (Real);  use Rotate;

  Zero : constant Real := +0.0;
  One  : constant Real := +1.0;
  Two  : constant Real := +2.0;

  ----------------
  -- U_Identity --
  ----------------

  function U_Identity return U_Matrix is
     U : U_Matrix := (others => (others => Zero));
  begin
     for r in Row_Index loop
        U(r, r) := One;
     end loop;
     return U;
  end U_Identity;
  
  ----------------
  -- V_Identity --
  ----------------

  function V_Identity return V_Matrix is
     V : V_Matrix := (others => (others => Zero));
  begin
     for c in Col_Index loop
        V(c, c) := One;
     end loop;
     return V;
  end V_Identity;
  
  --------------------
  -- Bi_Diagonalize --
  --------------------

  -- Operates only on square real blocks.
  -- 
  -- Want to use orthogonal transforms to make A into B,
  -- a bi-diagonal matrix with the 1st sub-diagonal non-zero.  
  -- Let Uj be a 2x2 givens rotation matrix,
  -- and let Uj' be its transpose (and inverse).  Then form
  --
  --   A = non-zero
  --
  --   (U1*...*Un) * (Un'*...*U1') * A * (V1*...*Vn) * (Vn'*...*V1') =
  --
  --    U * B * V, 
  --
  -- where B = U' * A * V'.
  -- and   A = U  * B * V.
  --
  --  A = U * B * V  where B is bi-diagonal, with the 1st lower off-diagonal 
  --  non-zero, and all upper off-diagonals zero.  (So B is lower triangular.)
  -- 

  procedure Bi_Diagonalize
    (A                : in out A_Matrix;  -- A becomes the B in A = U * B * V'
     V                :    out V_Matrix;  -- 
     U                :    out U_Matrix;  -- Initialized with Identity
     Initial_V_Matrix : in     V_Matrix;  -- Normally just use function V_Identity.
     Initial_U_Matrix : in     U_Matrix;  -- Normally just use function U_Identity.
     Starting_Col     : in     Col_Index  := Col_Index'First;
     Final_Col        : in     Col_Index  := Col_Index'Last;
     Final_Row        : in     Row_Index  := Row_Index'Last;
     Matrix_U_Desired : in     Boolean    := True;
     Matrix_V_Desired : in     Boolean    := True)
  is
     Starting_Row : constant Row_Index := Row_Index (Starting_Col);

     type Integer64 is range -2**63+1 .. 2**63-1;
     No_of_Rows : constant Integer64 := Integer64(Final_Row)-Integer64(Starting_Row)+1;
     No_of_Cols : constant Integer64 := Integer64(Final_Col)-Integer64(Starting_Col)+1;
     pragma Assert (No_of_Rows >= No_of_Cols);

     Pivot_Row : Row_Index;
     Final_Pivot_Col, Lo_Col : Col_Index;

     ---------------------------------------
     -- Rotate_to_Kill_Element_Lo_of_pCol --
     ---------------------------------------

     -- Zero out A(Lo_Row, pCol)

     procedure Rotate_to_Kill_Element_Lo_of_pCol
       (pCol   : in Col_Index; 
        Hi_Row : in Row_Index;
        Lo_Row : in Row_Index)
     is
        Pivot_Col : Col_Index renames pCol;
        Pivot_Row : Row_Index renames Hi_Row;
        Low_Row   : Row_Index renames Lo_Row;

        pragma Assert(Row_Index (Pivot_Col) = Pivot_Row); 
        -- R_index contains Col_Index

        sn, cs, cs_minus_1, sn_minus_1 : Real;
        hypot : Real;
        P_bigger_than_L    : Boolean;
        Skip_Rotation      : Boolean;
        P : constant Real := A(Pivot_Row, Pivot_Col); -- P is for Pivot
        L : constant Real := A(Low_Row, Pivot_Col);
        A_pvt, A_low, U_pvt, U_low : Real;
     begin
        Get_Rotation_That_Zeros_Out_Low
          (P, L, sn, cs, cs_minus_1, sn_minus_1, hypot, P_bigger_than_L, Skip_Rotation);

        if Skip_Rotation then  return;  end if;
  
        -- Rotate rows. 
        --
        -- Want    U B V' = A     (B = bi-diagonal.)
        --
        -- So generally,   I A I = A becomes,
        --                (I * G' * G) * A * (F * F' * I) = A
        --                          U  * B  * V'          = A
        --
        -- So the desired U will be the product of the G' matrices
        -- which we obtain by repeatedly multiplying I on the RHS by G'.
        --
        -- Start by multiplying A on LHS by givens rotation G.
 
        if P_bigger_than_L then   -- |s| < |c|
  
           for c in Pivot_Col .. Final_Col  loop
              A_pvt := A(Pivot_Row, c);
              A_low := A(Low_Row, c);
              A(Pivot_Row, c) := A_pvt + ( cs_minus_1*A_pvt + sn*A_low);
              A(Low_Row, c)   := A_low + (-sn*A_pvt + cs_minus_1*A_low);
           end loop;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1

           for c in Pivot_Col .. Final_Col  loop
              A_pvt := A(Pivot_Row, c);
              A_low := A(Low_Row, c);
              A(Pivot_Row, c) := A_low + ( cs*A_pvt + sn_minus_1*A_low);
              A(Low_Row, c)   :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_low);
           end loop;
 
        end if;
  
        A(Pivot_Row, Pivot_Col) := Real'Copy_Sign(hypot, A(Pivot_Row, Pivot_Col));

        -- Rotate corresponding columns of U. (Multiply on RHS by G', the
        -- transpose of above givens matrix to accumulate full U.)
        --
        -- U is (Row_Index x Row_Index).

        if Matrix_U_Desired then
        if P_bigger_than_L then   -- |s| < |c|
  
           for r in Starting_Row .. Final_Row  loop 
              U_pvt := U(r, Pivot_Row);
              U_low := U(r, Low_Row);
              U(r, Pivot_Row) := U_pvt + ( cs_minus_1*U_pvt + sn*U_low);
              U(r, Low_Row)   := U_low + (-sn*U_pvt + cs_minus_1*U_low);
           end loop;
 
        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1
  
           for r in Starting_Row .. Final_Row  loop 
              U_pvt := U(r, Pivot_Row);
              U_low := U(r, Low_Row);
              U(r, Pivot_Row) := U_low + ( cs*U_pvt + sn_minus_1*U_low);
              U(r, Low_Row)   :=-U_pvt + (-sn_minus_1*U_pvt + cs*U_low);
           end loop;
  
        end if;
        end if;
 
     end Rotate_to_Kill_Element_Lo_of_pCol;

     ---------------------------------------
     -- Rotate_to_Kill_Element_Hi_of_pRow --
     ---------------------------------------

     -- Zero out A(pRow, Hi_Col)

     procedure Rotate_to_Kill_Element_Hi_of_pRow
       (pRow   : in Row_Index; 
        Lo_Col : in Col_Index;
        Hi_Col : in Col_Index)
     is
        sn, cs, cs_minus_1, sn_minus_1 : Real;
        hypot : Real;
        P_bigger_than_H    : Boolean;
        Skip_Rotation      : Boolean;
        Pivot_Row : Row_Index renames pRow;
        Pivot_Col : Col_Index renames Lo_Col;
        A_pvt, A_hi , V_pvt, V_hi  : Real;
        P : constant Real := A(Pivot_Row, Pivot_Col); -- P is for Pivot
        H : constant Real := A(Pivot_Row, Hi_Col);

        pragma Assert(Row_Index (Pivot_Col) = Pivot_Row+1);
     begin

        Get_Rotation_That_Zeros_Out_Low
          (P, H, sn, cs, cs_minus_1, sn_minus_1, hypot, P_bigger_than_H, Skip_Rotation);

        if Skip_Rotation then  return;  end if;
  
        -- Rotate columns. Multiply on RHS by the above Givens matrix. (Hi_Col
        -- is to the rt. visually, and its index is higher than Pivot's.)
  
        if P_bigger_than_H then   -- |s| < |c|

           for r in Pivot_Row .. Final_Row  loop
              A_pvt := A(r, Pivot_Col);
              A_hi  := A(r, Hi_Col);
              A(r, Pivot_Col) := A_pvt + ( cs_minus_1*A_pvt + sn*A_hi);
              A(r, Hi_Col)    := A_hi  + (-sn*A_pvt + cs_minus_1*A_hi);
           end loop;
  
        else  -- Abs_P <= Abs_H, so     abs t := abs (P / H) <= 1
  
           for r in Pivot_Row .. Final_Row  loop
              A_pvt := A(r, Pivot_Col);
              A_hi  := A(r, Hi_Col);
              A(r, Pivot_Col) := A_hi  + ( cs*A_pvt + sn_minus_1*A_hi);
              A(r, Hi_Col)    :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_hi);
           end loop;
  
        end if;
 
        A(Pivot_Row, Pivot_Col) := Real'Copy_Sign(hypot, A(Pivot_Row, Pivot_Col));

        -- Rotate corresponding rows of V', hence the columns of V. 
        -- (Multiply on LHS of V' by transpose of above givens matrix to 
        -- accumulate full V.)
        --
        -- V is (Col_Index x Col_Index).
        --
        -- We are calculating  U * B * V' by inserting pairs of givens
        -- rotation matrices between the B and the V'.  When we rotate rows
        -- of V', we are rotating the 2 columns (Pivot_Col, Hi_Col) of V: 

        if Matrix_V_Desired then
        if P_bigger_than_H then   -- |s| < |c|

           for r in Starting_Col .. Final_Col  loop
              V_pvt := V(r, Pivot_Col);
              V_hi  := V(r, Hi_Col);
              V(r, Pivot_Col) := V_pvt + ( cs_minus_1*V_pvt + sn*V_hi);
              V(r, Hi_Col)    := V_hi  + (-sn*V_pvt + cs_minus_1*V_hi);
           end loop;
  
        else  -- Abs_P <= Abs_H, so     abs t := abs (P / H) <= 1
  
           for r in Starting_Col .. Final_Col  loop
              V_pvt := V(r, Pivot_Col);
              V_hi  := V(r, Hi_Col);
              V(r, Pivot_Col) := V_hi  + ( cs*V_pvt + sn_minus_1*V_hi);
              V(r, Hi_Col)    :=-V_pvt + (-sn_minus_1*V_pvt + cs*V_hi);
           end loop;
  
        end if;
        end if;
  
     end Rotate_to_Kill_Element_Hi_of_pRow;

     type Column is array(Row_Index) of Real;

     procedure Get_Sqrt_of_Sum_of_Sqrs_of_Col
       (Col_Id       : in     Col_Index;
        Starting_Row : in     Row_Index;
        Ending_Row   : in     Row_Index;
        Col_Sums     :    out Column)
     is
        Max : Real := Zero;
        Un_Scale, Scale : Real := Zero;
        Emin : constant Integer := Real'Machine_Emin;
        Emax : constant Integer := Real'Machine_Emax;
        Abs_A : Real := Zero;
        M_exp : Integer := 0;
     begin

        Max := Abs A(Starting_Row, Col_Id);
        if Ending_Row > Starting_Row then
           for i in Starting_Row+1 .. Ending_Row loop
              Abs_A := Abs A(i, Col_Id);
              if Abs_A > Max then Max := Abs_A; end if;
           end loop;
        end if;

        if Max < Two**Emin then Col_Sums := (others => Zero); return; end if;

        if Max < Two**(Emin / 2)  or else  Max > Two**(Emax / 4) then 
           M_Exp    := Real'Exponent (Max);
           Scale    := Two ** (-M_Exp);
           Un_Scale := Two ** ( M_Exp);
        else
           Scale    := One;
           Un_Scale := One;
        end if;

        Col_Sums(Starting_Row) := Abs A(Starting_Row, Col_Id);

        if Ending_Row > Starting_Row then

           Compensated_Sum:   
           declare
              val, Sum_tmp, Err, Sum : Real := Zero;
           begin
              for r in Starting_Row .. Ending_Row loop
                 Val := (Scale * A(r, Col_Id)) ** 2;
                 Val     := Val - Err;        -- correction to Val, next term in sum.
                 Sum_tmp := Sum + Val;        -- now increment Sum
                 Err     := (Sum_tmp - Sum) - Val;
                 Sum     := Sum_tmp;
                 Col_Sums(r) := Un_Scale * Sqrt (Sum);
              end loop;
           end Compensated_Sum;

        end if;

     end Get_Sqrt_of_Sum_of_Sqrs_of_Col;


     type Row is array(Col_Index) of Real;

     procedure Get_Sqrt_of_Sum_of_Sqrs_of_Row
       (Row_Id       : in     Row_Index;
        Starting_Col : in     Col_Index;
        Ending_Col   : in     Col_Index;
        Row_Sums     :    out Row)
     is
        Emin  : constant Integer := Real'Machine_Emin;
        Emax  : constant Integer := Real'Machine_Emax;
        Max : Real := Zero;
        Un_Scale, Scale : Real := Zero;
        Abs_A : Real := Zero;
        M_exp : Integer := 0;
     begin

        Max := Abs A(Row_Id, Starting_Col);
        if Ending_Col > Starting_Col then
           for c in Starting_Col+1 .. Ending_Col loop
              Abs_A := Abs A(Row_Id, c);
              if Abs_A > Max then Max := Abs_A; end if;
           end loop;
        end if;

        if Max < Two**Emin then Row_Sums := (others => Zero); return; end if;

        if Max < Two**(Emin / 2) or else  Max > Two**(Emax / 4) then 
           M_Exp    := Real'Exponent (Max);
           Scale    := Two ** (-M_Exp);
           Un_Scale := Two ** ( M_Exp);
        else
           Scale    := One;
           Un_Scale := One;
        end if;

        Row_Sums(Starting_Col) := Abs A(Row_Id, Starting_Col);

        if Ending_Col > Starting_Col then

           Compensated_Sum:   
           declare
              Term, Sum_tmp, Err, Sum : Real := Zero;
           begin
              for c in Starting_Col .. Ending_Col loop
                 Term    := (Scale * A(Row_Id, c)) ** 2;
                 Term    := Term - Err;        -- correction to Term, next term in sum.
                 Sum_tmp := Sum + Term;        -- now increment Sum
                 Err     := (Sum_tmp - Sum) - Term;
                 Sum     := Sum_tmp;
                 Row_Sums(c) := Un_Scale * Sqrt (Sum);
              end loop;
           end Compensated_Sum;

        end if;

     end Get_Sqrt_of_Sum_of_Sqrs_of_Row;

     Row_Sums : Row    := (others => Zero);
     Col_Sums : Column := (others => Zero);

  begin

     V := Initial_V_Matrix;  
     U := Initial_U_Matrix;  
     -- Always identity unless you preprocessed A with an L*Q transform etc.

     if  No_of_Cols < 2 then  return;  end if; -- can do a 2x2

     if Row_Index (Final_Col) < Final_Row then
        Final_Pivot_Col := Final_Col;
     else
        Final_Pivot_Col := Final_Col - 1;  -- A is square
     end if;

     for Pivot_Col in Starting_Col .. Final_Pivot_Col loop

        -- take lo elements in col = Pivot_Col, and rotate them to diagonal.
 
        Pivot_Row := Row_Index (Pivot_Col);
        -- rotate element from Low_Row to Pivot_Row

        Get_Sqrt_of_Sum_of_Sqrs_of_Col -- gd
          (Col_id       => Pivot_Col,
           Starting_Row => Pivot_Row,
           Ending_Row   => Final_Row,
           Col_Sums     => Col_Sums);

        for Low_Row in Pivot_Row+1 .. Final_Row loop -- order important.
        -- 50% slower because of M(r,c) c stuff

           Rotate_to_Kill_Element_Lo_of_pCol   -- zero out A(Lo_Row, pCol)
             (pCol   => Pivot_Col,
              Hi_Row => Pivot_Row, -- Hi = high to eye; its id is lower than Lo_Row.
              Lo_Row => Low_Row);

           A(Low_Row, Pivot_Col)   := Zero;
           A(Pivot_Row, Pivot_Col) := 
              Real'Copy_Sign (Col_Sums(Low_Row), A(Pivot_Row, Pivot_Col));

        end loop;

--        A(Pivot_Row, Pivot_Col) := 
--           Real'Copy_Sign (Col_Sums(Final_Row), A(Pivot_Row, Pivot_Col));

        -- Rotate columns.
        -- Take hi elements in row = Pivot_Col, and Zero them out by
        -- rotating them to 1st superdiagonal (column Pivot_Col + 1).

        if Pivot_Col < Final_Col-1 then    -- Lo_Col+1 below is Pivot_Col+2

           Pivot_Row := Row_Index (Pivot_Col);
           Lo_Col    := Pivot_Col + 1;  -- low to the eye
   
           Get_Sqrt_of_Sum_of_Sqrs_of_Row 
             (Row_Id       => Pivot_Row,
              Starting_Col => Lo_Col,
              Ending_Col   => Final_Col,
              Row_Sums     => Row_Sums);

           for Hi_Col in Lo_Col+1 .. Final_Col loop

              Rotate_to_Kill_Element_Hi_of_pRow   -- zero out A(pRow, Hi_Col)
                (pRow   => Pivot_Row,
                 Lo_Col => Lo_Col,
                 Hi_Col => Hi_Col);
     
              A(Pivot_Row, Hi_Col) := Zero;
              A(Pivot_Row, Lo_Col) := 
                 Real'Copy_Sign (Row_Sums(Hi_Col), A(Pivot_Row, Lo_Col));
     
           end loop;  --  over Hi_Col
--           A(Pivot_Row, Lo_Col) := 
--              Real'Copy_Sign (Row_Sums(Final_Col), A(Pivot_Row, Lo_Col));

        end if;

     end loop; -- over Pivot_Col

  end Bi_Diagonalize;

  -- Must input matrices V and U.  A must be bidiagonal - zero everywhere,
  -- except the diagonal and the superdiagonal. This isn't checked.
  -- Not optimized for speed.

  procedure Zero_Shift_Bidiagonal_QR
    (A                : in out A_Matrix;  -- A becomes the B in A = U * B * V'
     V                : in out V_Matrix;  -- 
     U                : in out U_Matrix;  -- Initialized with Identity
     Starting_Col     : in     Col_Index  := Col_Index'First;
     Final_Col        : in     Col_Index  := Col_Index'Last;
     Final_Row        : in     Row_Index  := Row_Index'Last;
     Matrix_U_Desired : in     Boolean    := True;
     Matrix_V_Desired : in     Boolean    := True)
  is
     Starting_Row : constant Row_Index := Row_Index (Starting_Col);

     type Integer64 is range -2**63+1 .. 2**63-1;
     No_of_Rows : constant Integer64 := Integer64(Final_Row)-Integer64(Starting_Row)+1;
     No_of_Cols : constant Integer64 := Integer64(Final_Col)-Integer64(Starting_Col)+1;
     pragma Assert (No_of_Rows >= No_of_Cols);

     Pivot_Row : Row_Index;
     Final_Pivot_Col : Col_Index;

     ---------------------------------------
     -- Rotate_to_Kill_Element_Lo_of_pCol --
     ---------------------------------------

     -- Zero out A(Lo_Row, pCol)

     procedure Rotate_to_Kill_Element_Lo_of_pCol
       (pCol   : in Col_Index; 
        Hi_Row : in Row_Index;
        Lo_Row : in Row_Index)
     is
        Pivot_Col : Col_Index renames pCol;
        Pivot_Row : Row_Index renames Hi_Row;
        Low_Row   : Row_Index renames Lo_Row;

        pragma Assert(Row_Index (Pivot_Col) = Pivot_Row); 

        sn, cs, cs_minus_1, sn_minus_1 : Real;
        hypot : Real;
        P_bigger_than_L    : Boolean;
        Skip_Rotation      : Boolean;
        P : constant Real := A(Pivot_Row, Pivot_Col); -- P is for Pivot
        L : constant Real := A(Low_Row, Pivot_Col);
        A_pvt, A_low, U_pvt, U_low : Real;
     begin
        Get_Rotation_That_Zeros_Out_Low
          (P, L, sn, cs, cs_minus_1, sn_minus_1, hypot, P_bigger_than_L, Skip_Rotation);

        if Skip_Rotation then  return;  end if;
  
        -- Rotate rows. 
        --
        -- Want    U B V' = A     (B = bi-diagonal.)
        --
        -- So generally,   I A I = A becomes,
        --                (I * G' * G) * A * (F * F' * I) = A
        --                          U  * B  * V'          = A
        --
        -- So the desired U will be the product of the G' matrices
        -- which we obtain by repeatedly multiplying I on the RHS by G'.
        --
        -- Start by multiplying A on LHS by givens rotation G.
 
        if P_bigger_than_L then   -- |s| < |c|
  
           for c in Pivot_Col .. Pivot_Col+1  loop
              A_pvt := A(Pivot_Row, c);
              A_low := A(Low_Row, c);
              A(Pivot_Row, c) := A_pvt + ( cs_minus_1*A_pvt + sn*A_low);
              A(Low_Row, c)   := A_low + (-sn*A_pvt + cs_minus_1*A_low);
           end loop;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1

           for c in Pivot_Col .. Pivot_Col+1  loop
              A_pvt := A(Pivot_Row, c);
              A_low := A(Low_Row, c);
              A(Pivot_Row, c) := A_low + ( cs*A_pvt + sn_minus_1*A_low);
              A(Low_Row, c)   :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_low);
           end loop;
 
        end if;
  
        A(Pivot_Row, Pivot_Col) := Real'Copy_Sign(hypot, A(Pivot_Row, Pivot_Col));

        -- Rotate corresponding columns of U. (Multiply on RHS by G', the
        -- transpose of above givens matrix to accumulate full U.)
        --
        -- U is (Row_Index x Row_Index).

        if Matrix_U_Desired then
        if P_bigger_than_L then   -- |s| < |c|
  
           for r in Starting_Row .. Final_Row  loop 
              U_pvt := U(r, Pivot_Row);
              U_low := U(r, Low_Row);
              U(r, Pivot_Row) := U_pvt + ( cs_minus_1*U_pvt + sn*U_low);
              U(r, Low_Row)   := U_low + (-sn*U_pvt + cs_minus_1*U_low);
           end loop;
 
        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1
  
           for r in Starting_Row .. Final_Row  loop 
              U_pvt := U(r, Pivot_Row);
              U_low := U(r, Low_Row);
              U(r, Pivot_Row) := U_low + ( cs*U_pvt + sn_minus_1*U_low);
              U(r, Low_Row)   :=-U_pvt + (-sn_minus_1*U_pvt + cs*U_low);
           end loop;
  
        end if;
        end if;
 
     end Rotate_to_Kill_Element_Lo_of_pCol;

     ---------------------------------------
     -- Rotate_to_Kill_Element_Hi_of_pRow --
     ---------------------------------------

     -- Zero out A(pRow, Hi_Col)

     procedure Rotate_to_Kill_Element_Hi_of_pRow
       (pRow   : in Row_Index; 
        Lo_Col : in Col_Index;
        Hi_Col : in Col_Index)
     is
        sn, cs, cs_minus_1, sn_minus_1 : Real;
        hypot : Real;
        P_bigger_than_H    : Boolean;
        Skip_Rotation      : Boolean;
        Pivot_Row : Row_Index renames pRow;
        Pivot_Col : Col_Index renames Lo_Col;
        A_pvt, A_hi , V_pvt, V_hi  : Real;
        P : constant Real := A(Pivot_Row, Pivot_Col); -- P is for Pivot
        H : constant Real := A(Pivot_Row, Hi_Col);

        pragma Assert(Row_Index (Pivot_Col) = Pivot_Row);
     begin

        Get_Rotation_That_Zeros_Out_Low
          (P, H, sn, cs, cs_minus_1, sn_minus_1, hypot, P_bigger_than_H, Skip_Rotation);

        if Skip_Rotation then  return;  end if;
  
        -- Rotate columns. Multiply on RHS by the above Givens matrix. (Hi_Col
        -- is to the rt. visually, and its index is higher than Pivot's.)
  
        if P_bigger_than_H then   -- |s| < |c|

           for r in Pivot_Row .. Pivot_Row+1  loop
              A_pvt := A(r, Pivot_Col);
              A_hi  := A(r, Hi_Col);
              A(r, Pivot_Col) := A_pvt + ( cs_minus_1*A_pvt + sn*A_hi);
              A(r, Hi_Col)    := A_hi  + (-sn*A_pvt + cs_minus_1*A_hi);
           end loop;
  
        else  -- Abs_P <= Abs_H, so     abs t := abs (P / H) <= 1
  
           for r in Pivot_Row .. Pivot_Row+1  loop
              A_pvt := A(r, Pivot_Col);
              A_hi  := A(r, Hi_Col);
              A(r, Pivot_Col) := A_hi  + ( cs*A_pvt + sn_minus_1*A_hi);
              A(r, Hi_Col)    :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_hi);
           end loop;
  
        end if;
 
        A(Pivot_Row, Pivot_Col) := Real'Copy_Sign(hypot, A(Pivot_Row, Pivot_Col));

        -- Rotate corresponding rows of V', hence the columns of V. 
        -- (Multiply on LHS of V' by transpose of above givens matrix to 
        -- accumulate full V.)
        --
        -- V is (Col_Index x Col_Index).
        --
        -- We are calculating  U * B * V' by inserting pairs of givens
        -- rotation matrices between the B and the V'.  When we rotate rows
        -- of V', we are rotating the 2 columns (Pivot_Col, Hi_Col) of V: 

        if Matrix_V_Desired then
        if P_bigger_than_H then   -- |s| < |c|

           for r in Starting_Col .. Final_Col  loop
              V_pvt := V(r, Pivot_Col);
              V_hi  := V(r, Hi_Col);
              V(r, Pivot_Col) := V_pvt + ( cs_minus_1*V_pvt + sn*V_hi);
              V(r, Hi_Col)    := V_hi  + (-sn*V_pvt + cs_minus_1*V_hi);
           end loop;
  
        else  -- Abs_P <= Abs_H, so     abs t := abs (P / H) <= 1
  
           for r in Starting_Col .. Final_Col  loop
              V_pvt := V(r, Pivot_Col);
              V_hi  := V(r, Hi_Col);
              V(r, Pivot_Col) := V_hi  + ( cs*V_pvt + sn_minus_1*V_hi);
              V(r, Hi_Col)    :=-V_pvt + (-sn_minus_1*V_pvt + cs*V_hi);
           end loop;
  
        end if;
        end if;
  
     end Rotate_to_Kill_Element_Hi_of_pRow;

  begin

     if  No_of_Cols < 2  then  return;  end if; -- can do a 2x2

     if Row_Index (Final_Col) < Final_Row then
        Final_Pivot_Col := Final_Col;
     else
        Final_Pivot_Col := Final_Col - 1;  -- A is square
     end if;

     -- Rotate away the superdiagonal by rotating Pivot_Col and Pivot_Col+1

     for Pivot_Col in Starting_Col .. Final_Col-1 loop

        -- Rotate columns.
        -- Take hi element in (row, col) = (Pivot_Col, Pivot_Col+1) and
        -- Zero it out by rotating to diagonal (column Pivot_Col).

        Pivot_Row := Row_Index (Pivot_Col);

        Rotate_to_Kill_Element_Hi_of_pRow   -- zero out A(pRow, Hi_Col)
          (pRow   => Pivot_Row,
           Lo_Col => Pivot_Col,
           Hi_Col => Pivot_Col + 1);

        A(Pivot_Row, Pivot_Col + 1) := Zero;
  
     end loop; -- over Pivot_Col

     -- Rotate away the subdiagonal by rotating Pivot_Row and Pivot_Row+1

     for Pivot_Col in Starting_Col .. Final_Pivot_Col loop

        -- Rotate rows.
        -- take lo element in col = Pivot_Col, and rotate it to diagonal.
 
        Pivot_Row := Row_Index (Pivot_Col);

        Rotate_to_Kill_Element_Lo_of_pCol   -- zero out A(Lo_Row, pCol)
          (pCol   => Pivot_Col,
           Hi_Row => Pivot_Row, -- Hi = high to eye; its id is lower than Lo_Row.
           Lo_Row => Pivot_Row + 1);

        A(Pivot_Row + 1, Pivot_Col) := Zero;

     end loop;

  end Zero_Shift_Bidiagonal_QR;

end Bidiagonal;


