
---------------------------------------------------------------------------
-- package body Hessenberg
-- Copyright (C) 2011-2018 Jonathan S. Parker
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
with Givens_Rotation;

package body Hessenberg is

  package math is new Ada.Numerics.Generic_Elementary_Functions (Real); use math;

  package Rotate is new Givens_Rotation (Real);  use Rotate;

  Zero : constant Real := +0.0;
  One  : constant Real := +1.0;
  Two  : constant Real := +2.0;

  function Identity return A_Matrix is
     Q : A_Matrix;
  begin
     Q := (others => (others => Zero));
     for c in C_Index loop
        Q(c, c) := One;
     end loop;
     return Q;
  end Identity;


  ----------------------
  -- Upper_Hessenberg --
  ----------------------

  -- Operates only on square real blocks.
  --
  -- Want to use similarity transforms to make A into H,
  -- upper Hessenberg.  Let Qj be a 2x2 givens rotation matrix,
  -- and let Qj' be its transpose (and inverse).  Then form
  --
  --   A =
  --
  --   (Q1*...*Qn) * (Qn'*...*Q1') * A * (Q1*...*Qn) * (Qn'*...*Q1') =
  --
  --    Q * H * Q',
  --
  -- where H = Q' * A * Q.
  --
  -- To complete the decomposition of A to H, insert Qj * Qj' into
  -- Q * H * Q' to get   Q * (Qj * Qj') * H * (Qj * Qj') * Q'.
  --
  -- So to develop Q, we rotate columns of Q by multiplying on RHS with
  -- Qj.  H gets multiplied on the LHS by Qj' (rotating rows to zero out
  -- the lower triangular region) and on the RHS by Qj.
  --
  -- Wind up with the eigenvalue equation  A =  Q * H * Q' which becomes
  --
  --   A * Q = Q * H.
  --
  -- If H were diagonal, then the column vectors of Q would be the eigenvectors
  -- and the diagonal elements of H would be the eigenvalues.

  procedure Upper_Hessenberg
    (A            : in out A_Matrix;
     Q            :    out A_Matrix;
     Starting_Col : in     C_Index  := C_Index'First;
     Final_Col    : in     C_Index  := C_Index'Last;
     Initial_Q    : in     A_Matrix := Identity)
  is
     Final_Row : constant C_Index := Final_Col;
     Pivot_Row : R_Index;

     ---------------------------------------
     -- Rotate_to_Kill_Element_Lo_of_pCol --
     ---------------------------------------

     -- Zero out A(Lo_Row, pCol) with a similarity transformation. In
     -- other words, multiply A on left by Q_tr and on right by Q:
     --    A_h = Q_transpose * A * Q

     procedure Rotate_to_Kill_Element_Lo_of_pCol
       (pCol   : in C_Index;
        Hi_Row : in R_Index;
        Lo_Row : in R_Index)
     is
        sn, cs : Real;
        cs_minus_1 : Real;
        sn_minus_1 : Real;
        hypot : Real;
        P_bigger_than_L    : Boolean;
        Skip_Rotation      : Boolean;
        Pivot_Col : C_Index renames pCol;
        Pivot_Row : R_Index renames Hi_Row;
        Low_Row   : R_Index renames Lo_Row;
        A_pvt, A_low, Q_pvt, Q_low : Real;
        P : constant Real := A(Pivot_Row, Pivot_Col); -- P is for Pivot
        L : constant Real := A(Low_Row, Pivot_Col);
     begin
        Get_Rotation_That_Zeros_Out_Low
          (P, L, sn, cs, cs_minus_1, sn_minus_1, hypot, P_bigger_than_L, Skip_Rotation);

        if Skip_Rotation then  return;  end if;

        -- Rotate rows. Multiply on LHS by givens rotation G.
        -- Want    Q' A Q = H = upper hessenberg.
        -- Each step is:    G A G' = partial H
        -- So the desired Q will be the product of the G' matrices
        -- which we obtain by repeatedly multiplying I on the RHS by G'.

        if P_bigger_than_L then   -- |s| < |c|

           for c in Pivot_Col .. Final_Col loop
              A_pvt := A(Pivot_Row, c);
              A_low := A(Low_Row, c);
              A(Pivot_Row, c) := A_pvt + ( cs_minus_1*A_pvt + sn*A_low);
              A(Low_Row, c)   := A_low + (-sn*A_pvt + cs_minus_1*A_low);
           end loop;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1

           for c in Pivot_Col .. Final_Col loop
              A_pvt := A(Pivot_Row, c);
              A_low := A(Low_Row, c);
              A(Pivot_Row, c) := A_low + ( cs*A_pvt + sn_minus_1*A_low);
              A(Low_Row, c)   :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_low);
           end loop;

        end if;

        A(Pivot_Row, Pivot_Col) := Real'Copy_Sign (Hypot, A(Pivot_Row, Pivot_Col));

        -- Rotate corresponding columns. Multiply on RHS by transpose
        -- of above givens matrix (second step of similarity transformation).
        -- (Low_Row is Lo visually, but its index is higher than Pivot's.)

        if P_bigger_than_L then   -- |s| < |c|

           for r in Starting_Col .. Final_Col loop
              A_pvt := A(r, Pivot_Row);
              A_low := A(r, Low_Row);
              A(r, Pivot_Row) := A_pvt + ( cs_minus_1*A_pvt + sn*A_low);
              A(r, Low_Row)   := A_low + (-sn*A_pvt + cs_minus_1*A_low);
           end loop;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1

           for r in Starting_Col .. Final_Col loop
              A_pvt := A(r, Pivot_Row);
              A_low := A(r, Low_Row);
              A(r, Pivot_Row) := A_low + ( cs*A_pvt + sn_minus_1*A_low);
              A(r, Low_Row)   :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_low);
           end loop;

        end if;

        -- Rotate corresponding columns of Q. (Multiply on RHS by transpose
        -- of above givens matrix to accumulate full Q.)

        if P_bigger_than_L then   -- |s| < |c|

           for r in Starting_Col .. Final_Col loop
              Q_pvt := Q(r, Pivot_Row);
              Q_low := Q(r, Low_Row);
              Q(r, Pivot_Row) := Q_pvt + ( cs_minus_1*Q_pvt + sn*Q_low);
              Q(r, Low_Row)   := Q_low + (-sn*Q_pvt + cs_minus_1*Q_low);
           end loop;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1

           for r in Starting_Col .. Final_Col loop
              Q_pvt := Q(r, Pivot_Row);
              Q_low := Q(r, Low_Row);
              Q(r, Pivot_Row) := Q_low + ( cs*Q_pvt + sn_minus_1*Q_low);
              Q(r, Low_Row)   :=-Q_pvt + (-sn_minus_1*Q_pvt + cs*Q_low);
           end loop;

        end if;

     end Rotate_to_Kill_Element_Lo_of_pCol;

     type Column is array(R_Index) of Real;

     procedure Get_Sqrt_of_Sum_of_Sqrs_of_Col
       (Col_Id       : in     C_Index;
        Starting_Row : in     R_Index;
        Ending_Row   : in     R_Index;
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

     Col_Sums : Column := (others => Zero);

  begin

     Q := Initial_Q;

     if (Final_Col - Starting_Col) < 2 then   return;  end if;

     for Pivot_Col in Starting_Col .. Final_Col-2 loop

        Pivot_Row := Pivot_Col + 1;

        Get_Sqrt_of_Sum_of_Sqrs_of_Col -- version _2 Only takes Sqrt of full column sum
          (Col_id       => Pivot_Col,
           Starting_Row => Pivot_Row,
           Ending_Row   => Final_Row,
           Col_Sums     => Col_Sums);

        for Low_Row in Pivot_Row+1 .. Final_Row loop

           Rotate_to_Kill_Element_Lo_of_pCol   -- zero out A(Lo_Row, pCol)
             (pCol   => Pivot_Col,
              Hi_Row => Pivot_Row, -- Hi = high to eye; its id is lower than Lo_Row.
              Lo_Row => Low_Row);

           A(Low_Row, Pivot_Col)   := Zero;
           A(Pivot_Row, Pivot_Col) :=
              Real'Copy_Sign (Col_Sums(Low_Row), A(Pivot_Row, Pivot_Col));

        end loop;  --  over Low_Row

     end loop; -- over Pivot_Col

  end Upper_Hessenberg;

  ----------------------
  -- Lower_Hessenberg --
  ----------------------

  procedure Lower_Hessenberg
    (A            : in out A_Matrix;
     Q            :    out A_Matrix;
     Starting_Col : in     C_Index := C_Index'First;
     Final_Col    : in     C_Index := C_Index'Last;
     Initial_Q    : in     A_Matrix := Identity)
  is
     Starting_Row : constant C_Index := Starting_Col;
     Final_Row    : constant C_Index := Final_Col;
     Pivot_Col : R_Index := Starting_Col;

     ---------------------------------------
     -- Rotate_to_Kill_Element_Hi_of_pRow --
     ---------------------------------------

     -- Zero out A(Lo_Row, pCol) with a similarity transformation. In
     -- other words, multiply A on left by Q_tr and on right by Q:
     --    A_h = Q_transpose * A * Q
     -- Use enhanced precision rotations here.

     procedure Rotate_to_Kill_Element_Hi_of_pRow
       (pRow   : in R_Index;
        Lo_Col : in R_Index;
        Hi_Col : in R_Index)
     is
        sn, cs : Real;
        cs_minus_1 : Real;
        sn_minus_1 : Real;
        hypot : Real;
        P_bigger_than_H    : Boolean;
        Skip_Rotation      : Boolean;
        Pivot_Row : R_Index renames pRow;
        Pivot_Col : C_Index renames Lo_Col;
        A_pvt, A_hi , Q_pvt, Q_hi  : Real;
        P : constant Real := A(Pivot_Row, Pivot_Col); -- P is for Pivot
        H : constant Real := A(Pivot_Row, Hi_Col);
        pragma Assert(Pivot_Col = Pivot_Row+1);
     begin

        Get_Rotation_That_Zeros_Out_Low
          (P, H, sn, cs, cs_minus_1, sn_minus_1, hypot, P_bigger_than_H, Skip_Rotation);

        if Skip_Rotation then  return;  end if;

        -- Rotate columns. Multiply on RHS by transpose
        -- of above givens matrix (second step of similarity transformation).
        -- (Hi_Col is to the rt. visually, and its index is higher than Pivot's.)

        if P_bigger_than_H then   -- |s| < |c|

           for r in Starting_Row .. Final_Row  loop
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

        -- Rotate corresponding columns of Q. (Multiply on RHS by transpose
        -- of above givens matrix to accumulate full Q.)

        if P_bigger_than_H then   -- |s| < |c|

           for r in Starting_Col .. Final_Col  loop
              Q_pvt := Q(r, Pivot_Col);
              Q_hi  := Q(r, Hi_Col);
              Q(r, Pivot_Col) := Q_pvt + ( cs_minus_1*Q_pvt + sn*Q_hi);
              Q(r, Hi_Col)    := Q_hi  + (-sn*Q_pvt + cs_minus_1*Q_hi);
           end loop;
  
        else  -- Abs_P <= Abs_H, so     abs t := abs (P / H) <= 1
  
           for r in Starting_Col .. Final_Col  loop
              Q_pvt := Q(r, Pivot_Col);
              Q_hi  := Q(r, Hi_Col);
              Q(r, Pivot_Col) := Q_hi  + ( cs*Q_pvt + sn_minus_1*Q_hi);
              Q(r, Hi_Col)    :=-Q_pvt + (-sn_minus_1*Q_pvt + cs*Q_hi);
           end loop;

        end if;

        -- Rotate rows. Multiply on LHS by givens rotation G.
        -- Want    Q' A Q = H = upper hessenberg.
        -- Each step is:    G A G' = partial H
        -- So the desired Q will be the product of the G' matrices
        -- which we obtain by repeatedly multiplying I on the RHS by G'.

        if P_bigger_than_H then   -- |s| < |c|

           for c in Starting_Col .. Final_Col  loop
              A_pvt := A(Pivot_Col, c);
              A_hi  := A(Hi_Col, c);
              A(Pivot_Col, c) := A_pvt + ( cs_minus_1*A_pvt + sn*A_hi);
              A(Hi_Col, c)    := A_hi  + (-sn*A_pvt + cs_minus_1*A_hi);
           end loop;

        else  -- Abs_P <= Abs_H, so     abs t := abs (P / H) <= 1

           for c in Starting_Col .. Final_Col  loop
              A_pvt := A(Pivot_Col, c);
              A_hi  := A(Hi_Col, c);
              A(Pivot_Col, c) := A_hi  + ( cs*A_pvt + sn_minus_1*A_hi);
              A(Hi_Col, c)    :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_hi);
           end loop;

        end if;

     end Rotate_to_Kill_Element_Hi_of_pRow;

  begin

     Q := Initial_Q;

     if (Final_Col - Starting_Col) < 2 then   return;  end if;

     for Pivot_Row in Starting_Row .. Final_Row-2 loop

        Pivot_Col := Pivot_Row + 1; -- the lower off-diagonal

        for Hi_Col in Pivot_Col+1 .. Final_Col loop

           Rotate_to_Kill_Element_Hi_of_pRow   -- zero out A(pRow, Hi_Col)
             (pRow   => Pivot_Row,
              Lo_Col => Pivot_Col,
              Hi_Col => Hi_Col);

           A(Pivot_Row, Hi_Col) := Zero;

        end loop;  --  over Low_Row

     end loop; -- over Pivot_Col

  end Lower_Hessenberg;

end Hessenberg;
