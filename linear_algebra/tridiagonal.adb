
---------------------------------------------------------------------------
-- package body Tridiagonal, symmetric matrix tridiagonalization
-- Copyright (C) 2018 Jonathan S. Parker
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

package body Tridiagonal is

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

  type Rotation is record
     sn              : Real    := Zero;
     cs              : Real    := One;
     cs_minus_1      : Real    := Zero;
     sn_minus_1      : Real    := Zero;
     hypotenuse      : Real    := Zero;
     P_bigger_than_L : Boolean := True;
     Skip_Rotation   : Boolean := True;
     Pivot_Col       : C_Index := C_Index'First;
     Hi_Row          : R_Index := R_Index'First;
     Lo_Row          : R_Index := R_Index'First;
  end record;

  type Row is array(C_Index) of Real;

  -- Sum = Small + Large.  Lost_Bits = Small - (Sum - Large)

  procedure Sum_with_Dropped_Bits
    (A, B : in Real;
     Sum  : out Real;
     Dropped_Bits : out Real)
  is
  begin
     Sum := A + B;
     if Abs A > Abs B then
        Dropped_Bits := B - (Sum - A);
     else
        Dropped_Bits := A - (Sum - B);
     end if;
  end Sum_with_Dropped_Bits;

  ---------------------------------------
  -- Arg_1_is_Negligible_Respect_Arg_2 --
  ---------------------------------------

  -- If Abs_x is negligible in comparison to Abs_y  then  return True.

  function Arg_1_is_Negligible_Respect_Arg_2 (x, y : Real) return Boolean is
     Abs_x      : constant Real := Abs x;
     Abs_y      : constant Real := Abs y;
     Min_Allowed_Real : constant Real := Two**(Real'Machine_Emin + 16);

     Added_Range  : constant      := 3;  -- 2, 3 seem bst; 3 gives 2**(-53)
     Eps_Factor   : constant Real := Real'Epsilon * Two**(-Added_Range);
  begin
     if Abs_x < Min_Allowed_Real and then Abs_x <= Abs_y then  -- eg, Abs_x = 0
        return True;
     elsif Abs_x < Abs_y * Eps_Factor then
        return True;
     else
        return False;
     end if;
  end Arg_1_is_Negligible_Respect_Arg_2;

  ---------------------------------
  -- Lower_Diagonal_QR_Iteration --
  ---------------------------------

  -- Operates only on square real blocks.

  procedure Lower_Diagonal_QR_Iteration
    (A                : in out A_Matrix;
     Q                : in out A_Matrix;
     Shift            : in     Real;
     Final_Shift_Col  : in     C_Index  := C_Index'Last;
     Starting_Col     : in     C_Index  := C_Index'First;
     Final_Col        : in     C_Index  := C_Index'Last;
     Q_Matrix_Desired : in     Boolean  := True)
  is
     Rotations : array (C_Index) of Rotation;

     -----------------------------------------
     -- Multiply_A_on_RHS_with_Transpose_of --
     -----------------------------------------

     -- multiply A on the right by R_transpose:
     --    A = A * R_transpose

     procedure Multiply_A_on_RHS_with_Transpose_of
       (R : in Rotation)
     is
        sn                 : Real    renames R.sn;
        cs                 : Real    renames R.cs;
        cs_minus_1         : Real    renames R.cs_minus_1;
        sn_minus_1         : Real    renames R.sn_minus_1;
        P_bigger_than_L    : Boolean renames R.P_bigger_than_L;
        Skip_Rotation      : Boolean renames R.Skip_Rotation;
        Pivot_Row          : R_Index renames R.Hi_Row;
        Low_Row            : R_Index renames R.Lo_Row;

        A_pvt, A_low : Real;
     begin

        if  Skip_Rotation  then  return; end if;

        -- Rotate corresponding columns. Multiply on RHS by transpose
        -- of above givens matrix (second step of similarity transformation).
        -- (Low_Row is Lo visually, but its index is higher than Pivot's.)

        if P_bigger_than_L then   -- |s| < |c|

           for r in Index'Base'Max(Starting_Col, Pivot_Row-1) .. Pivot_Row+1  loop
              A_pvt := A(r, Pivot_Row);
              A_low := A(r, Low_Row);
              A(r, Pivot_Row) := A_pvt + ( cs_minus_1*A_pvt + sn*A_low);
              A(r, Low_Row)   := A_low + (-sn*A_pvt + cs_minus_1*A_low);
           end loop;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1

           for r in Index'Base'Max(Starting_Col, Pivot_Row-1) .. Pivot_Row+1  loop
              A_pvt := A(r, Pivot_Row);
              A_low := A(r, Low_Row);
              A(r, Pivot_Row) := A_low + ( cs*A_pvt + sn_minus_1*A_low);
              A(r, Low_Row)   :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_low);
           end loop;

        end if;

     end Multiply_A_on_RHS_with_Transpose_of;

     ---------------------------------------
     -- Rotate_to_Kill_Element_Lo_of_pCol --
     ---------------------------------------

     -- Try to zero out A(Lo_Row, Pivot_Col) with a similarity transformation.
     -- In other words, multiply A on left by R:
     --
     --    A = R * A
     --
     -- and multiply Q on right by R_transpose:
     --
     --    Q = Q * R_transpose
     --

     procedure Rotate_to_Kill_Element_Lo_of_pCol
       (R : in Rotation)
     is
        sn              : constant Real := R.sn;
        cs              : constant Real := R.cs;
        cs_minus_1      : constant Real := R.cs_minus_1;
        sn_minus_1      : constant Real := R.sn_minus_1;
        P_bigger_than_L : constant Boolean := R.P_bigger_than_L;
        Skip_Rotation   : constant Boolean := R.Skip_Rotation;
        Pivot_Col       : constant C_Index := R.Pivot_Col;
        Pivot_Row       : constant R_Index := R.Hi_Row;
        Low_Row         : constant R_Index := R.Lo_Row;

        A_pvt, A_low, Q_pvt, Q_low : Real;

     begin

        if  Skip_Rotation  then  return; end if;

        if P_bigger_than_L then   -- |s| < |c|

         --for c in Starting_Col .. Final_Col loop
         --for c in Pivot_Col .. Final_Col loop  -- works only for upper hessenbergs
           for c in Pivot_Col .. Pivot_Col+2 loop  -- works only for Tridiagonals
              A_pvt := A(Pivot_Row, c);
              A_low := A(Low_Row, c);
              A(Pivot_Row, c) := A_pvt + ( cs_minus_1*A_pvt + sn*A_low);
              A(Low_Row, c)   := A_low + (-sn*A_pvt + cs_minus_1*A_low);
           end loop;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1

         --for c in Starting_Col .. Final_Col loop
         --for c in Pivot_Col .. Final_Col loop  -- works only for upper hessenbergs
           for c in Pivot_Col .. Pivot_Col+2 loop  -- works only for Tridiagonals
              A_pvt := A(Pivot_Row, c);
              A_low := A(Low_Row, c);
              A(Pivot_Row, c) := A_low + ( cs*A_pvt + sn_minus_1*A_low);
              A(Low_Row, c)   :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_low);
           end loop;

        end if;

        -- Rotate corresponding columns of Q. (Multiply on RHS by transpose
        -- of above givens matrix to accumulate full Q.)

        if Q_Matrix_Desired then
        if P_bigger_than_L then   -- |s| < |c|

           for r in Starting_Col .. Final_Col  loop
              Q_pvt := Q(r, Pivot_Row);
              Q_low := Q(r, Low_Row);
              Q(r, Pivot_Row) := Q_pvt + (cs_minus_1*Q_pvt +  sn*Q_low);
              Q(r, Low_Row)   := Q_low + (-sn*Q_pvt + cs_minus_1*Q_low);
           end loop;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1

           for r in Starting_Col .. Final_Col  loop
              Q_pvt := Q(r, Pivot_Row);
              Q_low := Q(r, Low_Row);
              Q(r, Pivot_Row) := Q_low + (cs*Q_pvt +  sn_minus_1*Q_low);
              Q(r, Low_Row)   :=-Q_pvt + (-sn_minus_1*Q_pvt + cs*Q_low);
           end loop;

        end if;  -- P_bigger_than_L
        end if;  -- Q_Matrix_Desired

     end Rotate_to_Kill_Element_Lo_of_pCol;

     type Diag_Storage is array(Index) of Real;

     Matrix_Size : constant Integer := Integer (Final_Col)-Integer (Starting_Col)+1;

     tst1, tst2 : Real;

     Lost_Bits : Diag_Storage;

     hypot, sn, cs, cs_minus_1, sn_minus_1 : Real;
     P_bigger_than_L    : Boolean;
     Skip_Rotation      : Boolean;

     Hi_Row, Lo_Row : R_Index;
     P, L : Real;

  begin

     if  Matrix_Size < 3  then  return;  end if;  -- can't QR it.

     Lost_Bits := (others => 0.0);

     -- Subtract 'Shift' from each diagonal element of A.
     -- Sum = A(c, c) + (-Shift) 
     -- Sum = Small + Large.  Lost_Bits = Small - (Sum - Large)
     declare
        Sum, Dropped_Bits : Real;
     begin
        for c in Starting_Col .. Final_Shift_Col loop
           Sum_with_Dropped_Bits (A(c,c), -Shift, Sum, Dropped_Bits);
           A(c, c)      := Sum;
           Lost_Bits(c) := Dropped_Bits;
        end loop;
     end;

     for Pivot_Col in Starting_Col .. Final_Shift_Col-1 loop

        Hi_Row := Pivot_Col;
        Lo_Row := Hi_Row + 1;

        P := A(Hi_Row, Pivot_Col);
        L := A(Lo_Row, Pivot_Col);

        tst1 := Abs A(Lo_Row, Pivot_Col);
        tst2 := Abs A(Pivot_Col, Pivot_Col) + Abs A(Lo_Row, Lo_Row);

        if Arg_1_is_Negligible_Respect_Arg_2 (tst1, tst2) then

           Rotations(Pivot_Col).Skip_Rotation := True; -- Rotations is initialized.
           A(Lo_Row, Pivot_Col) := Zero;
           A(Pivot_Col, Lo_Row) := Zero;

        else

           Get_Rotation_That_Zeros_Out_Low
             (P, L, sn, cs, cs_minus_1, sn_minus_1, hypot, P_bigger_than_L, Skip_Rotation);

           Rotations(Pivot_Col) :=
             (sn              => sn,
              cs              => cs,
              cs_minus_1      => cs_minus_1,
              sn_minus_1      => sn_minus_1,
              hypotenuse      => hypot,
              P_bigger_than_L => P_bigger_than_L,
              Skip_Rotation   => Skip_Rotation,
              Pivot_Col       => Pivot_Col,
              Hi_Row          => Hi_Row,
              Lo_Row          => Lo_Row);

           -- Zero out A(Lo_Row, Pivot_Col). Update A and Q as global memory.
           -- Apply rotation by multiplying Givens Matrix on LHS of A.
           -- Then multiply transpose of Givens Matrix on RHS of Q.

           Rotate_to_Kill_Element_Lo_of_pCol (Rotations(Pivot_Col));

           A(Hi_Row, Pivot_Col) := Real'Copy_Sign (hypot, A(Hi_Row, Pivot_Col));
           A(Lo_Row, Pivot_Col) := Zero;

        end if;

     end loop; -- over Pivot_Col

     -- These can be done inside the above loop (after a delay of 1 step):

     for Pivot_Col in Starting_Col .. Final_Shift_Col-1 loop
        Multiply_A_on_RHS_with_Transpose_of (Rotations(Pivot_Col));
     end loop;

     -- Add Shift back to A:

     for c in Starting_Col .. Final_Shift_Col loop
        A(c, c) := (A(c, c) + Shift) + Lost_Bits(c);  -- best default
     end loop;

     for c in Starting_Col+1 .. Final_Col loop
      --A(c-1, c) := Half * (A(c, c-1) + A(c-1, c));
        A(c, c-1) := A(c-1, c);
     end loop;

     for c in Starting_Col+2 .. Final_Col loop
        A(c-2, c) := Zero;
     end loop;

  end Lower_Diagonal_QR_Iteration;

  --------------------
  -- Tridiagonalize --
  --------------------

  -- Operates only on square real blocks.
  --
  -- Want to use similarity transforms to make A into T,
  -- Tridiagonal.  Let Qj be a 2x2 givens rotation matrix,
  -- and let Qj' be its transpose (and inverse).  Then form
  --
  --   A =
  --
  --   (Q1*...*Qn) * (Qn'*...*Q1') * A * (Q1*...*Qn) * (Qn'*...*Q1') =
  --
  --    Q * T * Q',
  --
  -- where T = Q' * A * Q.
  --
  -- To complete the decomposition of A to T, insert Qj * Qj' into
  -- Q * T * Q' to get   Q * (Qj * Qj') * T * (Qj * Qj') * Q'.
  --
  -- So to develop Q, we rotate columns of Q by multiplying on RHS with
  -- Qj.  H gets multiplied on the LHS by Qj' (rotating rows to zero out
  -- the lower triangular region) and on the RHS by Qj.
  --
  -- Wind up with the eigenvalue equation  A =  Q * T * Q' which becomes
  --
  --   A * Q = Q * T.
  --
  -- If T were diagonal, then the column vectors of Q would be the eigenvectors
  -- and the diagonal elements of T would be the eigenvalues.

  procedure Tridiagonalize
    (A                : in out A_Matrix;
     Q                :    out A_Matrix;
     Starting_Col     : in     C_Index  := C_Index'First;
     Final_Col        : in     C_Index  := C_Index'Last;
     Initial_Q        : in     A_Matrix := Identity;
     Q_Matrix_Desired : in     Boolean  := True)
  is
     -----------------------------------------
     -- Rotate_to_Kill_Element_Hi_of_pRow --
     -----------------------------------------

     -- Zero out A(Lo_Row, pCol) with a similarity transformation. In
     -- other words, multiply A on left by Q_tr and on right by Q:
     --    A_t = Q_transpose * A * Q
     -- Use enhanced precision rotations here.

     procedure Rotate_to_Kill_Element_Hi_of_pRow
       (pRow   : in R_Index;
        Hi_Col : in C_Index;
        Lo_Col : in C_Index)
     is
        sn, cs : Real;
        cs_minus_1 : Real;
        sn_minus_1 : Real;
        P_bigger_than_L    : Boolean;
        Skip_Rotation      : Boolean;
        Pivot_Row : R_Index renames pRow;
        Pivot_Col : C_Index renames Lo_Col;
        A_pvt, A_low, Q_pvt, Q_low : Real;
        P : constant Real := A(Pivot_Row, Pivot_Col); -- P is for Pivot
        L : constant Real := A(Pivot_Row, Hi_Col);
        hypot : Real;
     begin
        --if not A(Pivot_Row, Pivot_Col)'valid or else not A(Pivot_Row, Hi_Col)'valid then
        --   raise Constraint_Error with "Invalid input in Rotate_to_Kill..";
        --end if;

        Get_Rotation_That_Zeros_Out_Low
          (P, L, sn, cs, cs_minus_1, sn_minus_1, hypot, P_bigger_than_L, Skip_Rotation);

        if Skip_Rotation then  return;  end if;

        -- Rotate rows. Multiply on LHS by givens rotation G.
        -- Want    Q' A Q = H = upper hessenberg.
        -- Each step is:    G A G' = partial H
        -- So the desired Q will be the product of the G' matrices
        -- which we obtain by repeatedly multiplying I on the RHS by G'.

        if P_bigger_than_L then   -- |s| < |c|

           for r in Pivot_Row .. Final_Col  loop
              A_pvt := A(r, Pivot_Col);
              A_low := A(r, Hi_Col);
              A(r, Pivot_Col) := A_pvt + ( cs_minus_1*A_pvt + sn*A_low);
              A(r, Hi_Col)    := A_low + (-sn*A_pvt + cs_minus_1*A_low);
           end loop;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1

           for r in Pivot_Row .. Final_Col  loop
              A_pvt := A(r, Pivot_Col);
              A_low := A(r, Hi_Col);
              A(r, Pivot_Col) := A_low + ( cs*A_pvt + sn_minus_1*A_low);
              A(r, Hi_Col)    :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_low);
           end loop;

        end if;

        -- Rotate corresponding columns. Multiply on RHS by transpose
        -- of above givens matrix (second step of similarity transformation).
        -- (Hi_Col is Lo visually, but its index is higher than Pivot's.)

        -- Do the above 2 rotations in a single step (tho' it hardly matters):
        --
        --   A(Pivot_Col, Pivot_Col) := cs**2*x + cs*sn*2.0*z + sn**2*y;
        --   A(Hi_Col, Hi_Col)       := sn**2*x - cs*sn*2.0*z + cs**2*y;
        --   A(Pivot_Col, Hi_Col)    := cn**2*z - sn**2*z + sn*cs*(y - x);
        --
        --   x_new := A(Pivot_Col, Pivot_Col) := x + (cs*sn*2.0*z + sn**2*(y - x));
        --   y_new := A(Hi_Col, Hi_Col)       := y - (cs*sn*2.0*z + sn**2*(y - x));
        --   z_new := A(Pivot_Col, Hi_Col)    := z - sn**2*z2 + sn*cs*(y - x);
        --
        --   x_new := A(Pivot_Col, Pivot_Col) := x + sn*(cs*2.0*z + sn*(y - x));
        --   y_new := A(Hi_Col, Hi_Col)       := y - sn*(cs*2.0*z + sn*(y - x));
        --   z_new := A(Pivot_Col, Hi_Col)    := z - sn*(sn*2.0*z - cs*(y - x));

        if P_bigger_than_L then   -- |s| < |c|
  
           declare r : constant Index := Pivot_Col;  begin 
              A_pvt := A(Pivot_Col, r);
              A_low := A(Hi_Col, r);
              A(Pivot_Col, r) := A_pvt + ( cs_minus_1*A_pvt + sn*A_low);
              A(Hi_Col, r)    := A_low + (-sn*A_pvt + cs_minus_1*A_low);
           end;
           declare r : constant Index := Hi_Col;  begin
              A_pvt := A(Pivot_Col, r);
              A_low := A(Hi_Col, r);
            --A(Pivot_Col, r) := A_pvt + ( cs_minus_1*A_pvt + sn*A_low);
              A(Hi_Col, r)    := A_low + (-sn*A_pvt + cs_minus_1*A_low);
           end;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1
  
           declare r : constant Index := Pivot_Col;  begin 
              A_pvt := A(Pivot_Col, r);
              A_low := A(Hi_Col, r);
              A(Pivot_Col, r) := A_low + ( cs*A_pvt + sn_minus_1*A_low);
              A(Hi_Col, r)    :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_low);
           end;
           declare r : constant Index := Hi_Col;  begin 
              A_pvt := A(Pivot_Col, r);
              A_low := A(Hi_Col, r);
            --A(Pivot_Col, r) := A_low + ( cs*A_pvt + sn_minus_1*A_low);
              A(Hi_Col, r)    :=-A_pvt + (-sn_minus_1*A_pvt + cs*A_low);
           end;
        end if;

        -- Have modified 2 cols of entire (symmetric) matrix: Pivot_Col, Hi_Col.
        -- Now copy the cols to the 2 rows you get by transposing the 2 cols.

        for c in Starting_Col .. Final_Col  loop
           A(Pivot_Col, c) := A(c, Pivot_Col);
        end loop;

        for c in Starting_Col .. Final_Col  loop
           A(Hi_Col, c) := A(c, Hi_Col);
        end loop;

        -- Rotate corresponding columns of Q. (Multiply on RHS by transpose
        -- of above givens matrix to accumulate full Q.)

        if Q_Matrix_Desired then
        if P_bigger_than_L then   -- |s| < |c|

           for r in Starting_Col .. Final_Col  loop
              Q_pvt := Q(r, Pivot_Col);
              Q_low := Q(r, Hi_Col);
              Q(r, Pivot_Col) := Q_pvt + ( cs_minus_1*Q_pvt + sn*Q_low);
              Q(r, Hi_Col)    := Q_low + (-sn*Q_pvt + cs_minus_1*Q_low);
           end loop;

        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1

           for r in Starting_Col .. Final_Col  loop
              Q_pvt := Q(r, Pivot_Col);
              Q_low := Q(r, Hi_Col);
              Q(r, Pivot_Col) := Q_low + ( cs*Q_pvt + sn_minus_1*Q_low);
              Q(r, Hi_Col)    :=-Q_pvt + (-sn_minus_1*Q_pvt + cs*Q_low);
           end loop;

        end if;  -- P_bigger_than_L
        end if;  -- Q_Matrix_Desired

     end Rotate_to_Kill_Element_Hi_of_pRow;


     procedure Get_Sqrt_of_Sum_of_Sqrs_of_Row
       (Row_Id       : in     R_Index;
        Starting_Col : in     C_Index;
        Ending_Col   : in     C_Index;
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

     Row_Sums : Row := (others => Zero);

     type Int64 is range -2**63+1 .. 2**63-1;

     Matrix_Size : constant Int64 := Int64 (Final_Col)-Int64 (Starting_Col)+1;
     Pivot_Col : C_Index;
  begin

     if Matrix_Size < 4 then   return;  end if;  -- already tridiagonalized.

     Q := Initial_Q;

     for Pivot_Row in Starting_Col .. Final_Col-2 loop

        Pivot_Col := Pivot_Row + 1;

        Get_Sqrt_of_Sum_of_Sqrs_of_Row
          (Row_Id       => Pivot_Row,
           Starting_Col => Pivot_Col,
           Ending_Col   => Final_Col,
           Row_Sums     => Row_Sums);

        for Hi_Col in Pivot_Col+1 .. Final_Col loop

           Rotate_to_Kill_Element_Hi_of_pRow   -- zero out A(Pivot_Row, Hi_Col)
             (pRow     => Pivot_Row,
              Hi_Col   => Hi_Col, -- Hi = high to eye; its id is lower than Lo_Row.
              Lo_Col   => Pivot_Col);

           A(Pivot_Row, Hi_Col)    := Zero;
           A(Hi_Col, Pivot_Row)    := Zero;

           A(Pivot_Row, Pivot_Col) :=
              Real'Copy_Sign (Row_Sums(Hi_Col), A(Pivot_Row, Pivot_Col));

           A(Pivot_Col, Pivot_Row) := A(Pivot_Row, Pivot_Col);

        end loop;  --  over Hi_Col

     end loop; -- over Pivot_Col; max val of Pivot_Col is Final_Col-1

  end Tridiagonalize;

end Tridiagonal;

