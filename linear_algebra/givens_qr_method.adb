
---------------------------------------------------------------------------
-- package body Givens_QR_Method
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

package body Givens_QR_Method is

  package math is new Ada.Numerics.Generic_Elementary_Functions (Real);
  use math;

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
     cn              : Real    := One;
     cn_minus_1      : Real    := Zero;
     sn_minus_1      : Real    := Zero;
     P_bigger_than_L : Boolean := True;
     Skip_Rotation   : Boolean := True;
     Pivot_Col       : C_Index := C_Index'First;
     Hi_Row          : R_Index := R_Index'First;
     Lo_Row          : R_Index := R_Index'First;
  end record;

  ---------------------------------
  -- Lower_Diagonal_QR_Iteration --
  ---------------------------------

  -- Operates only on square real blocks.

  procedure Lower_Diagonal_QR_Iteration
    (A            : in out A_Matrix;
     Q            : in out A_Matrix;
     Shift        : in     Real;
     Starting_Col : in     C_Index  := C_Index'First;
     Final_Col    : in     C_Index  := C_Index'Last)
  is
     Hi_Row, Lo_Row : R_Index;
     P, L : Real;

     Rotations : array (C_Index) of Rotation;

     sn, cn          : Real;
     cn_minus_1      : Real;
     sn_minus_1      : Real;
     P_bigger_than_L : Boolean;
     Skip_Rotation   : Boolean;

     -------------------------------------------
     -- e_Multiply_A_on_RHS_with_Transpose_of --
     -------------------------------------------

     -- multiply A on the right by R_transpose:
     --    A = A * R_transpose
     --
     -- Use enhanced precision rotations here.

     procedure e_Multiply_A_on_RHS_with_Transpose_of
       (R : in Rotation)
     is
        sn                 : constant Real    := R.sn;
        cn                 : constant Real    := R.cn;
        cn_minus_1         : constant Real    := R.cn_minus_1;
        sn_minus_1         : constant Real    := R.sn_minus_1;
        P_bigger_than_L    : constant Boolean := R.P_bigger_than_L;
        Skip_Rotation      : constant Boolean := R.Skip_Rotation;
        Pivot_Row          : constant R_Index := R.Hi_Row;
        Low_Row            : constant R_Index := R.Lo_Row;

        A_pvt, A_low : Real;
     begin

        if  Skip_Rotation  then  return; end if;
 
        -- Rotate corresponding columns. Multiply on RHS by transpose
        -- of above givens matrix (second step of similarity transformation).
        -- (Low_Row is Lo visually, but its index is higher than Pivot's.)
  
        if P_bigger_than_L then   -- |s| < |c|
  
           for r in Starting_Col .. Final_Col  loop 
              A_pvt := A(r, Pivot_Row);
              A_low := A(r, Low_Row);
              A(r, Pivot_Row) := A_pvt + (cn_minus_1*A_pvt +  sn * A_low);
              A(r, Low_Row)   := A_low + (-sn * A_pvt + cn_minus_1*A_low);
           end loop;
  
        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1
  
           for r in Starting_Col .. Final_Col  loop 
              A_pvt := A(r, Pivot_Row);
              A_low := A(r, Low_Row);
              A(r, Pivot_Row) := A_low + (cn * A_pvt +  sn_minus_1*A_low);
              A(r, Low_Row)   :=-A_pvt + (-sn_minus_1*A_pvt + cn * A_low);
           end loop;
  
        end if;
  
     end e_Multiply_A_on_RHS_with_Transpose_of;

     -----------------------------------------
     -- Rotate_to_Kill_Element_Lo_of_pCol --
     -----------------------------------------

     -- Try to zero out A(Lo_Row, Pivot_Col) with a similarity transformation.
     -- In other words, multiply A on left by R:
     --
     --    A = R * A
     --
     -- and multiply Q on right by R_transpose:
     --
     --    Q = Q * R_transpose

     procedure Rotate_to_Kill_Element_Lo_of_pCol
       (R : in Rotation)
     is
        sn              : constant Real    := R.sn;
        cn              : constant Real    := R.cn;
        cn_minus_1      : constant Real    := R.cn_minus_1;
        sn_minus_1      : constant Real    := R.sn_minus_1;
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
           for c in Pivot_Col .. Final_Col loop   -- works only for upper hessenbergs
              A_pvt := A(Pivot_Row, c);
              A_low := A(Low_Row, c);
              A(Pivot_Row, c) := A_pvt + (cn_minus_1*A_pvt +  sn * A_low);
              A(Low_Row, c)   := A_low + (-sn * A_pvt + cn_minus_1*A_low);
           end loop;
 
        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1
  
       --for c in Starting_Col .. Final_Col loop 
         for c in Pivot_Col .. Final_Col loop   -- works only for upper hessenbergs
              A_pvt := A(Pivot_Row, c);
              A_low := A(Low_Row, c);
              A(Pivot_Row, c) := A_low + (cn * A_pvt +  sn_minus_1*A_low);
              A(Low_Row, c)   :=-A_pvt + (-sn_minus_1*A_pvt + cn * A_low);
           end loop;
 
        end if;
  
        -- Rotate corresponding columns of Q. (Multiply on RHS by transpose
        -- of above givens matrix to accumulate full Q.)

        if P_bigger_than_L then   -- |s| < |c|
  
           for r in Starting_Col .. Final_Col  loop 
              Q_pvt := Q(r, Pivot_Row);
              Q_low := Q(r, Low_Row);
              Q(r, Pivot_Row) := Q_pvt + (cn_minus_1*Q_pvt +  sn * Q_low);
              Q(r, Low_Row)   := Q_low + (-sn * Q_pvt + cn_minus_1*Q_low);
           end loop;
 
        else  -- Abs_P <= Abs_L, so     abs t := abs (P / L) <= 1
  
           for r in Starting_Col .. Final_Col  loop 
              Q_pvt := Q(r, Pivot_Row);
              Q_low := Q(r, Low_Row);
              Q(r, Pivot_Row) := Q_low + (cn * Q_pvt +  sn_minus_1*Q_low);
              Q(r, Low_Row)   :=-Q_pvt + (-sn_minus_1*Q_pvt + cn * Q_low);
           end loop;
  
        end if;
 
     end Rotate_to_Kill_Element_Lo_of_pCol;

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

     type Diag_Storage is array(C_Index) of Real;

     Lost_Bits : Diag_Storage;

     hypot : Real := Zero;

  begin

     if (Final_Col - Starting_Col) < 2 then   return;  end if;

     Lost_Bits := (others => 0.0);

     -- Subtract 'Shift' from each diagonal element of A.
     -- Sum = A(c, c) + (-Shift) 
     -- Sum = Small + Large.  Lost_Bits = Small - (Sum - Large)
     declare
        Sum, Dropped_Bits : Real;
     begin
        if Abs Shift > Zero then
           for c in Starting_Col .. Final_Col loop
              Sum_with_Dropped_Bits (A(c,c), -Shift, Sum, Dropped_Bits);
              A(c, c)      := Sum;
              Lost_Bits(c) := Dropped_Bits;
           end loop;
        else
           Lost_Bits := (others => 0.0);
        end if;
     end;

     for Pivot_Col in Starting_Col .. Final_Col-1 loop

        Hi_Row := Pivot_Col;
        Lo_Row := Hi_Row + 1;
 
        P := A(Hi_Row, Pivot_Col); 
        L := A(Lo_Row, Pivot_Col);

        Get_Rotation_That_Zeros_Out_Low
          (P, L, sn, cn, cn_minus_1, sn_minus_1, hypot, P_bigger_than_L, Skip_Rotation);

        Rotations(Pivot_Col).sn              := sn;
        Rotations(Pivot_Col).cn              := cn;
        Rotations(Pivot_Col).cn_minus_1      := cn_minus_1;
        Rotations(Pivot_Col).sn_minus_1      := sn_minus_1;
        Rotations(Pivot_Col).P_bigger_than_L := P_bigger_than_L;
        Rotations(Pivot_Col).Skip_Rotation   := Skip_Rotation;
        Rotations(Pivot_Col).Pivot_Col       := Pivot_Col;
        Rotations(Pivot_Col).Hi_Row          := Hi_Row;
        Rotations(Pivot_Col).Lo_Row          := Lo_Row;
 
      --Rotations(Pivot_Col).Skip_Rotation   := false; -- for testing

        Rotate_to_Kill_Element_Lo_of_pCol (Rotations(Pivot_Col));
        -- Zeroes out A(Lo_Row, Pivot_Col)
        -- Updates A and Q as global memory.
        -- Applies rotation by multiplying Givens Matrix on LHS of A.
        -- Then multiplies transpose of Givens Matrix on RHS of Q.

        A(Lo_Row, Pivot_Col) := Zero;
        A(Hi_Row, Pivot_Col) := Real'Copy_Sign(hypot, A(Hi_Row, Pivot_Col));

     end loop; -- over Pivot_Col

     -- These can be done inside the above loop (after a delay of 1 step):

     for Pivot_Col in Starting_Col .. Final_Col-1 loop
        e_Multiply_A_on_RHS_with_Transpose_of (Rotations(Pivot_Col));
     end loop;

     -- Add Shift back to A:

     for c in Starting_Col .. Final_Col loop
        A(c, c) := (A(c, c) + Shift) + Lost_Bits(c);  -- best default
     end loop;
 
  end Lower_Diagonal_QR_Iteration;

end Givens_QR_Method;
