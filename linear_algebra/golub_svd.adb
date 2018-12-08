
---------------------------------------------------------------------------
-- package body Golub_SVD, Singular Value Decomposition
-- Copyright (C) 2018 Jonathan S. Parker.
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

with Hypot;
with Ada.numerics.generic_elementary_functions;

package body Golub_SVD is

  package Hypo is new Hypot (Real); use Hypo;

  package math is new Ada.numerics.generic_elementary_functions (Real); use Math;

  Half : constant Real := +0.5;
  Zero : constant Real := +0.0;
  One  : constant Real := +1.0;
  Two  : constant Real := +2.0;
  --  Convenient for cases in which package is instantiated with
  --  extended precision floating point: just make sure a "+" operator
  --  exists that translates the numbers into abstract type Real.
 
  type Int64 is range -2**63+1 .. 2**63-1;

  ----------------------------------------
  --  Arg_1_Is_Negligible_Respect_Arg_2 --
  ----------------------------------------

  function Arg_1_Is_Negligible_Respect_Arg_2 (x, y : Real) return Boolean is
     Abs_x      : constant Real := Abs x;
     Abs_y      : constant Real := Abs y;
     Eps_Factor : constant Real := Real'Epsilon * Two**(-2); -- 6 and 7 fail
     --  Stnd setting for Eps_Factor is:  Real'Epsilon * Two**(-2). -1 is gd.
     --  Convergence failure at 2.0**(-6): random_anti_sym, and compan_0.
     --
     --  Usually Real'Epsilon = 2**(-50) if Real is 15 digits, but
     --  2**(-53) is closer to the actual working precision of IEEE 754.
     --  So using 2.0**(-3) here for 15 digit arith would be like linpack's
     --  original threshold policy.
     --  2.0**(-1) might be best.
     --
     Min_Exp          : constant Integer := Real'Machine_Emin;
     Min_Allowed_Real : constant Real    := Two**(Min_Exp - Min_Exp/4);

  begin
     if Abs_x = Zero then
        return True;
     elsif Abs_x < Min_Allowed_Real and then Abs_x <= Abs_y then       -- policy D
        return True;
     elsif Abs_x <= Abs_y*Eps_Factor then --Abs_x is negligible in comparison to Abs_y
        return True;
     else
        return False;
     end if;
  end Arg_1_Is_Negligible_Respect_Arg_2;

  ----------------------
  -- Rotation_Factors --
  ----------------------

  --  sn = a / sqrt(a**2 + b**2)
  --  cs = b / sqrt(a**2 + b**2)
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
  --  hypot = sqrt(a**2 + b**2)
  --  hypot = |a| * sqrt(1+t*t) = |a| * (1  + t*t/(1+sqrt(1+t*t)) = hi + lo if t<1
  --
  --  hypot := Abs_a * sq;
  --  hypot := Abs_a + Abs_b * (Abs (t) / (sq + One));
  --  hypot := Abs_a + Abs_b * (Abs_b / (Abs_a * (sq + One)));
  --  hypot := Abs_a + Abs_b * (Abs_b / (hypot + Abs_a)); -- iterate
  --

  procedure Rotation_Factors
    (a, b         : in     Real;
     sn_lo, sn_hi :    out Real;  -- a / Sqrt(a*a + b*b)
     cs_lo, cs_hi :    out Real;  -- b / Sqrt(a*a + b*b)
     hypot        :    out Real;
     cs_hi_is_one :    out Boolean)
  is
     Abs_a, Abs_b : Real;
   --min_arg_over_hypot, max_arg_over_hypot_minus_1 : Real;
     Emin : constant Integer := Real'Machine_Emin;
     Min_Allowed_Real : constant Real := Two**(Emin - Emin / 32);
  begin

     if (not a'Valid) or else (not b'Valid) then
        raise Constraint_Error with "a, b invalid in Rotation_Factors";
     end if;

     Abs_a := Abs a;
     Abs_b := Abs b;

     sn_lo := Zero;
     sn_hi := Zero;
     cs_lo := Zero;
     cs_hi := One;
     hypot := Abs_b;
     cs_hi_is_one := True;  --  default: no rotation

     if Abs_a >= Abs_b then

       if Abs_a > Min_Allowed_Real then
     --if Abs_a > Zero then

          cs_hi_is_one := False;
          hypot        := Hypotenuse (Abs_b, Abs_a);
        --Get_Hypotenuse (Abs_b, Abs_a, hypot, min_arg_over_hypot, max_arg_over_hypot_minus_1);

          cs_hi := Zero;
          cs_lo := b / hypot;
        --cs_lo := Real'Copy_Sign (min_arg_over_hypot, b);

          sn_hi := One;
          sn_lo := - cs_lo * (b / (hypot + Abs_a));
        --sn_lo := - |b|**2 / (hypot * (hypot + |a|)) = |a| / hypot - 1;
        --sn_lo := max_arg_over_hypot_minus_1;

          if a < Zero then sn_hi := -sn_hi; sn_lo := -sn_lo; end if;

       end if;

     else  -- Abs_a < Abs_b

       if Abs_b > Min_Allowed_Real then
     --if Abs_b > Zero then

          cs_hi_is_one := True;
          hypot        := Hypotenuse (Abs_b, Abs_a);
        --Get_Hypotenuse (Abs_b, Abs_a, hypot, min_arg_over_hypot, max_arg_over_hypot_minus_1);

          sn_hi := Zero;
          sn_lo := a / hypot;
        --sn_lo := Real'Copy_Sign (min_arg_over_hypot, a);

          cs_hi := One;
          cs_lo := - sn_lo * (a / (hypot + Abs_b));
        --cs_lo := - |a|**2 / (hypot * (hypot + max)) = |b| / hypot - 1;
        --cs_lo := max_arg_over_hypot_minus_1;

          if b < Zero then cs_lo := -cs_lo; cs_hi := -cs_hi; end if;
       end if;
     end if;

  end Rotation_Factors;

  procedure Householder_Bidiagonalize
   (A                : in out A_Matrix;                    --  m x n
    U                :    out U_matrix;                    --  m x m
    V                :    out V_matrix;                    --  n x n
    E, S             :    out Singular_Vector;             --  1 x n
    Bidiag_Final_Col :    out Col_Index;
    Final_Row        : in     Row_Index;
    Starting_Col     : in     Col_Index;         -- Starting_Row is set to Starting_Col
    Final_Col        : in     Col_Index;
    Matrix_U_Desired : in     Boolean   := True;
    Matrix_V_Desired : in     Boolean   := True)
  is

    Min_Exp             : constant Integer := Real'Machine_Emin;
    Min_Allowed_Pivot   : constant Real    := Two**(Min_Exp - Min_Exp/8);
    Min_Allowed_Real_tr : constant Real    := Two**(Min_Exp - Min_Exp/32);
    Min_Allowed_Denom   : constant Real    := Two**(Min_Exp - Min_Exp/32);

    Starting_Row : constant Row_Index := Row_Index (Integer (Starting_Col));

    Work : array (Starting_Row .. Final_Row) of Real;

    Working_Final_Col : Col_Index;

    t, t_Scaler, A_val, E_val, U_val, V_val : Real;
    Scaling_Factor : Real;

    nct, nrt : Col_Index;

    procedure Get_Scale_Factor_of_Col (
       Scale_Factor         : out Real;
       Inverse_Scale_Factor : out Real;
       Desired_Col          : in  Col_Index;
       Start_of_Row         : in  Row_Index)
    is
       Max_A : Real := Zero;
       Max_Exp : Integer;
       Min_Real : constant Real := Two ** (Min_Exp + 4);
    begin
       for r in Start_of_Row .. Final_Row loop
          if Abs A(r, Desired_Col) > Max_A then Max_A := Abs A(r, Desired_Col); end if;
       end loop;

       Max_Exp              := Real'Exponent (Max_A + Min_Real);
       Scale_Factor         := Two ** (-Max_Exp);
       Inverse_Scale_Factor := One / Scale_Factor;

    end Get_Scale_Factor_of_Col;

    procedure Get_Scale_Factor_of_Row (
       Scale_Factor         : out Real;
       Inverse_Scale_Factor : out Real;
       Desired_Row          : in  Singular_Vector;
       Start_Col            : in  Col_Index)
    is
       Max_A : Real := Zero;
       Max_Exp : Integer;
       Min_Real : constant Real := Two ** (Min_Exp + 4);
    begin
       for c in Start_Col .. Final_Col loop
          if Abs Desired_Row(c) > Max_A then Max_A := Abs Desired_Row(c); end if;
       end loop;

       Max_Exp              := Real'Exponent (Max_A + Min_Real);
       Scale_Factor         := Two ** (-Max_Exp);
       Inverse_Scale_Factor := One / Scale_Factor;

    end Get_Scale_Factor_of_Row;

    Scale_Factor, Inverse_Scale_Factor, Sum, tmp : Real;

  begin

    --  If these mats are undesirable, don't fill up memory with Zero's:

    if Matrix_U_Desired then
       U := (others => (others => Zero));
    end if;
    if Matrix_V_Desired then
       V := (others => (others => Zero));
    end if;

    E := (others => Zero);
    S := (others => Zero);

    --  Can use Min to do the calculation only because Starting_Row = Starting_Col:

    nct := Col_Index (Int64'Min (Int64 (Final_Row-1), Int64 (Final_Col)));
    nrt := Col_Index (Int64'Min (Int64 (Final_Row),   Int64 (Final_Col-2)));

    Working_Final_Col := Col_Index'Max (nct, nrt);

    --  No_Rows = m,  No_Cols = n                (Bidiagonal_Order == mn)
    --  m = n-3,      Working_Final_Col = n-3     Bidiagonal_Order = n-2
    --  m = n-2,      Working_Final_Col = n-2     Bidiagonal_Order = n-1
    --  m = n-1,      Working_Final_Col = n-2     Bidiagonal_Order = n-0
    --  m = n-0,      Working_Final_Col = n-1     Bidiagonal_Order = n-0
    --  m = n+1,      Working_Final_Col = n+0     Bidiagonal_Order = n+0
    --  m = n+2,      Working_Final_Col = n+0     Bidiagonal_Order = n+0
    --
    --  so notice that  Working_Final_Col  is truncated by  m=No_Rows.

    Bidiag_Final_Col := Col_Index(Int64'Min(Int64 (Final_Row+1),Int64 (Final_Col)));
    --  The final bidiagonal matrix has order mn = (Bidiag_Final_Col - Starting_Col + 1).

    --  Reduce A to bidiagonal form, storing the diagonal elements
    --  in S and the super-diagonal elements in E.

    Householder:
    declare
      Pivot_Row   : Row_Index;
      U_Pivot_Col : U_Col_Index;
    begin

    for Pivot_Col in Starting_Col .. Working_Final_Col loop

      Pivot_Row   := Row_Index (Pivot_Col); -- No of Rows >= No of Cols, and row 1st = col 1st
      U_Pivot_Col := U_Col_Index (Pivot_Col);

      --  Compute the transformation for the L-th column and
      --  place the L-th diagonal in S(L).

      if  Pivot_Col <= nct  then

        Get_Scale_Factor_of_Col
          (Scale_Factor, Inverse_Scale_Factor, Pivot_Col, Pivot_Row);

        Sum := Zero;
        for r in Pivot_Row .. Final_Row loop
          tmp := A(r, Pivot_Col) * Scale_Factor;  -- scale so max A is near 1
          Sum := Sum + tmp * tmp;
        end loop;
        S(Pivot_Col) := Inverse_Scale_Factor * Sqrt (Sum);

        if  Abs S(Pivot_Col) > Min_Allowed_Pivot then

          if  A(Pivot_Row, Pivot_Col) < Zero  then
            S(Pivot_Col) := -S(Pivot_Col);
          end if;
          Scaling_Factor := One / S(Pivot_Col);

          for r in Pivot_Row .. Final_Row loop
             A(r, Pivot_Col) := A(r, Pivot_Col) * Scaling_Factor;
          end loop;
          A(Pivot_Row, Pivot_Col) := One + A(Pivot_Row, Pivot_Col);

        end if;

        S(Pivot_Col) := -S(Pivot_Col);

      end if;

      --  Apply the transformation.

      if Pivot_Col <= nct and then Abs (S(Pivot_Col)) > Min_Allowed_Real_tr then
         A_val    := A(Pivot_Row, Pivot_Col);
         t_scaler := One / (A_val + Real'Copy_Sign (Min_Allowed_Denom, A_val));
 
         declare
            t, t2 : Real;
         begin
            t := Zero;

            if Pivot_Col < Final_Col then

               for r in Pivot_Row .. Final_Row loop
                  t := t + A(r, Pivot_Col) * A(r, Pivot_Col+1);
               end loop;
               t := -t * t_scaler;
    
               for c in Pivot_Col+1 .. Final_Col-1 loop
                  t2  := Zero; 
                  for r in Pivot_Row .. Final_Row loop
                     A(r, c) := A(r, c) + t * A(r, Pivot_Col);
                     t2      := t2 + A(r, c+1) * A(r, Pivot_Col);
                  end loop;
                  t := -t2 * t_scaler;
               end loop;
 
               for r in Pivot_Row .. Final_Row loop
                  A(r, Final_Col) := A(r, Final_Col) + t * A(r, Pivot_Col);
               end loop;

            end if;

          end;

      end if;

      --  Place the L-th row of A into E for the
      --  subsequent calculation of the row transformation.

      for c in Pivot_Col+1 .. Final_Col loop
        E(c) := A(Pivot_Row, c);
      end loop;

      --  Place the transformation in U for subsequent back multiplication.

      if  Matrix_U_Desired and Pivot_Col <= nct  then
        for r in Pivot_Row .. Final_Row loop
          U(r, U_Pivot_Col) := A(r, Pivot_Col);
        end loop;
      end if;

      --  Compute the L-th row transformation and place the
      --  L-th superdiagonal in E(L).

      if  Pivot_Col <= nrt  then

         E(Pivot_Col) := Zero;
 
         Get_Scale_Factor_of_Row
           (Scale_Factor, Inverse_Scale_Factor, E, Pivot_Col);
 
         Sum := Zero;
         for c in Pivot_Col+1 .. Final_Col loop
            tmp := E(c) * Scale_Factor; --so largest E is near 1
            Sum := Sum + tmp * tmp;
         end loop;
         E(Pivot_Col) := Inverse_Scale_Factor * Sqrt (Sum);
 
         if  Abs E(Pivot_Col) > Min_Allowed_Pivot  then
 
           if  E(Pivot_Col+1) < Zero  then
              E(Pivot_Col) := -E(Pivot_Col);  --  so has same sign now as E(Pivot_Col+1)
           end if;
           Scaling_Factor := One / E(Pivot_Col);
 
           for c in Pivot_Col+1 .. Final_Col loop
              E(c) := E(c) * Scaling_Factor;
           end loop;
           E(Pivot_Col+1) := One + E(Pivot_Col+1);
 
        end if;

        E(Pivot_Col) := -E(Pivot_Col);

        --  Apply the transformation.

        if Pivot_Row < Final_Row and then Abs E(Pivot_Col) > Min_Allowed_Real_tr then

           E_val    := E(Pivot_Col+1);
           t_scaler := One / (E_val + Real'Copy_Sign (Min_Allowed_Denom, E_val));
 
           for r in Pivot_Row+1 .. Final_Row loop
              Work(r) := Zero;
           end loop;
 
           for c in Pivot_Col+1 .. Final_Col loop
              for r in Pivot_Row+1 .. Final_Row loop
                 Work(r) := Work(r) + E(c) * A(r, c);
              end loop;
           end loop;
 
           for c in Pivot_Col+1 .. Final_Col loop
              t := -E(c) * t_scaler;
              for r in Pivot_Row+1 .. Final_Row loop
                 A(r, c) := A(r, c) + t * Work(r);
              end loop;
           end loop;

        end if;

        --  Place the transformation in V for subsequent back multiplication.

        if  Matrix_V_Desired  then
           for r in Pivot_Col+1 .. Final_Col loop  -- V is n x n = Col_Index x Col_Index
              V(r, Pivot_Col) := E(r);
           end loop;
        end if;

      end if;

    end loop; --  Pivot_Col

    end Householder;

    --  Set up the bidiag matrix with col range: Starting_Col .. Bidiag_Final_Col:

    if  nct < Final_Col  then
       S(nct+1) := A(Row_Index (nct+1), nct+1);
    end if;

    if  Int64 (Bidiag_Final_Col) > Int64 (Final_Row)  then
       S(Bidiag_Final_Col) := Zero;
    end if;

    if  nrt+1 < Bidiag_Final_Col  then
       E(nrt+1) := A(Row_Index (nrt+1), Bidiag_Final_Col);
    end if;

    E(Bidiag_Final_Col) := Zero;

    --  If so desired, generate U.

    if  Matrix_U_Desired  then       --  U is    Row_Index x U_Col_Index

    Make_Matrix_U:
    declare
       U_Pivot_Row : Row_Index;
       U_Pivot_Col : U_Col_Index;
    begin

       -- (nct+1) <= U_Col_Index'Last = Row_Index'Last. Always true
       
       for r in Starting_Row .. Final_Row loop
       for c in U_Col_Index (nct+1) .. U_Col_Index'Last loop
          U(r, c) := Zero;
       end loop;
       end loop;
 
       for j in U_Col_Index (nct+1) .. U_Col_Index'Last loop
          U(Row_Index(j), j) := One;
       end loop;
 
       for Pivot_Col in reverse Starting_Col .. nct loop  -- nct < Final_Row
 
          U_Pivot_Row :=   Row_Index (Pivot_Col);
          U_Pivot_Col := U_Col_Index (Pivot_Col);
  
          if  Abs S(Pivot_Col) > Zero  then
  
             U_val    := U(U_Pivot_Row, U_Pivot_Col);
             t_Scaler := One / (U_val + Real'Copy_Sign(Min_Allowed_Denom, U_val));
   
             for c in U_Pivot_Col+1 .. U_Col_Index'Last loop
   
                t := Zero;
                for r in U_Pivot_Row .. Final_Row loop
                   t := t +  U(r, c) * U(r, U_Pivot_Col);
                end loop;
   
                t := -t * t_Scaler;
                for r in U_Pivot_Row .. Final_Row loop
                   U(r, c) := U(r, c) +  t * U(r, U_Pivot_Col);
                end loop;
   
             end loop; --  in c
   
             for r in U_Pivot_Row .. Final_Row loop
                U(r, U_Pivot_Col) := -U(r, U_Pivot_Col);
             end loop;
             U(U_Pivot_Row, U_Pivot_Col) := One + U(U_Pivot_Row, U_Pivot_Col);
  
             for r in Starting_Row .. U_Pivot_Row-1 loop
                U(r, U_Pivot_Col) := Zero;
             end loop;
 
          else

             for r in Starting_Row .. Final_Row loop
                U(r, U_Pivot_Col) := Zero;
             end loop;
             U(U_Pivot_Row, U_Pivot_Col) := One;
 
          end if;  --  Abs S(Pivot_Col) > Zero

       end loop; --  in Pivot_Col

    end Make_Matrix_U;

    end if; --  Matrix_U_Desired

    --  If so required, generate V.

    if  Matrix_V_Desired  then  --  V is defined on   Col_Index x Col_Index

       for Pivot_Col in reverse Starting_Col .. Final_Col loop
 
          if  Pivot_Col <= nrt and Abs E(Pivot_Col) > Zero  then
  
             V_val := V(Pivot_Col+1, Pivot_Col);
             t_Scaler := One / (V_val + Real'Copy_Sign (Min_Allowed_Denom, V_val));
   
             for c in Pivot_Col+1 .. Final_Col loop
                t := Zero;
                for r in Pivot_Col+1 .. Final_Col loop
                   t := t + V(r, c) * V(r, Pivot_Col);
                end loop;
    
                t := -t * t_Scaler;
                for r in Pivot_Col+1 .. Final_Col loop
                   V(r, c) := V(r, c) + t * V(r, Pivot_Col);
                end loop;
             end loop;
   
          end if;
 
          for r in Starting_Col .. Final_Col loop
             V(r, Pivot_Col) := Zero;
          end loop;
          V(Pivot_Col, Pivot_Col) := One;
 
       end loop; --  in Pivot_Col
 
    end if; --  V desired

  end Householder_Bidiagonalize;

  --------------------------
  -- QR_For_Singular_Vals --
  --------------------------

  -- Translated from the LINPACK SVD.

  procedure QR_For_Singular_Vals
    (E, S             : in out Singular_Vector;
     U                : in out U_Matrix;
     V                : in out V_Matrix;
     Id_of_1st_S_Val  : in out Col_Index;
     Starting_Col     : in     Col_Index;
     Final_Col        : in     Col_Index;
     Starting_Row     : in     Row_Index;
     Final_Row        : in     Row_Index;
     Matrix_U_Desired : in     Boolean;
     Matrix_V_Desired : in     Boolean)
  is
     No_of_Iterations : Integer := 0;
     mn : Col_Index := Final_Col; -- mn will be decremented to Starting_Col
     l  : Col_Index'base := mn - 1;
     ls : Col_Index'base := mn;
     kase : Integer;
     g, f : Real;
     Test, tmp : Real;
     Max_measured_Iterations : Integer := 0;
     Min_Exp          : constant Integer := Real'Machine_Emin;
     Min_Allowed_Real : constant Real := 2.0**(Min_Exp + Min_Exp/32);
  begin

     QR_Iteration: loop

     if  No_of_Iterations > Max_measured_Iterations then
        Max_measured_Iterations := No_of_Iterations;
     end if;
     if  No_of_Iterations > Max_Allowed_No_of_Iterations then
        --  We failed to converge while working on Col = mn.  The correct
        --  Singular vals are in S(mn+1), S(mn+2), ...
        if mn = Final_Col then
           raise SVD_Convergence_Error;  -- First attempt failed to converge.
        else
           Id_of_1st_S_Val := mn+1;
        end if;
        return;
     end if;

     --  This section of the program inspects for
     --  negligible elements in the S and E arrays.
     --
     --  On completion the variables KASE and L are set as follows:
     --
     --  KASE := 1     if S(MN) and E(L-1) are negligible and L < MN;
     --  KASE := 2     if S(L) is negligible and L < MN;
     --  KASE := 3     if E(L-1) is negligible, L < MN, and;
     --                  S(L), ..., S(MN) are not negligible (QR step).
     --  KASE := 4     if E(MN-1) is negligible (convergence).;

     -- initialize l:

     for ll in Starting_Col .. mn loop

        l := mn - ll + Starting_Col - 1;

        if ll = mn then
           exit;   --  l = Starting_Col-1  persists below
        end if;

        if Arg_1_Is_Negligible_Respect_Arg_2 (Abs E(l), Abs S(l) + Abs S(l+1)) then
           E(l) := Zero;
           exit;
        end if;

     end loop;

     if  l = mn - 1  then

        kase := 4;

     else

        for lls in l+1 .. mn+1 loop

           ls := mn - lls + l + 1;

           if  ls = l  then
              exit;
           end if;

           Test := Zero;
           if  ls /= mn  then
              Test := Test + Abs E(ls);
           end if;

           if  ls /= l + 1  then
              Test := Test + Abs E(ls-1);
           end if;

           if  Arg_1_Is_Negligible_Respect_Arg_2 (Abs S(ls), Test)  then
              S(ls) := Zero;
              exit;
           end if;

        end loop;

        if  ls = l  then
           kase := 3;
        elsif  ls = mn  then
           kase := 1;
        else
           kase := 2;
           l := ls;
        end if;

     end if;

     l := l + 1;

     if  kase = 1  then

        Deflate_Negligible_S_Vals:
        declare
           sn_lo, sn_hi, cs_lo, cs_hi, hypot : Real;
           v0, v1 : Real;
           cs_hi_is_one : Boolean;
        begin

           f       := E(mn-1);
           E(mn-1) := Zero;

           for k in reverse l .. mn-1 loop

              --  hypot :=        Sqrt(f*f + S(k)*S(k));
              --  sn    := f    / Sqrt(f*f + S(k)*S(k));
              --  cs    := S(k) / Sqrt(f*f + S(k)*S(k));

              Rotation_Factors
                (f, S(k), sn_lo, sn_hi, cs_lo, cs_hi, hypot, cs_hi_is_one);

              S(k) := hypot;

              if  k /= l  then
                 f      := -sn_hi * E(k-1) - sn_lo * E(k-1);
                 E(k-1) :=  cs_hi * E(k-1) + cs_lo * E(k-1);
              end if;

              if  Matrix_V_Desired  then
              if cs_hi_is_one then                 -- cs_hi = +/-1, sn_hi = 0
                 for r in Starting_Col .. Final_Col loop
                    v0 := V(r, k);
                    v1 := V(r, mn);
                    V(r, k)  :=  cs_hi*v0 + ( cs_lo*v0 + sn_lo*v1);
                    V(r, mn) :=  cs_hi*v1 + (-sn_lo*v0 + cs_lo*v1);
                 end loop;
              else                                 -- cs_hi = 0, |sn_hi| = 1
                 for r in Starting_Col .. Final_Col loop
                    v0 := V(r, k);
                    v1 := V(r, mn);
                    V(r, k)  :=  sn_hi*v1 + ( cs_lo*v0 + sn_lo*v1);
                    V(r, mn) := -sn_hi*v0 + (-sn_lo*v0 + cs_lo*v1);
                 end loop;
              end if;
              end if;

           end loop;
        end Deflate_Negligible_S_Vals;

     elsif  kase = 2  then

        Split_at_Negligible_S_of_l:
        declare
           sn_lo, sn_hi, cs_lo, cs_hi, hypot : Real;
           u0, u1 : Real;
           cs_hi_is_one : Boolean;
        begin

           f      := E(l-1);
           E(l-1) := Zero;

           for k in l .. mn loop

              --  hypot :=        Sqrt(f*f + S(k)*S(k));
              --  sn    := f    / Sqrt(f*f + S(k)*S(k));
              --  cs    := S(k) / Sqrt(f*f + S(k)*S(k));

              Rotation_Factors
                (f, S(k), sn_lo, sn_hi, cs_lo, cs_hi, hypot, cs_hi_is_one);

              S(k) := hypot;

              f    := -sn_lo * E(k) - sn_hi * E(k);
              E(k) :=  cs_lo * E(k) + cs_hi * E(k);

              if  Matrix_U_Desired  then
              declare
                 kc : constant Row_Index := Row_Index(k);
                 lc : constant Row_Index := Row_Index(l-1);
              begin
                 if cs_hi_is_one then              -- cs_hi = +/-1, sn_hi = 0
                    for r in Starting_Row .. Final_Row loop -- rotate cols k, l-1 in U
                       u0 := U(r, kc);
                       u1 := U(r, lc);
                       U(r, kc) :=  cs_hi*u0 + ( cs_lo*u0 + sn_lo*u1);
                       U(r, lc) :=  cs_hi*u1 + (-sn_lo*u0 + cs_lo*u1);
                    end loop;
                 else                              -- cs_hi = 0, |sn_hi| = 1
                    for r in Starting_Row .. Final_Row loop -- rotate cols k, l-1 in U
                       u0 := U(r, kc);
                       u1 := U(r, lc);
                       U(r, kc) :=  sn_hi*u1 + ( cs_lo*u0 + sn_lo*u1);
                       U(r, lc) := -sn_hi*u0 + (-sn_lo*u0 + cs_lo*u1);
                    end loop;
                 end if;
              end;
              end if;

           end loop;
        end Split_at_Negligible_S_of_l;

     --  Case 3, Perform one QR step.

     elsif  kase = 3  then

        Get_Shift:
        declare
           sm, smm1, emm1, Scale, sl, el, Beta, b, Hyp, Shift : Real;
        begin

           Scale := Real'Max (Abs S(mn),   Real'Max (Abs S(mn-1),
                    Real'Max (Abs E(mn-1), Real'Max (Abs S(l), Abs E(l)))));

           Scale := (Scale + Min_Allowed_Real); --  Scale >= 0

           Scale := Two ** (-Real'Exponent (Scale));

           sl   := S(l)    * Scale;
           el   := E(l)    * Scale;
           sm   := S(mn)   * Scale;
           smm1 := S(mn-1) * Scale;
           emm1 := E(mn-1) * Scale;
           Beta := ((smm1 + sm) * (smm1 - sm) + emm1*emm1) * Half;
           b    := sm * emm1;

           Shift := Zero;

           if  Abs Beta > Zero  or  Abs b > Zero  then
              Hyp := Hypo.Hypotenuse (Beta, b);
              if  Beta < Zero  then  Hyp := -Hyp;  end if;
              Shift := b * (b / (Beta + Hyp));
           end if;

         --Eig_1 := smm1**2 + emm1**2 + Shift;
         --Eig_2 := sm**2 - Shift;             -- e(mn) = negligable
         --f := sl**2 - Eig_2;

           f := (sl + sm) * (sl - sm) + Shift;
           g := sl * el;

        end Get_Shift;


        Chase_Zeros:
        declare
           sn_lo, sn_hi, cs_lo, cs_hi, hypot : Real;
           v0, v1, u0, u1, tmp : Real;
           cs_hi_is_one : Boolean;
        begin
           for k in l .. mn-1 loop

              --  sn := g / sqrt (g*g + f*f) = g / hypot
              --  cs := f / sqrt (g*g + f*f) = f / hypot

              Rotation_Factors
                (g, f, sn_lo, sn_hi, cs_lo, cs_hi, hypot, cs_hi_is_one);

              if  k /= l  then
                 E(k-1) := hypot;
              end if;

              tmp    :=       (cs_lo*S(k) + sn_lo*E(k));
              f      := tmp  + cs_hi*S(k) + sn_hi*E(k);
              tmp    :=       (cs_lo*E(k) - sn_lo*S(k));
              E(k)   := tmp  + cs_hi*E(k) - sn_hi*S(k);

              g      := sn_lo * S(k+1) + sn_hi * S(k+1);
              S(k+1) := cs_lo * S(k+1) + cs_hi * S(k+1);

              if  Matrix_V_Desired  then
              if cs_hi_is_one then              -- cs_hi = (+/-)1, sn_hi = 0
                 for r in Starting_Col .. Final_Col loop
                    v0 := V(r, k);
                    v1 := V(r, k+1);
                    V(r, k)   :=  cs_hi*v0 + ( cs_lo*v0 + sn_lo*v1);
                    V(r, k+1) :=  cs_hi*v1 + (-sn_lo*v0 + cs_lo*v1);
                 end loop;
              else                              -- cs_hi = 0, sn_hi = (+/-)1
                 for r in Starting_Col .. Final_Col loop
                    v0 := V(r, k);
                    v1 := V(r, k+1);
                    V(r, k)   :=  sn_hi*v1 + ( cs_lo*v0 + sn_lo*v1);
                    V(r, k+1) := -sn_hi*v0 + (-sn_lo*v0 + cs_lo*v1);
                 end loop;
              end if;
              end if;

              Rotation_Factors
                (g, f, sn_lo, sn_hi, cs_lo, cs_hi, hypot, cs_hi_is_one);

              S(k) := hypot;

              tmp    :=         cs_lo*E(k) + sn_lo*S(k+1);
              f      := tmp   + cs_hi*E(k) + sn_hi*S(k+1);
              tmp    :=        -sn_lo*E(k) + cs_lo*S(k+1);
              S(k+1) := tmp    -sn_hi*E(k) + cs_hi*S(k+1);

              g      := sn_lo * E(k+1) + sn_hi * E(k+1);
              E(k+1) := cs_lo * E(k+1) + cs_hi * E(k+1);

              if  Matrix_U_Desired and Row_Index (k) < Row_Index'Last then
              declare
                 k0 : constant Row_Index := Row_Index(k);
                 k1 : constant Row_Index := Row_Index(k+1);
              begin
                 if cs_hi_is_one then            -- |cs_hi| = (+/-)1, sn_hi = 0
                    for r in Starting_Row .. Final_Row loop
                       u0 := U(r, k0);
                       u1 := U(r, k1);
                       U(r, k0) :=  cs_hi*u0 + ( cs_lo*u0 + sn_lo*u1);
                       U(r, k1) :=  cs_hi*u1 + (-sn_lo*u0 + cs_lo*u1);
                    end loop;
                 else                                 -- cs_hi = 0, |sn_hi| = 1
                    for r in Starting_Row .. Final_Row loop
                       u0 := U(r, k0);
                       u1 := U(r, k1);
                       U(r, k0) :=  sn_hi*u1 + ( cs_lo*u0 + sn_lo*u1);
                       U(r, k1) := -sn_hi*u0 + (-sn_lo*u0 + cs_lo*u1);
                    end loop;
                 end if;
              end;
              end if;

           end loop; --  k in l .. mn-1

           E(mn-1) := f;
           No_of_Iterations := No_of_Iterations + 1;

        end Chase_Zeros; -- Completed 1 QR step.

     --  Case 4, Convergence:

     elsif  kase = 4  then

        --  Make the singular value nonnegative.

        if  S(l) < Zero  then
           S(l) := -S(l);
           if  Matrix_V_Desired  then
              for r in Starting_Col .. Final_Col loop
                 V(r,l) := -V(r,l);
              end loop;
           end if;
        end if;

        --  Order singular values largest 1st to smallest last in array S:

        Sort_Singular_Vals: loop

           if  l = Final_Col  then
             exit Sort_Singular_Vals;
           end if;

           if  S(l+1) <= S(l)  then
             exit Sort_Singular_Vals;
           end if;

           tmp    := S(l);
           S(l)   := S(l+1);
           S(l+1) := tmp;

           if  Matrix_V_Desired and l < Final_Col  then
              for r in Starting_Col .. Final_Col loop
                 tmp       := V(r, l);
                 V(r, l)   := V(r, l+1);
                 V(r, l+1) := tmp;
              end loop;
           end if;

           if  Matrix_U_Desired and Row_Index (l) < Row_Index'Last  then
           declare
              l0 : constant Row_Index := Row_Index(l);
              l1 : constant Row_Index := Row_Index(l+1);
           begin
              for r in Starting_Row .. Final_Row loop
                 tmp      := U(r, l0);
                 U(r, l0) := U(r, l1);
                 U(r, l1) := tmp;
              end loop;
           end;
           end if;

           l := l + 1;

        end loop Sort_Singular_Vals;

        No_of_Iterations := 0;

        if mn = Starting_Col then  exit QR_Iteration;  else  mn := mn - 1; end if;

     end if;  --  kase 1 through 4

     end loop QR_Iteration; --  while mn > Starting_Col

  end QR_For_Singular_Vals;

  --------------------
  --  SVD_Decompose --
  --------------------

  procedure SVD_Decompose
    (A                : in out A_Matrix;                    --  m x n
     U                :    out U_matrix;                    --  m x m
     V                :    out V_matrix;                    --  n x n
     S                :    out Singular_Vector;
     Id_of_1st_S_Val  :    out Col_Index;
     Starting_Col     : in     Col_Index;        -- Starting_Row is set to Starting_Col
     Final_Col        : in     Col_Index;
     Final_Row        : in     Row_Index;
     Matrix_U_Desired : in     Boolean   := True;
     Matrix_V_Desired : in     Boolean   := True)
  is
     Starting_Row : constant Row_Index := Row_Index (Starting_Col);
     E : Singular_Vector;
     Bidiag_Final_Col : Col_Index;
  begin

     --  Init out parameters:

     Id_of_1st_S_Val := Starting_Col;  --  modified only if convergence fails.

     --  If U and V are undesirable, then don't fill up memory with Zero's.
     --  Sometimes that frees up memory.

     if Matrix_U_Desired then
        U := (others => (others => Zero));
        for r in Row_Index loop
           U(r, r) := One;
        end loop;
     end if;

     if Matrix_V_Desired then
        V := (others => (others => Zero));
        for c in Col_Index loop
           V(c, c) := One;
        end loop;
     end if;

     E := (others => Zero);
     S := (others => Zero);

     if Int64 (Row_Index'Last) - Int64 (Row_Index'First) < 1 then
        raise Constraint_Error with "Matrix too small";
     end if;

     if Int64 (Col_Index'Last) - Int64 (Col_Index'First) < 1 then
        raise Constraint_Error with "Matrix too small";
     end if;

     -- Must have Starting_Col = Starting_Row. Just to make things less error
     -- prone, let's require the same of the Matrix Index types:

     if Int64 (Col_Index'First) /= Int64 (Row_Index'First) then
        raise Constraint_Error with "Must have Col_Index'First = Row_Index'First";
     end if;

     -- Can SVD square matrices only, or matrices with more Rows than Columns:

     if Int64 (Final_Row) < Int64 (Final_Col) then
        raise Constraint_Error with "Can't have Final_Row < Final_Col";
     end if;

     -- Can SVD matrices with 3 or more Columns:

     if Int64 (Final_Col) < Int64 (Starting_Col) + 2 then
        raise Constraint_Error with "Can't have Final_Col < Starting_Col + 2";
     end if;

     --  Reduce A to bidiagonal form, storing the diagonal elements
     --  in S and the super-diagonal elements in E.

     Householder_Bidiagonalize
       (A,                                -- A becomes the B in A = U * B * V'
        U, V,                             -- out params
        E, S,                             -- out params
        Bidiag_Final_Col,
        Final_Row,
        Starting_Col,
        Final_Col,
        Matrix_U_Desired,
        Matrix_V_Desired);

     --  Iterate for the singular values.

     QR_For_Singular_Vals
       (E, S,
        U, V,
        Id_of_1st_S_Val,
        Starting_Col, Bidiag_Final_Col,
        Starting_Row, Final_Row,
        Matrix_U_Desired, Matrix_V_Desired);

  end SVD_Decompose;

end Golub_SVD;
