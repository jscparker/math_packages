
---------------------------------------------------------------------------
-- package body Test_Matrices
-- Copyright (C) 2008-2018 Jonathan S. Parker.
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

with Ada.Numerics.Discrete_Random;
--with Ada.Numerics.Generic_Elementary_Functions;
--with Ada.Numerics;

package body Test_Matrices is

  Zero : constant Real := +0.0;
  Half : constant Real := +0.5;
  One  : constant Real := +1.0;
  Two  : constant Real := +2.0;
  Gamma : constant Real := 
    +0.57721_56649_01532_86060_65120_90082_40243_10421_59335_93992;

  Max_Allowed_Real : constant Real := Two ** (Real'Machine_Emax - 20);

  --package Math is new Ada.Numerics.Generic_Elementary_Functions(Real);
  --use Math;

  --  for Random nums:
 
  type Unsigned32 is mod 2**32;
 
  package Discrete_32_bit is new Ada.Numerics.Discrete_Random (Unsigned32);
 
  Random_Stream_id : Discrete_32_bit.Generator;
 
  type Bits_per_Ran is range 1 .. 32;

  ----------------
  -- Symmetrize --
  ----------------

  procedure Symmetrize
    (A : in out Matrix)
  is
     Result : Matrix;
  begin
     for r in Index loop
     for c in Index loop
        Result(r,c) := 0.5 * (A(r,c) + A(c,r));
     end loop;
     end loop;
     A := Result;
  end Symmetrize;

  ---------------
  -- Transpose --
  ---------------

  procedure Transpose 
    (A              : in out Matrix;
     Starting_Index : in     Index     := Index'First;
     Max_Index      : in     Index     := Index'Last)
  is
     tmp : Real;
  begin
     for r in Starting_Index .. Max_Index-1 loop
     for c in r+1 .. Max_Index loop
        tmp     := A(r, c);
        A(r, c) := A(c, r);
        A(c, r) := tmp;
     end loop;
     end loop;
  end Transpose;

  -----------------
  -- Real_Random --
  -----------------
 
  -- if you use  No_of_Bits = 1, then get 0.5 and 0.0.
  -- if you use  No_of_Bits = 2, then get 0.75, 0.5, 0.25, 0.0.
 
  function Real_Random (No_of_Bits : Bits_per_Ran) return Real is
    N : constant Integer := Integer (No_of_Bits);
  begin
 
    return Two**(-N) * Real (Discrete_32_bit.Random (Random_Stream_id) / 2**(32-N));
 
  end Real_Random;
 
  function "+"(M : Matrix; Right : Real) return Matrix is
     Result : Matrix;
  begin
     for r in Index loop
     for c in Index loop
        Result(r, c) :=  M(r, c) + Right;
     end loop;
     end loop;
     return Result;
  end "+";
 
  ----------------
  -- Get_Zielke --
  ----------------
 
  procedure Get_Zielke
    (M              : out Matrix;
     Starting_Index : in Index;
     Max_Index      : in Index;
     Z              : in Real)
   is
     x : constant Real := One;
     y : constant Real := Two;
     Start   : constant Index'Base := Starting_Index;
     Finish  : constant Index'Base := Max_Index;
     M_Order : constant Real := Real (Finish) - Real (Start) + One;
     Test : Real;
   begin
 
     for i in Start .. Finish loop
     for j in Start .. Finish loop
        Test := (Real (i) - Real (Start) + One) + (Real (j) - Real (Start) + One);
        if i = j then
           if Test <= M_Order then
              M(i, j) := x + y + z;
           elsif Test < Two*M_Order then
              M(i, j) := x + z;
           else
              M(i, j) := x - y + z;
           end if;
        else
           if Test <= M_Order then
              M(i, j) := x + y;
           else
              M(i, j) := x;
           end if;
        end if;
     end loop;
     end loop;
 
  end Get_Zielke;
 
  ---------------------
  -- Get_Companion_B --
  ---------------------
 
  -- Eigvals are zeros of poly w/ coefficients in 1st Col:
  --
  --    P(x) = x**n + C(1)*x^(n-1) +  .. . + C(n-1)*x^1 + C(n)*x^0
  --
  -- where  
  --
  --    M(i, Starting_Index) = -C(1+i-Starting_Index)
  --
  -- with
  --
  --    i in Starting_Index+0 .. Starting_Index+(n-1)

  procedure Get_Companion_B
    (M              : out Matrix;
     Starting_Index : in Index;
     Max_Index      : in Index;
     B              : in Real := 1.1892071150027210667175) -- Sqrt_of_Sqrt_of_2
   is
     Start_I : constant Integer := Integer (Starting_Index);
     Exp : Integer;
     Pow : constant := 1; -- 1 is stronger than 3.
   --Sqrt_of_2 : constant := 1.4142135623730950488017;
  begin
     M := (others => (others => Zero));

     for Row in Starting_Index .. Max_Index-1 loop
        M(Row, Row+1) := One;
     end loop;

     for Row in Starting_Index .. Max_Index loop
        Exp := Integer'Min (Integer(Row) - Start_I, Real'Machine_Emax/Pow-9);
        M(Row, Starting_Index) := Two - B ** Exp;
     end loop;

  end Get_Companion_B;
 
  --------------------------
  -- Get_Pas_Fib_Rescaled --
  --------------------------
 
  procedure Get_Pas_Fib_Rescaled
    (M              :     out Matrix;
     Starting_Index : in      Index := Index'First;
     Max_Index      : in      Index := Index'Last)
  is
     M_max, Scale_Factor : Real;
   --Exp : constant Integer := Real'Machine_Emin;  -- roughly -1010 usually
   --Smallest_Allowed_Real : constant Real := Two**(Exp/2 - Exp/8);
  begin
     M := (others => (others => One));
 
     if Integer (Max_Index) - Integer (Starting_Index) > 0 then
     for Row in Starting_Index+1 .. Max_Index loop
 
        M_max := Zero;
        for Col in Starting_Index+1 .. Max_Index loop
           M(Row, Col) := M(Row-1, Col-1) + M(Row-1, Col);
           if abs M(Row, Col) > M_max then M_max := abs M(Row, Col); end if;
        end loop;
        -- M_max gets larger each step.
 
        -- Scale the Matrix to peak of about 1:
 
        if M_max > Max_Allowed_Real then
           Scale_Factor := One / M_max;
           for r in Starting_Index .. Max_Index loop
           for c in Starting_Index .. Max_Index loop
              M(r, c) := Scale_Factor * M(r, c);
           end loop;
           end loop;
           M_max := One;
        end if;
 
     end loop;
     end if;
 
--     for r in M'Range(1) loop
--     for c in M'Range(2) loop
--        if Abs M(r,c) < Smallest_Allowed_Real * M_max then
--           M(r,c) := Smallest_Allowed_Real * M_max;
--        end if;
--     end loop;
--     end loop;
 
  end Get_Pas_Fib_Rescaled;
 
  -------------------------
  -- Get_Pascal_Rescaled --
  -------------------------
 
  procedure Get_Pascal_Rescaled
    (M              :     out Matrix;
     Starting_Index : in      Index := Index'First;
     Max_Index      : in      Index := Index'Last)
  is
     M_max, Scale_Factor : Real;
   --Exp : constant Integer := Real'Machine_Emin;  -- roughly -1010 usually
   --Smallest_Allowed_Real : constant Real := Two**(Exp/2 - Exp/8);
  begin
     M := (others => (others => Zero));
 
     for Col in Index loop
       M(Col, Col) := One;
     end loop;
 
     for Row in Starting_Index .. Max_Index loop
       M(Row, Starting_Index) := One;
     end loop;
 
     if Integer (Max_Index) - Integer (Starting_Index) > 1 then
     for Row in Starting_Index+2 .. Max_Index loop
 
        M_max := Zero;
        for Col in Starting_Index+1 .. Row-1 loop
           M(Row, Col) := (M(Row-1, Col-1) + M(Row-1, Col));
           if abs M(Row, Col) > M_max then M_max := abs M(Row, Col); end if;
        end loop;
 
        -- rescale the whole Matrix:
 
        if M_max > Max_Allowed_Real then
           Scale_Factor := One / M_max;
           for r in Starting_Index .. Max_Index loop
           for c in Starting_Index .. Max_Index loop
              M(r, c) := Scale_Factor * M(r, c);
           end loop;
           end loop;
           M_max := One;
        end if;
 
     end loop;
     end if;
 
--     for r in M'Range(1) loop
--     for c in M'Range(2) loop
--        if Abs M(r,c) < Smallest_Allowed_Real * M_max then
--           M(r,c) := Smallest_Allowed_Real * M_max;
--        end if;
--     end loop;
--     end loop;
 
  end Get_Pascal_Rescaled;

  -----------------
  -- Get_Redheff --
  -----------------
 
  procedure Get_Redheff
    (M              :     out Matrix;
     Starting_Index : in      Index := Index'First;
     Max_Index      : in      Index := Index'Last)
  is
     function I_divides_J (I, J : Integer) return Boolean is
     begin
        if (J / I) * I = J then
           return True;
        else
           return False;
        end if;
     end I_divides_J;
     i, j : Integer;
  begin
     M := (others => (others => Zero));

     for r in Starting_Index .. Max_Index loop
     for c in Starting_Index .. Max_Index loop
        i := Integer (c) - Integer (Starting_Index) + 1;
        j := Integer (r) - Integer (Starting_Index) + 1;
        if I_divides_J (i, j) or j = 1 then
           M(r, c) := One;
        else
           M(r, c) := Zero;
        end if;
     end loop;
     end loop;

  end Get_Redheff;

  ------------------
  -- Get_Sampling --
  ------------------
 
  procedure Get_Sampling
    (M              :     out Matrix;
     Starting_Index : in      Index := Index'First;
     Max_Index      : in      Index := Index'Last;
     Basic          : in      Boolean := True)
  is
     Denom, Sum : Real;
     X : array (Index) of Real := (others => Zero);
     Half_Way : constant Integer := (Integer(Max_Index) - Integer(Starting_Index)) / 2;

     Emin             : constant Integer := Real'Machine_Emin;
     Min_Allowed_Real : constant Real := Two ** (Emin - Emin/16);
  begin
     M := (others => (others => Zero));

     if Basic then
        for i in Starting_Index .. Max_Index loop
           X(i) := Real (i) - Real (Starting_Index) + One;
        end loop;
     else
        X(Starting_Index) := Two ** Integer'Min (Half_Way, Abs (Emin) - 2);
        for i in Starting_Index+1 .. Max_Index loop
           if X(i-1) > Min_Allowed_Real then
              X(i) := X(i-1) * Half;
           else
              X(i) := X(i-1) * 701.0; -- prevent 0's in extreme limit
           end if;
        end loop;
     end if;

     for r in Starting_Index .. Max_Index loop
     for c in Starting_Index .. Max_Index loop
        Denom := X(r) - X(c);
        if Abs Denom < Min_Allowed_Real then Denom := Min_Allowed_Real; end if;
        if r /= c then
           M(r, c) := X(r) / Denom;
        end if;
     end loop;
     end loop;

     for c in Starting_Index .. Max_Index loop
        Sum := Zero;
        for r in Starting_Index .. Max_Index loop
           if r /= c then
              Sum := Sum + M(r, c);
           end if;
        end loop;
        M(c, c) := Sum;
     end loop;

  end Get_Sampling;

  ------------------
  -- Get_Laguerre --
  ------------------
 
  -- Eig(i) = (-1)**(i-1) / (i-1)!
  -- Upper triangular(???), so Eigs are just diagonal elements.
  -- Notice the TRANSPOSE at end

  procedure Get_Laguerre
    (M              :     out Matrix;
     Starting_Index : in      Index := Index'First;
     Max_Index      : in      Index := Index'Last)
  is
   --N : Index'Base := Max_Index - Starting_Index + 1;
   --Scale : Real :=  Two ** (- Integer(N) / 16);
     i, j : Integer;
  begin
     M := (others => (others => Zero));

     M(Starting_Index, Starting_Index) := One;
   
     if Max_Index = Starting_Index then
        return;
     end if;
   
     M(Starting_Index+1, Starting_Index)   :=  One;
     M(Starting_Index+1, Starting_Index+1) := -One;
   
     if Max_Index = Starting_Index+1 then
        return;
     end if;
   
     for r in Starting_Index+2 .. Max_Index loop
     for c in Starting_Index   .. Max_Index loop

        i := Integer (r) - Integer (Starting_Index) + 1;
        j := Integer (c) - Integer (Starting_Index) + 1;
   
        if j = 1 then
          M(r,c) := (Real (2*i - 3) * M(r-1,c)
                   - Real (  i - 2) * M(r-2,c))
                   / Real (  i - 1);
        else
          M(r,c) := (Real (2*i - 3) * M(r-1,c) - M(r-1,c-1)
                   - Real (  i - 2) * M(r-2,c))
                   / Real (  i - 1);
        end if;

        --if Abs M(r,c) <  Two ** (Real'Machine_Emin - Real'Machine_Emin / 8) then
        --   M(r,c) := Zero;
        --end if;

     end loop;
     end loop;
 
     Transpose (M);

  end Get_Laguerre; 

  ----------------
  -- Get_Lotkin --
  ----------------
 
  procedure Get_Lotkin
    (M              :    out Matrix;
     Starting_Index : in     Index := Index'First;
     Max_Index      : in     Index := Index'Last)
  is
     i, j : Integer;
     Denominator : Real;
     Prime_Factors : constant Real := (+166966608033225.0)*(+580027.0) * Two**(-68);
  begin
     M := (others => (others => Zero));

     -- Prime_Factors=3.0*3.0*3.0*5.0*5.0*7.0*11.0*13.0*17.0*19.0*23.0*29.0*31.0*37.0;
     -- so  Prime_Factors / D  is exactly represented in 15 digit floating point
     -- up to D = 39 (allowing an 20x20 matrix).  Prime_Factors = 166966608033225.0
     -- Prime_Factors := +166966608033225.0 * Two**(-47); -- bring it near to One
  
     for r in Starting_Index .. Index'Last loop
     for c in Starting_Index .. Index'Last loop

        i := Integer (r) - Integer (Starting_Index) + 1; 
        j := Integer (c) - Integer (Starting_Index) + 1; 

        Denominator :=  Real(i + j) - One;
        M(r, c)     :=  Prime_Factors / Denominator;
 
     end loop;
     end loop;

     for r in Starting_Index .. Max_Index loop
        M(r, Starting_Index) := One;
     end loop;
 
  end Get_Lotkin;

  -----------------
  -- Get_Hilbert --
  -----------------
 
  procedure Get_Hilbert
    (M              :    out Matrix;
     Starting_Index : in     Index := Index'First)
  is
     i, j : Integer;
     Denominator : Real;
     Prime_Factors : constant Real := (+166966608033225.0)*(+580027.0) * Two**(-68);
  begin
     M := (others => (others => Zero));

     -- Prime_Factors=3.0*3.0*3.0*5.0*5.0*7.0*11.0*13.0*17.0*19.0*23.0*29.0*31.0*37.0;
     -- so  Prime_Factors / D  is exactly represented in 15 digit floating point
     -- up to D = 39 (allowing an 20x20 matrix).  Prime_Factors = 166966608033225.0
     -- Prime_Factors := +166966608033225.0 * Two**(-47); -- bring it near to One
  
     for r in Starting_Index .. Index'Last loop
     for c in Starting_Index .. Index'Last loop

        i := Integer (r) - Integer (Starting_Index) + 1; 
        j := Integer (c) - Integer (Starting_Index) + 1; 

        Denominator :=  Real(i + j) - One;
        M(r, c)     :=  Prime_Factors / Denominator;
 
     end loop;
     end loop;

  end Get_Hilbert;

  -------------------
  -- Get_Ding_Dong --
  -------------------
 
  procedure Get_Ding_Dong
    (M              :    out Matrix;
     Starting_Index : in     Index := Index'First;
     Max_Index      : in     Index := Index'Last)
  is
     i, j : Integer;
     Denominator : Real;
  begin
     M := (others => (others => Zero));

     for r in Starting_Index .. Index'Last loop
     for c in Starting_Index .. Index'Last loop

        i := Integer (r) - Integer (Starting_Index) + 1; 
        j := Integer (c) - Integer (Starting_Index) + 1; 

        -- Can use Prime_Factors := 166966608033225.0 * Two**(-47) to
        -- bring it near exact integer elements, but lose the Pi valued eigs.
  
        Denominator := (Real (Max_Index) - Real (Starting_Index) + One) - Real (i + j) + 1.5;
        M(r, c) := One / Denominator;
 
     end loop;
     end loop;

  end Get_Ding_Dong;

  -----------------
  -- Get_Gregory --
  -----------------
 
  procedure Get_Gregory
    (M              :    out Matrix;
     Starting_Index : in     Index := Index'First;
     Max_Index      : in     Index := Index'Last)
  is
     i, j : Integer;
  begin
     M := (others => (others => Zero));

     for r in Starting_Index .. Max_Index loop
     for c in Starting_Index .. Max_Index loop

        i := Integer (r) - Integer (Starting_Index) + 1; 
        j := Integer (c) - Integer (Starting_Index) + 1; 
        if  r = Max_Index  then
           M(r,c) := Real (j);
        elsif  c = Max_Index  then
           M(r,c) := Real (i);
        elsif  r = c  then
           M(r,c) := One;
        end if;
 
     end loop;
     end loop;

  end Get_Gregory;

  ---------------
  -- Get_Chow3 --
  ---------------

  -- full rather than lower Hessenberg. Non-Symmetric.

  procedure Get_Chow3
    (M              :     out Matrix;
     Starting_Index : in      Index := Index'First;
     Max_Index      : in      Index := Index'Last;
     Alpha          : in      Real  := Gamma)        -- 1.0 is harder
  is
     i_r, j_c : Integer;
  begin
     M := (others => (others => Zero));

     -- lower Hessenberg   (r+1 >= c)

     for r in Starting_Index .. Max_Index loop
     for c in Starting_Index .. Max_Index loop
        i_r := Integer (r) - Integer (Starting_Index) + 1;
        j_c := Integer (c) - Integer (Starting_Index) + 1;
        M(r, c) := Alpha ** (i_r + 1 - j_c);
        if Abs M(r,c) <  Two ** (Real'Machine_Emin - Real'Machine_Emin / 8) then
           M(r, c) := Zero;
        end if;
     end loop;
     end loop;


     declare 
        Beta : constant Real := Zero;
     begin
        for r in Starting_Index .. Max_Index loop
           M(r, r) := M(r, r) + Beta;
        end loop;
     end;

  end Get_Chow3;

  --------------
  -- Get_Chow --
  --------------
 
  -- Described by T S Chow; has eigs 4*alpha*cos(k*pi/(n+2))^2
  -- and Floor [N / 2] zeros (if no diagonal addends).

  procedure Get_Chow
    (M              :     out Matrix;
     Starting_Index : in      Index := Index'First;
     Max_Index      : in      Index := Index'Last;
     Alpha          : in      Real  := One;
     Beta           : in      Real  := Zero)
  is
     i, j : Integer;
  begin
     M := (others => (others => Zero));

     -- lower Hessenberg   (r+1 >= c)

     for c in Starting_Index .. Max_Index loop
     for r in Starting_Index .. Max_Index loop
        i := Integer (r) - Integer (Starting_Index) + 1;
        j := Integer (c) - Integer (Starting_Index) + 1;
        if i + 1 >= j then
           M(r, c) := Alpha ** (i + 1 - j);
        end if;
        if Abs M(r,c) <  Two ** (Real'Machine_Emin - Real'Machine_Emin / 8) then
           M(r, c) := Zero;
        end if;
     end loop;
     end loop;

     for r in Starting_Index .. Max_Index loop
        M(r, r) := M(r, r) + Beta;
      --M(r, r) := M(r, r) + 2.0**32;
     end loop;

  end Get_Chow;

  ----------------
  -- Get_Lehmer --
  ----------------
 
  procedure Get_Lehmer
    (M              :     out Matrix;
     Starting_Index : in      Index := Index'First;
     Max_Index      : in      Index := Index'Last)
  is
     i, j : Integer;
  begin
     M := (others => (others => Zero));

     for r in Starting_Index .. Max_Index loop
     for c in Starting_Index .. Max_Index loop
        i := Integer (r) - Integer (Starting_Index) + 1;
        j := Integer (c) - Integer (Starting_Index) + 1;
        M(r, c) := Real (Integer'Min (i, j)) / Real (Integer'Max (i, j));
     end loop;
     end loop;

  end Get_Lehmer;

  --------------
  -- Get_Lesp --
  --------------

  --  From John Burkardt:
  --
  --    The eigenvalues are real, and smoothly distributed in [-2*N-3.5, -4.5].
  --
  --    The eigenvalues are sensitive.
  --
  --    The matrix is similar to the symmetric tridiagonal matrix with
  --    the same diagonal entries and with off-diagonal entries 1,
  --    via a similarity transformation using the diagonal matrix
  --    D = diagonal ( 1!, 2!,  .. ., N! ).
  --
  procedure Get_Lesp
    (M              :    out Matrix;
     Starting_Index : in     Index := Index'First;
     Max_Index      : in     Index := Index'Last)
  is
     i, j : Integer;
  begin
     M := (others => (others => Zero));

     for r in Starting_Index .. Max_Index loop
     for c in Starting_Index .. Max_Index loop
        i := Integer (r) - Integer (Starting_Index) + 1;
        j := Integer (c) - Integer (Starting_Index) + 1;
        if (i - j) = 1 then
           M(r, c) := One / Real (i);
        elsif i = j then
           M(r, c) := -Real (2*i + 3);
        elsif (i - j) = -1 then
           M(r, c) := Real (j);
        else
           M(r, c) := Zero;
        end if;
     end loop;
     end loop;

     Transpose (M);

  end Get_Lesp;

  ----------------
  -- Get_Trench --
  ----------------

  procedure Get_Trench
    (M              :    out Matrix;
     Starting_Index : in     Index := Index'First;
     Max_Index      : in     Index := Index'Last;
     Alpha          : in     Real)
  is
     i, j : Integer;
  begin
     M := (others => (others => Zero));

     for r in Starting_Index .. Max_Index loop
     for c in Starting_Index .. Max_Index loop
        i := Integer (r) - Integer (Starting_Index) + 1;
        j := Integer (c) - Integer (Starting_Index) + 1;
        if i = j then
           M(r, c) := Alpha;
        else
           M(r, c) := Two ** (-Abs (i - j) - 1);
        end if;
     end loop;
     end loop;
  end Get_Trench;

  -----------------
  -- Init_Matrix --
  -----------------

  procedure Init_Matrix
    (M              : out Matrix;
     Desired_Matrix : in  Matrix_Id := Easy_Matrix;
     Starting_Index : in  Index     := Index'First;
     Max_Index      : in  Index     := Index'Last)
  is
     A, OffFactor : Real;

  begin
 
    if Max_Index <= Starting_Index then
       raise Constraint_Error with "Can't have Max_Index <= Starting_Index";
    end if;

    M := (others => (others => Zero));     --essential init.
    --the 0.0 is essential.
 
    case Desired_Matrix is

    when Laguerre =>

       Get_Laguerre (M, Starting_Index, Max_Index); 

    when Lesp =>

       Get_Lesp (M, Starting_Index, Max_Index); 

    when Lehmer =>

       Get_Lehmer (M, Starting_Index, Max_Index); 

    when Chow =>

       Get_Chow (M, Starting_Index, Max_Index, One, Zero); 

    when Chow1 =>

       Get_Chow (M, Starting_Index, Max_Index, -1.05, Zero); 

    when Chow2 =>

       Get_Chow (M, Starting_Index, Max_Index, Gamma, -One);

    when Chow3 =>

       Get_Chow3 (M, Starting_Index, Max_Index, 1.05); 

    when Redheff =>

       Get_Redheff (M, Starting_Index, Max_Index); 

    when Sampling =>

       Get_Sampling (M, Starting_Index, Max_Index); 

    when Sampling_1 =>

       Get_Sampling (M, Starting_Index, Max_Index, False); 

    when Gregory =>

       Get_Gregory (M, Starting_Index, Max_Index); 
 
    when Anti_Hadamard_Upper_Tri =>
 
       M := (others => (others => Zero));     --essential init.
 
       for Col in Index'First+1 .. Index'Last loop
       for Row in Index'First .. Col-1 loop
          if ((Integer(Col) + Integer(Row)) mod 2) = 1 then
             M(Row, Col) := One;
          else
             M(Row, Col) := Zero;
          end if;
       end loop;
       end loop;
 
       for Col in Index loop
          M(Col, Col) :=  One;
       end loop;
 
       M := M + Matrix_Addend;
 
    when Anti_Hadamard_Lower_Tri =>
 
       M := (others => (others => Zero));     --essential init.
 
       for Col in Index'First .. Index'Last-1 loop
       for Row in Col+1 .. Index'Last loop
          if ((Integer(Col) + Integer(Row)) mod 2) = 1 then
             M(Row, Col) := One;
          else
             M(Row, Col) := Zero;
          end if;
       end loop;
       end loop;
 
       for Col in Index loop
          M(Col, Col) :=  One;
       end loop;
 
       M := M + Matrix_Addend;
 
    when Symmetric_Banded =>
 
       OffFactor := 2.101010101;
 
       for I in Index loop
         M(I, I) := One;
       end loop;
 
       declare Col : Index; begin

       for BottomDiagonal in Starting_Index+1 .. Index'Last loop
       for Row in BottomDiagonal .. Index'Last loop
         Col := Index (Integer(Row) - Integer(BottomDiagonal) + Integer(Starting_Index));
         M(Row, Col) :=
            OffFactor * (Real(BottomDiagonal) - Real(Starting_Index) + 1.0);
       end loop;
       end loop;

       end;
 
       for Row in Starting_Index+1 .. Index'Last loop
       for Col in Starting_Index .. Row-1 loop
          M(Col, Row) :=  M(Row, Col);
       end loop;
       end loop;
 
    when Easy_Matrix =>
 
       OffFactor := 0.10101010101;
 
       for I in Index loop
         M(I, I) := 1.010101010101;
       end loop;

       declare Col : Index; begin

       for BottomDiagonal in Starting_Index+1 .. Index'Last loop
       for Row in BottomDiagonal .. Index'Last loop
         Col := Index (Integer(Row) - Integer(BottomDiagonal) + Integer(Starting_Index));
         M(Row, Col)
          := 0.031 * (Real(Row) - Real(Starting_Index) + One) / Real(Max_Index)
             + OffFactor / (Real(BottomDiagonal) - Real(Starting_Index));
       end loop;
       end loop;

       end;
 
       for Row in Starting_Index+1 .. Index'Last loop
       for Col in Starting_Index .. Row-1 loop
          M(Col, Row) :=  M(Row, Col) + 0.333;
       end loop;
       end loop;
 
    when Small_Diagonal  =>
 
       OffFactor := 101010.01;
 
       for I in Index loop
         M(I, I) := 0.01010101010101010101;
       end loop;
 
       declare Col : Index; begin

       for BottomDiagonal in Starting_Index+1 .. Index'Last loop
       for Row in BottomDiagonal .. Index'Last loop
         Col := Index (Integer(Row) - Integer(BottomDiagonal) + Integer(Starting_Index));
         M(Row, Col)
          := 0.013 * (Real(Row) - Real(Starting_Index) + 1.0) / Real(Max_Index)
              + OffFactor / (Real(BottomDiagonal) - Real (Starting_Index));
       end loop;
       end loop;

       end;
 
       for Row in Starting_Index+1 .. Index'Last loop
       for Col in Starting_Index .. Row-1 loop
          M(Col, Row) :=  M(Row, Col) + 0.333; -- the 333 illconditions it
       end loop;
       end loop;
 
    when Pascal_Col_Scaled  =>
 
       Get_Pascal_Rescaled (M, Starting_Index, Max_Index);
 
       declare
          Max, Scaling : Real;
          Exp : Integer;
          Col_Norm : array(Index) of Real;
       begin
          for Col in Starting_Index .. Max_Index loop
             Max := Zero;
             for Row in Starting_Index .. Max_Index loop
                if not M(Row, Col)'Valid then raise constraint_error; end if;
                if Abs M(Row, Col) > Max then Max := Abs M(Row, Col); end if;
             end loop;
             Col_Norm(Col) := Max;
          end loop;
 
          for Col in Starting_Index .. Max_Index loop
             Exp := Real'Exponent (Col_Norm(Col));
             if Exp < Real'Machine_Emin then
                Exp := Real'Machine_Emin + 16;
             elsif Exp > Real'Machine_Emax then
                Exp := Real'Machine_Emax - 16;
             end if;
             Scaling := Two ** (-Exp);
             for Row in Starting_Index .. Max_Index loop
                M(Row, Col) :=  M(Row, Col) * Scaling;
             end loop;
          end loop;
       end;
 
       M := M + Matrix_Addend;
 
    when Pascal_Row_Scaled  =>
 
       Get_Pascal_Rescaled (M, Starting_Index, Max_Index);
 
       declare
           Max, Scaling : Real;
           Exp : Integer;
           Row_Norm : array(Index) of Real := (others => Zero);
       begin
         for Row in Starting_Index .. Max_Index loop
            Max := Zero;
            for Col in Starting_Index .. Max_Index loop
               if not M(Row, Col)'Valid then raise constraint_error; end if;
               if Abs M(Row, Col) > Max then Max := Abs M(Row, Col); end if;
            end loop;
            Row_Norm(Row) := Max;
         end loop;
 
         for Row in Starting_Index .. Max_Index loop
            Exp := Real'Exponent (Row_Norm(Row));
            if Exp < Real'Machine_Emin then
               Exp := Real'Machine_Emin + 16;
            elsif Exp > Real'Machine_Emax then
               Exp := Real'Machine_Emax - 16;
            end if;
            Scaling := Two ** (-Exp);
            for Col in Starting_Index .. Max_Index loop
               M(Row, Col) :=  M(Row, Col) * Scaling;
            end loop;
         end loop;
       end;
 
       M := M + Matrix_Addend;
 
    when Pas_Fib =>
 
       Get_Pas_Fib_Rescaled (M, Starting_Index, Max_Index);
 
    when Pascal_Symmetric  =>
 
       Get_Pascal_Rescaled (M, Starting_Index, Max_Index);
 
       for Row in Starting_Index .. Max_Index-1 loop
         for Col in Row+1 .. Max_Index loop
            M(Row, Col) :=  M(Col, Row);
         end loop;
       end loop;
 
    when Pascal  =>
 
       Get_Pascal_Rescaled (M, Starting_Index, Max_Index);
 
       M := M + Matrix_Addend;
 
    when Forsythe_Symmetric  =>
 
       M := (others => (others => Zero));
 
       A := Two**(-8);
 
       for I in Index loop
          M(I, I) := A;
       end loop;
 
       for Col in Starting_Index+1 .. Index'Last loop
          M(Col, Col-1) := One;
       end loop;
 
       for Col in Starting_Index+1 .. Index'Last loop
          M(Col-1, Col) := One;
       end loop;
 
       M(Index'First, Max_Index) := One;
       M(Max_Index, Index'First) := One;
 
    when Forsythe_0 =>
 
       M := (others => (others => Zero));
 
       A := Zero;
 
       for I in Index loop
          M(I, I) := A;
       end loop;
 
       for Col in Starting_Index+1 .. Index'Last loop
         M(Col-1, Col) := One;
       end loop;
 
       M(Max_Index, Starting_Index) := One;
 
    when Forsythe_1 =>
 
       M := (others => (others => Zero));
 
       A := Two**(0);
 
       for I in Index loop
          M(I, I) := A;
       end loop;
 
       for Col in Starting_Index+1 .. Index'Last loop
         M(Col-1, Col) := One;
       end loop;
 
       M(Max_Index, Starting_Index) := One;
 
    when Zero_Diagonal  =>
 
       M := (others => (others => Zero));
 
      --for I in Index loop
      --   M(I, I) := Two**(-63);
      --end loop;
 
       for Col in Starting_Index+1 .. Index'Last loop
         M(Col-1, Col) := One;
       end loop;
 
       for Col in Starting_Index .. Index'Last-1 loop
         M(Col+1, Col) := One;
       end loop;
 
    when Ring_Adjacency_0  =>
 
       M := (others => (others => Zero));
 
       M(Starting_Index, Max_Index) := One;
       M(Max_Index, Starting_Index) := One;
 
       for Col in Starting_Index+1 .. Max_Index loop
          M(Col-1, Col) := One;
       end loop;
 
       for Col in Starting_Index .. Max_Index-1 loop
          M(Col+1, Col) := One;
       end loop;
 
    when Ring_Adjacency_1  =>
 
       M := (others => (others => Zero));
 
       M(Starting_Index, Max_Index) :=  One;
       M(Max_Index, Starting_Index) := -One;
 
       for Col in Starting_Index+1 .. Max_Index loop
          M(Col-1, Col) := One;
       end loop;
 
       for Col in Starting_Index .. Max_Index-1 loop
          M(Col+1, Col) := One;
       end loop;
 
    when Wilkinson_Plus_2I  =>
 
       declare
          Mid : constant Real := Half * (Real (Max_Index) - Real (Starting_Index) + One);
       begin 
          for I in Index loop
             M(I, I) := Abs (Mid - (Real (I) - Real (Starting_Index))) + Two;
          end loop;
       end;
 
       for Col in Index'First+1 .. Index'Last loop
         M(Col-1, Col) := One;
       end loop;
 
       for Col in Index'First .. Index'Last-1 loop
         M(Col+1, Col) := One;
       end loop;
 
    when Wilkinson_Plus  =>
 
       declare
          Mid : constant Real := Half * (Real (Max_Index) - Real (Starting_Index) + One);
       begin 
          for I in Index loop
             M(I, I) := Abs (Mid - (Real (I) - Real (Starting_Index)));
          end loop;
       end;
 
       for Col in Index'First+1 .. Index'Last loop
          M(Col-1, Col) := One;
       end loop;
 
       for Col in Index'First .. Index'Last-1 loop
          M(Col+1, Col) := One;
       end loop;
 
    when Wilkinson_Minus  =>

       declare
          Mid : constant Real := Half * (Real (Max_Index) - Real (Starting_Index) + One);
       begin 
          for I in Index loop
             M(I, I) := Mid - (Real (I) - Real (Starting_Index));
          end loop;
       end;
 
       for Col in Index'First+1 .. Index'Last loop
          M(Col-1, Col) := One;
       end loop;
 
       for Col in Index'First .. Index'Last-1 loop
          M(Col+1, Col) := One;
       end loop;
 
    when QR_Test =>
 
       for Col in Starting_Index+1 .. Index'Last loop
          M(Col-1, Col) := One;
       end loop;
 
       for Col in Starting_Index .. Index'Last-1 loop
          M(Col+1, Col) := One;
       end loop;

       if Index'Last > Starting_Index+1 then
          for i in Starting_Index .. Index'Last-1 loop
             if (Integer (i) - Integer (Starting_Index)) mod 4 = 1 then
                M(i, i+1) := -320.0 * Real'Epsilon;
                M(i+1, i) :=  320.0 * Real'Epsilon;
             end if;
          end loop;
       end if;
 
    when Upper_Ones  =>
 
       for Row in Starting_Index .. Index'Last loop
       for Col in Row .. Index'Last loop
          M(Row, Col) := One;
       end loop;
       end loop;
 
       for Col in Starting_Index .. Index'Last loop
          M(Col, Col) := One;
       end loop;
 
       M := M + Matrix_Addend;
 
    when Lower_Ones  =>
 
       for Row in Starting_Index .. Index'Last loop
       for Col in Row .. Index'Last loop
          M(Col, Row) := One;
       end loop;
       end loop;
 
       for Col in Starting_Index .. Index'Last loop
          M(Col, Col) := One;
       end loop;
 
       M := M + Matrix_Addend;
 
    when Lower_Integers =>
 
       for Row in Starting_Index .. Index'Last loop
       for Col in Row .. Index'Last loop
          M(Col, Row) := One;
       end loop;
       end loop;
 
       for Col in Starting_Index .. Index'Last loop
          M(Col, Col) := Real (Col) - Real (Starting_Index) + One;
       end loop;
 
       M := M + Matrix_Addend;
 
    when Zero_Cols_and_Rows  =>
 
       M := (others => (others => Zero));
 
       if Real (Starting_Index) + Two <= Real (Max_Index) then
          for Row in Starting_Index+2 .. Max_Index loop
          for Col in Row .. Max_Index loop
             M(Col, Row) := One;
          end loop;
          end loop;
       end if;
 
       M := M + Matrix_Addend;

    when Frank_0  =>   -- upper Hessenberg

       -- order > 2:
       -- Max_Index - Starting_Index + 1 > 2
 
       for Row in Starting_Index .. Max_Index loop
       for Col in Starting_Index .. Max_Index loop
          M(Row, Col) := Real (Max_Index) - Real (Starting_Index) + One
                      - (Real (Col) - Real (Starting_Index));
       end loop;
       end loop;
 
       for Col in Starting_Index .. Max_Index-2 loop
       for Row in Col+2 .. Max_Index loop
          M(Row, Col) := Zero;
       end loop;
       end loop;
 
       for Row in Starting_Index+1 .. Max_Index loop
          M(Row, Row-1) := M(Row-1, Row);
       end loop;

       M := M + Matrix_Addend;
 
    when Frank_1  =>   -- upper Hessenberg
 
       declare
          i, j : Real;
       begin
          for Row in Starting_Index .. Max_Index loop
          for Col in Starting_Index .. Max_Index loop
             i := Real (Row) - Real (Starting_Index) + 1.0;
             j := Real (Col) - Real (Starting_Index) + 1.0;
             M(Row, Col) := Real (Max_Index) - Real (Starting_Index) + 1.0
                         - (Real'Min (i, j) - 1.0);
          end loop;
          end loop;
    
          for Col in Index'Base (Starting_Index) .. Index'Base (Max_Index)-2 loop
          for Row in Col+2 .. Max_Index loop
             M(Row, Col) := Zero;
          end loop;
          end loop;
       end;
    
       M := M + Matrix_Addend;
 
    when Frank_2  =>
 
       for Row in Starting_Index .. Index'Last loop
       for Col in Starting_Index .. Index'Last loop
          M(Row, Col) := Real(Index'Min(Col,Row)) - Real(Starting_Index) + One;
       end loop;
       end loop;
 
    when Moler  =>
 
       for Row in Starting_Index .. Index'Last loop
       for Col in Starting_Index .. Index'Last loop
          M(Row, Col) :=  (Real(Index'Min(Col,Row)) -  Real(Starting_Index) - One);
       end loop;
       end loop;
 
       for Col in Starting_Index .. Index'Last loop
          M(Col, Col) :=  (Real(Col) -  Real(Starting_Index) + One);
       end loop;
 
    when Random_32_bit  =>
 
       Discrete_32_bit.Reset (Random_Stream_id);
 
       for Row in Starting_Index .. Max_Index loop
       for Col in Starting_Index .. Max_Index loop
          M(Row, Col) := Real_Random(32);
       end loop;
       end loop;
 
    when Random_1_bit_anti  =>
 
       Discrete_32_bit.Reset (Random_Stream_id);
 
       for Row in Starting_Index .. Index'Last loop
       for Col in Starting_Index .. Row loop
          M(Row, Col) := Real_Random(1) - Half;
       end loop;
       end loop;
 
       for Row in Starting_Index .. Index'Last loop
       for Col in Row .. Index'Last loop
          M(Row, Col) := -M(Col, Row);
       end loop;
       end loop;
 
       for Row in Index loop
          M(Row, Row) := Zero;
       end loop;
 
    when Random_1_bit  =>
 
       Discrete_32_bit.Reset (Random_Stream_id);
 
       for Row in Starting_Index .. Max_Index loop
       for Col in Starting_Index .. Max_Index loop
          M(Row, Col) := Real_Random(1) - Half;
       end loop;
       end loop;
 
    when Ding_Dong  =>

       Get_Ding_Dong (M,  Starting_Index, Max_Index);
 
    when Clustered  =>
 
       Clustered_Eigs:
       declare
          Value : Real := One;
       begin
          for Row in Starting_Index .. Max_Index loop
          for Col in Starting_Index .. Max_Index loop
             M(Row, Col) := Value;
             Value := Value + One;
          end loop;
          end loop;
       end Clustered_Eigs;
 
    when Hilbert  =>

       Get_Hilbert (M, Starting_Index);
 
    when Lotkin  =>
 
       Get_Lotkin (M, Starting_Index, Max_Index); 
 
    when Fiedler_0  =>
 
       declare
          C : array(Index) of Real;
          Length : constant Real := (Real (Max_Index) - Real (Starting_Index) + One);
          Scale  : constant Real := Two**(-Real'Exponent (Length) + 8);
       begin
 
          for i in Starting_Index .. Max_Index loop
             C(i) := Scale / (Real(i) - Real(Starting_Index) + One)**2;
          end loop;
 
          for Row in Starting_Index .. Max_Index loop
          for Col in Starting_Index .. Max_Index loop
             M(Col, Row) := Abs (C(Row) - C(Col));
          end loop;
          end loop;
 
       end;
 
    when Fiedler_1  =>   -- not Fiedler's matrix
 
       M := (others => (others => Zero));

       declare
          C : array(Index) of Real;
          Alpha : constant Real := One;
       begin
 
          for i in Starting_Index .. Max_Index loop
             C(i) := Real(i) - Real(Starting_Index) + One;
          end loop;
 
          for Row in Starting_Index .. Max_Index loop
          for Col in Starting_Index .. Row loop
             M(Row, Col) := Abs (C(Row) - C(Col));
             M(Col, Row) := Alpha;
          end loop;
          end loop;

--          for Col in Starting_Index .. Max_Index loop
--             M(Col, Col) := 1.0e-5;
--          end loop;
       end;
 
    when U_Hard  =>
 
       -- hard on calculation of  U, V
 
       M := (others => (others => One));
 
       for Col in Index loop
        --M(Starting_Index, Col) := Zero;
          M(Col, Starting_Index) := Zero;
       end loop;
 
    when Companion_2  =>
 
       -- Eigvals are zeros of poly w/ coefficients in Final Row of M:
       --   P(x) = x**n - M(n-1)*x**(n-1) -  .. . - M(0)
       -- where M(i) = M(Max_Index, Starting_Index + i).
       -- Eig_Vectors are (1, lambda, lambda^2,  .. . lambda^(n-1)).
       -- M is not diagonizable if there exist multiple roots.

       Get_Companion_2:
       declare
          Start_I : constant Integer := Integer (Starting_Index);
          Exp, Exp0 : Integer;
          Pow : constant := 3;
       begin
          M := (others => (others => Zero));
 
          for Row in Starting_Index .. Max_Index-1 loop
             M(Row, Row+1) := One;
          end loop;

          for Col in Starting_Index .. Max_Index loop
             Exp0 := Integer (Max_Index) - Integer (Col) + Integer (Starting_Index);
             Exp := Integer'Min (Exp0 - Start_I, Real'Machine_Emax/Pow-9);
             M(Max_Index, Col) := Two - 1.01**Exp;
          end loop;

       end Get_Companion_2;
 
    when Companion_1  =>

       --Sqrt_of_2         : constant := 1.4142135623730950488017;
       --Sqrt_of_Sqrt_of_2 : constant := 1.1892071150027210667175;
       --
       -- Uses default: B := Sqrt_of_Sqrt_of_2

       Get_Companion_B (M, Starting_Index, Max_Index, 1.0625);
 
    when Companion_0  =>
 
       -- Eigvals are zeros of poly w/ coefficients in 1st Col:
       --
       --    P(x) = x**n + C(1)*x^(n-1) +  .. . + C(n-1)*x^1 + C(n)*x^0
       --
       -- where  
       --
       --    M(i, Starting_Index) = -C(1+i-Starting_Index)
       --
       -- with
       --
       --    i in Starting_Index+0 .. Starting_Index+(n-1)

       Get_Companion_0:
       declare
          C : Array (Index) of Real := (others => Zero);
       begin
          M := (others => (others => Zero));

          for i in Index loop
             C(i) :=  Real (i) - Real (Starting_Index) - One;
          end loop;
 
          for r in Index loop
             M(r, Starting_Index) := -C(r);
          end loop;

          for r in Index'First .. Index'Last-1 loop
             M(r, r+1) := One;
          end loop;

       end Get_Companion_0;
 
    when Companion  =>
 
       -- eigvals are zeros of poly w/ coefficients in 1st Col:
       --   P(x) = x**n + M(1)*x**(n-1) +  .. . + M(n)
       -- here M = n, n-1,  .. .

       Get_Companion:
       begin
       M := (others => (others => Zero));
 
       for r in Starting_Index .. Max_Index-1 loop
          M(r, r+1) := One;
       end loop;

       for r in Starting_Index .. Max_Index loop
          M(r, Starting_Index) := One/Real(-Integer(r)+Integer(Starting_Index)-1)**5;
       end loop;

       end Get_Companion;
 
    when Gear_0  =>
 
       M := (others => (others => Zero));
 
       for Col in Index'First .. Index'Last-1 loop
          M(Col+1, Col) := One;
       end loop;
 
       for Col in Index'First+1 .. Index'Last loop
          M(Col-1, Col) := One;
       end loop;
 
       M(Max_Index, Starting_Index) := -(One + Two**(-8));
       M(Starting_Index, Max_Index) :=  (One + Two**(-8));
       -- the "-" makes it almost singular if M = One
 
    when Gear_1  =>
 
       M := (others => (others => Zero));
 
       for Col in Index'First .. Index'Last-1 loop
          M(Col+1, Col) := One;
       end loop;
 
       for Col in Index'First+1 .. Index'Last loop
          M(Col-1, Col) := One;
       end loop;
 
       M(Max_Index, Starting_Index) := -One;
       M(Starting_Index, Max_Index) :=  One;
       -- the "-" makes it almost singular if M = One
 
    when Fibonacci  =>
 
       -- golden mean eigval p = (1+sqrt(5))/2 (all same val); maximally defective;
       -- Has eigenvalue p with eigenvector (1, p, p^2,  .. ., p^(N-1)).
 
       M := (others => (others => Zero));
 
       for Col in Index loop
         M(Col, Col) := One;
       end loop;
 
       for Col in Index'First .. Index'Last-1 loop
         M(Col+1, Col) := One;
       end loop;
 
       M(Starting_Index, Starting_Index)   := Zero;
       M(Starting_Index, Starting_Index+1) := One;
 
       M := M + Matrix_Addend;
 
    when Peters  =>  -- one small eig., lower triangular
 
       for Col in Starting_Index .. Max_Index loop
       for Row in Starting_Index .. Col loop
          M(Col, Row) := -One;
       end loop;
       end loop;
 
       for Col in Index loop
          M(Col, Col) := One;
       end loop;
 
       M := M + Matrix_Addend;
 
    when Peters_0  =>  -- one small eig., upper triangular, defeats Gauss. Elim.
 
       for Col in Starting_Index .. Max_Index loop
       for Row in Starting_Index .. Col loop
          M(Row, Col) := -One;
       end loop;
       end loop;
 
       for Col in Index loop
          M(Col, Col) := One;
       end loop;
 
       M := M + Matrix_Addend;
 
    when Combin  =>
 
       M := (others => (others => Half));
 
       for Col in Index loop
          M(Col, Col) := Half + Half**11;
       end loop;
 
    when Peters_1  =>  -- transpose of peters_2 matrix. row pivoting likes this:
 
       for Row in Starting_Index .. Max_Index loop
       for Col in Row .. Max_Index loop
          M(Row, Col) := -One;
       end loop;
       end loop;
 
       for Col in Starting_Index .. Max_Index loop
          M(Max_Index, Col) := One;
       end loop;
 
       for Col in Starting_Index .. Max_Index loop
          M(Col, Col) := One;
       end loop;
 
    when Peters_2  =>  -- Real peters matrix. row pivoting hates this:
 
       declare
       -- Eigval : constant Real := 2.0;
          Eigval : constant Real := Zero;
       begin
          for Row in Starting_Index .. Index'Last loop
          for Col in Starting_Index .. Row        loop
             M(Row, Col) := -One + Eigval; -- Eigval is one of the eigs
          end loop;
          end loop;
       end;
 
       for Row in Starting_Index .. Index'Last loop
        --M(Row, Max_Index) := One - 1.0E-9; -- kills the error growth
          M(Row, Max_Index) := One;
       end loop;
 
       for Col in Starting_Index .. Index'Last loop
          M(Col, Col) := One;
       end loop;
 
    when Diag_Test =>
 
       M := (others => (others => Zero));  -- Essential init. The 0.0 is essential.
 
       for Col in Starting_Index+1 .. Max_Index loop
          M(Col-1, Col) := One;
       end loop;
 
       for Col in Starting_Index .. Max_Index loop
          M(Col, Col) := 0.5;
       end loop;
 
       M := M + Matrix_Addend;
 
    when Upper_Tri_K  =>
 
       M := (others => (others => Zero));
       --the 0.0 is essential.
 
       declare
          c : constant Real := Two;  -- higher condition num;
        --c := constant Half;
       begin
 
          for Row in Starting_Index .. Max_Index loop
          for Col in  Row .. Max_Index      loop
             M(Row, Col) := -c;
          end loop;
          end loop;
 
          for Col in Starting_Index .. Max_Index loop
             M(Col, Col) := One;
          end loop;
 
       end;
 
       M := M + Matrix_Addend;
 
    when Lower_Tri_K  =>
 
       M := (others => (others => Zero));
       --the 0.0 is essential.
 
       declare
          c : constant Real := Two;  -- higher condition num;
        --c := constant Half;
       begin
 
          for Row in Starting_Index .. Max_Index loop
          for Col in Starting_Index .. Row       loop
             M(Row, Col) := -c;
          end loop;
          end loop;
 
          for Col in Starting_Index .. Max_Index loop
             M(Col, Col) := One;
          end loop;
 
       end;
 
       M := M + Matrix_Addend;
 
    when Kahan_Row_Scaled  =>  --  make max element of Kahan Row = 1
 
       M := (others => (others => Zero));
       --the 0.0 is essential.
 
       declare
          c : Real;
       begin
 
          c := 0.70710678118654752440;
 
          for Row in Starting_Index .. Max_Index loop
          for Col in  Row .. Max_Index      loop
             M(Row, Col) := -c;
          end loop;
          end loop;
 
          for Col in Starting_Index .. Max_Index loop
             M(Col, Col) := One;
          end loop;
 
       end;
 
       M := M + Matrix_Addend;
 
    when Kahan  =>  --  need s**2 + c**2 = 1
 
       M := (others => (others => Zero));     --essential init.
       --the 0.0 is essential.
 
       declare
          s, c, s_i : Real;
          i : Integer;
       begin
 
          s := 0.70710678118654752440;
          c := 0.70710678118654752440;
 
          for Row in Starting_Index .. Max_Index loop
          for Col in  Row .. Max_Index      loop
             M(Row, Col) := -c;
          end loop;
          end loop;
 
          for Col in Starting_Index .. Max_Index loop
             M(Col, Col) := One;
          end loop;
 
          i := 0;
          for Row in Starting_Index .. Max_Index loop
             s_i := s**i;
             i   := i + 1;
             for Col in Starting_Index .. Max_Index loop
                M(Row, Col) := s_i * M(Row, Col);
             end loop;
          end loop;
 
       end;
 
    when Kahan_Col_Scaled  =>  --  make length of Col = 1
 
       M := (others => (others => Zero));     --essential init.
       --the 0.0 is essential.
 
       declare
          s, c, s_i : Real;
          i : Integer;
          Col_Norm : array(Index) of Real;
       begin
 
          s := 0.70710678118654752440;
          c := 0.70710678118654752440;
 
          for Row in Starting_Index .. Max_Index loop
          for Col in  Row .. Max_Index      loop
             M(Row, Col) := -c;
          end loop;
          end loop;
          for Col in Starting_Index .. Max_Index loop
            M(Col, Col) := One;
          end loop;
 
          i := 0;
          for Row in Starting_Index .. Max_Index loop
             s_i := s**i;
             i   := i + 1;
             for Col in Starting_Index .. Max_Index loop
                M(Row, Col) := s_i * M(Row, Col);
             end loop;
          end loop;
 
          -- Got Kahan, now normalize the cols:
 
          for Col in Starting_Index .. Max_Index loop
             s_i := Zero;
             for Row in Starting_Index .. Max_Index loop
                s_i := s_i + Abs M(Row, Col);
             end loop;
             Col_Norm(Col) := S_i;
          end loop;
 
          for Col in Starting_Index .. Max_Index loop
            for Row in Starting_Index .. Max_Index loop
               M(Row, Col) :=  M(Row, Col) / Col_Norm(Col);
            end loop;
          end loop;
 
       end;
 
       M := M + Matrix_Addend;
 
    when Kahan_2  =>  --  need s**2 + c**2 = 1; make max element of Col = 1
 
       M := (others => (others => Zero));     --essential init.
       --the 0.0 is essential.
 
       declare
          s, c, s_i : Real;
          i : Integer;
       begin
 
          s := 0.25;
          c := 0.96824583655185422129;
 
          for Row in Starting_Index .. Max_Index loop
          for Col in  Row .. Max_Index      loop
             M(Row, Col) := -c;
          end loop;
          end loop;
 
          for Col in Starting_Index .. Max_Index loop
             M(Col, Col) := One;
          end loop;
 
          i := 0;
          for Row in Starting_Index .. Max_Index loop
             s_i := s**i;
             i   := i + 1;
             for Col in Starting_Index .. Max_Index loop
                M(Row, Col) := s_i * M(Row, Col);
             end loop;
          end loop;
 
          for Col in Starting_Index+1 .. Max_Index loop
             for Row in Starting_Index .. Max_Index loop
                M(Row, Col) := M(Row, Col) / (c);
             end loop;
          end loop;
 
       end;
 
       M := M + Matrix_Addend;
 
    when Vandermonde  =>
 
       -- small elements underflow; max element near ~1.
 
       Vandermondes_Matrix:
       declare
          Exp, X_Count : Integer;
          X : Real;
          Half_No_Of_Rows : constant Integer
                     := (Integer(Max_Index) - Integer(Starting_Index) + 1)/2;
          B : constant Real := One / Real (Half_No_Of_Rows);
          A : constant Real := Two ** (Real'Exponent(B) - 1);
          --Smallest_Allowed_Real : constant Real := Two ** (Real'Machine_Emin + 512);
       begin
          for Row in Starting_Index .. Max_Index loop
          for Col in Starting_Index .. Max_Index loop
             Exp         := Integer(Col) - Integer(Starting_Index);
             X_Count     := Integer(Row) - Integer(Starting_Index);
             X           := A * (Real (X_Count - Half_No_Of_Rows));
             M(Row, Col) := X ** Exp;
          end loop;
          end loop;
 
--          for r in M'Range(1) loop
--          for c in M'Range(2) loop
--             if not M(r,c)'Valid  or
--                Abs M(r,c) < Smallest_Allowed_Real
--             then
--                M(r,c) := Smallest_Allowed_Real;
--             end if;
--          end loop;
--          end loop;

       end Vandermondes_Matrix;
 
    when Trench_1 =>
 
       -- Alpha = -One; -- singular submatrices

       Get_Trench (M, Starting_Index, Max_Index, Alpha => -One);
  
    when Trench =>
 
       -- Alpha = One; -- highest condition number

       Get_Trench (M, Starting_Index, Max_Index, Alpha => One);
  
    when Zielke_0 =>
 
       -- Z = 1.0; -- Singular
 
       Get_Zielke (M, Starting_Index, Max_Index, Z => One);
 
    when Zielke_1 =>
 
       -- Z = -1.0; -- ill conditioned
 
       Get_Zielke (M, Starting_Index, Max_Index, Z => -One);
 
    when Zielke_2 =>
 
       Get_Zielke (M, Starting_Index, Max_Index, (+2.0));
 
    when All_Zeros  =>
 
       M := (others => (others => Zero));
 
    when All_Ones  =>
 
       M := (others => (others => One));
 
    end case;
 
  end Init_Matrix;
 
end Test_Matrices;
