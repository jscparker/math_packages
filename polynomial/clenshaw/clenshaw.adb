
-----------------------------------------------------------------------
-- package body Clenshaw. Generates functions from recurrance relations.
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

package body Clenshaw is

   Zero : constant Real := +0.0;

   -----------------
   -- Evaluate_Qs --
   -----------------

   --  The recurrance relation for the Q's at X is easily written as matrix
   --  equation.  In the following, f(X) is the given function for Q_0:
   --
   --  Q(0)  = Q_0(p,X);
   --  Q(1)  = 0 + Alpha_1*Q_0;
   --  Q(2)  = 0 + Alpha_2*Q_1 + Beta_2*Q_0;
   --  Q(3)  = 0 + Alpha_3*Q_2 + Beta_3*Q_1;
   --  ...
   --  Q(N)  = 0 + Alpha_N*Q_N-1 + Beta_N*Q_N-2
   --
   -- In matrix form, M*Q = (f(X), 0, 0, ...) , this becomes:
   --
   --   | 1        0        0        0 |  |Q(0)|     | Q_0(p,X)|    | C_0 |
   --   | E_1      1        0        0 |  |Q(1)|  =  | 0       |  = | C_1 |
   --   | B_2      E_2      1        0 |  |Q(2)|     | 0       |    | C_2 |
   --   | 0        B_3      E_3      1 |  |Q(3)|     | 0       |    | C_3 |
   --
   --  where E_m = -Alpha_m, B_m = -Beta_m.
   --
   --  So Q = M_inverse * C is the desired solution, but there may be numerical
   --  error in the calculation of M_inverse by back-substitution.  The
   --  solution vector Q can be improved numerically by iterative refinement
   --  via Newton's method:
   --
   --             Q_new = Q_old + M_inverse * (C - M*Q_old)
   --
   --  where Q = M_inverse * C is the calculation of Q given at the top.
   --
   procedure Evaluate_Qs
     (X                : in     Real;
      Q                : in out Poly_Values;
      Max_Poly_ID      : in     Poly_ID_Type;
      P                : in     Real          := 0.0;
      No_Of_Iterations : in     Positive      := 1)
   is
      Product, Del : Poly_Values;
      m            : Poly_ID_Type;
   begin
      --
      -- Step 0. Initialize excess values of Q to Zero.
      --
      if Max_Poly_ID < Poly_ID_Type'Last then
         for m in Max_Poly_ID+1 .. Poly_ID_Type'Last loop
            Q(m) := Zero;
         end loop;
      end if;

      --
      -- Step 0. Want zeroth order poly  Q_0(p,X).  No work to do.
      --
      if Max_Poly_ID = Poly_ID_Type'First then
         m      := Poly_ID_Type'First;
         Q(m)   := Q_0(p,X);
      end if;

      --
      -- Step 0b. Poly is 1st order.  Almost no work to do.
      -- Don't do any iteration.
      --
      if Max_Poly_ID = Poly_ID_Type'First + 1 then
         m      := Poly_ID_Type'First;
         Q(m)   := Q_0(p,X);
         m      := Poly_ID_Type'First+1;
         Q(m)   := Alpha(m,p,X) * Q(m-1);
      end if;

      --
      -- Step 1.  We now know that Max_Poly_ID > 1.
      -- Start by getting starting value of Q by solving M*Q = f.
      -- Use recurrence relation to get Q at X.
      -- Start with special formulas for the 1st two Q's:
      --
      if Max_Poly_ID > Poly_ID_Type'First + 1 then

         m      := Poly_ID_Type'First;
         Q(m)   := Q_0(p,X);
         m      := Poly_ID_Type'First+1;
         Q(m)   := Alpha(m,p,X) * Q(m-1);

         for m in Poly_ID_Type'First+2..Max_Poly_ID loop
	    Q(m) :=   Alpha(m,p,X) * Q(m-1)
	            +  Beta(m,p,X) * Q(m-2);
         end loop;

         --
         -- Step 2. Improve Q numerically through Newton iteration.
         --             Q_new = Q_old + M_inverse * (C - M*Q_old)
         --
         Iterate: for Iter in 2..No_Of_Iterations loop

           --  Get Product = M*Q_old:

           m            := Poly_ID_Type'First;
           Product(m)   := Q(m);
           m            := Poly_ID_Type'First+1;
           Product(m)   := Q(m) - Alpha(m,p,X)*Q(m-1);

           for m in Poly_ID_Type'First+2..Max_Poly_ID loop
              Product(m) :=  Q(m) - Alpha(m,p,X)*Q(m-1)
                                  -  Beta(m,p,X)*Q(m-2);
           end loop;

           --  Get Residual = (Q_0(p,X), 0, ... , 0) - M*D_old.  Reuse the
           --  array Product to hold the value of Residual:

           Product(Poly_ID_Type'First) := Zero;
           --  Residual is always exactly 0.0 here

           for m in Poly_ID_Type'First+1 .. Max_Poly_ID loop
              Product(m) := - Product(m);
           end loop;

           --  Get Del = M_inverse * (C - M*Q_old) = M_inverse * Product:

           m      := Max_Poly_ID;
           Del(m) := Product(m);

           m      := Max_Poly_ID - 1;
           Del(m) := Product(m) + Alpha(m,p,X)*Del(m-1);

           for m in Poly_ID_Type'First+2 .. Max_Poly_ID loop
              Del(m) := Product(m) + Alpha(m,p,X)*Del(m-1)
                                   +  Beta(m,p,X)*Del(m-2);
           end loop;

           --  Get Q_new = Q_old + Del;

           for m in Poly_ID_Type'First..Max_Poly_ID loop
              Q(m) := Q(m) + Del(m);
           end loop;

        end loop Iterate;

      end if;

   end Evaluate_Qs;

   -------------
   -- M_times --
   -------------

   --  M is Upper-Triangular.
   --  The elements of M are 1 down the diagonal, and -Alpha(m,p,X) and
   --  -Beta(m,p,X) down the the off-diagonals.
   --
   function M_times
     (D         : in Coefficients;
      X         : in Real;
      P         : in Real;
      Sum_Limit : in Poly_ID_Type)
      return Coefficients
   is
      Product : Coefficients;
      m       : Poly_ID_Type;
   begin

      --  These inits are amazingly slow!
      for m in Sum_Limit .. Poly_ID_Type'Last loop
          Product(m) := Zero;
      end loop;

      --  Get Product = M*D:

      m          := Sum_Limit;
      Product(m) := D(m);

      if Sum_Limit > Poly_ID_Type'First then
         m          := Sum_Limit - 1;
         Product(m) := D(m) - Alpha(m+1,p,X)*D(m+1);
      end if;

      if Sum_Limit > Poly_ID_Type'First+1 then
         for m in Poly_ID_Type'First .. Sum_Limit-2 loop
            Product(m) := D(m) - Alpha(m+1,p,X)*D(m+1)
                               -  Beta(m+2,p,X)*D(m+2);
         end loop;
      end if;

      return Product;

   end M_times;

   pragma Inline (M_times);

   ---------------------
   -- M_inverse_times --
   ---------------------

   --  M is Upper-Triangular so solution is by back-substitution.
   --  The elements of M are 1 down the diagonal, and -Alpha and
   --  -Beta down the off-diagonals.
   --
   function M_inverse_times
     (C         : in Coefficients;
      X         : in Real;
      P         : in Real;
      Sum_Limit : in Poly_ID_Type)
      return         Coefficients
   is
      Result : Coefficients;
      m      : Poly_ID_Type;
   begin

      --  These inits are amazingly slow!
      for m in Sum_Limit .. Poly_ID_Type'Last loop
          Result(m) := Zero;
      end loop;

      m         := Sum_Limit;
      Result(m) := C(m);

      if Sum_Limit > Poly_ID_Type'First then
         m         := Sum_Limit - 1;
         Result(m) := C(m) + Alpha(m+1,p,X) * Result(m+1);
      end if;

      if Sum_Limit > Poly_ID_Type'First+1 then
         for m in reverse Poly_ID_Type'First .. Sum_Limit-2 loop
            Result(m) := C(m) + Alpha(m+1,p,X) * Result(m+1)
	                      +  Beta(m+2,p,X) * Result(m+2);
         end loop;
      end if;

      return Result;

   end M_inverse_times;

   pragma Inline (M_inverse_times);

   ---------
   -- Sum --
   ---------

   --  This is easily written as matrix equation, with Sum = D(0):
   --
   --  D_n   = C_n;
   --  D_n-1 = C_n-1 + Alpha_n*D_n;
   --  D_n-2 = C_n-2 + Alpha_n-1*D_n-1 + Beta_n-2*D_n-2;
   --  ...
   --  D_1   = C_1   + Alpha_2*D_2     + Beta_3*D_3
   --  D_0   = C_0   + Alpha_1*D_1     + Beta_2*D_2
   --
   -- In matrix form, M*D = C, this becomes:
   --
   --
   --   | 1        E_1      B_2      0 |  |D(0)  |         | C(0)   |
   --   | 0        1        E_2    B_3 |  |D(1)  |         | C(1)   |
   --
   --    ...
   --
   --   | 1        E_n-2    B_n-1    0 |  |D(n-3)|    =    | C(n-3) |
   --   | 0        1        E_n-1  B_n |  |D(n-2)|         | C(n-2) |
   --   | 0        0        1      E_n |  |D(n-1)|         | C(n-1) |
   --   | 0        0        0        1 |  |D(n)  |         | C(n)   |
   --
   --  where E_m = -Alpha_m, B_m = -Beta_m.
   --
   --  Can attemp iterative refinement of D with Newton's
   --  method:
   --             D_new = D_old + M_inverse * (C - M*D_old)
   --
   --  where D = M_inverse * C is the calculation of D given at the top.  if the
   --  said calculation of D is numerically imperfect, then the iteration above
   --  will produce improved values of D.  Of course, if the Coefficients of
   --  the polynomials C are numerically poor, then this effort may be wasted.
   --
   function  Sum
     (X                : in Real;
      C                : in Coefficients;
      Sum_Limit        : in Poly_ID_Type;
      P                : in Real          := 0.0;
      No_Of_Iterations : in Positive      := 1)
      return                Real
   is
      Product, Del : Coefficients; -- initialized by  M_inverse_times and M_times.
      D : Coefficients;            -- initialized by  M_inverse_times.
      Result : Real := Zero;
   begin
      --
      -- Step 1. Getting starting value of D (D_old) by solving M*D = C.
      --
      D := M_inverse_times (C, X, p, Sum_Limit);

      --
      -- Step 2. Improve D numerically through Newton iteration.
      --             D_new = D_old + M_inverse * (C - M*D_old)
      --
      Iterate: for k in 2..No_Of_Iterations loop

         --  Get Product = M*D_old:

         Product := M_times (D, X, p, Sum_Limit);

         --  Get Residual = C - M*D_old.  Reuse the array Product
         --  to hold the value of Residual:

         for m in Poly_ID_Type'First..Sum_Limit loop
            Product(m) := C(m) - Product(m);
         end loop;

         --  Get Del = M_inverse * (A - M*D_old) = M_inverse * Product:

         Del :=  M_inverse_times (Product, X, p, Sum_Limit);

         --  Get D_new = D_old + Del;

         for m in Poly_ID_Type'First..Sum_Limit loop
            D(m) := D(m) + Del(m);
         end loop;

      end loop Iterate;

      Result := D(0) * Q_0 (p, X);

      return Result;

   end Sum;

end Clenshaw;

