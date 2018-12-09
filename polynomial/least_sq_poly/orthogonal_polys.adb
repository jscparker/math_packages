
-------------------------------------------------------------------------------
-- package body Orthogonal_Polys, Gram-Schmidt orthogonal polynomials on discrete grid points.
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
-------------------------------------------------------------------------------

with Text_IO; use Text_IO;

package body Orthogonal_Polys is

  -- Global range for operations on data vectors:

  Data_First : Points_Index := Points_index'First;
  Data_Last  : Points_Index := Points_index'Last;

  -------------
  -- Make_Re --
  -------------

  --  Converts Coeff_Index to Real in a way that simplifies things
  --  when Real is a private extended precision floating point type.
  --  Its slow, but it doesn't slow down any critical inner loops.

  function Make_Re (N : Coeff_Index) return Real is
    Result : Real := Zero;
  begin
    for i in 1 .. N loop
      Result := Result + One;
    end loop;
    return Result;
  end Make_Re;

  ------------------------------
  -- Set_Limits_On_Vector_Ops --
  ------------------------------

  procedure Set_Limits_On_Vector_Ops (First, Last : Points_Index) is
  begin
     Data_First := First;
     Data_Last  := Last;
  end;

  --------------------------------
  --  Max_Permissable_Degree_Of --
  --------------------------------

  function Max_Permissable_Degree_Of (P : Polynomials) return Coeff_Index is
  begin
     return P.Max_Permissable_Degree_Of_Poly;
  end Max_Permissable_Degree_Of;

   -------------------
   -- Inner product --
   -------------------

   function Inner_Product
     (X, Y        : in Data;
      First, Last : in Points_Index;
      Weights     : in Data)
      return Real
   is
      Sum : Real := Zero;
   begin
      for i in First .. Last loop
         Sum := Sum + Weights(i) * Y(i) * X(i);
      end loop;
      return Sum;
   end Inner_Product;

   ---------
   -- "-" --
   ---------

   function "-"
     (Left  : in Data;
      Right : in Data)
      return Data
   is
      Result : Data;
   begin
      for i in Data_First .. Data_Last loop
         Result(i) := Left(i) - Right(i);
      end loop;
      return Result;
   end "-";

   ---------
   -- "+" --
   ---------

   function "+"
     (Left  : in Data;
      Right : in Data)
      return Data
   is
      Result : Data;
   begin
      for i in Data_First .. Data_Last loop
         Result(i) := Left(i) + Right(i);
      end loop;
      return Result;
   end "+";

   ---------
   -- "-" --
   ---------

   function "-"
     (Left  : in Data;
      Right : in Real)
      return Data
   is
      Result : Data;
   begin
      for i in Data_First .. Data_Last loop
         Result(i) := Left(i) - Right;
      end loop;
      return Result;
   end "-";

   ---------
   -- "*" --
   ---------

   function "*"
     (Left  : in Real;
      Right : in Data)
      return Data
   is
      Result : Data;
   begin
      for i in Data_First .. Data_Last loop
         Result(i) := Left * Right(i);
      end loop;
      return Result;
   end "*";

   ---------
   -- "*" --
   ---------

   -- fortran 90 Array * Array multiplication.

   function "*"
     (Left  : in Data;
      Right : in Data)
      return Data
   is
      Result : Data;
   begin
      for i in Data_First .. Data_Last loop
         Result(i) := Left(i) * Right(i);
      end loop;
      return Result;
   end "*";

   ----------------------------------------
   -- Start_Gram_Schmidt_Poly_Recursion  --
   ----------------------------------------

   --  Given a set of discrete points X
   --        X = (X1, X2, ...)
   --  and the weights
   --        W = (W1, W2, ...)
   --  the Gram-Schmidt recurrance method yield a unique set of
   --  orthogonal polynomials (functions defined at the grid points Xj).
   --  It is assumed that the the zeroth order poly is constant 1,
   --  and the 1st order poly is X.  (Both have normalization factors
   --  that are not applied to the polys here...instead the norm factor
   --  is calculated and returned as field in the poly data structure.)

   procedure Start_Gram_Schmidt_Poly_Recursion
     (X_axis, Weights : in     Data;
      First, Last     : in     Points_Index;
      Poly_0, Poly_1  : in out Poly_Data;
      Poly_Set        :    out Polynomials)
   is
      X_max, X_min, Max_Delta_X : Real;
      Slope, Const : Real;
      X_Poly_0 : Data;
      Alpha, Beta : Real;
      X_Scaled    : Data;

      type Local_Float is digits 9;
      Data_Length, No_Of_Weightless_Data_Points : Local_Float := 0.0;
   begin
      Set_Limits_On_Vector_Ops (First, Last);

      -- Step 0.  Make sure we have enough data points to calculate
      -- a polynomial of the desired degree. For example we don't want
      -- to try to fit a parabola to just two data points.  So we start by
      -- calculating the number of data points to be used. if Weight(i) is
      -- less than Smallest_Weight then don't count this as a data
      -- point. (However the data point will be used
      -- in the calculation below.)  We require a minimum of 2 data points.

      No_Of_Weightless_Data_Points := 0.0;
      for i in First .. Last loop
         if Weights(i) < Smallest_Weight then
            No_Of_Weightless_Data_Points := No_Of_Weightless_Data_Points - 1.0;
         end if;
      end loop;

      Data_Length := Local_Float (Last) - Local_Float (First) + 1.0;
      Data_length := Data_Length - No_Of_Weightless_Data_Points;
      if Data_Length < 2.0 then
         put_line ("Need at last 2 data points with positive weight.");
         raise Constraint_Error;
      end if;
      --  Because we make a first order poly in this procedure.

      if Data_Length - 1.0 <= Local_Float (Max_Order_Of_Poly) then
         Poly_Set.Max_Permissable_Degree_Of_Poly := Coeff_Index(Data_Length-1.0);
      else
         Poly_Set.Max_Permissable_Degree_Of_Poly := Max_Order_Of_Poly;
      end if;

      Poly_Set.Degree_Of_Poly := 1; -- Below we make 0 and 1.


      -- Step 1.  Make sure that the DeltaX is > Smallest_Delta_X for all
      -- all X.  So no two points can have the same X, and we can't
      -- have Data_X(i+1) < Data_X(i) for any i.

      for i in First+1 .. Last loop
         if X_axis (i) - X_axis (i-1) < Smallest_Delta_X then
            put_line("Data Points Must Be Ordered In X and Distinct.");
            raise Constraint_Error;
         end if;
      end loop;


      -- Step 2.  Accuracy can be much improved if X is in the interval
      -- [-2,2], so the X axis is scaled such that the X values lie in
      -- the interval [-2,2].  We
      -- scale X with the following equation: X_new = X * Slope + Const,
      -- where
      --            Slope =  4.0 / (X_max - X_min)
      -- and
      --            Const = -2.0 * (X_max + X_min) / (X_max - X_min).
      --
      -- The final poly for Y will be correct, but the coefficients
      -- of powers of X in the power form of Poly calculated below
      -- must be scaled.
      --
      -- The results will be stored in Poly_Set, for subsequent use in
      -- generating polynomials, and unscaling the calculations done on
      -- these polynomials.

      X_max := X_axis (First);
      X_min := X_axis (First);
      for i in First+1 .. Last loop
         if not (X_axis (i) < X_max) then
            X_max := X_axis (i);
         end if;
         if X_axis(i) < X_min then
            X_min := X_axis (i);
         end if;
      end loop;

      Max_Delta_X := (X_max - X_min);
      if Max_Delta_X < Smallest_Delta_X then
         put_line ("Data Points Too Close Together In X");
         raise Constraint_Error;
      end if;

      Slope :=  Four / Max_Delta_X;
      Const := -Two * ((X_max + X_min) / Max_Delta_X);

      X_scaled := Slope * X_axis - (-Const);  --  Vector operations.

      --  Store the results in Poly_Set:
      Poly_Set.X_scaled      := X_scaled;
      Poly_Set.Scale_Factors := X_Axis_Scale'(Slope, Const);


      -- Step 3. Get the Polynomials.  Vector op limits have been set above.

      -- The zeroth order poly (unnormalized) is just 1.0:

      Poly_Set.Alpha(0) := Zero;
      Poly_Set.Beta(0)  := Zero;

      Poly_0.Points  := (others => One);
      Poly_0.First   := First;
      Poly_0.Last    := Last;
      Poly_0.Squared :=
         Inner_Product (Poly_0.Points, Poly_0.Points, First, Last, Weights);
      Poly_0.Degree  := 0;

      -- Get the 1st order Polynomial.  Unnormalized, its just X - alpha;

      X_Poly_0 := X_Scaled * Poly_0.Points;
      Alpha := 
         Inner_Product(X_Poly_0, Poly_0.Points, First, Last, Weights) / Poly_0.Squared;
      Beta  := Zero;
      Poly_Set.Alpha(1) := Alpha;
      Poly_Set.Beta(1)  := Beta;

      Poly_1.Points  := X_scaled - Alpha;
      Poly_1.Squared := 
         Inner_Product (Poly_1.Points, Poly_1.Points, First, Last, Weights);
      Poly_1.First   := First;
      Poly_1.Last    := Last;
      Poly_1.Degree  := 1;

   end Start_Gram_Schmidt_Poly_Recursion;

   -------------------
   -- Get_Next_Poly --
   -------------------

   -- We want Q_m, Alpha_m, and Beta_m, given the previous values.
   --
   --                Q_0  =  1
   --                Q_1  = (X - Alpha_1)
   --                Q_m  = (X - Alpha_m) * Q_m-1  -  Beta_m*Q_m-2
   -- where
   --                Alpha_m  =  (X*Q_m-1, Q_m-1) / (Q_m-1, Q_m-1)
   --                Beta_m   =  (X*Q_m-1, Q_m-2) / (Q_m-2, Q_m-2)
   --
   -- Can be shown:  Beta_m   =    (Q_m-1, Q_m-1) / (Q_m-2, Q_m-2) which is
   -- the form used below.

   procedure Get_Next_Poly 
     (Poly_0, Poly_1 : in     Poly_Data;
      Weights        : in     Data;
      Poly_2         : in out Poly_Data;
      Poly_Set       : in out Polynomials) 
   is
      X_Poly_1    : Data;
      Alpha, Beta : Real;
      X_scaled    : Data renames Poly_Set.X_scaled;
      Degree_2    : Coeff_Index;
      First : Points_Index renames Poly_1.First;
      Last  : Points_Index renames Poly_1.Last;
   begin

      Set_Limits_On_Vector_Ops (First, Last);

      --  Have to tell the vector ops that the data goes from First..Last.
      --  Not really necessary, because Start_Gram_.. has already done this.
      --  But we do it anyway.  Next some checks:
   
      if Poly_0.First /= Poly_1.First or Poly_0.Last /= Poly_1.Last then
         put_line ("Some error in input polys for Get_Next_Poly.");
         raise Constraint_Error;
      end if;
   
      -- Must have Poly_Set.Degree_Of_Poly = Poly_1.Degree = Poly_0.Degree+1
   
      if Poly_Set.Degree_Of_Poly /= Poly_0.Degree + 1 or
         Poly_Set.Degree_Of_Poly /= Poly_1.Degree then
         put_line ("Some error in input polys for Get_Next_Poly.");
         raise Constraint_Error;
      end if;
   
      --  The purpose of this is to raise the degree of poly_set by one.
      --  Can we do that?
   
      if Poly_Set.Degree_Of_Poly >= Poly_Set.Max_Permissable_Degree_Of_Poly then
         put_line ("Cannot make a poly of that order with so few points.");
         raise Constraint_Error;
      end if;
   
      -- Now we can construct the next higher order polynomial: Poly_2
   
      X_Poly_1 := X_Scaled * Poly_1.Points;
      Alpha := Inner_Product(X_Poly_1, Poly_1.Points,
             First, Last, Weights) / Poly_1.Squared;
      Beta  := Poly_1.Squared / Poly_0.Squared;
   
      Degree_2                  := Poly_Set.Degree_Of_Poly + 1;
      Poly_Set.Degree_Of_Poly   := Degree_2;
      Poly_Set.Beta (Degree_2)  := Beta;
      Poly_Set.Alpha (Degree_2) := Alpha;
   
      Poly_2.Points  := (X_scaled - Alpha) * Poly_1.Points
          - Beta * Poly_0.Points;
      Poly_2.Squared := Inner_Product (Poly_2.Points, Poly_2.Points,
            First, Last, Weights);
      Poly_2.First   := Poly_1.First;
      Poly_2.Last    := Poly_1.Last;
      Poly_2.Degree  := Poly_1.Degree + 1;
   
   end Get_Next_Poly;

   -------------------------------
   -- Get_Coeffs_Of_Powers_Of_X --
   -------------------------------

   -- Calculate the Coefficients of powers of X in the best fit poly, using
   -- Alpha, Beta, C as calculated above and the following formula for
   -- the orthogonal polynomials:
   --             Q_0  =  1
   --             Q_1  = (X - Alpha_1)
   --             Q_k  = (X - Alpha_k) * Q_k-1  -  Beta_k*Q_k-2
   -- and the best fit poly is SUM {C_k * Q_k}.  The Coefficients of X**k
   -- are put in array Poly_Coefficients(k), which is set to 0.0 initially,
   -- since the formula may assign values to only a small subset of it.
   -- the E_k's in SUM {E_k * X**k} go into Poly_Coefficients(k).
   -- Clenshaw's formula is used to get Poly_Coefficients.  The
   -- coefficients of the following D polynomials are put into arrays
   -- D_0, D_1, and D_2, and advanced recursively until the final
   -- D_0 = SUM {C_k * Q_k} is found. This will be named Poly_Coefficients.
   -- The recursion formula for the Coefficients follows from the formula for D(X):
   --
   --     D_n+2(X) = 0
   --     D_n+1(X) = 0
   --     D_m(X)   = C_m + (X - Alpha_m+1)*D_m+1(X) - Beta_m+2*D_m+2(X)
   --
   -- where n = Desired_Poly_Degree and m is in 0..n.
   --
   -- Now suppose we want the coefficients of powers of X for the above D polys.
   -- (In the end that will give us the coeffs of the actual poly, D_0.)
   -- Well, the first poly D_n is 0-th order and equals C(n).  The second poly, D_n-1,
   -- gets a contribution to its coefficient of X**1 from the X*D_n term.  That
   -- contribution is the coefficient of X**0 in D_n.  Its X**0 coeff. gets a
   -- contribution from C(n-1) and one from -Alpha(n)*D_n at X**0, or -Alpha(n)*D_n(0).
   -- Now we re-use the D
   -- polynomial arrays to store these coefficients in the obvious place.
   -- The arrays D_0, D_1, and D_2 are initialized to 0.0.
   -- D_0, D_1, and D_2 contain be the coeff's of powers of X in the orthogonal
   -- polynomials D_m, D_m+1, D_m+2, respectively.  The formulas above
   -- imply:  for m in Desired_Poly_Degree .. 0:
   --
   --   D_0(0) := C(m);
   --   for k in 1 .. Desired_Poly_Degree loop
   --     D_0(k) := D_1(k-1);
   --   end loop;
   --   --  The above initalizes D_0.
   --
   --   for k in 0 .. Desired_Poly_Degree loop
   --     D_0(k) := D_0(k) - Alpha(m+1)*D_1(k) - Beta(m+2)*D_2(k);
   --   end loop;
   --
   -- So if we define shl_1 = Shift_Array_Left_by_1_and_Put_0_at_Index_0, the
   -- formula in vector notation is:
   --
   -- D_0 = Y(m) + shl_1 (D_1) - Alpha(m+1)*D_1 - Beta(m+2)*D_2
   --
   -- where Y(m) = (C(m),0,0,....).
   -- the above step is repeated recursivly using D_2 = D_1, and D_1 = D_0.
   -- Notice that the above formula must be modified at m = n and m = n-1.
   -- In matrix notation, using vector D, and vector Y this becomes:
   --
   --
   --   | 1        0        0        0 |  |D(n)  |         | Y(n)   |
   --   | A_n      1        0        0 |  |D(n-1)|    =    | Y(n-1) |
   --   | B_n      A_n-1    1        0 |  |D(n-2)|         | Y(n-2) |
   --   | 0        B_n-1    A_n-2    1 |  |D(n-3)|         | Y(n-3) |
   --
   -- where A_m = -(shl_1 - Alpha_m), B_m = Beta_m, and Y(m) = (C(m),0,0,....).
   -- In the end, D(0) should be an array that contains the Coefficients of Powers
   -- of X.  The operator is not a standard matrix operator, but its linear and we
   -- know its inverse and forward operation, so Newton's method gives the
   -- iterative refinement.  The array D(k)(m) is possibly very large!

   procedure Get_Coeffs_Of_Powers_Of_X
     (Coeffs : out Powers_Of_X_Coeffs;
      C : in  Poly_Sum_Coeffs;
      Poly_Set : in  Polynomials) 
   is
      D_1, D_2 : Recursion_Coeffs := (others => Zero);
      D_0      : Recursion_Coeffs := (others => Zero);

      m : Coeff_Index;
      Const2   : Real := Zero;

      A : Recursion_Coeffs renames Poly_Set.Alpha;
      B : Recursion_Coeffs renames Poly_Set.Beta;
      Scale_Factors : X_Axis_Scale renames Poly_Set.Scale_Factors;
      Poly_Degree   : Coeff_Index renames Poly_Set.Degree_Of_Poly;
   begin
      Coeffs := (others => Zero);

       -- Special Case.  Calculate D_2 (i.e. the D_m+2 poly):

      m := Poly_Degree;
      D_2(0) := C(m);

      if Poly_Degree = 0 then
         D_0 := D_2;
      end if;

       -- Special Case.  Calculate D_1 (i.e. the D_m+1 poly):

      if Poly_Degree > 0 then

         m := Poly_Degree - 1;

         D_1(0) := C(m);
    
         for k in Coeff_Index range 1 .. Poly_Degree-m loop
           D_1(k) := D_2(k-1);
         end loop;
    
         -- The previous 2 assigments have initialized D_1. Now:
    
         for k in Coeff_Index'First .. Poly_Degree-m-1 loop
           D_1(k) := D_1(k) - A(m+1) * D_2(k);
         end loop;
    
      end if;
    
      if Poly_Degree = 1 then
         D_0 := D_1;
      end if;
    
      -- Calculate D's for D_n-2 and lower:
    
      if Poly_Degree > 1 then
    
         for m in reverse Coeff_Index'First .. Poly_Degree-2 loop
       
            D_0(0) := C(m);
       
            for k in Coeff_Index range 1 .. Poly_Degree-m loop
              D_0(k) := D_1(k-1);
            end loop;
       
            -- The previous 2 assigments have initialized D_0. Now:
       
            for k in Coeff_Index'First .. Poly_Degree-m-1 loop
              D_0(k) := D_0(k) - A(m+1) * D_1(k);
            end loop;
       
            for k in Coeff_Index'First .. Poly_Degree-m-2 loop
              D_0(k) := D_0(k) - B(m+2) * D_2(k);
            end loop;
       
            D_2 := D_1;
            D_1 := D_0;
        
         end loop;
        
      end if;

      -- Now we have the coefficients of powers of X for the poly P1 (Z(X))
      -- whose Z is in the range [-2,2].  How do we get the coeffs of
      -- poly P2 (X) = P1 (Z(X)) whose X is in the range [a, b].  The
      -- relation between Z and X is Z = 2*(2*X - (a + b)) / (a - b).
      -- or Z = Slope * (X - Const2) where Slope = 4 / (a - b) and
      -- Const2 = (a + b) / 2.  We have P1 (Z).  The first step in getting
      -- P2 (X) is to get P1 (X - Const2) by multiplying the Coeffs of
      -- of P1 by powers of Slope.
      -- This is a common source of overflow.
      -- The following method is slower but slightly more overflow resistant
      -- than the more obvious method.
     
      for j in 1 .. Poly_Degree loop
         for k in Coeff_Index'First+j .. Poly_Degree loop
            D_0(k) := D_0(k) * Scale_Factors.Slope;
         end loop;
      end loop;
     
      -- Next we want coefficients of powers of X in P2 where P2 (X) is
      -- defined P2 (X) = P1 (X - Const2).  In other words we want the
      -- coefficients E_n in
      --
      -- P2 (X) = E_n*X**n + E_n-1*X**n-1 .. + E_1*X + E_0.
      --
      -- We know that
      --
      -- P1 (X) = F_n*X**n + F_n-1*X**n-1 .. + F_1*X + F_0,
      --
      -- where the F's are given by Coeff_0 above,
      -- and P2 (X + Const2) = P1 (X).  Use synthetic division to
      -- shift P1 in X as follows.  (See Mathew and Walker).
      -- P2 (X + Const2) = E_n*(X + Const2)**n + .. + E_0.
      -- So if we divide P2 (X + Const2) by (X + Const2) the remainder
      -- is E_0.  if we repeat the division, the remainder is E_1.
      -- So we use synthetic division to divide P1 (X) by (X + Const2).
     
      -- Synthetic division: multiply (X + Const2) by
      -- F_n*X**(n-1) = D_0(n)*X**(n-1) and subtract from P1.
      -- Repeat as required.
     
      -- What is Const2?
      -- Slope :=  4.0 / Max_Delta_X;
      -- Const := -2.0 * (X_max + X_min) / Max_Delta_X; -- 2 (a + b)/(b - a)
      -- X_scaled = Z = X * Slope + Const.
      -- Want Z = Slope * (X - Const2).  Therefore Const2 = - Const/Slope
     
      Const2 := -Scale_Factors.Const / Scale_Factors.Slope;
     
      for m in Coeff_Index range 1 .. Poly_Degree loop
         for k in reverse Coeff_Index range m .. Poly_Degree loop
            D_0 (k-1) := D_0 (k-1) - D_0 (k) * Const2;
         end loop;
         Coeffs (m-1) := D_0 (m-1);
      end loop;
      Coeffs (Poly_Degree) := D_0 (Poly_Degree);
  
   end Get_Coeffs_Of_Powers_Of_X;

   ----------------
   -- Horner_Sum --
   ----------------

   --  Want   Sum  =  a_0 + a_1*X + a_2*X**2 + ... + a_n*X**n.
   --  or in Horner's form: Sum  =  a_0 + X*(a_1 + ... + X*(a_n-1 + X*a_n)))))).
   --  This is easily written as matrix equation, with Sum = S_0:
   --
   --  S_n = a_n;  S_n-1 = a_n-1 + X*S_n;  S_1 = a_1 + X*S_2; S_0 = a_0 + X*S_1;
   --
   --  In matrix form, vector S is the solution to matrix equation M*S = A,
   --  where A = (a_0,...,a_n), S = (S_0,...,S_n) and matrix M is equal to
   --  the unit matrix i minus X*O1, where O1 is all 1's on the 1st lower off-
   --  diagonal.  The reason this form is chosen is that the solution vector
   --  S can be improved numerically by iterative refinement with Newton's
   --  method:
   --             S(k+1) = S(k) + M_inverse * (A - M*S(k))
   --
   --  where S = M_inverse * A is the calculation of S given above.  if the
   --  said calculation of S is numerically imperfect, then the iteration above
   --  will produce improved values of S.  Of course, if the Coefficients of
   --  the polynomial A are numerically poor, then this effort may be wasted.
   --
   function Horner_Sum
     (A                : in Recursion_Coeffs;
      Coeff_Last       : in Coeff_Index;
      X                : in Real;
      No_Of_Iterations : in Natural)
      return Real
   is
      S : Recursion_Coeffs := (others => Zero);
      Del, Product : Recursion_Coeffs;
   begin

      if Coeff_Last = Coeff_Index'First then
         return A(Coeff_Index'First);
      end if;
      --  Poly is zeroth order = A(Index'First).  No work to do.  Go home.

      --  Now solve for S in the matrix equation M*S = A. first iteration:

      S(Coeff_Last) := A(Coeff_Last);
      for n in reverse Coeff_Index'First .. Coeff_Last-1 loop
         S(n) := A(n) + X * S(n+1);
      end loop;

      --  Now iterate as required.  We have the first S, S(1), now get S(2) from
      --     S(k+1) = S(k) + M_inverse * (A - M*S(k))

      Iterate: for k in 1..No_Of_Iterations loop

         --  Get Product = M*S(k):

         Product(Coeff_Last) := S(Coeff_Last);
         for n in reverse Coeff_Index'First..Coeff_Last-1 loop
            Product(n) := S(n) - X * S(n+1);
         end loop;

         --  Get  Product = Residual = A - M*S(k):

         for n in Coeff_Index'First .. Coeff_Last loop
            Product(n) := A(n) - Product(n);
         end loop;

         --  Get Del = M_inverse * (A - M*S(k)) = M_inverse * Product:

         Del(Coeff_Last) := Product(Coeff_Last);
         for n in reverse Coeff_Index'First .. Coeff_Last-1 loop
            Del(n) := Product(n) + X * Del(n+1);
         end loop;

         --  Get S(k+1) = S(k) + Del;

         for n in Coeff_Index'First .. Coeff_Last loop
            S(n) := S(n) + Del(n);
         end loop;

       end loop Iterate;

       return S(Coeff_Index'First);

   end Horner_Sum;

   -------------------
   -- Poly_Integral --
   -------------------

   --  Use Clenshaw summation to get coefficients of powers of X,
   --  then sum analytically integrated polynomial using Horner's rule.
   --  Poly_Integral returns the indefinite integral.  Integral on an 
   --  interval [A, B] is Poly_Integral(B) - Poly_Integral(A).
   --
   function Poly_Integral
     (X                    : in Real;
      C                    : in Poly_Sum_Coeffs;
      Poly_Set             : in Polynomials;
      Order_Of_Integration : in Coeff_Index   := 1)
      return Real
   is
      D_1, D_2 : Recursion_Coeffs := (others => Zero);
      D_0      : Recursion_Coeffs := (others => Zero);

      m : Coeff_Index;
      Result   : Real := Zero;
      Denom, X_scaled : Real := Zero;

      A : Recursion_Coeffs renames Poly_Set.Alpha;
      B : Recursion_Coeffs renames Poly_Set.Beta;
      Poly_Degree   : Coeff_Index renames Poly_Set.Degree_Of_Poly;
   begin

      -- Special Case.  Calculate D_2 (i.e. the D_m+2 poly):

      m := Poly_Degree;
      D_2(0) := C(m);

      if Poly_Degree = 0 then
         D_0 := D_2;
      end if;

       -- Special Case.  Calculate D_1 (i.e. the D_m+1 poly):

      if Poly_Degree > 0 then

         m := Poly_Degree - 1;

         D_1(0) := C(m);

         for k in Coeff_Index range 1 .. Poly_Degree-m loop
            D_1(k) := D_2(k-1);
         end loop;

         -- The previous 2 assigments have initialized D_1. Now:

         for k in Coeff_Index'First .. Poly_Degree-m-1 loop
            D_1(k) := D_1(k) - A(m+1) * D_2(k);
         end loop;

      end if;

      if Poly_Degree = 1 then
         D_0 := D_1;
      end if;

      -- Calculate D's for D_n-2 and lower:

      if Poly_Degree > 1 then

         for m in reverse Coeff_Index'First .. Poly_Degree-2 loop

            D_0(0) := C(m);

            for k in Coeff_Index range 1 .. Poly_Degree-m loop
               D_0(k) := D_1(k-1);
            end loop;

            -- The previous 2 assigments have initialized D_0. Now:

            for k in Coeff_Index'First .. Poly_Degree-m-1 loop
               D_0(k) := D_0(k) - A(m+1) * D_1(k);
            end loop;

            for k in Coeff_Index'First .. Poly_Degree-m-2 loop
               D_0(k) := D_0(k) - B(m+2) * D_2(k);
            end loop;

            D_2 := D_1;
            D_1 := D_0;

         end loop;

      end if;

      --  The unscaled Coeffs of X**m are D_0(m).  Integrate once,
      --  (brings a 1/(m+1) down) then use Horner's rule for poly sum:

      --  First scale X from [a, b] to [-2, 2]:
      X_Scaled := X * Poly_Set.Scale_Factors.Slope + Poly_Set.Scale_Factors.Const;

      for m in reverse Coeff_Index'First .. Poly_Degree loop

         Denom := One;
         for i in 1 .. Order_Of_Integration loop
            Denom := Denom * (Make_Re (m + i));
         end loop;

         D_0(m) := D_0(m) / Denom;
      end loop;

      Result :=
         Horner_Sum
           (A                => D_0,
            Coeff_Last       => Poly_Degree,
            X                => X_scaled,
            No_Of_Iterations => 1);

      Result := Result * X_scaled ** Integer(Order_Of_Integration);
      --  This X was neglected above in Horner_Sum.

      --  The integral was on a scaled interval [-2, X]. Unscale the result:

      Result :=
         Result / Poly_Set.Scale_Factors.Slope ** Integer(Order_Of_Integration);

      return Result;

   end Poly_Integral;

   ----------------------
   -- Poly_Derivatives --
   ----------------------

   -- How do we get the derivatives of the best-fit polynomial?  Just
   -- take the derivative of the Clenshaw recurrence formula given above.
   -- In the special case of orthogonal polynomials it is particularly
   -- easy.  By differentiating the formula given above p times its easy
   -- to see that the p-th derivative of the D_m functions of X satisfy:
   --
   -- p = order of derivative = 0:
   --
   --   D_n+2(0,X) = 0
   --   D_n+1(0,X) = 0
   --   D_m(0,X)  = C_m + (X - Alpha(m+1)) * D_m+1(0,X) + Beta(m+2) * D_m+2(0,X)
   --
   -- p = order of derivative > 0:
   --
   --   D_n+2(p,X) = 0
   --   D_n+1(p,X) = 0
   --   D_m(p,X)
   --    = p*D_m+1(p-1,X) + (X - Alpha(m+1))*D_m+1(p,X) - Beta(m+2)*D_m+2(p,X)
   --
   -- for m in 0..n,
   -- where D(p,X) is the pth derivative of D(X).  It follows that the
   -- p-th derivative of the sum over m of C_m*Q_m(X) equals D_0(p,X).
   --
   -- We still aren't finished.  What we really want is the derivative
   -- respect the UNSCALED variable, Y.  Here X = X_Scaled is in the range
   -- [-2,2] and X_scaled = Slope * Y + Constant.  So d/dY = Slope * d/dX.
   -- Usually Y is in (say) 1..100, and X is in -2..2, so Slope is << 1.
   -- It follows that the recurrence relation for the p-th derivative
   -- of the D polynomials respect Y is
   --
   --   D_n+2(p,X) = 0
   --   D_n+1(p,X) = 0
   --   D_m(p,X) = p * Slope * D_m+1(p-1,X)
   --                      + (X - Alpha(m+1))*D_m+1(p,X) - Beta(m+2)*D_m+2(p,X)
   --
   -- for m in 0..n,
   -- where D(p,X) is the p-th derivative of D(X).  It follows that the
   -- p-th derivative the sum over m of C_m*Q_m(X) equals D_0(p,X).
   --
   -- To perform the calculation, the 0th derivative (p=0) D is calculated
   -- first, then used as a constant in the recursion relation to get the
   -- p=1 D.  These steps are repeated recursively.
   --
   -- In the code that follows D is an array only in "m".  X is a constant,
   -- input by the user of the procedure, and p is reduced to 2 values,
   -- "Hi" and "Low".  So we calculate D_low(m) where derivative order
   -- p = Low, which starts at 0, and use D_low(m) to get D_hi(m), where
   -- Hi = Low + 1.
   --
   -- p = order of derivative == Low = 0:
   --
   --   D_low(n+2) = 0
   --   D_low(n+1) = 0
   --   D_low(m)  = C_m + (X - Alpha(m+1)) * D_low(m+1) + Beta(m+2) * D_low(m+2)
   --
   -- p = order of derivative == hi > 0
   --
   --   D_hi(n+2) = 0
   --   D_hi(n+1) = 0
   --   D_hi(m) = p * Slope * D_low(m+1)
   --                      + (X - Alpha(m+1))*D_hi(m+1) - Beta(m+2)*D_hi(m+2)
   --
   -- Next iterative refinement is optionally performed.  For each value of
   -- of p, the following matrix equation represents the recursive equations
   -- above.  Remember, in the following, D_low is a previously calculated
   -- constant:
   --
   --   | 1        0        0        0 |  |D_hi(n)  |         | Y(n)   |
   --   | A_n      1        0        0 |  |D_hi(n-1)|    =    | Y(n-1) |
   --   | B_n      A_n-1    1        0 |  |D_hi(n-2)|         | Y(n-2) |
   --   | 0        B_n-1    A_n-2    1 |  |D_hi(n-3)|         | Y(n-3) |
   --
   -- where A_m = -(X - Alpha_m), B_m = Beta_m, and Y(m) = C(m) if p=0, and
   -- Y(m) = p * Slope * D_low(m+1) if p > 0.  (Remember, D_any(n+1) = 0.0).
   -- So the refinement is in the m iteration not the p iteration.
   -- The iteration can actually be done in both p and m, but in that case D
   -- must be stored as a 2-d array.  In that case the matrix equation
   -- is a block Lower triangular matrix.  We do it the less sophisticated way
   -- here.

   procedure Poly_Derivatives
     (Derivatives    : in out Derivative_List;
      X              : in     Real;
      Order_Of_Deriv : in     Derivatives_Index;
      C              : in     Poly_Sum_Coeffs;
      Poly_Set       : in     Polynomials)
   is
      Order             : Real;
      Order_times_Slope : Real;
      X_Scaled          : Real;

      D_Hi, D_Low : Recursion_Coeffs; 
      -- D_Hi is the higher deriv. in the recurrence relation.

      Local_Order_Of_Deriv : Derivatives_Index;

      A : Recursion_Coeffs renames Poly_Set.Alpha;
      B : Recursion_Coeffs renames Poly_Set.Beta;
      Scale_Factors  : X_Axis_Scale renames Poly_Set.Scale_Factors;
      Poly_Degree    : Coeff_Index renames Poly_Set.Degree_Of_Poly;

   begin
      --  The derivatives of a polynomial are zero if their order is
      --  greater than the degree of the polynomial, so in that case
      --  don't bother to get them:
      Derivatives := (others => Zero);

      if Order_Of_Deriv > Poly_Degree then
         Local_Order_Of_Deriv := Poly_Degree;
      else
         Local_Order_Of_Deriv := Order_Of_Deriv;
      end if;

      -- Scale X to the interval [-2,2].
      X_Scaled := X * Scale_Factors.Slope + Scale_Factors.Const;

      -- Step 1.  We need a 0th order poly to start off the recurrence
      -- relation.  Start by getting undifferentiated D's (p=0).
      -- Store them in array D_Low.  The Low is for "Lower Order".
      -- Use recurrence relation to get Poly at X.

      -- Start with special formulas for the 1st 2 D-Polys:

      D_Low(Poly_Degree)   := C(Poly_Degree);

      if Poly_Degree > Coeff_Index'First then
         D_low(Poly_Degree-1) := C(Poly_Degree-1) +
           (X_Scaled - A(Poly_Degree)) * D_low(Poly_Degree);
      end if;

      for m in reverse Coeff_Index'First+2 .. Poly_Degree loop
         D_Low(m-2) := C(m-2) +
           (X_Scaled - A(m-1))*D_Low(m-1) - B(m)*D_low(m);
      end loop;

      Derivatives (Derivatives_Index'First) := D_Low(Coeff_Index'First);

      -- Step 2.  Use the recurrence relation to get next higher
      -- higher derivative.  Store it in array D_Hi.

      for p in Derivatives_Index'First+1 .. Local_Order_Of_Deriv loop

         Order             := Make_Re (p);
         D_Hi(Poly_Degree) := Zero;
         Order_times_Slope := Order * Scale_Factors.Slope;

         if Poly_Degree > Coeff_Index'First then
            D_Hi(Poly_Degree-1) := Order_times_Slope * D_Low(Poly_Degree) +
               (X_Scaled - A(Poly_Degree)) * D_Hi(Poly_Degree);
         end if;

         for m in reverse Coeff_Index'First+2 .. Poly_Degree loop
            D_Hi(m-2) := Order_times_Slope * D_low(m-1) +
               (X_Scaled - A(m-1))*D_Hi(m-1) - B(m)*D_Hi(m);
         end loop;

         Derivatives (p) := D_Hi(Coeff_Index'First);
         D_Low := D_Hi;

      end loop;

   end Poly_Derivatives;

   ----------------
   -- Poly_Value --
   ----------------

   --  This is easily written as matrix equation, with Sum = S_0:
   --
   --  D_n   = C_n;
   --  D_n-1 = C_n-1 + (X - A_n)*D_n;
   --  D_n-2 = C_n-2 + (X - A_n-1)*D_n-1 - B_n-2*D_n-2;
   --  ...
   --  D_0   = C_0 +   (X - A_1)*D_1     - B_2*D_2
   --
   -- In matrix form, M*D = C, this becomes:
   --
   --   | 1        0        0        0 |  |D(n)  |         | C(n)   |
   --   | E_n      1        0        0 |  |D(n-1)|    =    | C(n-1) |
   --   | B_n      E_n-1    1        0 |  |D(n-2)|         | C(n-2) |
   --   | 0        B_n-1    E_n-2    1 |  |D(n-3)|         | C(n-3) |
   --
   --  where E_m = (A_m - X), B_m = B_m.
   --
   --  D can be improved numerically by iterative refinement with Newton's
   --  method:
   --             D_new = D_old + M_inverse * (C - M*D_old)
   --
   --  where D = M_inverse * C is the calculation of D given at the top.  if the
   --  said calculation of D is numerically imperfect, then the iteration above
   --  will produce improved values of D.  Of course, if the Coefficients of
   --  the polynomials C are numerically poor, then this effort may be wasted.

   function Poly_Value
     (X        : in Real;
      C        : in Poly_Sum_Coeffs;
      Poly_Set : in Polynomials)
      return Real
   is
      D, Product, Del : Recursion_Coeffs := (others => Zero);
      X_Scaled      : Real;
      m             : Coeff_Index;

      A : Recursion_Coeffs renames Poly_Set.Alpha;
      B : Recursion_Coeffs renames Poly_Set.Beta;
      Scale_Factors : X_Axis_Scale renames Poly_Set.Scale_Factors;
      Poly_Degree   : Coeff_Index  renames Poly_Set.Degree_Of_Poly;

      No_Of_Iterations : constant Natural := 0;
   begin
      -- Scale X to the interval [-2,2].
      X_Scaled := X * Scale_Factors.Slope + Scale_Factors.Const;

      -- Step 0. Poly is zeroth order = C(Index'First).  No work to do.

      if Poly_Degree = Coeff_Index'First then
         m      := Poly_Degree;
         D(m)   := C(m);
         return D(Coeff_Index'First);
      end if;

      -- Step 0b. Poly is 1st order.  Almost no work to do.
      -- Don't do any iteration.

      if Poly_Degree = Coeff_Index'First + 1 then
         m      := Poly_Degree;
         D(m)   := C(m);
         m      := Poly_Degree - 1;
         D(m)   := C(m) - (A(m+1) - X_Scaled)*D(m+1);
         return D(Coeff_Index'First);
      end if;


      -- Step 1.  We now know henceforth that Poly_Degree > 1.
      -- Start by getting starting value of D by solving M*D = C.
      -- Use recurrence relation to get Poly at X.
      -- Start with special formulas for the 1st two Polys:

      m      := Poly_Degree;
      D(m)   := C(m);

      m      := Poly_Degree - 1;
      D(m)   := C(m) - (A(m+1) - X_Scaled)*D(m+1);

      for m in reverse Coeff_Index'First .. Poly_Degree-2 loop
         D(m) := C(m) - (A(m+1) - X_Scaled)*D(m+1) - B(m+2)*D(m+2);
      end loop;

      -- Step 2. Improve D numerically through Newton iteration.
      --             D_new = D_old + M_inverse * (C - M*D_old)

      Iterate: for k in 1 .. No_Of_Iterations loop

         --  Get Product = M*D(k):

         m          := Poly_Degree;
         Product(m) := D(m);

         m          := Poly_Degree - 1;
         Product(m) := D(m) + (A(m+1) - X_Scaled)*D(m+1);

         for m in reverse Coeff_Index'First .. Poly_Degree-2 loop
            Product(m) := D(m) + (A(m+1) - X_Scaled)*D(m+1) + B(m+2)*D(m+2);
         end loop;

         --  Get Residual = C - M*D(k) and set it equal to Product:

         for m in Coeff_Index'First .. Poly_Degree loop
            Product(m) := C(m) - Product(m);
         end loop;

         --  Get Del = M_inverse * (A - M*S(k)) = M_inverse * Product:

         m      := Poly_Degree;
         Del(m) := Product(m);

         m      := Poly_Degree - 1;
         Del(m) := Product(m) - (A(m+1) - X_Scaled)*Del(m+1);

         for m in reverse Coeff_Index'First .. Poly_Degree-2 loop
            Del(m) := Product(m) - (A(m+1) - X_Scaled)*Del(m+1) - B(m+2)*Del(m+2);
         end loop;

         --  Get D(k+1) = D(k) + Del;

         for m in Coeff_Index'First .. Poly_Degree loop
            D(m) := D(m) + Del(m);
         end loop;

      end loop Iterate;

      return D(Coeff_Index'First);

   end Poly_Value;

   --------------
   -- Poly_Fit --
   --------------

   --  Generate orthogonal polys and project them onto the data
   --  with the Inner_Product function in order to calculate
   --  C_k, the Best_Fit_Coefficients.

   procedure Poly_Fit
     (Data_To_Fit         : in     Data;
      X_axis, Weights     : in     Data;
      First, Last         : in     Points_Index;
      Desired_Poly_Degree : in     Coeff_Index;
      Best_Fit_Poly       : in out Poly_Data;
      Best_Fit_Coeffs     : in out Poly_Sum_Coeffs;
      Poly_Set            : in out Polynomials;
      Mean_Square_Error   :    out Real)
   is
      Poly_0, Poly_1, Poly_2 : Poly_Data;
      Local_Poly_Degree : Coeff_Index := Desired_Poly_Degree;
      Data_Length : Real;

      C : Poly_Sum_Coeffs renames Best_Fit_Coeffs;
      Dat      : Data renames Data_To_Fit;
      Best_Fit : Data renames Best_Fit_Poly.Points;
   begin
      Start_Gram_Schmidt_Poly_Recursion
        (X_axis, Weights, First, Last, Poly_0, Poly_1, Poly_Set);

       -- Get Degree of poly: may be less than the desired degree
       -- if too few points exist in the subset defined by X.

       if Local_Poly_Degree > Poly_Set.Max_Permissable_Degree_Of_Poly then
          Local_Poly_Degree := Poly_Set.Max_Permissable_Degree_Of_Poly;
       end if;

       Data_Length := Make_Re (Local_Poly_Degree) + (One);
       --  By definition of Local_Poly_Degree.

       Set_Limits_On_Vector_Ops (First, Last);

       --  Get C(0) = Coefficient of 0th orthog. poly.

       C(0) :=
          Inner_Product (Poly_0.Points, Dat, First, Last, Weights) / Poly_0.Squared;

       Best_Fit := C(0) * Poly_0.Points;

       --  Get C(1) = Coefficient of 1st orthog. poly.

       if Local_Poly_Degree > 0 then

          C(1) := Inner_Product (Poly_1.Points, Dat - Best_Fit,
          First, Last, Weights) / Poly_1.Squared;

          Best_Fit := Best_Fit + C(1) * Poly_1.Points;

       end if;

       --  Get C(2) = Coefficient of 2nd orthog. poly.

       if Local_Poly_Degree > 1 then

          Get_Next_Poly (Poly_0, Poly_1, Weights, Poly_2, Poly_Set);

          C(2) :=
             Inner_Product
               (Poly_2.Points, Dat - Best_Fit, First, Last, Weights) / Poly_2.Squared;

          Best_Fit := Best_Fit + C(2) * Poly_2.Points;

       end if;

       -- Higher order Polynomials: get C(i) which is the coefficient of
       -- the Ith orthogonal polynomial.  Also get the Best_Fit_Polynomial.
       -- Notice that the formula used to get C(i) is a little more complicated
       -- than the the one written in the prologue above.  The formula used
       -- below gives substantially better numerical results, and is
       -- mathematically identical to formula given in the prologue.

       for N in Coeff_Index range 3 .. Local_Poly_Degree loop

          Poly_0   := Poly_1;
          Poly_1   := Poly_2;
          Get_Next_Poly (Poly_0, Poly_1, Weights, Poly_2, Poly_Set);

          C(N) := 
             Inner_Product 
               (Poly_2.Points, Dat - Best_Fit, First, Last, Weights) / Poly_2.Squared;

          Best_Fit := Best_Fit + C(N) * Poly_2.Points;

       end loop;

       -- Reuse Poly_0 to calculate Error**2 per Point = Mean_Square_Error:

       Poly_0.Points     := Dat - Best_Fit;
       Mean_Square_Error := Inner_Product (Poly_0.Points, Poly_0.Points,
              First, Last, Weights) / Data_Length;

       --  Finish filling Best_Fit_Poly and Poly_Set:

       Poly_Set.Degree_Of_Poly := Local_Poly_Degree;

     --Best_Fit_Poly.Points  := Best_Fit;  -- through renaming.
       Best_Fit_Poly.First   := First;
       Best_Fit_Poly.Last    := Last;
       Best_Fit_Poly.Degree  := Local_Poly_Degree;
       Best_Fit_Poly.Squared :=
          Inner_Product (Best_Fit, Best_Fit, First, Last, Weights);

   end Poly_Fit;

begin

   if  Max_Order_Of_Poly > Max_No_Of_Data_Points - 1  then
      put_line("Max poly order must be less than max number of data points.");
      raise Constraint_Error;
   end if;

end Orthogonal_Polys;

