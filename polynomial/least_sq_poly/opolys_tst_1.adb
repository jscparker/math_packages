
with Orthogonal_Polys;
with text_io; use text_io;
with Ada.Numerics.Generic_Elementary_Functions;

procedure opolys_tst_1 is  -- fit data with a nth order polynomial

   type Real is digits 15;

   package Math is new Ada.Numerics.Generic_Elementary_Functions (Real);

   function Sin (X : Real) return Real renames Math.Sin;
   function Cos (X : Real) return Real renames Math.Cos;

   type Points_Index is range 0..81;
   --  This determines the maximum storage size of the Data array.
   --  the Data array will start at Points_Index'First and end
   --  at Points_Index'Last.
   --  Least squares fitting can be done on any sub interval of
   --  this array by setting the obvious parameters in the procedure
   --  below.

    type Data is array (Points_Index) of Real;

    package LS is new Orthogonal_Polys (Real, Points_Index, Data, 0.0, 1.0);
    --                          "*", "+", "-", "/", "<", "**", "-");
    use LS;

    X_axis  : Data;
    Y_Data  : Data;
    Weights : constant Data := (others => 1.0);

    F : constant Points_Index := Data'First;
    L : constant Points_Index := Data'Last;
    Poly_A, Poly_0, Poly_1, Poly_2 : Poly_Data; -- already initialized.
    Norm01, Norm12, Norm02 : Real;
    NormA2 : Real;
    Delta_X : Real;
    Order_Of_Deriv : Coeff_Index;

    Q1, Q2, Q3, Q4, Max1, Max2, Max3, d1, d2, d3 : Real;

     -- Variables for IO:
    Answer            : character;

    X, Y : Real;
    X0, Y0: Real;
    X1, Y1 : Real;
    Integral_At_0 : Real := 0.0;

    Factorial      : Real := 1.0;
    E              : Real;
    Order          : Real;
    Degree_Of_Poly : Coeff_Index;
    Best_Fit       : Poly_Data;
    Poly_Set       : Polynomials;
    C              : Poly_Sum_Coeffs;
    Derivatives    : Derivative_List;
    Deg_Of_Test_Poly_Plus : Coeff_Index;

    X_Coeffs, X_Coeffs_2 : Powers_Of_X_Coeffs;

    First   : constant Points_Index := Points_Index'First;
    Last    : constant Points_Index := Points_Index'Last;
    X_Start : constant Real := -Real (2.0);
    X_end   : constant Real := 2.0;
    Degree_Of_Test_Poly : constant Coeff_Index := Coeff_Index'(5);

    package rio is new Text_IO.FLoat_IO(Real);
    use rio;

    ---------
    -- Sum --    
    ---------
    
    function Sum (X : Real; 
                 Coeffs : Powers_Of_X_Coeffs) return Real is
      Result : Real := 0.0;
    begin
      Result := Coeffs(Coeff_Index'Last);
      for I in reverse Coeff_Index'First..Coeff_Index'Last-1 loop
         Result := Result * X;
         Result := Result + Coeffs(I);
      end loop;

      return Result;
   end Sum;

   -----------
   -- Pause --    
   -----------

   procedure Pause (s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13 : string := "") is
      Continue : Character := ' ';
   begin
      new_line;
      if S0 /= "" then put_line (S0); end if;
      if S1 /= "" then put_line (S1); end if;
      if S2 /= "" then put_line (S2); end if;
      if S3 /= "" then put_line (S3); end if;
      if S4 /= "" then put_line (S4); end if;
      if S5 /= "" then put_line (S5); end if;
      if S6 /= "" then put_line (S6); end if;
      if S7 /= "" then put_line (S7); end if;
      if S8 /= "" then put_line (S8); end if;
      if S9 /= "" then put_line (S9); end if;
      if S10 /= "" then put_line (S10); end if;
      if S11 /= "" then put_line (S11); end if;
      if S12 /= "" then put_line (S12); end if;
      if S13 /= "" then put_line (S13); end if;
      new_line;
      begin
         Put ("Press a key to continue: ");
         Get_Immediate (Continue); New_Line;
      exception
         when others => null;
      end;
   end pause;

    --------------------
    -- Make_True_Poly --
    --------------------

    -- Make a polynomial of order Poly_Order.  All of the
    -- the coefficients of powers of X are 1.0.
    --     Y = 15.0 + 14.0*X + 13.0*X**2 + 12.0*X**3 + ...
    -- Below we fit least-squares polynomials to this data
    -- in order to test package Least_Squares_Poly.

    procedure Make_True_Poly
      (Poly_Order    : in  Coeff_Index;
       Coeffs        : in  Powers_Of_X_Coeffs;
       X_Axis, Data1 : out Data;
       X_Start,X_end : in  Real)
    is
       X, Y    : Real;
       Delta_X : constant Real :=
         (X_end - X_Start) / Real(Points_Index'Last - Points_Index'First);

    begin

    for I in Points_Index loop

       X := X_Start + Delta_X * Real (I - Points_Index'First);
       Y := Coeffs (Poly_Order);    -- highest order coefficient

       if Poly_Order > 0 then
       for Power_Of_X in reverse Coeff_Index range 1..Poly_Order loop
           Y     := Y * X;
           Y     := Y + Coeffs (Power_Of_X - 1);
       end loop;
       end if;

       X_axis(I) := X;
       Data1 (I) := Y;

    end loop;

   end Make_True_Poly;

 begin

    Delta_X := (X_end - X_Start) / Real (X_axis'Length - 1);
    for I in Points_Index loop
       X_axis(I) := X_Start + Real (I - Points_Index'First) * Delta_X;
    end loop;

    -------------
    -- Test 1. --
    -------------

    Pause ("Test 1.  Check orthogonality of each polynomial with the three",
           "preceding polys.  Print out inner product... should be about 0.");

      begin

      Start_Gram_Schmidt_Poly_Recursion
        (X_axis, Weights, Data'First, Data'Last, Poly_0, Poly_1, Poly_Set);

       Norm01 := Inner_Product (Poly_0.Points, Poly_1.Points, F, L, Weights);
       new_line;
       put (Norm01);
       new_line;

       for I in Coeff_Index range 2..Coeff_Index'Last loop

         Get_Next_Poly (Poly_0, Poly_1, Weights, poly_2, Poly_Set);

         NormA2 := Inner_Product (Poly_A.Points, Poly_2.Points, F, L, Weights);
         Norm02 := Inner_Product (Poly_0.Points, Poly_2.Points, F, L, Weights);
         Norm12 := Inner_Product (Poly_1.Points, Poly_2.Points, F, L, Weights);
         new_line;
         --put (Norm2 / Poly_2.Squared);
         put (NormA2); put(" "); put (Norm12); put (" "); put (Norm02);

         Poly_A := Poly_0;
         Poly_0 := Poly_1;
         Poly_1 := Poly_2;

       end loop;

      exception
         when others => 
         put_line ("Some error.");
      end;

    Pause ("");

    -------------
    -- Test 2. --
    -------------

    Pause ("Test 2 compares various calculations of poly values at the X-axis",
      "grid points.  Four different methods are used.",
      "0: Use the original polynomial - the values stored in Poly_Data.",
      "1: Calculate the 0th order derivative using the derivative function.", 
      "2: Call Poly_Value (which also uses the Clenshaw summation formula).",
      "3: Calculate the coefficients of powers of X, and then sum the polynomial",
      "   using these coefficients. The differences will be printed below.",
      " ",
      "The first column will print the differences between method 0 and method",
      "1 above.  The second column prints differences between 0 and 2 above.",
      " ",
      "The 3rd column prints the differences between method 0 above, and the",
      "value calculated by summing over coefficients of powers of X.",
      "You'll notice that method 3 is a terrible way to evaluate polynomials.");

      Start_Gram_Schmidt_Poly_Recursion
        (X_axis, Weights, Data'First, Data'Last, Poly_0, Poly_1, Poly_Set);

      for I in Coeff_Index range 2..Coeff_Index'Last loop

         Get_Next_Poly (Poly_0, Poly_1, Weights, poly_2, Poly_Set);

         Poly_0 := Poly_1;
         Poly_1 := Poly_2;

         -- In what follows look at poly value at each grid point X_axis(I)

         C              := (others => +0.0);
         C(I)           := +1.0;
         Get_Coeffs_Of_Powers_Of_X (X_Coeffs, C, Poly_Set);

         Max1 := 0.0; Max2 := 0.0; Max3 := 0.0;
         for J in Points_Index loop

              X              := X_axis(J);
              Order_Of_Deriv := 0;

              Poly_Derivatives (Derivatives, X, Order_Of_Deriv, C, Poly_Set);
              Q1 := Derivatives(0);

              Q2 := Poly_Value (X, C, Poly_Set); 

              Q3 := Sum (X, X_Coeffs);

              Q4 := Poly_2.Points(J);

              d1 := Q4 - Q1; d2 := Q4 - Q2; d3 := Q4 - Q3;
              if Abs(d1) > Max1 then Max1 := Abs(d1); end if;
              if Abs(d2) > Max2 then Max2 := Abs(d2); end if;
              if Abs(d3) > Max3 then Max3 := Abs(d3); end if;
              --new_line;
              --put(d1); put(" "); put(d2); put(" "); put(d3);

         end loop;

         new_line;
         put(Max1); put(" "); put(Max2); put(" "); put(Max3);

      end loop;

    -------------
    -- Test 1  --
    -------------

    New_Line;
    Pause ("New series of tests: create a test polynomial, calculate",
           "series of least squares polynomial fits to it, and print the comparison.");

    Pause ("Test 1: we create a polynomial of order 14, with coefficients",
       "equal to 16.0, 15.0, 14.0 etc.: Y = 16.0 + 15.0*X + 14.0*X**2 + ...",
       "In what follows, the user fits a",
       "least squares polynomial of some given order to the data.  As the user",
       "raises the order of the polynomial fit, the mean square error between",
       "the least-squares polynomial fit and the data should fall.  It",
       "should suddenly fall by a factor 10**23 at order 5 and higher.");

    loop
      begin

      put ("Input order of polynomial fit (try 6 first, and exit):");
      New_Line; Get (Order);
      Degree_Of_Poly  :=  Coeff_Index (Order);

      X_Coeffs_2 := (others => 0.0);
      for I in Coeff_Index range 0..Degree_Of_Test_Poly loop
          X_Coeffs_2 (I) := 16.0 - Real(I);
      end loop;

      Make_True_Poly
        (Poly_Order => Degree_Of_Test_Poly,
         Coeffs     => X_Coeffs_2,
         X_axis     => X_axis,
         Data1      => Y_Data,
         X_Start    => X_Start,
         X_end      => X_end);

      Poly_Fit (Y_Data, X_axis, Weights, First, Last, Degree_Of_Poly,
                                              Best_Fit, C, Poly_Set, E);

      put("Mean square error = "); put(E); New_Line;
      put_line ("Exit loop? if Yes, enter a 'y', else a 'n' ");
      get (Answer);
      if Answer = 'y' then
          Exit;
      end if;

      exception
         when others =>
         put_line ("Some error.  Try again: ");
      end;
    end loop;

    -------------
    -- Test 2  --
    -------------

    pause ("Test 2: print the differences between the LS fit and the data.",
           "This test compares the original data with the least squares fit",
           "polynomial returned by procedure Poly_Fit.");

    --  Print differences:

    for I in Points_Index loop
       New_Line; Put (Best_Fit.Points(I) - Y_Data(I)); put ("  ");
    end loop;

    -------------
    -- Test 3  --
    -------------

    pause;
    pause("Test 3: print the differences between the Poly and the Data.",
          "This is a test of the function Poly_Value, which uses Clenshaw's",
          "formula to sum the least-squares polynomial.  This is calculated",
          "in a different way from the polynomial test previously.  Also,",
          "we redo the least-squares fit to make sure that the order of the",
          "least-squares poly is greater than the order of the polynomial",
          "used to generate the test data.");

    if Degree_Of_Test_Poly < Coeff_Index'Last then
       Deg_Of_Test_Poly_Plus :=  Degree_Of_Test_Poly + 1;
    else
       Deg_Of_Test_Poly_Plus :=  Degree_Of_Test_Poly;
    end if;

    Poly_Fit (Y_Data, X_axis, Weights, First, Last, Deg_Of_Test_Poly_Plus,
             Best_Fit, C, Poly_Set, E);

    for I in Points_Index loop

       X := X_axis(I);
       Y := Poly_Value (X, C, Poly_Set);
       New_Line; Put (Y - Y_Data(I)); put ("  ");

    end loop;

    ------------
    -- Test 4 --
    ------------

    pause;
    pause("Test 4: attempt to get the coefficients of the powers of X",
          "in the least-squares fit polynomial.  This is a process that",
          "can result in a large loss of precision, underflows, overflows",
          "etc. if the degree of the polynomial is high.  In the present",
          "case, a 14th order poly, things work out OK if you sample the",
          "poly over a sufficiently large interval in X.");

    Get_Coeffs_Of_Powers_Of_X (X_Coeffs, C, Poly_Set);

    New_Line;
    for I in Coeff_Index loop
      New_Line; put (X_Coeffs(I));
    end loop;

    -------------
    -- Test 5  --
    -------------

    pause;
    pause("Test 5: take the first derivative of the polynomial,",
          "and compare with the predictions of the analytical",
          "derivative.");

    X_Coeffs_2 := (others => 0.0);
    Y_Data     := (others => 0.0);

    if Degree_Of_Test_Poly > 0 then
       for I in Coeff_Index range 0..Degree_Of_Test_Poly-1 loop
         X_Coeffs_2 (I) := Real(I+1) * (16.0 - Real(I+1));
       end loop;

       Make_True_Poly (Poly_Order   => Degree_Of_Test_Poly-1,
                           Coeffs   => X_Coeffs_2,
                             X_axis => X_axis,
                            Data1   => Y_Data,
                            X_Start => X_Start,
                            X_end   => X_end);
    end if;

    for I in Points_Index loop

      X := X_axis(I);
      Poly_Derivatives (Derivatives, X, 1, C, Poly_Set);
      New_Line; Put (Derivatives(1) - Y_Data(I)); put ("  ");

    end loop;

    -------------
    -- Test 6  --
    -------------

    pause;
    pause("Test 6: take the second derivative of the polynomial,",
          "and compare with the predictions of the analytical",
          "derivative.");

    X_Coeffs_2 := (others => 0.0);
    Y_Data := (others => 0.0);

    if Degree_Of_Test_Poly > 1 then
       for I in Coeff_Index range 0..Degree_Of_Test_Poly-2 loop
          X_Coeffs_2 (I) := Real(I+2) * Real(I+1) * (16.0 - Real(I+2));
       end loop;

       Make_True_Poly 
         (Poly_Order => Degree_Of_Test_Poly-2,
          Coeffs     => X_Coeffs_2,
          X_axis     => X_axis,
          Data1      => Y_Data,
          X_Start    => X_Start,
          X_end      => X_end);

    end if;

    for I in Points_Index loop

      X := X_axis(I);
      Poly_Derivatives (Derivatives, X, 2, C, Poly_Set);
      New_Line; Put (Derivatives(2) - Y_Data(I)); put ("  ");

    end loop;

    -------------
    -- Test 7  --
    -------------

    pause;
    pause("Test 7: take the third derivative of the polynomial,",
          "and compare with the predictions of the analytical",
          "derivative.");

    X_Coeffs_2  := (others => 0.0);
    Y_Data := (others => 0.0);

    if Degree_Of_Test_Poly > 2 then
       for I in Coeff_Index range 0..Degree_Of_Test_Poly-3 loop
         X_Coeffs_2 (I) :=
            Real(I+3) * Real(I+2) * Real(I+1) * (16.0 - Real(I+3));
       end loop;

       Make_True_Poly
         (Poly_Order => Degree_Of_Test_Poly-3,
          Coeffs     => X_Coeffs_2,
          X_axis     => X_axis,
          Data1      => Y_Data,
          X_Start    => X_Start,
          X_end      => X_end);

    end if;

    for I in Points_Index loop

      X := X_axis(I);
      Poly_Derivatives (Derivatives, X, 3, C, Poly_Set);
      New_Line; Put (Derivatives(3) - Y_Data(I)); put ("  ");

    end loop;

    -------------
    -- Test 8  --
    -------------

    pause;
    pause("Test 8: check that the Derivative routine also gets the",
          "correct zeroth order derivative:",
          "Below are the differences between the predictions of procedure",
          "Poly_Derivatives, and the least-squares fit polynomial.");
    for I in Points_Index loop

      X := X_axis(I);
      Poly_Derivatives (Derivatives, X, 2, C, Poly_Set);
      New_Line; Put (Derivatives(0) - Best_Fit.Points(I)); put ("  ");

    end loop;

    -------------
    -- Test 9  --
    -------------

    pause;
    pause("Test 9: check that the Derivative routine is able to correctly",
          "calculate the coefficients of powers of X of the polynomial.  We",
          "use the fact that the coefficient of the m-th power of X equals",
          "m-th derivative of the least squares fit polynomial evaluated at",
          "X = 0, and divided by m!.");

    Factorial := 1.0;
    Poly_Derivatives (Derivatives, 0.0, Coeff_Index'Last, C, Poly_Set);
    for m in Coeff_Index loop
      if m > Coeff_Index'First then
         Factorial := Factorial * Real(m);
      end if;
      New_Line; put (Derivatives(m) / Factorial);
    end loop;

    -------------
    -- Test 10 --
    -------------

    pause;
    pause("Test 10: take the indefinite Integral of the polynomial,",
          "and compare with the predictions of the analytical",
          "integral.");

    X_Coeffs_2 := (others => 0.0);
    Y_Data     := (others => 0.0);

    if Degree_Of_Test_Poly > 0 then
       for I in reverse Coeff_Index range 0..Degree_Of_Test_Poly loop
          X_Coeffs_2 (I+1) := (16.0 - Real(I)) / Real(I+1);
       end loop;
       X_Coeffs_2 (0) := (0.0);

       Make_True_Poly (Poly_Order   => Degree_Of_Test_Poly+1,
                           Coeffs   => X_Coeffs_2,
                             X_axis => X_axis,
                            Data1   => Y_Data,
                            X_Start => X_Start,
                            X_end   => X_end);
    end if;

    Integral_At_0 := Poly_Integral (0.0, C, Poly_Set);

    for I in Points_Index loop

      X := X_axis(I);
      Y := Poly_Integral (X, C, Poly_Set) - Integral_At_0;
      New_Line; Put (Y - Y_Data(I)); put ("  ");

    end loop;

    -------------
    -- Test 1  --
    -------------
 
    New_Line;
    Pause ("New series of tests: the previous series was",
         "designed to discover gross errors, but it didn't say much about",
         "the numerical accuracy of the fits, because the error may have",
         "been in the calculation of the original data. In the following",
         "series of tests we generate data with Sin's and Cos's so that we",
         "can reliably measure numerical error in the least squares fits");

    Delta_X := 16.0 / Real (X_axis'Length - 1);
    for I in Points_Index loop
      X_axis(I) := Real (I - Points_Index'First) * Delta_X;
    end loop;

    for I in Points_Index loop
       Y_Data(I) := Sin (X_axis(I));
    end loop;

    Poly_Fit (Y_Data, X_axis, Weights, First, Last, 36,
             Best_Fit, C, Poly_Set, E);

    pause("Test 1: Make sinusoidal data, fit a high order poly to it,",
          "and print the difference between the fit and the original data:");

     for I in Points_Index loop

        X := X_axis(I);
        Y := Poly_Value (X, C, Poly_Set);
        new_line; Put (Y - Y_Data(I)); put ("  ");

     end loop;

    pause;
    pause("Test 2: Make sinusoidal data, fit a high order poly to it,",
          "and print the difference between the 1st derivative of the fit",
          "and the 1st derivative of the original data:");

     for I in Points_Index loop

        X := X_axis(I);
        Poly_Derivatives (Derivatives, X, 1, C, Poly_Set);
        new_line; Put (Derivatives(1) - Cos(X_axis(I))); put ("  ");

     end loop;

    pause;
    pause("Test 3: Make sinusoidal data, fit a high order poly to it,",
          "and print the difference between the 2nd derivative of the fit",
          "and the 2nd derivative of the original data:");

     for I in Points_Index loop

        X := X_axis(I);
        Poly_Derivatives (Derivatives, X, 2, C, Poly_Set);
        new_line; Put (Derivatives(2) + Sin(X_axis(I))); put ("  ");

     end loop;

    pause;
    pause("Test 4: Make sinusoidal data, fit a high order poly to it,",
          "and print the difference between the 3rd derivative of the fit",
          "and the 3rd derivative of the original data:");

     for I in Points_Index loop

        X := X_axis(I);
        Poly_Derivatives (Derivatives, X, 3, C, Poly_Set);
        new_line; Put (Derivatives(3) + Cos(X_axis(I))); put ("  ");

     end loop;

    pause;
    pause("Test 5: Make sinusoidal data, fit a high order poly to it,",
          "and print the difference between the integral of the fit",
          "and the integral of the original data:");

    Integral_At_0 := Poly_Integral (0.0, C, Poly_Set);

     for I in Points_Index'first+1 .. Points_Index'last loop

        X0 := X_axis(I-1);
        Y0 := Poly_Integral (X0, C, Poly_Set);
        X1 := X_axis(I);
        Y1 := Poly_Integral (X1, C, Poly_Set);
        new_line; put (Y1 - Y0 - (-Cos (X1) + Cos (X0)) ); put ("  ");

     end loop;

end opolys_tst_1;
