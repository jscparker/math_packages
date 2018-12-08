
-- Test QR least squares equation solving real valued *square* matrices.

with Ada.Numerics.Generic_elementary_functions;
with Givens_QR;
with Test_Matrices;
With Text_IO; use Text_IO;

procedure givens_qr_tst_2 is

   type Real is digits 15;

   subtype Index is Integer range 1..191;

   Start_Index : constant Index := Index'First + 0;
   Max_Index   : constant Index := Index'Last - 0;


   subtype Row_Index is Index;
   subtype Col_Index is Index;

   Starting_Row : constant Row_Index := Start_Index;
   Starting_Col : constant Col_Index := Start_Index;

   Final_Row : constant Row_Index := Max_Index;
   Final_Col : constant Col_Index := Max_Index;

   type Matrix is array(Row_Index, Col_Index) of Real;

   type Matrix_inv is array(Col_Index, Row_Index) of Real;
   --  For inverses of A : Matrix; has shape of A_transpose.

 --pragma Convention (Fortran, Matrix); --No! This QR prefers Ada convention.

   package math is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use math;

   package QR is new Givens_QR
     (Real     => Real, 
      R_Index  => Index, 
      C_Index  => Index, 
      A_Matrix => Matrix);
   use QR;

   -- QR exports Row_Vector and Col_Vector

   package Make_Square_Matrix is new Test_Matrices (Real, Index, Matrix);
   use Make_Square_Matrix;

   package rio is new Float_IO(Real);
   use rio;

 --subtype Longer_Real  is Real;     -- general case, and for best speed
 --type Longer_Real  is digits 18;   -- 18 ok on intel, rarely useful

   Min_Real :constant Real := 2.0 ** (Real'Machine_Emin/2 + Real'Machine_Emin/4);

   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;

   --  Used by Frobenius_Norm and by Get_Err in reassembled A

   --------------------
   -- Frobenius_Norm --
   --------------------
  
   function Frobenius_Norm 
     (A            : in     Matrix)
    --Final_Row    : in     Index; 
    --Final_Col    : in     Index;
    --Starting_Row : in     Index; 
    --Starting_Col : in     Index)
      return Real
    is
      Max_A_Val : Real := Zero;
      Sum, Scaling, tmp : Real := Zero;
    begin
 
      Max_A_Val := Zero;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         if Max_A_Val < Abs A(Row, Col) then Max_A_Val := Abs A(Row, Col); end if;
      end loop;
      end loop;

      Max_A_Val := Max_A_Val + Two ** (Real'Machine_Emin + 4);
      Scaling := One / Max_A_Val;

      Sum := Zero;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * A(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;


   ------------
   -- Invert --
   ------------
  
   --  Get Inverse of the Matrix:
   
   procedure Get_Inverse 
     (R              : in     Matrix;
      Q              : in     Q_Matrix;
      Scale          : in     Col_Vector;
      Permute        : in     Permutation;
      A_inverse      :    out Matrix_inv)     -- shape is A_transpose
   is
      Solution_Row_Vector : Row_Vector; --  X
      Unit_Col_Vector     : Col_Vector := (others => 0.0); --  B
      --  We are going to solve for  X  in  A*X = B

      Final_Row    : Row_Index renames Q.Final_Row;
      Final_Col    : Col_Index renames Q.Final_Col;
      Starting_Row : Row_Index renames Q.Starting_Row;
      Starting_Col : Col_Index renames Q.Starting_Col;
   begin
   
      A_inverse := (others => (others => Zero));

      --new_line;
      --for i in Starting_Row_Col..Max_Row_Col loop
         --put(R(i,i));
      --end loop;
   
      for i in Starting_Row .. Final_Row loop

         if i > Starting_Row then 
            Unit_Col_Vector(i-1) := 0.0;
         end if;
         Unit_Col_Vector(i) := 1.0;
         --  Make all possible unit Col vectors:  (Final_Row-Starting_Row+1) of them


         QR_Solve 
           (X                  => Solution_Row_Vector,
            B                  => Unit_Col_Vector,
            R                  => R,
            Q                  => Q,
            Row_Scalings       => Scale,
            Col_Permutation    => Permute);
         --  Solve equ. A*Solution_Row_Vector = Unit_Col_Vector (for Solution_Row...).
	 --  Q contains the Starting_Row, Final_Row, Starting_Col, etc.

         for j in Starting_Col .. Final_Col  loop
            A_Inverse (j,i) := Solution_Row_Vector(j);
         end loop;
	 --  All vectors you multiply by A must be Row_Vectors, 
	 --  So the cols of A_inv are Row_Vectors.  

      end loop;
   
   end Get_Inverse;

   procedure Get_Err_in_Least_Squares_A
     (A              : in     Matrix;
      R              : in     Matrix;
      Q              : in     Q_Matrix;
      Scale          : in     Col_Vector;
      Permute        : in     Permutation;
      Err_in_Least_Squ_A        : out Real;
      Condition_Number_of_Least_Squ_A : out Real;
      Condition_Number_of_A           : out Real)
   is
      Max_Diag_Val_of_R : constant Real := Abs R(Starting_Col, Starting_Row);
      Min_Allowed_Diag_Val_of_R : constant Real := 
                              Max_Diag_Val_of_R * Default_Singularity_Cutoff;
      Min_Diag_Val_of_Least_Squ_R : Real;

      Least_Squares_Truncated_Final_Row : Row_Index := Final_Row; -- essential init

      Product_Vector, Col_of_R : Col_Vector := (others => Zero);
      Row : Row_Index;
      Err_Matrix : Matrix := (others => (others => Zero));

      Final_Row    : Row_Index renames Q.Final_Row;
      Final_Col    : Col_Index renames Q.Final_Col;
      Starting_Row : Row_Index renames Q.Starting_Row;
      Starting_Col : Col_Index renames Q.Starting_Col;
   begin

      -- The Columns of R have been permuted; unpermute before comparison of A with Q*R
      -- The Columns of R have been scaled. Before comparison of A with Q*R
      -- Must unscale each col of R by multiplying them with 1/Scale(Row).

      Row := Starting_Row;
      Find_Truncated_Final_Row:
      for Col in Starting_Col .. Final_Col loop
         Min_Diag_Val_of_Least_Squ_R := Abs R(Row, Col);
         if Min_Diag_Val_of_Least_Squ_R < Min_Allowed_Diag_Val_of_R then
	    Least_Squares_Truncated_Final_Row := Row - 1;
            Min_Diag_Val_of_Least_Squ_R := Abs R(Row-1, Col-1);
	    exit Find_Truncated_Final_Row;
	 end if;
         if Row < Final_Row then Row := Row + 1; end if;
      end loop Find_Truncated_Final_Row;

      for Col in Starting_Col .. Final_Col loop

         Col_of_R := (others => Zero);

         for Row in Starting_Row .. Least_Squares_Truncated_Final_Row loop
            Col_of_R(Row) := R(Row, Col);
	 end loop;

         Product_Vector := Q_x_Col_Vector (Q, Col_of_R);

         for Row in Starting_Row .. Final_Row loop
            Err_Matrix(Row, Col) := Abs (A(Row, Permute(Col)) - 
                --Product_Vector(Row) / (Scale(Permute(Col)) + Min_Real));
                  Product_Vector(Row) / (Scale(Row) + Min_Real));
         end loop;
      end loop;

      --  Froebenius norm fractional error = ||Err_Matrix|| / ||A||

      Err_in_Least_Squ_A := 
         Frobenius_Norm (Err_Matrix) / (Frobenius_Norm (A) + Min_Real); 

      Condition_Number_of_Least_Squ_A := 
         Max_Diag_Val_of_R / (Min_Diag_Val_of_Least_Squ_R + Min_Real);

      Condition_Number_of_A := 
         Max_Diag_Val_of_R / (Abs R(Final_Row, Final_Col) + Min_Real);

   end Get_Err_in_Least_Squares_A;

   -----------
   -- Pause --
   -----------

   procedure Pause (s0,s1,s2,s3,s4,s5,s6,s7,s8,s9 : String := "") is
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
     new_line;
     begin
	put ("Type a character to continue: ");
	get_immediate (Continue);
     exception
	when others => null;
     end;
  end pause;


  -------------------------
  --  Print_Err_in_LSQ_A --
  -------------------------

  procedure Print_Err_in_LSQ_A 
    (Chosen_Matrix : Matrix_id)
  is
     A_inv : Matrix_inv;
     A, R : Matrix;
     Q : Q_Matrix;
     Scale   : Col_Vector;
     Permute : Permutation;
     Err_in_Least_Squ_A, Condition_Number_of_Least_Squ_A, Condition_Number_of_A : Real;
  begin
     Init_Matrix (A, Chosen_Matrix, Start_Index, Max_Index);

     R := A; -- leave A unmolested; R will be returned as upper rt. triangular
 
     QR_Decompose 
       (A                  => R,  -- input matrix A, but call it R.
        Q                  => Q,
        Row_Scalings       => Scale,
        Col_Permutation    => Permute,
        Final_Row          => Final_Row,
        Final_Col          => Final_Col,
        Starting_Row       => Starting_Row,
        Starting_Col       => Starting_Col);

     --new_line;
     --for Row in Starting_Row .. Final_Row loop
       --put (R(Row, Col_Index (Starting_Col + Row - Starting_Row)));
     --end loop;

     Get_Inverse 
       (R, Q, Scale, Permute, A_inv);

     Get_Err_in_Least_Squares_A
       (A, R, Q, Scale, Permute,
        Err_in_Least_Squ_A, Condition_Number_of_Least_Squ_A, Condition_Number_of_A);

     new_line;
     put ("For matrix A of type  "); put(Matrix_id'Image(Chosen_Matrix));  put(":");

     new_line;
     put(" Estimate (usually low) of condition number of original A       ="); 
     put(Condition_Number_of_A);
     new_line;
     put(" Estimate (usually low) of condition number of least squares A  ="); 
     put(Condition_Number_of_Least_Squ_A);
     new_line;
     put(" How much A changed in taking least squares: ||A - A_ls||/||A|| ="); 
     put(Err_in_Least_Squ_A);
     new_line;

  end Print_Err_in_LSQ_A;

begin

   Pause( 
      "The following test demonstrates least squares equation solving using QR",
      "decomposition of A.  Least squares happens naturally in QR, because in",
      "obtaining A = Q*R we have written the columns of A as linear combinations",
      "of a complete set of orthonormal vectors (the Cols of Q). If you zero out",
      "one or more of the Cols of Q, and recalculate A = Q*R, then you have a",
      "least squares approximation of the columns of A in the basis spanned by the",
      "remaining cols of Q. For lack of a better term, let's call this the least",
      "squares A. We do this below explicitly. The idea is to show that you get a",
      "least squares version of A that differs from A by only 1 part in (say) 10**11,",
      "but has vastly better properties in equation solving (lower condition number)."
     );

   Pause( 
      "Below we print condition numbers of the original A, and the least squares A.",
      "The condition numbers are estimates obtained from the diagonal elements of R.",
      "They are rough estimates (you need SVD for accurate values), but are ok for",
      "the present demonstration.  The condition number of the least squares A is", 
      "that of the reduced rank matrix used in least squares equation solving.",
      "(And they are condition numbers of column scaled matrices, not the raw",
      "matrices emitted by package Matrix_Sampler.  SVD uses the raw matrices.)"
      );


   Print_Err_in_LSQ_A (Chosen_Matrix => KAHAN);

   Pause( 
      "We started with an ideal case: A is nearly singular (condition number > 10**18).",
      "You would need 23 digit arithmetic to get solutions correct to 4 digits. The",
      "least squares A is well-conditioned (condition number ~ 8), and notice we",
      "only had to change A by 2 parts in 10**16 to get a well-conditioned matrix."
      );

   Print_Err_in_LSQ_A (Chosen_Matrix => RANDOM_32_bit);

   Pause( 
      "All the elements in Random are generated by a random num generator. Its",
      "condition number should be considered typical. It is not ill-conditioned,",
      "and it was not modified by the least squares process.  Least squares fitting,",
      "only occurs when condition number exceeds some threshold (here about 10**11)."
      );


   Print_Err_in_LSQ_A (Chosen_Matrix => Hilbert);

   Pause( 
      "Hilbert's matrix is a worst case matrix. The condition number exceeds 10**11",
      "so least squares is performed.  The distance between A and least squares A",
      "is now significant, (the maximum likely to be observed).  The decrease in",
      "condition number is as predicted, but not as good as some other cases."
      );

   Print_Err_in_LSQ_A (Chosen_Matrix => Pascal);

   Pause( 
      "Pascal's matrix is another worst case. The condition number exceeds 10**11",
      "so the matrix is modified.  The difference between A and the least squares A",
      "is now significant, (the maximum likely to be observed).  The decrease in",
      "condition number is as predicted, but not as good as some other cases."
      );

   Pause( 
      "We now run through all the matrices in the test suite.  Most are singular,",
      "ill-conditioned, badly scaled, or several of the above. The least squares",
      "process described above is the default, because when the matrix is singular",
      "(infinite condition number, or 2**231 here) then it is the only practical",
      "choice. Badly behaved matrices are likely to occur sporadically in most",
      "problem domains, but in some problem domains they are the norm."
      );


   for Chosen_Matrix in Matrix_id loop

      Print_Err_in_LSQ_A (Chosen_Matrix);

   end loop;

end;
