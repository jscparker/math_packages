
with Tridiagonal_LU;
With Text_IO; use Text_IO;

procedure tridiag_tst_1 is

   type Real is digits 15;

   package rio is new Text_IO.Float_IO(Real);
   use rio;
 
   Type Index is range 1..40;
 
   Package Tri is new Tridiagonal_LU (Real, Index);
   use Tri;
 
   SolutionVector : Column;
   D            : Matrix := (others => (others => 0.0));
   A            : Matrix;
   UnitVector   : Column;
   ZeroVector   : constant Column := (others => 0.0);
   Index_First  : Index := Index'First;
   Index_Last   : Index := Index'Last;
   RealMaxIndex : Real;
   Test         : Real;
 
begin

    Put("Input First Index Of Matrix To Invert. (e.g. 2)"); New_Line;
    get(RealMaxIndex);
    Index_First := Index (RealMaxIndex);

    Put("Input Last Index Of Matrix To Invert. (e.g. 8)"); New_Line;
    get(RealMaxIndex);
    Index_Last := Index (RealMaxIndex);

    -- if Banded_Matrix_Desired then
    -- Construct a banded matrix:
    for I in Index loop
       D(0)(I) := 0.3;
    end loop;
    for Row in Index'First..Index'Last loop
       D(1)(Row) := 0.9;
    end loop;
    for Row in Index'First..Index'Last loop
       D(-1)(Row) := 0.5;
    end loop;

    A := D;
    LU_Decompose (A, Index_First, Index_Last);

    -- Get Nth column of the inverse matrix
    -- Col       := N;

    Put("Output should be the Identity matrix."); new_Line;

    for Col in Index range Index_First..Index_Last loop

       UnitVector      := ZeroVector;
       UnitVector(Col) := 1.0;
       Solve (SolutionVector, A, UnitVector, Index_First, Index_Last);
 
       New_Line;

       -- Multiply SolutionVector times tridiagonal D:

       for Row in Index_First..Index_Last loop
 
          Test := 0.0;
          if Row > Index_First then
             Test :=  Test + D(-1)(Row) * SolutionVector(Row-1);
          end if;
 
          Test :=  Test + D(0)(Row) * SolutionVector(Row);
          if Row < Index_Last then
             Test :=  Test + D(1)(Row) * SolutionVector(Row+1);
          end if;
 
          Put (Test, 2, 3, 3); Put (" ");
       end loop;

    end loop;

end;
