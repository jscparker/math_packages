
-- PACKAGE Fourier8
--
-- Standard Radix 8 (decimation in frequency) Cooley-Tukey FFT.
-- Procedure FFT does a Discrete Fourier Transform on data sets
-- that are powers-of-2 in length. The package is pure.
--
-- Radix 2 FFT's read and write the entire data set to the data array
-- log(N) times per transform (where N is the data length).
-- Radix 2**k FFT's read and write the entire data set to the data array
-- log(N) / k times per transform (where N is the data length).
-- These reads and writes slow down the calculation significantly,
-- especially as arrays get large and spill out of the machine's
-- cache ram. So the higher Radix FFT's are usually quite a bit faster
-- than the simpler Radix 2 FFT's when the data sets get large.
-- (The code, on the other hand, is not pretty!)
--
--
-- 2. Notes on use.
--
-- If length of the data set to be transformed is not a power
-- of 2, then it is padded with Zeros to make it a power of 2.
--
-- Data is stored as 2 arrays of Real numbers, one for the Real part,
-- one for the complex part.
--
-- Only data in the range 0..N is transformed.
-- N is input as: Input_Data_Last.
-- If N+1 is not a power of 2, then the Data array will be padded
-- with zeros out to the nearest power of 2.
-- The FFT is then performed on the enlarged set of data.
--
-- More precisely if Input_Data_Last is not 2**N-1 for some integer N
-- then all data points in the range Input_Data_Last+1..2**M-1
-- *will be set to zero*,  i.e. data is padded with Zeros to the nearest
-- power of two.  The FFT will be performed on the interval 0..2**M-1,
-- where 2**M-1 = Transformed_Data_Last = smallest power of 2 (minus
-- one) that is greater than or equal to Input_Data_Last.
-- (You input the number of data points *minus one*, so that this number
-- can have the same type as the array index, which has range 0..2**N-1).
--
-- The user sets the maximum size of the storage array by
-- setting the generic parameter Log_Of_Max_Data_Length.
--
-- The user has the option of Normalizing the data by dividing
-- the results by SQRT (N), where N is the number of data points.
-- If you want normalization by 1/N instead, then choose the no
-- normalization option and perform the normalization yourself.
--
--
-- 3. Notes on the Discrete Fourier Transform (DFT)
--
-- The N point Discrete Fourier Transform F(W) of a function of time G(T):
--
--                        N-1
--   F(W_n) = DFT(N){G} = SUM { G(T_m) * Exp (-i * W_n * T_m)/Sqrt (N) },
--                        m=0
--
-- where the sum is over m and goes from 0 to N-1. (N is the
-- number of data points and is a power of 2.)  Notice that this
-- is the discretized version of the integral by which one calculates
-- the Fourier coefficients of a Fourier series.  In other words,
-- F(W_n) is the projection of the function G (which we wish to
-- write as a sum of the Fourier basis functions Exp) onto the
-- appropriate basis function.
--
-- Dividing by SQRT (N) normalizes the Fourier basis function Exp.
-- If this is done, then the Fourier coefficients F(W) give the
-- power spectrum: Power(W) = F(W) * Conjugate (F(W)).
--
-- The discrete components of time T_m are given by (Tmax/N)*m
-- where Tmax is the total time interval of the transform, and m
-- goes from 0 to N-1.  Similarly, the dicrete components of
-- frequency are W_n = (2*Pi/Tmax)*n where n goes from 0 to N-1.
-- Consequently the Dicrete Fourier Transform is
--
--             N-1
--    F(W_n) = SUM { G(T_m) * Exp (-i * (2*Pi/N) * n * m)/Sqrt(N) }.
--             m=0
--
-- In calculating F(W) we are finding the coefficients F(W_n) of
-- plane waves Exp (i * W_n * T_m) such that a sum over these plane
-- waves times the appropriate coefficient F(W_n) reproduces (approx.)
-- the funtion of time G(T). The plane waves Exp (i * W_n * T_m) are
-- said to form a complete set of states for functions G(T) that
-- are restricted to an interval [0,Tmax).  They satisfy periodic
-- boundary conditions: State (n, T) = State (n, T + Tmax).  It
-- follows from the periodic boundary conditions and from the
-- definition State (n, T) =  Exp (i * W_n * T) that the discrete
-- frequencies W_n have the form W_n = (2*Pi/Tmax)*n. Because the
-- the sum is truncated at N-1 the Fourier series won't exactly
-- approximate the function of time G(T).  Finally, notice that
-- the (un-normalized) states State (n, T) look like:
--
--  { 1, Exp (i*(2*Pi/Tmax)*T), Exp (i*(2*Pi/Tmax)*2*T), ... }
--
-- If we take the FFT of State (n, T) we should a get back a
-- delta function F(W) peaked at the nth data point W_n.  A test
-- procedure demonstrates this.
--
-- It is important to note from these definitions that we are
-- assuming that the function of time G(T) is restricted to an
-- interval of time [0, Tmax).  If we want to use the FFT to approximate
-- the integral of G(T)*Exp(i*W*T) on another time interval, say
-- [-Tmax/2, Tmax/2], then must transform coordinates of the integral
-- to the time interval [0,Tmax).  The
-- result of this transformation is to put an W (hence n) dependent phase
-- factor out front of the prediction of the FFT.  In fact, in the case
-- just described, the phase factor is
--    Exp (i * W_n * Tmax/2) = Exp (i * W_n * 2 * Pi * n / (W_n * 2)).
-- This equals Exp (i * Pi * n), which is +1 or -1 depending on whether
-- n is even of odd.
--
--
-- 4. Notes on the radix 8 fast fourier transform (FFT)
--    see below
--
--*****************************************************************
generic

   type Real is digits <>;
   type Array_Index is range <>;
   type Data_Array is array (Array_Index) of Real;

   Log_Of_Max_Data_Length : Positive;
   --  The FFT can only operate on range 0..2**N-1, where the maximum
   --  value N can have is N_max = Log_Of_Max_Data_Length.
   --
   --  The generic formal Data_Array has arbitrary limits.  This is
   --  useful, but the FFT algorithm works on a restricted range. It
   --  only operates on data in power-of-2 ranges: 0..2**N-1 for some N.
   --  The user can find out about this restriction at runtime or at
   --  compile time.  To make sure he finds out about it at compile time
   --  he is required to enter the Maximum value N is allowed to have.
   --  That value is N_max = Log_Of_Max_Data_Length.
   --  The package checks that 0..2**N_max-1 is in in the range of
   --  Array_Index, during compilation, and also exports a
   --  subtype of Array_Index that has range 0..2**N_max-1.  This
   --  subtype is called Data_Index, and all computation is done on it.

package Fourier8 is

   pragma Pure (Fourier8);

   Data_Index_Last : constant Array_Index
                          := Array_Index (2**Log_Of_Max_Data_Length-1);

   subtype Data_Index is Array_Index range 0..Data_Index_Last;

   --  Data_Index:  the index on which all computation is performed.  One
   --  of the reasons this subset is defined is to make sure that the generic
   --  formal type Array_Index contains the range 0..2**N-1
   --  where N is the log of the max number of data points.  The FFT
   --  is always done on some sub range of this: 0..M where M <= 2**N-1.

   type Exp_Storage is private;
   --  Must declare an object of this type and pass it to FFT as a parameter.
   --  That is all you ever have to do. The program does eveything else
   --  for you. (It contains the Sin's for the FFT ... keeps the package Pure.)
   --  So all you do is declare:
   --
   --     My_Exp_Table : Exp_Storage;
   --
   -- and then in the call to FFT, you add the parameter:
   --
   --     Exp_Table => My_Exp_Table,
   --

   procedure FFT
     (Data_Re, Data_Im        : in out Data_Array; --destroys original data
      Transformed_Data_Last   :    out Data_Index;
      Input_Data_Last         : in     Data_Index;
      Exp_Table               : in out Exp_Storage;
      Inverse_FFT_Desired     : in     Boolean     := False;
      Normalized_Data_Desired : in     Boolean     := False;
      Bit_Reversal_Desired    : in     Boolean     := True);


   --  FFT's the arrays Data_Re and Data_Im, and puts the results back
   --  into the same arrays.  The original data is destroyed.
   --
   --  The user inputs:  Input_Data_Last
   --
   --  The procedure performs an FFT on data that lies in the interval
   --
   --           0 .. Input_Data_Last.
   --
   --  Data beyond this point will be set to Complex_Zero out to the
   --  nearest power of 2: Transformed_Data_Last.
   --  The transformed data is returned in the array Data, in the range
   --
   --           0 .. Transformed_Data_Last.
   --
   --  Transformed_Data_Last+1 is a power of 2
   --  Input_Data_Last+1 need not be a power of 2

   pragma Inline (FFT);

private

   -- Types for a table of SIN's and COS's called Exp_Table. Makes the
   -- package pure by moving declaration of the table to the client program.
   -- The Table is passed to the FFT as an in/out parameter.

   Half_Max_Data_Length : constant Array_Index := Data_Index_Last / 2;

   subtype Exp_Mode_Index is Data_Index range 0..Half_Max_Data_Length+29;
   --  The +29 sometimes helps. (the usual -1 is ok).

   type Sinusoid_Storage is Array(Exp_Mode_Index) of Real;

   type Exp_Storage is
   record
     Re : Sinusoid_Storage;
     Im : Sinusoid_Storage;
     Current_Size_Of_Exp_Table : Data_Index := 0;
   end record;

   --  Current_Size_Of_Exp_Table
   --  tells the FFT routine whether to reconstruct Exp_Table.
   --  Current_Size_Of_Exp_Table is initialized to 0 so that the table
   --  is correctly initialized on first call to Make_Exp_Table.

end Fourier8;

--*****************************************************************
-- Notes on the radix 8 fast fourier transform (FFT)
--
-- The FFT is an algorithm for computing the discrete Fourier transform
-- with K*N*Log(N) arithmetic operation (where K is about 4 or 5).  A
-- direct calulation of the DFT uses O(N**2) operations.  (N is the
-- number of Data points; always a power of 2 in what follows.)  Here's
-- the decimation in frequency way.  We break the N point DFT (DFT(N))
-- into two DFT's of N/2 points.  Recursive application of this gives
-- the Cooley_Tukey FFT.  One of the 2 DFTs will give the even points
-- of the final result; the other gives the odd points of the final
-- result. Here is the derivation of first of the 2 DFTs.  We start by
-- defining F as the N point discrete Fourier transform of G: DFT(N){G}.
--
--                          N-1
--    r in [0,N-1]:  F(r) = SUM { W(N)**(rn) * G(n) }
--                          n=0
--
-- where W(N) = exp(-i*2*Pi/N), for the direct DFT, and exp(i*2*Pi/N) for
-- the inverse DFT.  For the even-indexed data points this is:
--
--                             N/2-1
--    r in [0,N/2-1]:  F(2r) =  SUM  { W(N)**(2rn) * G(n) }
--                              n=0
--
--                             N/2-1
--                      +       SUM  { W(N)**(2r(n+N/2)) * G(n+N/2) }
--                              n=0
--
-- Now use the fact that W(N)**(2rn) = W(N/2)**(rn) and that W(N)**(rN)
-- equals 1 to get:
--
--  r in [0,N/2-1]:    F(2r)   = DFT(N/2) { G(n) + G(N+N/2) }.
--
-- Now get the odd-indexed data points by the same process.  The result:
--
--  r in [0,N/2-1]:
--
--     F(2r)   = DFT(N/2) { (G(n) + G(N+N/2)) }
--     F(2r+1) = DFT(N/2) { (G(n) - G(N+N/2)) * W(N)**n }.
--
-- So now if you wrote a recursive procedure that performs the above
-- process, making N a variable (a parameter of the recursive procedure)
-- then you are done.  (i.e. next step is to break the two N/2 point
-- DFTs into 4 N/4 point DFTs, etc.)  In practice we just do it directly.
-- Each stage you form two new data sequences of length N/2 as defined
-- above: G(n) + G(N+N/2) and [G(n) + G(N+N/2)] * W(N)**n .
-- These operation are called the Butterflies.  Then you
-- perform the N/2 point DFT on both by doing the same thing all over
-- again but with N -> N/2. For example the W(N)**n term goes to W(N/2)**n.
-- This process is repeated Log2(N) times.  The variable in the code is
-- called Stage: Stage is in 0..Log2(N)-1.  Actually we have just described
-- the radix 2 FFT.  Below is the radix 4 FFT, in which you directly break
-- DFT into 4 DFT's of length N/4.  This gives you a 20-30% reduction
-- in floating point ops.  (And you use half as many stages, which is
-- very important in the out-of-core FFT in this set, which writes
-- the entire data set to the hard disk each stage.)  Here is the radix
-- 4 FFT.  To derive it, you just follow the above procedure.
--
-- for r in [0,N/4-1]:
--
-- F(4r)   = DFT(N/4){ [G(n) +  G(n+N/4) + G(n+2N/4) +  G(n+3N/4)] * W(N)**0n}
-- F(4r+1) = DFT(N/4){ [G(n) -i*G(n+N/4) - G(n+2N/4) +i*G(n+3N/4)] * W(N)**1n}
-- F(4r+2) = DFT(N/4){ [G(n) -  G(n+N/4) + G(n+2N/4) -  G(n+3N/4)] * W(N)**2n}
-- F(4r+3) = DFT(N/4){ [G(n) +i*G(n+N/4) - G(n+2N/4) -i*G(n+3N/4)] * W(N)**3n}
--
-- The improvement in computation comes from the fact that one need not
-- perform any operations in multiplying by i, and from the fact that
-- several quantities, like G(n) + G(N+N/2), appear more than once in
-- the expessions.  The 4 expressions in brackets [] form a 4 point DFT.
-- Notice that the factors of G are given by exp (-i k m (2Pi/N)) where
-- N = 4 and k, m go from 0..3. i.e. exp (-i k m (2Pi/N)) =
-- m=      0    1    2    3
--
-- k=0     1    1    1    1
-- k=1     1   -i   -1    i
-- k=2     1   -1    1   -1
-- k=3     1    i   -1   -i
--
-- This is the kernel to which all higher DFT's are reduced.  The
-- W(N)**kn terms are called twiddle factors, I think.  If the number
-- of Radix 2 stages (log(N)) is not divisible by 2, then a final radix
-- 2 stage must be performed to complete the FFT.  Finally, be aware that
-- to get the inverse DFT, all of the i's and exp's above must be
-- replaced with their complex conjugates.
-- The factors of G for Radix 8 are given by exp (-i k m (2Pi/N)) where
-- N = 8 and k, m go from 0..7. i.e. exp (-i k m (2Pi/N)) =
-- m=      0       1       2       3       4       5       6       7
--
-- k=0     1       1       1       1       1       1       1       1
-- k=1     1    a(1-i)    -i    a(-1-i)   -1   -a(1-i)     i    -a(-1-i)
-- k=2     1      -i      -1       i       1      -i      -1       i
-- k=3     1    a(-1-i)    i    a(1-i)    -1   -a(-1-i)   -i    -a(1-i)
-- k=4     1      -1       1      -1       1      -1       1      -1
-- k=5     1    a(-1+i)   -i    a(1+i)    -1   -a(-1+i)    i    -a(1+i)
-- k=6     1       i      -1      -i       1       i      -1      -i
-- k=7     1    a(1+i)     i    a(-1+i)   -1   -a(1+i)    -i    -a(-1+i)
--
-- In the above, a = 1/Sqrt(2.0).  If you call the above matrix Exp_km,
-- the FFT becomes:
--
-- for each r in [0, N/8-1] we have:
--
--  for k in 0..7:
--   F(8r+k) = DFT(N/8){ [Exp_km * G(n+mN/8)] * W(N)**(k*n)}
--
-- where W(N) = exp(-i*2*Pi/N), for the direct DFT, and exp(i*2*Pi/N) for
-- the inverse DFT.
-- Each stage of the FFT you transform the data according to the formula in
-- in the { }.  Then the idea is to perform eight N/8 point FFT's.  But each
-- of these is again performed according the the above formula, giving the
-- next stage in the FFT.  The above turns an N point FFT into 8 N/8 point
-- FFT's.  The next stage turns those eight N/8 pt. FFT's into 64 N/64 pt
-- FFT's.  To apply the above formula to that latter step, N must be replaced
-- by N/8 in the above formula.
-- You don't get a big reduction in floating point going from Radix 4 to 8,
-- but you read the array fewer times from memory,
-- and you can perform more floating point per expression which often helps.
-- Accuracy is improved because fewer complex mults are done.
-- Below let Dm = G(n+mN/8).  (D stands for Data.)
--
-- D0new := ((D0+D4) + (D2+D6) +   ( (D1+D5) + (D3+D7))) * Exp0;
-- D2new := ((D0+D4) - (D2+D6) - i*( (D1+D5) - (D3+D7))) * Exp2;
-- D4new := ((D0+D4) + (D2+D6) -   ( (D1+D5) + (D3+D7))) * Exp4;
-- D6new := ((D0+D4) - (D2+D6) + i*( (D1+D5) - (D3+D7))) * Exp6;
--
-- D1new := ((D0-D4) - i(D2-D6) + a( 1-i)*(D1-D5) + a(-1-i)*(D3-D7))*Exp1;
-- D3new := ((D0-D4) + i(D2-D6) + a(-1-i)*(D1-D5) + a( 1-i)*(D3-D7))*Exp3;
-- D5new := ((D0-D4) - i(D2-D6) + a(-1+i)*(D1-D5) + a( 1+i)*(D3-D7))*Exp5;
-- D7new := ((D0-D4) + i(D2-D6) + a( 1+i)*(D1-D5) + a(-1+i)*(D3-D7))*Exp7;
--
-- The last set can be written:
--
--D1new = ((D0-D4) - i(D2-D6) - a(i((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7))))*Exp1;
--D3new = ((D0-D4) + i(D2-D6) - a(i((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7))))*Exp3;
--D5new = ((D0-D4) - i(D2-D6) + a(i((D1-D5)+(D3-D7)) - ((D1-D5)-(D3-D7))))*Exp5;
--D7new = ((D0-D4) + i(D2-D6) + a(i((D1-D5)+(D3-D7)) + ((D1-D5)-(D3-D7))))*Exp7;
--
--*****************************************************************
