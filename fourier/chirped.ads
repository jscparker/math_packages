
-- PACKAGE Chirped
--
-- The Chirped FFT performs a true Discrete Fourier Transform (DFT).
-- Unlike standard FFT's, the Chirped FFT performs a Discrete Fourier
-- Transform on data sets of any length. Data set lengths are not restricted 
-- to powers of 2 or any other special number. 
-- Like the FFT, the Chirped FFT's running time scales as O(N*Log(N)), 
-- where N is the Data length. But N can be a prime number like 700_001, 
-- and the true Discrete Fourier Transform will be performed on 0..N-1.
--
-- Just a reminder of why this matters: the simplest implementation of a
-- true Discrete Fourier Transform (DFT) scales as O(N**2). So on a big
-- data set, the ratio of the DFT to the Chirped FFT is O(N / Log(N)).
-- If N = 700_001, then that ratio is in the 10_000 range; the DFT would
-- take hours, the Chirped FFT seconds.
--
-- Ordinary FFT's are restricted to certain special values of N.
-- (N = 2**k is typical.) If the data is some other length, then 
-- you usually pad the array with 0's out the required length. You might 
-- also have to artificially damp the data to 0 at its end point with a
-- windowing function.  This process usually alters the original data. 
-- The chirped FFT requires no padding, which reduces the need to alter 
-- the data. That's why the chirped FFT is better suited to getting least
-- squares fits of Sines and Cosines to data sets of arbitrary length. 
-- If instead the standard FFT were used for this, you'ld be doing least
-- squares fits to data that in most cases is padded with artificial data
-- out to some length. The answers would be completely different.
--
-- The standard FFT is often essential in real-time signal processing
-- where speed is important, and where N can be set to an optimal length 
-- for run-time efficiency. The Chirped FFT would be preferred more
-- generally in data analysis problems.  In these problems N is rarely
-- chosen for run-time efficiency, and speed is rarely critical to
-- success.
--
-- The Chirp algorithm uses FFT's to do a convolution on the data
-- set.  That means that 2 or 3 FFT's are performed per call to 
-- FFT_Chirped and that they operate on arrays that are potentially 
-- twice as long as the original data set. That's why FFT_Chirped is
-- quite a bit slower than a standard FFT. (Expect a factor of 10
-- times slower.)
--
-- If the call to FFT_Chirped is in an inner loop, and the length of
-- data being transformed remains constant, and the choice of inverse
-- remains constant, then only 2 calls to FFT are performed per call to
-- FFT_Chirped, and certain tables of exponentials are reused between
-- calls and need not be calculated.
--
-- 1. Notes on use
--
-- The procedure operates only on data in the interval 0..N. 
-- N is input by the user. 
-- N is called Input_Data_Last.
-- The data to be transformed always starts at 0 for efficiency reasons.
--
generic

   type Real is digits <>;
   type Array_Index is range <>;
   type Data_Array is array (Array_Index) of Real;
   --  Must be at least twice as long as the data set to be FFT'ed.  In
   --  other words must include range 0..2**(Log_Of_Max_Data_Length+1)-1.
   --  Data will only be fft'ed on range 0..Some_Number, where
   --  Some_Number <= 2**(Log_Of_Max_Data_Length)-1.
     
   Log_Of_Max_Data_Length : Positive;
   --  The generic formal Data_Array has arbitrary limits.  This is
   --  useful, but the FFT algorithm works on a restricted range.
   --  It only operates on data in the range 0..2**N-1 for some N.
   --  The user can find out about this restriction at runtime or at
   --  compile time.  To make sure he finds out about it at compile time
   --  he is required to enter the Maximum value that he wishes N to have.
   --  That value is N_max = Log_Of_Max_Data_Length.

package Chirped is

   Work_Array_Last : constant Array_Index
                      := Array_Index (2**(Log_Of_Max_Data_Length+1)-1);
   
   --  The array Data_Array has to be twice as large as the data to be FFT'ed.
   --  These declaration perform the check.
                      
   Data_Index_Last : constant Array_Index := Work_Array_Last / 2;

   subtype Data_Index is Array_Index range 0..Data_Index_Last;
   
   --  Only data in range 0..Some_Number, 
   --  Some_Number <= Data_Index'Last will be FFT'ed.

   procedure FFT_Chirped 
     (Data_Re, Data_Im        : in out Data_Array;
      Input_Data_Last         : in     Data_Index;
      Inverse_FFT_Desired     : in     Boolean := False;
      Normalized_Data_Desired : in     Boolean := True);
   --  The procedure will perform an FFT on data in the interval
   --           0..Input_Data_Last.
   --  Data outside this range will be ignored, and may be destroyed.
   --  The transformed data is returned in the array Data, in the range
   --           0..Input_Data_Last.
   --  Data beyond this point will be set to Complex_Zero.
   --  See additional notes on use in the text above.

end Chirped;

