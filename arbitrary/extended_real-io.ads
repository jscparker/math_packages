
generic
package Extended_Real.IO is

   --  E_real Text translation.

   function e_Real_Image 
     (X   : in e_Real; 
      Aft : in Positive := Positive'Last) 
      return String;
   --  e_Real to Text translation.    
   --  Minimum length of Aft is 12.  If less, it gets set to 12.
   --  The number of digits in the result is the minimum of the 3 quantities:
   --  Aft+1, Max_No_Of_Digits, and Max_Practical_String_Length.  Aft is
   --  input by the user: it's the number of decimal digits beyond the
   --  decimal point.  Max_No_Of_Digits is slightly greater than the
   --  Desired_Decimal_Digit_Precision, the generic parameter. (see below.)


   --  Text to E_real translation.

   --  The following Sets define white space. In translation from ASCII to
   --  e_Real only white space may separate the numbers in ASCII format. 
   --  The present setting is standard white space on UNIX workstations and
   --  Crays.  The parenthesis are used by complex numbers.

   type Set is Array(Character) of Boolean;
   pragma Pack (Set);

   Is_White_Space : constant Set := 
      (' ' | ',' | '(' | ')' => True, others => False);

   Is_Exp_Symbol : constant Set  :=
      ('E' | 'e' | 'D' | 'd' => True, others => False);

   
   procedure e_Real_Val 
     (X    : in     String; 
      Y    :    out e_Real;
      Last :    out Natural);
   --  Accepts many non Ada standard formats of the sort found in ascii output
   --  on common machines.  Also accepts formats of the sort one is likely
   --  to type in at the key board.  Accepts the following formats:
   --  
   --  INTEGER     : 1234567 or +1234567 or -1234567
   --  DECIMAL     : 12.34567 or -.1234567 or +.1234567 or 1234567.
   --  EXPONENTIAL : 1234.567E+002 or .1234567E002 or 123467.E-03
   --  NON_DECIMAL_EXPONENTIAL : -1234567E-003
   --
   --  Notice that
   --  0) Both the 'E' and the '.' are optional.
   --  1) If an 'E' exists, then the '.' may be anywhere to the left of the 'E'.
   --  2) The Leading sign is optional if its '+'.
   --  3) The sign of the exponent is optional if its '+'.
   --  4) The set Is_Exp_Symbol determines which Exp symbols are acceptible.
   --
   --  "Last" is the index of the last character of the Number.  The purpose
   --  is to allow one to read numbers sequentially from a string by inputting
   --  string X(Last+1..X'Length) once the value of Last is determined.
   --
   --  We leave it up to the extended arithm. package to raise constraint_error       
   --  when the Exponent is out of range.
   
   E_Format_Error : Exception;

private
   
   --  Parameters that determine the maximum size of the output string:

   No_Of_Decimal_Digits_Per_Chunk : constant := 6;
   -- 8 is good if Radix is big enough (see next Assert).

   pragma Assert (10.0**No_Of_Decimal_Digits_Per_Chunk < Radix_Minus_1);
   --  The Decimal value of e_Real is actually calculated in 
   --  Radix 10.0**No_Of_Decimal_Digits_Per_Chunk. i.e. in chunks of 5-8 digits.

   Bits_Per_e_Real : constant := No_Of_Bits_In_Radix * Mantissa'Length;

   Decimal_Digits_Per_e_Real : constant := 1 + (Bits_Per_e_Real*1000 - 1) / 3322;
   -- Ceiling (Bits_Per_e_Real / 3.322).

   Chunks_Per_e_Real : constant :=
          1 + (Decimal_Digits_Per_e_Real - 1) / No_Of_Decimal_Digits_Per_Chunk;
   --  Ceiling (Decimal_Digits_Per_e_Real / No_Of_Decimal_Digits_Per_Chunk)

   Max_No_Of_Digits : constant Positive :=
                   No_Of_Decimal_Digits_Per_Chunk * Chunks_Per_e_Real - 4; 
   --
   --  Notice: it keeps some of the digits beyond the number supposedly held
   --  in the input e_Real (Item). By keeping these extra digits at the end,
   --  you reduce the Max error in translating from binary to ascii back to 
   --  binary again.
   --
   --  Subtract 4 because they always seem to be noise anyway.

   pragma Assert (Max_Exponent <= 2 ** (e_Integer'Size - 5));
   --  Just making sure.  Doesn't work for excessively large Max_Exponent.

   pragma Assert (e_Integer'Size <= Integer'Size);
   --  Essentially requires 32 bit Integer or greater.

   Max_Practical_String_Length : constant Positive := Max_No_Of_Digits + 19;
   --  Maximum size of the output strings.
   --  Set to whatever you think is appropriate. 
   --  Strings this size are actually created, so memory must be sufficient.
        
end Extended_Real.IO;
