
-----------------------------------------------------------------------
-- package body Extended_Real.IO, translations between extended precision and text
-- Copyright (C) 2008-2018 Jonathan S. Parker
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

-- Test for E_Real_Machine_Emax etc

package body Extended_Real.IO is

   --  The following are some global constant sets for ASCII to e_Real
   --  Translation.

   Is_Numeral   : constant Set 
                    := Set'('0'..'9' => True, others => False);

   Is_Sign : constant Set 
                    := Set'('-' | '+'  => True, others => False);
                 
   Is_Decimal_Pt : constant Set 
                    := Set'('.' => True, others => False);
                 

   -- A Chunk is a digit in Radix 10**6, which is the radix of decimal
   -- number we're translating the binary number into..6 decimal digits
   -- at a time.  The following are used by both the ASCII to e_Real
   -- and e_Real to ASCII routines.

   type Chunk_Integer is range -2**31+1 .. 2**31-1;
   Chunk : Chunk_Integer;
   Chunk_Width : constant := No_Of_Decimal_Digits_Per_Chunk;  -- up in spec.  
   Ten_To_The_Chunk_Width : constant E_Digit
                                 := Make_E_Digit (10.0**Chunk_Width); 
                     

   --  The following are for translating the exponent to a string.
   --  Exp_String is used by both "e_Real_Image" and "Exp_Image"

   Exp_String_Width : constant Positive := E_Integer'Width;
   --  Should include space for the sign.
     
   subtype Exp_String is String (1..Exp_String_Width);

   Exp_Full_Text : Exp_String := (others => ' ');

   Log_Base_10_Of_2 : constant := 0.3010299956639812;
   
   Zero : constant e_Real := +0.0;
   Ten  : constant e_Real := +10.0;    
   
   
   ---------------
   -- Exp_Image --
   ---------------
   
   --  Also returns the sign of Exp_Val: either '+' or '-'.
   --  Uses fact that first char of E_Integer'Image is ' ' or '-'.
   --  Left justifies the digits of the Exp in a field of blanks.
   --  Another option is to right justify in a field of 0's.
   
   function Exp_Image (Exp_Val : E_Integer) return Exp_String is
     Result : Exp_String := (others => ' ');
     E_String : constant String := E_Integer'Image (Exp_Val);
   begin
      
     Result (1..E_String'Length) :=  E_String;
     if Result(1) = ' ' then
         Result(1) := '+';
     end if;
     return Result;
      
   end Exp_Image;
   
   ------------------------------
   -- Count_of_Trailing_Blanks --
   ------------------------------
   
   function Count_of_Trailing_Blanks (Exp : Exp_String) return Integer is
     cnt : Integer := 0;
   begin
     for i in reverse 1..Exp_String'Length loop
        if Exp(i) = ' ' then
           cnt := cnt + 1;
        else
           exit;
        end if;
     end loop;
     return cnt;
   end Count_of_Trailing_Blanks;
   
   ------------------
   -- e_Real_Image --
   ------------------
   
   -- Extended Real (e_Real) to Text translation

   function e_Real_Image
     (X   : in e_Real; 
      Aft : in Positive := Positive'Last)
      return String
   is

     --  Make strings that are large enough to return any desired
     --  result that is inside the bounds stated in the spec.  The 
     --  extra chunk is there to fill in for lost digits due to leading
     --  zeros.

     Max_Result_Length : constant Positive 
       := Max_Practical_String_Length + 2*Chunk_Width + Exp_String_Width + 3;

     subtype Result_String is String(1..Max_Result_Length);                        
     Result          : Result_String := (others => ' ');
     Mantissa_String : Result_String := (others => '0');
    
     Sign  : Character := ' '; -- init important.
     Mantissa_Length : Positive;
     Result_Length   : Positive;
     Exp_Stripped_Length : Integer;
    

     --  Types for translating from extended e_Real to Decimal.

     Exp_Base_2 : Real;
     Exp_Base_10_Shift : Real;
     Exp_Shift, Exp_Val : E_Integer;
     I_Exp_Shift : Integer;
     Leading_Zeros : Natural;
    
     Y : e_Real := Abs (X);
     Leading_Chunk, Trailing_Chunks : e_Real;
     
     Stage_Last : Positive; -- Initialized in Step 0. 
    

     --  Read global memory Mantissa_Length.
     --  Update global memory "Mantissa_String".  Right justify 
     --  Chunk_Integer'Image(Chunk).  Use fact that Mantissa_String is
     --  initialized all '0', so no need to add '0's to the left of the 
     --  right justified Chunk_Integer'Image(Chunk).   Recall  Mantissa_..
     --  allows an extra chunk at the end to adjust for removal of leading
     --  zeros.  

     procedure Add_To_Mantissa_String (Chunk : Chunk_Integer; 
                                       Stage : Positive) is
       Ascii_Chunk : constant String  := Chunk_Integer'Image (Chunk);
       Len         : constant Integer := Ascii_Chunk'Length;
       Start, Ending : Integer;
     begin
    
      Start  := (Stage - Positive'First + 1)*Chunk_Width - Len + 2;
      Ending := (Stage - Positive'First + 1)*Chunk_Width;
      --  Right justifies digits in a field of '0's.
      
      if Start > Mantissa_Length + 2*Chunk_Width then
         return;
      end if;
      if Ending > Mantissa_Length + 2*Chunk_Width  then
         Ending := Mantissa_Length + 2*Chunk_Width;
      end if;
      
      Mantissa_String (Start..Ending) := Ascii_Chunk(2..Ending-Start+2);
      --  Ascii_Chunk always starts with ' ' or '-'.  Here it's always positive.
      --  Its min length (Len) is always 2, even if chunk = '0'.
    
     end Add_To_Mantissa_String;
                     
  begin


    -- STEP 0-. Special cases. Zero and Infinity.

    if Are_Equal (X, Zero) then
      return "0.0";
    end if;
    
    if Are_Equal (X, Positive_Infinity) then
      return "+inf";
    end if;
    
    if Are_Equal (X, Negative_Infinity) then
      return "-inf";
    end if;

    -- STEP 0. Determine number of decimal digits in Mantissa.
    -- Its Min (Aft+1, Max_No_Of_Digits, and Max_Practical_String_Length).
    -- It may not be less than 12.

    if Aft = Positive'Last then  
       Mantissa_Length := Aft;  --  Positive'Last is the upper limit.
    else     
       Mantissa_Length := Aft + 1;  -- the usual case.
    end if;
    
    if Mantissa_Length < 13 then Mantissa_Length := 13; end if;
    
    if Mantissa_Length > Max_No_Of_Digits then
       Mantissa_Length := Max_No_Of_Digits;
    end if;
    
    if Mantissa_Length > Max_Practical_String_Length then
       Mantissa_Length := Max_Practical_String_Length;
    end if;
 
    Stage_Last := (Mantissa_Length - 1) / Chunk_Width + 3;
    --  This is ceiling (Mantissa_Length / Chunk_Width) + 1. (+ 2 really)
    --  This is the number of steps required to strip digits of number
    --  Mantissa_Length from the e_Real in chunks of Chunk_Width.  The
    --  extra 2 chunk2 are usually required to fill in for digits lost due
    --  to leading zeros in the first chunk. (Normalization of the decimal
    --  representation.)
        

    -- STEP 1. Multiply Item by a power of 10.0 to make it less than 1.0.
    -- Why a power of ten?  Cause then we can shift the final Radix 10 exp
    -- by exactly the value of the exponent.  The formula is 
    --       Shift  =  -Ceiling (Log_Base_10_Of_2 * Exp_Base_2)
    -- Exp_Base_2 is the normalized exp:  i.e. 2**(-Exp_Base_2) * Item
    -- is slightly less than 1.  But we want a power of 10 with this property,
    -- not a power of 2.  10**(-N) = 2**(-Exp_Base_2) implies
    -- N = Log_Base_10_Of_2 * Exp_Base_2, except that we must make N the
    -- appropriate nearest integer.  If we do that by rounding N up (taking 
    -- the ceiling) then 10**(-N) < 2**(-Exp_Base_2) which implies that
    -- 10**(-N) * X < 2**(-Exp_Base_2) * X, so 10**(-N) * X < 1.0,
    -- as required.

    Exp_Base_2 := Real (No_Of_Bits_In_Radix) * Real (Exponent (X));
    
    Exp_Base_10_Shift := -Real'Ceiling (Log_Base_10_Of_2 * Exp_Base_2);    
    I_Exp_Shift :=   Integer (Exp_Base_10_Shift);
    Exp_Shift   := E_Integer (Exp_Base_10_Shift);
  

    -- STEP 2. Multiplying X by 10**Shift will make it less than 1.
    -- Want to multiply Item by 10**(6 + Shift) to get no more than 6 decimal
    -- digits sticking out to the left of the decimal point (and often 0),
    -- because we're stripping 6 decimal digits at a time from the left of the
    -- decimal point with the truncate function.  Remember Y := Abs(X).
    -- Loop translates Y into Radix 10**6.
    -- Each Chunk is a digit in Radix 10**6.

    Y := Ten ** (No_Of_Decimal_Digits_Per_Chunk + I_Exp_Shift) * Y;
    --  Y := Machine(Y);
    --  Machine would round away all guard digits.
    
    for Stage in 1..Stage_Last loop 
    
      Leading_Chunk   := Truncation (Y);
      Trailing_Chunks := Y - Leading_Chunk;
      
      Chunk := Chunk_Integer (Make_Real (Leading_Chunk));
      Add_To_Mantissa_String (Chunk, Stage);
      
      Y  := Ten_To_The_Chunk_Width * Trailing_Chunks;
      --  Shift another 6 decimal digits to the left of the decimal point.
      
    end loop;


    -- STEP 3. Construct the string.  Get leading sign.  Strip away leading
    -- Zeros.  Exp is the amount we shifted by above, adjusted for stripped
    -- leading zeros.  (ie, move decimal point to right of leading zeros
    -- and 1st non-zero digit; decrement Exp by 1 for each leading zero etc.)

    if X < Zero then    
       Sign := '-';
    end if;
    --  Set the sign.  Sign is initialized ' '.  Recall Y = Abs(X).
    
    --  Count leading zeros: 
    Leading_Zeros := 0;
    for I in 1..20 loop -- Should be no more than 8 if digit 1..10^9.
       if Mantissa_String (I) /= '0' then
          exit;
       else
          Leading_Zeros := Leading_Zeros + 1;
       end if;
    end loop;
    --  Right now the virtual decimal point sits to the left of every digit
    --  in Mantissa_String, and the Exponent is given by -Shift.  We want
    --  to shift the decimal point to the right of each leading zero, and to
    --  the right of the first non-zero digit.  Then DECREASE the exponent
    --  by that shift.

    Exp_Val := -Exp_Shift - (1 + E_Integer(Leading_Zeros)); 
    --  the 1 accounts for the 1st non-zero digit to the left of decimal point.

    Exp_Full_Text       := Exp_Image (Exp_Val);
    Exp_Stripped_Length := Exp_String_Width - Count_of_Trailing_Blanks (Exp_Full_Text);

    Result_Length := Mantissa_Length + Exp_Stripped_Length + 3;

    Result(1..Result_Length) :=
       Sign & Mantissa_String (1+Leading_Zeros) & '.' 
            & Mantissa_String (2+Leading_Zeros..Mantissa_Length+Leading_Zeros)
            & 'E' & Exp_Full_Text (1..Exp_Stripped_Length);
         
    return Result(1..Result_Length);
    
   end e_Real_Image;
   
   ----------------------
   -- Integer_Value_Of --
   ----------------------
   
   --  Assumes that a decimal point is immediately to the right of the last digit.
   --  Translate each chunk (Digit_n) of 8 decimal digits into e_Real, and sum
   --  the polynomial in powers of 10**8.  (Do it this way because we can then
   --  use E_Digit*e_Real operations for efficiency.)  Horner's rule is
   --
   --  Digit_0 + 10**8*(Digit_1 + 10**8*(Digit_2 + ... + 10**8*(Digit_n)))))
   --
   --  Below, the 10**8 is called Ten_To_The_Chunk_Width.  It's a Global constant.
   --  IMPORTANT to return ZERO if string has 0 length.  e_Real_Val requires it.

   function Integer_Value_Of (F : String) return e_Real is
     Result, Digit_i : e_Real; -- initialized to zero. (Important).
     No_Of_Full_Chunk_Iterations : Natural;
     First_Partial_Chunk_Width   : Natural;
     Start, Ending : Integer;
   begin
      if F'Length = 0 then
         return Zero;
      end if;
      
      No_Of_Full_Chunk_Iterations := F'Length  /  Chunk_Width;
      First_Partial_Chunk_Width   := F'Length MOD Chunk_Width;

      --  Special case for highest order Digit_i, the one of width 
      --  First_Partial_Chunk_Width.  If this partial chunk is the only
      --  chunk (ie., No_Of_Full_Chunk_Iterations = 0) then there is no
      --  multiplication of the result by 10**8.  (See formula above.)
      
      if First_Partial_Chunk_Width > 0 then     
         Start  := F'First;
         Ending := F'First + First_Partial_Chunk_Width - 1;
         
         Digit_i := +Real (Chunk_Integer'Value (F(Start..Ending)));
         
         if No_Of_Full_Chunk_Iterations > 0 then
            Result := Ten_To_The_Chunk_Width * Digit_i;
         else
            Result := Digit_i;
         end if;   
      end if;
      
      --  Do the lower order Digits_i's.  The lowest order Digit_i has no
      --  10**8 multiplied by it.
      
      for i in 0..No_Of_Full_Chunk_Iterations-1 loop
      
         Start  := F'First + First_Partial_Chunk_Width + i * Chunk_Width;
         Ending := Start + Chunk_Width - 1;
         
         Digit_i := +Real (Chunk_Integer'Value (F(Start..Ending)));
         
         if i < No_Of_Full_Chunk_Iterations-1 then
            Result := Ten_To_The_Chunk_Width * (Result + Digit_i);
         else
            Result := Result + Digit_i;
         end if;   
         
      end loop;
      
      return Result;
      
   end Integer_Value_Of;

   -------------------------
   -- Fractional_Value_Of --
   -------------------------
   
   --  Assumes that a decimal point is immediately to the left of the 1st digit.
   --  IMPORTANT to return ZERO if string has 0 length.  e_Real_Val requires it.
   --
   function Fractional_Value_Of (F : String) return e_Real is
      Result : e_Real;
   begin
   
      if F'Length = 0 then
         return Zero;
      end if;
      
      Result := Integer_Value_Of(F);
      --  We have Val as though decimal point were to the right of all digits.
      --  Want it to the left of all digits: just multiply by Ten**(-F'Length)
      --  We do it this way so can use efficiency of E_Digit*e_Real in routine
      --  Integer_Value_Of.  But many possible optimizations are neglected.
      --  If performance is an issue, then maybe try e_Real / E_Digit method
      --  where E_Digit = 10**8.
      
      Result := Result * Ten**(-F'Length);
   
      return Result;
   
   end Fractional_Value_Of;
   
   ----------------
   -- e_Real_Val --
   ----------------

   --  Accepts the following formats:
   --  INTEGER : 1234567
   --  DECIMAL : 12.34567 or -.1234567 or 1234567.
   --  EXPONENTIAL : 1234.567E+002 or .1234567E2 or 123467.E-03
   --  NON_DECIMAL_EXPONENTIAL : -1234567E-003
   --
   --  Notice that 
   --  1)  The Decimal point may be anywhere to the left of the E.
   --  2)  The Leading sign is optional if it's '+'.
   --  3)  The sign of the exponent is optional if it's '+'.
   --     
   --  Start_Of_Num is the first non-white space.  If the first char of the
   --  string is not white-space, then this will be Start_Of_Num. End_Of_Num
   --  is the non-white char just before the start of more white-space.  BUT if
   --  the string comes to an end before any white-space re-appears, then the
   --  end of the string is taken as End_Of_Num.
   --
   procedure e_Real_Val 
     (X    : in  String;
      Y    : out e_Real;
      Last : out Natural)
   is
      No_Of_Decimal_Pts, No_Of_Exp_Symbols, No_Of_Exp_Signs : Natural := 0;
      Decimal_Pt_Exists, Exp_Symbol_Exists, Exp_Sign_Exists : Boolean := False;
      Leading_Sign_Exists : Boolean := False;
      Decimal_Pt_Pos, Exp_Symbol_Pos, Exp_Sign_Pos : Positive;
      Start_Of_Aft : Positive;
      Num_Is_Positive : Boolean;
      Char : Character;
      
      type Format is 
         (Int, Decimal_Pt_Only, Exponential, No_Decimal_Pt_Exponential);
      The_Format : Format;
      
      Start_Of_Num, End_Of_Num, Start_Of_Exp : Natural := 0;
      Fore_Width, Aft_Width, Exp_Width : Natural := 0;
      Exp_Str : Exp_String := (others => ' ');
      Exp_Val : Integer := 0;
      
      Fore, Aft, Result : e_Real;

   begin
   
      --  Handle null strings:
      
      if X'Length = 0 then
         raise E_Format_Error;
      end if;
   

      -- STEP 1. Strip away leading whitespace and sign. Start_Of_Num is the
      -- first non-white space character.  If it's '-' or '+', then increment 
      -- Start_Of_Num by 1. 

      for I in X'Range loop
         if not Is_White_Space (X(I)) then
            Start_Of_Num := I;
            exit;
         end if;
      end loop;
   
      if Start_Of_Num = 0 then -- The string is all white-space, so:
         raise E_Format_Error;
      end if;
   
      --  If there is a leading sign, make a note of it, and then bypass it.
      
      Num_Is_Positive := True; -- No sign means positive
      Char := X(Start_Of_Num);   
      if Is_Sign (Char) then
        Start_Of_Num := Start_Of_Num + 1;
        Leading_Sign_Exists := True;
        if Char = '+' then
           Num_Is_Positive := True;
        elsif Char = '-' then
           Num_Is_Positive := False;
        end if;
      end if;
      
      if Leading_Sign_Exists and then Start_Of_Num > X'Last then
         raise E_Format_Error;
      end if;
      --  when we incremented Start_Of_Num, we went beyond X'Last.  So only char
      --  is the Sign.
      

      -- STEP 2. Scan everything beyond the leading sign.  End_Of_Num is
      -- initialized to 0.  Here is where we update it. 

      for I in Start_Of_Num..X'Last loop
   
         Char := X(I);
         
         if Is_White_Space (Char) then -- we know Start.. points to non-whitespace
            End_Of_Num := I-1;
            exit;
         end if;
         
         if Is_Decimal_Pt (Char) then
            Decimal_Pt_Exists := True;
            Decimal_Pt_Pos    := I;
            No_Of_Decimal_Pts := No_Of_Decimal_Pts + 1;
         elsif Is_Exp_Symbol (Char) then
            Exp_Symbol_Exists := True;
            Exp_Symbol_Pos    := I;
            No_Of_Exp_Symbols := No_Of_Exp_Symbols + 1;
         elsif Is_Sign (Char) then
            Exp_Sign_Exists := True;
            Exp_Sign_Pos    := I;
            No_Of_Exp_Signs := No_Of_Exp_Signs + 1;
         elsif not Is_Numeral (Char) then
            raise E_Format_Error;
         end if;
   
      end loop;
      
      if End_Of_Num = 0 then --  Reached the end of string w/o detecting whitespace
         End_Of_Num := X'Last;
      end if;
   
      --  Do some error checking:
      
      if No_Of_Decimal_Pts > 1 then
         raise E_Format_Error;
      end if;
      
      if No_Of_Exp_Signs > 1 then
         raise E_Format_Error;
      end if;
      
      if No_Of_Exp_Symbols > 1 then
         raise E_Format_Error;
      end if;
   
      if Decimal_Pt_Exists and Exp_Symbol_Exists then
          if Decimal_Pt_Pos >= Exp_Symbol_Pos then
             raise E_Format_Error;
          end if;
      end if;
   
      --  if there's an 'E' and a '+', then the sign must directly follow the
      --  the 'E'. Neither can be at the beginning or end of the num.  if the
      --  sign exists, then 'E' must exist, but not vice-versa.  
      if Exp_Sign_Exists then
         if Exp_Sign_Pos = End_Of_Num or else Exp_Sign_Pos = Start_Of_Num then
            raise E_Format_Error;
         end if;
      end if;
     
      if Exp_Symbol_Exists then
         if Exp_Symbol_Pos = End_Of_Num or else Exp_Symbol_Pos = Start_Of_Num then
            raise E_Format_Error;
         end if;
      end if;
      
      if Exp_Sign_Exists and Exp_Symbol_Exists then
         if not (Exp_Sign_Pos = Exp_Symbol_Pos + 1) then
            raise E_Format_Error;
         end if;
      end if;
      
      if Exp_Sign_Exists and not Exp_Symbol_Exists then
         raise E_Format_Error; 
      end if;
      
      -- Diagnose the format of the number:
   
      if (not Decimal_Pt_Exists) and (not Exp_Symbol_Exists) then
          The_Format := Int;
      elsif Decimal_Pt_Exists and (not Exp_Symbol_Exists) then
          The_Format := Decimal_Pt_Only;
      elsif Decimal_Pt_Exists and Exp_Symbol_Exists then
          The_Format := Exponential;
      else
          The_Format := No_Decimal_Pt_Exponential;
      end if;
   
      -- STEP 3. Do the arithmetic.  The string goes like Fore . Aft E Exp,
      -- in the most general case.  The following is not really optimized
      -- for speed.

      case The_Format is
      when Int =>
 
         Fore_Width := End_Of_Num - Start_Of_Num + 1;
         if Fore_Width = 0 then
            raise E_Format_Error;
         end if;
         
         Result := Integer_Value_Of(X (Start_Of_Num..End_Of_Num));
         
      when Decimal_Pt_Only =>
 
         Fore_Width   := Decimal_Pt_Pos - Start_Of_Num;
         Aft_Width    := End_Of_Num - Decimal_Pt_Pos;
         Start_Of_Aft := Decimal_Pt_Pos + 1;
         
         if Fore_Width = 0 and Aft_Width = 0 then
            raise E_Format_Error;
         end if;
         --  Notice that Fore or Aft may be 0.  The funcs below return
         --  0.0 in that case.
         
         Fore := Integer_Value_Of (X (Start_Of_Num..Start_Of_Num + Fore_Width - 1));
         Aft  := Fractional_Value_Of (X (Start_Of_Aft..Start_Of_Aft + Aft_Width - 1));
         
         Result := Fore + Aft;
 
      when Exponential =>
 
         if Exp_Sign_Exists then
            Start_Of_Exp := Exp_Sign_Pos + 1;
         else
            Start_Of_Exp := Exp_Symbol_Pos + 1;
         end if;
 
         Start_Of_Aft := Decimal_Pt_Pos + 1;
         Fore_Width   := Decimal_Pt_Pos - Start_Of_Num;
         Aft_Width    := Exp_Symbol_Pos - Start_Of_Aft;
         Exp_Width    := End_Of_Num - Start_Of_Exp + 1;
         
         if Exp_Width = 0 then
            raise E_Format_Error;
         end if;
         if Exp_Width > Integer'Width then 
         -- E_Integer is derived Integer.  Both Integer and E_Integer are much
         -- wider than the allowed range of Exponents, so Contraint-Error will
         -- be raised in the arithmetic package, not here.
            raise E_Format_Error;
         end if;
         if Fore_Width = 0 and Aft_Width = 0 then
            raise E_Format_Error;
         end if;
         --  Notice that Fore or Aft may be 0.  The funcs below return
         --  0.0 in that case.
 
         Fore := Integer_Value_Of (X (Start_Of_Num..Start_Of_Num + Fore_Width - 1));
         Aft  := Fractional_Value_Of (X (Start_Of_Aft..Start_Of_Aft + Aft_Width - 1));
         Exp_Str(1..Exp_Width) := X (Start_Of_Exp..Start_Of_Exp + Exp_Width - 1);
         Exp_Val := Integer'Value(Exp_Str);
         if Exp_Sign_Exists and then X(Exp_Sign_Pos) = '-' then
            Exp_Val := -Exp_Val;
         end if;
                 
         Result := (Fore + Aft) * Ten**Exp_Val;
 
      when No_Decimal_Pt_Exponential =>
      
         --  unusual case.  Say there's no Aft, only Fore
 
         if Exp_Sign_Exists then
            Start_Of_Exp := Exp_Sign_Pos + 1;
         else
            Start_Of_Exp := Exp_Symbol_Pos + 1;
         end if;
 
         Fore_Width := Exp_Symbol_Pos - Start_Of_Num;
         Exp_Width  := End_Of_Num - Start_Of_Exp + 1;
         if Fore_Width = 0 then
            raise E_Format_Error;
         end if;
         if Exp_Width = 0 then
            raise E_Format_Error;
         end if;
         if Exp_Width > Integer'Width then
            raise E_Format_Error;
         end if;
 
         Fore := Integer_Value_Of (X (Start_Of_Num..Start_Of_Num + Fore_Width - 1));
         Exp_Str(1..Exp_Width) := X (Start_Of_Exp..Start_Of_Exp + Exp_Width - 1);
         Exp_Val := Integer'Value(Exp_Str);
         if Exp_Sign_Exists and then X(Exp_Sign_Pos) = '-' then
            Exp_Val := -Exp_Val;
         end if;
 
         Result := Fore * Ten**Exp_Val;
 
      end case;
 
      --  update the 2 out parameters:
      
      if Num_Is_Positive then
         Y := Result;
      else
         Y := -Result;     
      end if;
        
      Last := End_Of_Num;
        
   end e_Real_Val;
  
end Extended_Real.IO;

