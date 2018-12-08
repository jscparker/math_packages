
generic

   No_of_Bits_per_Segment : Integer;
   Segments_per_Vector    : Integer;

package Binary_Rank is

  type Unsigned_64 is mod 2**64;

  Bits_per_Segment : constant Integer := No_of_Bits_per_Segment;
  Bits_per_Vector  : constant Integer := Bits_per_Segment*Segments_per_Vector;

  subtype Unsigned_Segment is Unsigned_64 range 0 .. 2**Bits_per_Segment-1;

  subtype Matrix_Index is Integer range 1 .. Bits_per_Vector;

  subtype Segments is Matrix_Index range 1 .. Segments_Per_Vector;

  type Vector is array (Segments) of Unsigned_Segment;

  type Binary_Matrix is array(Matrix_Index) of Vector;

  procedure  Get_Rank
    (r       : in out Binary_Matrix;
     Max_Row : in Matrix_Index;
     Max_Col : in Matrix_Index;
     Rank    : out Integer);

  subtype Bit_Index is Integer range Matrix_Index'First .. Matrix_Index'First+63;
  type Binary_Matrix_Block is array(Bit_Index) of Unsigned_64;

  Bit_is_1_at : constant Binary_Matrix_Block := 
    (2**00, 2**01, 2**02, 2**03, 2**04, 2**05, 2**06, 2**07, 2**08, 2**09,
     2**10, 2**11, 2**12, 2**13, 2**14, 2**15, 2**16, 2**17, 2**18, 2**19,
     2**20, 2**21, 2**22, 2**23, 2**24, 2**25, 2**26, 2**27, 2**28, 2**29,
     2**30, 2**31, 2**32, 2**33, 2**34, 2**35, 2**36, 2**37, 2**38, 2**39,
     2**40, 2**41, 2**42, 2**43, 2**44, 2**45, 2**46, 2**47, 2**48, 2**49,
     2**50, 2**51, 2**52, 2**53, 2**54, 2**55, 2**56, 2**57, 2**58, 2**59,
     2**60, 2**61, 2**62, 2**63);
 
end Binary_Rank;

