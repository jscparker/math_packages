-----------------------------------------------------------------------
-- package body Binary_Rank, routines for Marsaglia's Rank Test.
-- Copyright (C) 1995-2018 Jonathan S. Parker
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

package body Binary_Rank is

 --------------
 -- Get_Rank --
 --------------

 procedure  Get_Rank
   (r       : in out Binary_Matrix;
    Max_Row : in Matrix_Index;
    Max_Col : in Matrix_Index;
    Rank    : out Integer)
 is
    tmp : Vector;
    j   : Matrix_Index := Max_Col;             -- essential init.
    i   : Matrix_Index := Matrix_Index'First;  -- essential init.

    -----------------------------
    -- Vector_has_a_1_at_Index --
    -----------------------------

    function Vector_has_a_1_at_Index
      (j : Matrix_Index;
       V : Vector)
       return Boolean
    is
       Result     : Boolean := False;
       Segment_id : Segments;
       Bit_id     : Matrix_Index;
    begin
       Segment_id := 1 + (j-1) / Bits_per_Segment;
       Bit_id     := 1 + (j-1) - (Segment_id-1) * Bits_per_Segment;
       if (V(Segment_id) AND Bit_is_1_at(Bit_id))  >  0 then
         Result := True;
       end if;
       return Result;
    end;

 begin

    Rank := 0;
 
    <<Reduce_Matrix_Rank_by_1>>
 
    -- Start at row i, the current row.
    -- Find row ii = i, i+1, ... that has a 1 in column j. 
    -- Then move this row to position i.  (Swap ii with i.)
    -- Reduce rank of matrix by 1 (by xor'ing row(i) with succeeding rows).
    -- Increment i, Decrement j.
    -- Repeat until no more i's or j's left.
 
    Row_Finder: for ii in i .. Max_Row loop
 
      if  Vector_has_a_1_at_Index (j, r(ii))  then
 
         Rank := Rank + 1;
  
         if i /= ii then
            tmp   := r(ii);
            r(ii) := r(i);
            r(i)  := tmp;
         end if;
 
         for k in i+1 .. Max_Row loop
            if  Vector_has_a_1_at_Index (j, r(k))  then
               for Seg_id in Segments loop
                  r(k)(Seg_id) := r(k)(Seg_id) XOR r(i)(Seg_id);
               end loop;
            end if;
         end loop;
 
         if  i = Max_Row  then return; end if;
         i := i + 1;
         if  j = Matrix_Index'First  then return; end if;
         j :=  j - 1;
         goto Reduce_Matrix_Rank_by_1;
 
      end if;
 
    end loop Row_Finder; -- ii loop

    if  j = Matrix_Index'First  then return; end if;
    j :=  j - 1;
    goto Reduce_Matrix_Rank_by_1;

 end Get_Rank;

end Binary_Rank;
