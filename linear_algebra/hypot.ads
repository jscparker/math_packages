
---------------------------------------------------------------------------
-- package Hypot
-- Copyright (C) 2018 Jonathan S. Parker.
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

--  Package Hypot provides routines for calculating
--
--     Hypot                      = Sqrt (Max**2 + Min**2)
--     Min_Arg_Over_Hypot         = Min / Sqrt (Max**2 + Min**2)
--     Max_Arg_Over_Hypot_Minus_1 = Max / Sqrt (Max**2 + Min**2) - One
--
--  The idea, mostly for fun, was to improve numerical accuracy on
--  Intel SSE hardware. The routines are slightly more accurate than
--  the simplest methods, and quite a bit slower.
--  Typically, they add an additional 1 or 2 calls to "/".

generic

  type Real is digits <>;

package Hypot is

  --  If Abs b > Abs a then
  --
  --     Max = Abs b,  Min = Abs a 
  --
  --  Hypot                      = Sqrt (Max**2 + Min**2)
  --  Min_Arg_Over_Hypot         = Min / Sqrt (Max**2 + Min**2)
  --  Max_Arg_Over_Hypot_Minus_1 = Max / Sqrt (Max**2 + Min**2) - One

  procedure Get_Hypotenuse 
    (a, b                       : in  Real;
     Hypot                      : out Real; 
     Min_Arg_Over_Hypot         : out Real;
     Max_Arg_Over_Hypot_Minus_1 : out Real);

  function Hypotenuse
    (a, b : in Real)
    return Real;

end Hypot;
