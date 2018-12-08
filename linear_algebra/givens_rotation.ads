

-- package Givens_Rotation
--
-- cos = P/r,   sin = L/r,   r = sqrt(P*P + L*L)
--
-- (P is for Pivot, L is for Low.
--
-- clockwise rotation: notice all rotations are centered on the diagonal
--
--   1  0  0  0  0       0       0
--   0  c  s  0  0       P       r
--   0 -s  c  0  0   x   L   =   0
--   0  0  0  1  0       0       0
--   0  0  0  0  1       0       0
--
-- hypot = r = Sqrt(P*P + L*L)
--

generic

   type Real is digits <>;

package Givens_Rotation is

  procedure Get_Rotation_That_Zeros_Out_Low
    (Pivot, Low      : in     Real;
     sn, cs          :    out Real;
     cs_minus_1      :    out Real;
     sn_minus_1      :    out Real;
     hypot           :    out Real;
     P_bigger_than_L :    out Boolean;
     Skip_Rotation   :    out Boolean);

end Givens_Rotation;

