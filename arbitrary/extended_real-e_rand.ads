
-- PACKAGE E_Rand
--
--  Minimal random number generator, mainly for basic test of 
--  extended arithmetic.
--
--  Outputs extended precision pseudo-random numbers on [0.0 .. 1.0)

generic
package Extended_Real.E_Rand is

   function Random return E_Real;  
   --  returns E_Reals in [0, 1)
  
   procedure Reset (Initiator : in Positive := 7777777);
   --  Resets both   function Random   and   function Next_Random_Int.


   type Random_Int is mod 2**64;

   function Next_Random_Int return Random_Int;
   --  61 bit generator, period 2**61-1.

end Extended_Real.E_Rand;
