
-------------------------------------------------------------------------------
-- package Disorderly, Linear Random Number Generator
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

-- PACKAGE Disorderly
--
-- 1. Procedure  Disorderly.Random.Get_Random
--
--   is a 61-bit non-linear pseudo random number generator optimized for
--   high statistical quality, and designed to be task-safe.  
--
--   procedure  Disorderly.Random.Get_Random  meets several design goals:
--
--    1. Uniform output in the range 0 .. 2**61-1.  
--       (ie, 61 random bits per call.)  
--       (2**61-1 is a Mersenne Prime; explanation below.) 
--    2. Period > 2^246. (Actually, period is about 1.813381 * 2^246.)
--    3. Period is the product of 4 large (8-byte) prime numbers. 
--    4. Generator is non-linear.
--    5. Generator is full-period, and 2 of the 3 component generators are
--       full-period. (The 3rd, the non-linear one, is optionally full-period.)
--    6. Generator is pure (the package is stateless) for convenient use
--       in multi-tasking simulations.
--    7. CPU time per call is constant (again for multi-tasking simulations).
--    8. Size of state per generator is 4 x 64 bits.
--
--   Items 1-5 are general characteristics of good generators.
--   Items 6-8 are attributes that are desirable in a language with built-in 
--   concurrency.
--
-- 2. procedure  Disorderly.Basic_Rand.Get_Random
--
--   is a stripped down linear version of Disorderly.Random generator - for
--   cases in which speed of Random number generation is the most important
--   consideration. It's more than twice as fast as Disorderly.Random. I give
--   it the unappealing name to encourage use of the Disorderly.Random generator
--   instead.
-- 
-- 3. package  Disorderly.Random.Deviates
--
--   is a package of random deviates: floating point random variables
--   satisfying the following distributions: 
--
--      Uniform, Normal (Gaussian), Exponential, Lorentzian (Cauchy),
--      Poissonian, Binomial, Negative Binomial, Weibull, Rayleigh, 
--      Student_t, Beta, Gamma, Chi_Squared, Log_Normal, Multivariate_Normal.
--
-- Rationale
--
-- I'd summarize the history of pseudo random number generation in 
-- five disgraceful chapters:
--
--  1. Failure of numerical analysts to provide algorithms that work.
--  2. Failure of programmers to choose correct algorithms to implement.
--  3. Failure of programmers to correctly implement the algorithms that do work.
--  4. Failure of implementors to communicate limitations in the algorithms
--     to users.
--  5. User preference for bad algorithms that don't work, on the grounds
--     that they are faster than good algorithms that do work.
--
-- By the early 1990's, most of the technical problems had been solved. It's now
-- straightforward to make high integrity random number generators with good
-- efficiency, but the less technical problems (2 - 5) remain unsolved.

package Disorderly is

   pragma Pure (Disorderly);

end Disorderly;

