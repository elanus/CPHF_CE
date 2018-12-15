# SCPHF_CE
Builds (Sherwood) Covering Perfect Hash Families by Conditional Expectation for general t

Parameters:
k -- the number of columns in the SCPHF
q -- the symbol set for the resulting covering array; the SCPHF is built on q^(t-1) symbols
t -- the strength of the SCPHF
N -- the final number of rows of the SCPHF

Parameter restrictions are based on the following:
1. Arithmetic tables are computed on the fly for primes but precomputed for prime powers (due to choosing polynomial). Prime power arithmetic tables have been provided for values up to and including 25.
2. Noncovering tuple list is precomputed to use noncovering check as a simple lookup. These are stored in a reduced canonical form where: a) only tuples of distinct symbols are stored, b) the first symbol is always 0, and c) the symbols appear in increasing order. As any noncovering tuple can either be immediately identified as noncovering due to a repeated symbol or can be put in this form by subtraction and sorting, this is a compact representation of a large list. The noncovering lists are produced in folder NCReduced for t=3, q <= 7 and t=4, q <= 9. Additional lists need to be computed for larger t, q.
3. Count of noncovering tuples is used for computation of the density. This allows for choosing a threshold value greater than the density to guarantee meeting bounds. For some values, this is an exact count; for others, an upper bound is used. Counts have been computed for t=3, q <= 7 and t=4, q <= 9 and coded into the program. Additional counts need to be computed for larger t, q.

The program lists the following options:
Debug -- this prints information about the inner calculations and running of the program. Set to 0 for faster execution.
Verbose -- prints more information about what the parameters are 
Add on -- allows you to enter a partial SCPHF to build from. The partial SCPHF is read from input.txt
Serial first row -- starts the first row with all the symbols from the set in order, repeated as needed. (E.g. 012345678012...) When k is large relative to q^(t-1)
k -- the columns
q -- symbols from the set q^(t-1)
t -- the strength
threshhold -- when choosing the best symbol, the percent of the row density to meet in order to early exit. (E.g. if the expectation is x and the threshhold is y, if some symbol has expectation xy, select that symbol and stop checking)
