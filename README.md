# Isotonic
Interior point method for computing the Isotonic regression (under Euclidean norm) on arbitrary directed acyclic graphs. The code accompanies the paper 'Fast, Provable Algorithms for Isotonic Regression in all l_p-norms' at NIPS 2015

This repository presently contains code for computing isotonic regression in directed acyclic graphs (DAG). The code is written in Matlab, and is meant to be called from Matlab. 

The main routine in Matlab is isotonicIPM. By default, the code uses cmg solver, and assumes that the user has installed cmg solver from http://www.cs.cmu.edu/~jkoutis/cmg.html. isotonicIPM takes as input the adjacency matrix  ‘a’ of a DAG, followed by a vector ‘v’ giving the initial values on the vertices. v should be a vector of length equal to the number of vertices.  
On input ‘a’ and ‘v’, the code computes a vector ‘x’ of same length as v, such that it minimizes the quantity:
sum_i (x(i) - v(i))^2
subject to x(i) <= x(j) if (i,j) is an edge in ‘a’.

Here is an example with 200*200 grid graph.

>> a = grid2(200,200);
>> a = triu(a);
>> v = randn(length(a),1);
>> [x, accuracy] = isotonicIPM(a,v);

The output contains a vector ‘x’ and an accuracy measure ‘accuracy’. x is a vector containing the values after l2-isotonic regression, while accuracy is equal to the duality gap of the computed value, which is a measure of closeness to the optimum value.

Note that the cmg solver has a few bugs. It occasionally fails to solve linear equations and runs into error. We suggest rerunning the program a few times and/or use a different solver. An alternate solver is provided with the code which is based on incomplete Cholesky factorization. You can use it by passing a third argument to isotonicIPM as follows:

>> [x, accuracy] = isotonicIPM(a,v,1);