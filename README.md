# Isotonic

	Interior point method for computing Isotonic Regression
    Part of the code for computing Isotonic regresstion
    Original code downloaded from https://github.com/sachdevasushant/Isotonic
	Copyright (C) 2015 Rasmus Kyng, Anup Rao, Sushant Sachdeva

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
Interior point method for computing the Isotonic regression (under Euclidean norm) on arbitrary directed acyclic graphs. The code accompanies the paper 'Fast, Provable Algorithms for Isotonic Regression in all l_p-norms' at NIPS 2015

This repository presently contains code for computing isotonic regression in directed acyclic graphs (DAG). The code is written in Matlab, and is meant to be called from Matlab. 

The main routine in Matlab is isotonicIPM. By default, the code uses cmg solver, and assumes that the user has installed cmg solver from http://www.cs.cmu.edu/~jkoutis/cmg.html. isotonicIPM takes as input the adjacency matrix  ‘a’ of a DAG, followed by a vector ‘v’ giving the initial values on the vertices. v should be a vector of length equal to the number of vertices.  
On input ‘a’ and ‘v’, the code computes a vector ‘x’ of same length as v, such that it minimizes the quantity:
sum_i (x(i) - v(i))^2
subject to x(i) <= x(j) if (i,j) is an edge in ‘a’.

Here is an example with 200*200 grid graph.
```
>> a = grid2(200,200);
>> a = triu(a);
>> v = randn(length(a),1);
>> [x, accuracy] = isotonicIPM(a,v);
```

The output contains a vector ‘x’ and an accuracy measure ‘accuracy’. x is a vector containing the values after l2-isotonic regression, while accuracy gives an upper bound on the error in the objective value at x, compared to the optimum. The upper bound is computed by constructing a dual certificate.

Note that the cmg solver has a few bugs. It occasionally fails to
solve linear equations and runs into error. We suggest rerunning the
program a few times and/or using a different solver. An alternate
solver based on incomplete Cholesky factorization is provided with the
code.
You can use it by passing a third argument to isotonicIPM as follows:
```
>> [x, accuracy] = isotonicIPM(a,v,1);
```
