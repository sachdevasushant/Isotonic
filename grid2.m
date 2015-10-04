function [a,bdry,xy] = grid2(k1,k2,ani,ani2)
% function [a,bdry,xy] = grid2(k1,k2,ani,ani2)
%
% a k1 x k2 grid graph
% if k2 is not specified, k2 = k1;
%
% ani is degree of anisotropy
% ani is 1 by default

%    Constructs the adjacency matrix of a 2d grid
%    Part of the code for computing Isotonic regresstion
%    Original code downloaded from https://github.com/sachdevasushant/Isotonic    
%    Copyright (C) 2015 Rasmus Kyng, Anup Rao, Sushant Sachdeva

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/ ...
%        licenses/>.

if (nargin < 4)
  ani2 = 1;
end

if (nargin < 3)
  ani = 1;
end

if (nargin < 2)
  k2 = k1;
end


n = k1*k2;

x = kron(ones(1,k2),[ones(1,k1-1), 0])';
a = ani2*spdiags(x, -1, n,n);

a = a + ani*spdiags(ones(1,n)', -k1, n, n);

a = a + a';

tm = zeros(k1,k2);
tm(1,:) = 1;
tm(end,:) = 1;
tm(:,1) = 1;
tm(:,end) = 1;
bdry = (tm(:) ==1);


x = kron(ones(1,k2),[1:k1]);
y = kron([1:k2],ones(1,k1));
xy = [x;y]';
