function [H, xNewton] = newtonStep(B, x, v, mu0, solver)
%[H, xNewton] = newtonStep( B,x,v,mu0 )
%   finds the Newton step for mu0*c'*s + log-barrier stuff.
%
%    Newton's step iteration for interior point methods
%    Part of the code for computing Isotonic regression
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
    
m = size(B, 1);
n = size(B, 2);

d1 = 2 * mu0 * speye(n);
t = B * x;
d2 = sparse(1:m, 1:m, t .^ -2, m, m);
H1 = B'*d2*B;

H = H1 + d1;

grad = 2 * mu0 * (x - v) + B'*(-1./t);

%xNewton = -pinv(full(H))*grad;
minEntryH = max(abs(t)) ^ 2;
if (solver==0)
    xNewton = - cmgSolver(minEntryH * H, grad) * minEntryH;
else
    xNewton = -iccSolver(minEntryH *H,grad)* minEntryH;
end

end