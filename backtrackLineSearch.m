function [xt] = backtrackLineSearch(F, gradF, B, xNewton, x)
%[xt] = backtrackLineSearch(F,gradF,B,xNewton ) takes in a function handle F,gradF,xNewton, x  and does a back
%track line search as in Boyd's book
    
%    Backtracking line search for Netwon's iteration
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
    
A0 = 0.01;
B0 = 0.5;

t = 1;
iter = 0;

xt = x + t * xNewton;

while ((notFeasible(xt, B)) || (F(xt) > F(x) + A0 * t * gradF'*xNewton))
    t = B0 * t;
    iter = iter + 1;
    if (iter >= 100)
        error('Number of iterations in line search exceeded the limit')
    end
 
    xt = x + t * xNewton;
end

end

function [y] = notFeasible(xt, B)

voltDiff = (B * xt);
y = (max(voltDiff) >= 0);

end