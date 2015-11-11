function [x,accuracy] = isotonicIPM(a,v,solver)
% [x,accuracy] = isotonicIPM(a,v) takes as input the adjacency matrix of a
% DAG a, and a vector v of real numbers corresponding to the vertices.
%
%    Computes the Isotonic Regression of v for the DAG given by a
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
    
default('solver',0);

n= length(a);
[ai,aj]  = find((a));
m = length(ai);

B = e2m([ai,aj],-1,n);
B = B';

Theta = m;
eps0 = 0.1;
beta2 = 10;



%initialization
ord = graphtopoorder(sparse(a));
x(ord) = 1:n;
x = x(:);


t = B*x;
if (max(t)>0)
    error('Not a DAG!')
end


eTa1 = -(x-v)'*B'*(-1./t)/norm(x-v);

if (sign((-(x)'*B'*(-1./t)/norm(x-v))) == 1)
    while eTa1<0
       x = 2*x; 
       eTa1 = (-(x-v)'*B'*(-1./t)/norm(x-v));
    end
else
    while eTa1<0
       x = x/2; 
       eTa1 = (-(x-v)'*B'*(-1./t)/norm(x-v));
    end
end

if(eTa1 <0)
    error('bad initialization!')
end
eTa2 = Theta/(eps0);
mu0 = eTa1;

[H,xNewton] = newtonStep( B,x,v,mu0,solver);
centMeasure = xNewton'*H*xNewton;


while mu0 < eTa2
    
    while(centMeasure>10^-2)
        
        F = @(y)(mu0*norm(y-v)^2 - sum(log(-B*y)));
        gradF = 2*mu0*(x-v) + B'*(-1./t);
        [x] = backtrackLineSearch(F,gradF,B,xNewton,x );
        [H,xNewton] = newtonStep( B,x,v,mu0, solver );
        centMeasure = xNewton'*H*xNewton;
        
    end
    
    dGap = dualGap(B,x,v,mu0);
    if ((sum((x-v).^2) - dGap)/sum((x-v).^2)<eps0*10^-1)
        break;
    end
    
    mu0 = mu0*beta2;
    [H,xNewton] = newtonStep( B,x,v,mu0,solver );
    centMeasure = xNewton'*H*xNewton;
end

accuracy = (sum((x-v).^2) - dGap)/sum((x-v).^2);


end


function dGap = dualGap(B,x,v,mu0)

y = -mu0*B*x; y = 1./y;
dGap = sum((B'*y).^2)/4 + y'*B*v - y'*(B*B')*y/2;

end