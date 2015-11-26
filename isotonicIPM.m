function [x,accuracy,i] = isotonicIPM(a,v,solver)
% [x,accuracy] = isotonicIPM(a,v) takes as input the adjacency matrix of a
% DAG a, and a vector v or real numbers correspdonding on the vertices.

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

% Note: this initial eTa1 choice is not as recommended by Boyd,
% but has worked well for us in practice
eTa1 = (-(x-v)'*B'*(-1./t)/norm(x-v))

while eTa1<0
    x = 2*x; 
    eTa1 = (-(x-v)'*B'*(-1./t)/norm(x-v));
end

if(eTa1 <0)
    error('bad initialization!')
end
eTa2 = Theta/(eps0);
mu0 = eTa1;

[H,xNewton] = newtonStep( B,x,v,mu0,solver);
centMeasure = xNewton'*H*xNewton;

i=0;

while mu0 < eTa2
    while(centMeasure>10^-2)
        t = B*x;
        F = @(y)(mu0*norm(y-v)^2 - sum(log(-B*y)));
        gradF = 2*mu0*(x-v) + B'*(-1./t);
        [x] = backtrackLineSearch(F,gradF,B,xNewton,x );
        [H,xNewton] = newtonStep( B,x,v,mu0, solver );
        centMeasure = xNewton'*H*xNewton;
        i = i+1;
    end
    
    dGap = dualGap(B,x,v,mu0);
    if ((sum((x-v).^2) - dGap)/sum((x-v).^2)<10^-2)
        break;
    end
    
    mu0 = mu0*beta2;
    [H,xNewton] = newtonStep( B,x,v,mu0,solver );
    centMeasure = xNewton'*H*xNewton;
    i=i+1;
end

accuracy = (sum((x-v).^2) - dGap)/sum((x-v).^2);

display(i);
end


function dGap = dualGap(B,x,v,mu0)

y = -mu0*B*x; y = 1./y;
dGap = y'*B*v - y'*(B*B')*y/4;

end