function [H, xNewton] = newtonStep(B, x, v, mu0, solver)
%[H, xNewton] = newtonStep( B,x,v,mu0 )
%   finds the Newton step for mu0*c'*s + log-barrier stuff.

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