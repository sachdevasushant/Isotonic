function [xt] = backtrackLineSearch(F, gradF, B, xNewton, x)
%[xt] = backtrackLineSearch(F,gradF,B,xNewton ) takes in a function handle F,gradF,xNewton, x  and does a back
%track line search as in Boyd's book

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