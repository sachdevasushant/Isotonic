function [x] = cmgSolver(la,b)
% function [x] = cmgSolver(la,b)
%
% this calls cmg to solve the system
% is designed to be instantiated as f = @(la,b)(cmgSolver(la,b))
% for use in general testing routines

pfun = cmg_sdd(la);
%[x,flag] = pcg(la, b, 1e-8, 100, pfun);
[x] = pcg(la, b, 1e-8, 100, pfun);

