function [x] = iccSolver(la,b,opts)
% function [x] = iccSolver(la,b,opts)
% function [f] = iccSolver(la,[],opts)
%
% this calls pcg with incomplete cholesky preconditioner,
% opts is passed to ichol, defaults are:
%   'nofill' and michol 'off'
%   tol 1e-6
%   maxit 100
%   if opts.L is supplied, it skips the call to ichol
%
% puts matrix in rcm order
%
%
%    Incomplete Cholesky Solver
%    Part of the code for computing Isotonic regresstion
%    Original code downloaded from https://github.com/sachdevasushant/Isotonic    
%    Copyright (C) 2015 Daniel Spielman, Yale University

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


default('b',[]);
default('opts','type','nofill');
default('opts','tol',1e-6);
default('opts','maxit',100);
default('opts','L',[]);

%opts.michol = 'on';

p = symrcm(la);
laperm = la(p,p);

icholOpts.type = opts.type;

if (isempty(opts.L))
  L2 = ichol(laperm,icholOpts);
else
  L2 = opts.L;
end

if isempty(b)
    f = @(b)(internal(laperm,p,L2,b,opts));
    x = f;
else
    x = internal(laperm,p,L2,b,opts);
end

end % main function

function x = internal(laperm,p,L2,b,opts)

  bperm = b(p);
  [xperm,flag,relres,iter] = pcg(laperm,bperm,opts.tol,opts.maxit,L2,L2');
  x(p) = xperm;
  x = x(:);

end
