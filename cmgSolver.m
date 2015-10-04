function [x] = cmgSolver(la,b)
% function [x] = cmgSolver(la,b)
%
% this calls cmg to solve the system
% is designed to be instantiated as f = @(la,b)(cmgSolver(la,b))
% for use in general testing routines

%    Wrapper code for cmg solver
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

pfun = cmg_sdd(la);
%[x,flag] = pcg(la, b, 1e-8, 100, pfun);
[x] = pcg(la, b, 1e-8, 100, pfun);

