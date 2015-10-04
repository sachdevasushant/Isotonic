function m = e2m(e,s,nv,w)
% function m = e2m(e,s,nv,w)
%
% converts edge list e to vertex-edge transfer matrix,
% so that vert-vec = m * edge-vec
%
% nv is number of vertices, default value is highest vertex in edge list
%
% if s is 1, or not there, all entries are 1.
% if s is -1, get signed matrix m s.t. m*m' is the laplacian
%
% if w is defined, entries are weighted by edge weights w
%
%    Converts an edge list to a vertex-edge matrix
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

default('s',1);

default('nv',max(max(e)));

default('w',1);

ne = size(e,1);

m = sparse(e(:,1),1:ne,w,nv,ne);
m = m + s*sparse(e(:,2),1:ne,w,nv,ne);
