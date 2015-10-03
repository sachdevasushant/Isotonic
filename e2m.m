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


default('s',1);

default('nv',max(max(e)));

default('w',1);

ne = size(e,1);

m = sparse(e(:,1),1:ne,w,nv,ne);
m = m + s*sparse(e(:,2),1:ne,w,nv,ne);
