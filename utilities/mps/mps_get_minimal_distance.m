function d_min = mps_get_minimal_distance(Xk, Xv)
% MPS_GET_MINIMAL_DISTANCE Returns minimal half-distance between centers.
%   d_min = MPS_GET_MINIMAL_DISTANCE(Xk, Xv) returns minimal distance
%   between points stored in inputs Xk and Xk. The input Xk stores points
%   that correspond to centers of valid period tuples. The input Xv is a
%   cell array which stores points that corresponds to vertices on the
%   edges of the hypercube (less the vertices at 0 and 2pi) and their
%   correspondent wrapped period order vectors and indices of true period
%   order vector. The first column of Xv stores points, the second column
%   of Xv stores wrapped period order vectors, and the third stores index
%   of true period order vector.
%
%   Returned distance d_min may be used to speed up classification: if the
%   test point is closer to the period tuple center point than d we know it
%   certainly belongs to that period tuple as all remaining centers are
%   further away.
%
%   See also MPS_GET_PROJECTION_MATRIX_AND_CENTERS.

% $Revision: 1.2 $  $Date: 2017/03/07 $
% $Author(s): Tomislav Petkovic $

narginchk(2, 2);

% Sanity check.
N = size(Xk, 1);

if ~iscell(Xv)
    % To retain backward compatibility also accept non-cell array input.
    K = size(Xv, 2);
    assert( all(size(Xv) == [N K]) );
    Xvtmp = Xv;
    Xv = cell(1, 2);
    Xv{1,1} = Xvtmp;
else
    K = size(Xv, 1);
    assert( all(size(Xv) == [K 3]) );
    for i = 1 : K
        assert( size(Xv{i,1}, 1) == N );        
    end
end

% Assemble points.
X = Xk;
for i = 1 : size(Xv, 1)
    X = cat(2, X, Xv{i,1});
end

% Compute pairwise distances and retain minimal one.
M = size(X, 2);
d_min = Inf;
if 1 == N
    
    for i = 1 : M
        dst = min( abs( X(i+1:end) - X(i) ) );
        if dst < d_min; d_min = dst; end;
    end
    d_min = 0.5 * d_min;
    
else
    
    for i = 1 : M
        dst = sum( bsxfun(@minus, X(:,i+1:end), X(:, i)).^2, 1 );
        dst = min(dst);
        if dst < d_min; d_min = dst; end;
    end
    d_min = 0.5 * sqrt(d_min);
    
end

% Return lower real number as the limit.
d_min = d_min - eps(10*d_min);
assert(0 < d_min);