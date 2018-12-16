function [P, idxk, d] = mps_unwrap_phase_NN(WP, k, O, Xk, Xv, n, wgt)
% MPS_UNWRAP_PHASE_NN Unwrap phase using the nearest neighbour approach.
%   P = MPS_UNWRAP_PHASE_NN(WP, k, O, Xk, Xv) returns normalized unwrapped
%   phase P from the set of wrapped phases stored in WP. Inputs are cell
%   array WP of wrapped phases and decoding parameters as follows:
%   k is the array holding valid period tuples as returned by
%   MPS_GET_PERIOD_TUPLES and O, Xk, and Xv are decoding parameters as
%   returned by MPS_GET_PROJECTION_MATRIX_AND_CENTERS.
%
%   This function implements optimal decoding based on minimising the true
%   Euclidean distance to the projected constellation of frequencies. It is
%   superior to L1 and pure integer methods.
%
%   [P, idxk, d] = MPS_UNWRAP_PHASE_NN(WP, k, O, Xk, Xv) also returns the
%   index of the best period-order vector in idxk and the minimal distance
%   d in the wrapped phase space.
%
%   P = MPS_UNWRAP_PHASE_NN(WP, k, O, Xk, Xv, n) uses n as fringe count
%   over the whole screen. If n is empty or omitted then n is computed from
%   k using n = max(k) + 1.
%
%   P = MPS_UNWRAP_PHASE_NN(WP, k, O, Xk, Xv, ~, wgt) applies weights wgt
%   when computing the unwraped phase P.
%
%   See also MPS_GET_PROJECTION_MATRIX_AND_CENTERS, MPS_GET_PERIOD_TUPLES,
%   MPS_GET_WEIGHTS.

% $Revision: 1.1 $  $Date: 2016/06/13 $
% $Author(s): Tomislav Petkovic $

narginchk(5, 7);
nargoutchk(1, 3);

% Check input wrapped phases for consistency and assemble wrapped phase matrix.
assert( iscell(WP) );
N = numel(WP);
sz = size(WP{1});
M = prod(sz);
WPall = zeros(N, M);
for i = 1 : N
    assert( ismatrix(WP{i}) && all(sz == size(WP{i})) );
    WPall(i, :) = reshape(WP{i}, 1, M);
end

% Check decoding parameters.
assert( ismatrix(k) );
K = size(k, 1);
assert( size(k, 2) == N );
assert( all(size(O) == [N-1 N]) );
assert( all(size(Xk) == [N-1 K]) );
assert( iscell(Xv) && (3 == size(Xv,2)) );

% Extend the period order tuples in Xk by adding data from Xv.
Xk_ext = Xk;
Xv_first = zeros(N - 1, 0);
k_ext = k;
k_map = int32(1 : K);
j = K + 1;
for i = 1 : size(Xv, 1)
    if 1 == size(Xv{i,1}, 2)
        Xk_ext(:, j) = Xv{i, 1};
        k_ext(j, :) = Xv{i, 2};
        k_map(j) = Xv{i, 3};
        j = j + 1;
    else
        assert( isempty(Xv_first) );
        Xv_first = Xv{i, 1};
        assert(1 == Xv{i, 3});
    end
end
assert( ( (j == K + size(Xv, 1)) && ~isempty(Xv_first) ) ~= ( (j - 1 == K + size(Xv, 1)) && isempty(Xv_first) ) );
K_ext = size(k_ext, 1);

% Ortographically project wrapped phases to the space of co-dimension 1.
WPO = O * WPall;

% Compute distances to all center points.
assert(K_ext < pow2(31) - 1);
idxk_ext = zeros(1, M, 'int32');
d = Inf(1, M);
if 1 == size(Xk_ext, 1)
    
    % For two frequency MPS decoding is especially simple as it reduces to
    % distances in dimension 1, i.e. Euclidean distance reduces to simple
    % absolute value.
    
    % Decode centers.
    for i = 1 : K_ext
        dst = abs(WPO - Xk_ext(i));
        better = dst < d;
        idxk_ext(better) = i;
        d(better) = dst(better);
    end
    
    % Clip vertices.
    for i = 1 : size(Xv_first, 2)
        dst = abs(WPO - Xv_first(i));
        better = dst < d;
        idxk_ext(better) = 1;
        d(better) = dst(better);
    end
    
else
    
    % For three or more frequency MPS decoding must compute true Euclidean
    % distances. This may be done using k-d trees or sequentially. Default
    % choice is k-d tree (if available).    
    if (2 == exist('knnsearch')) && (2 == exist('KDTreeSearcher'))
        
        % Always construct full k-d tree as we have many more points to
        % test than there is centers in Xk_ext.
        ns = KDTreeSearcher(Xk_ext.', 'Distance', 'euclidean', 'BucketSize', 1);
        [idxk_ext, d] = knnsearch(ns, WPO.');
        
    else
        
        % To speed up sequential search due to fixed constellation we may
        % remove any points that are closer than the minimal constellation
        % distance from the further distance computations.
        dst_min = mps_get_minimal_distance(Xk_ext, Xv_first).^2;
        %assert( dst_min <= mps_get_minimal_distance(Xk, Xv).^2 );
        unclassified = int32(1 : M);
        
        % Sequentially decode centers.
        for i = 1 : K_ext
            dst = sum( bsxfun(@minus, WPO(:, unclassified), Xk_ext(:, i)).^2, 1 );
            better = dst < d(unclassified);
            idx = unclassified(better);
            if 0 < numel(idx)
                idxk_ext(idx) = i;
                d(idx) = dst(better);
            end
            not_optimal = dst > dst_min;
            unclassified = unclassified(not_optimal);
        end
    end
    
    % Clip vertices.
    for i = 1 : size(Xv_first, 2)
        dst = sum( bsxfun(@minus, WPO, Xv_first(:, i)).^2, 1 );
        better = dst < d;
        idxk_ext(better) = 1;
        d(better) = dst(better);
    end
    
end

% Define fringe order numbers if not given.
if (6 > nargin) || isempty(n); n = double(max(k) + 1); end;
assert( isnumeric(n) && (N == numel(n)) && all(0 < n) );

% Unwrap phase. Note that we combine all unwrapped phases using the same
% weights. If noise is known for each of wrapped phases then weighted
% combination may be used to further improve wrapped phase estimate.
kpi = 2 * pi * double(k_ext);
valid = (1 <= idxk_ext) & (idxk_ext <= K_ext);
idxk_ext_valid = idxk_ext(valid);
accP = 0;
if (7 > nargin) || isempty(wgt)
    scl = 1 ./ ( 2 * pi * N * n(:) );
    for i = 1 : N
        tmpkpi = kpi(:, i);
        accP = accP + scl(i) * ( WP{i}(valid) + reshape(tmpkpi(idxk_ext_valid), size(idxk_ext_valid)) );
    end
else
    assert( isnumeric(wgt) && (numel(wgt) == N) && all(0 <= wgt) );
    scl = wgt(:) ./ ( 2 * pi * n(:) ) ./ sum(wgt(:));
    for i = 1 : N
        tmpkpi = kpi(:, i);
        accP = accP + scl(i) * ( WP{i}(valid) + reshape(tmpkpi(idxk_ext_valid), size(idxk_ext_valid)) );
    end
end

% Assemble output.
P = NaN(sz);
P(valid) = accP;
P(P < 0) = P(P < 0) + 1;

if (1 < nargout); idxk = reshape(k_map(idxk_ext), sz); end;
if (2 < nargout); d = reshape(sqrt(d), sz); end;