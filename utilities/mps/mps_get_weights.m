function wgt = mps_get_weights(lambda, sigma)
% MPS_GET_WEIGHTS Function returns weights for ideal phase unwrapping.
%   wgt = MPS_GET_WEIGHTS(lambda) returns vector of weights wgt for the
%   selected periods lambda.
%
%   wgt = MPS_GET_WEIGHTS(lambda, sigma) also takes into consideration the
%   phase variance sigma for each measurement. The vector sigma must be of
%   the same size as the vector lambda.
%
%   Returned weights are defining a convex combination of unwrapped phases
%   according to Eq. (29) from the article T. Petkoviæ, T. Pribaniæ, M.
%   Ðonliæ "Temporal phase unwrapping using orthographic projection"
%   (doi:10.1016/j.optlaseng.2016.09.006) and are optimal if unwrapped
%   phases have zero mean Gaussian noise with known variances sigma. Note
%   that if wrapped phase noise variance is the same for all phases then
%   weights do not depend on it.
%
%   Another popular choice of weights is to simply select the phase of the
%   smallest wavelength as the true solution. That may be achieved using
%   the following code snippet:
%      wgt = zeros(size(lambda));
%      [~, idx_best] = min(lambda);
%      wgt(idx_best) = 1;
%
%   See also MPS_UNWRAP_PHASE_NN.

% $Revision: 1.0 $  $Date: 2016/06/21 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);

assert( isnumeric(lambda) );
N = numel(lambda);

% Compute weights under the assumption of the same deviation sigma.
wgt = ones(N, 1);
wgt = wgt ./ ( double(lambda(:)).^2 );

% If phase deviations are supplied rebalance weights.
if 2 == nargin
    assert( numel(lambda) == numel(sigma) );
    wgt = wgt ./ ( double(sigma(:)).^2 );
end

% Normalize weights.
wgt = wgt ./ sum(wgt(:));