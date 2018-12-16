function K = mp_get_baseband_steps(N, k)
% MP_GET_BASEBAND_STEPS Return baseband phase steps.
%   K = MP_GET_BASEBAND_STEPS(N, k) returns baseband phase steps K which
%   correspond to phase steps in array k for the given number of phase
%   shifts N.
%
%   See also MP_REGULAR_PHASE_SHIFTS, MP_ARE_PHASE_STEPS_UNIQUE.

% $Revision: 1.0 $  $Date: 2017/03/08 $
% $Author(s): Tomislav Petkovic $

narginchk(2, 2);

% Sanity check.
assert( isnumeric(N) && (1 == numel(N)) && (0 < N) && (N == round(N)) );
assert( isnumeric(k) && all(round(k) == k) );

% Move all phase steps to primary interval.
K = mod(k, N);

% Mirror positive frequencies higher than half sampling rate.
idx = find(K > N/2);
if ~isempty(idx)
    K(idx) = N - K(idx);
end