function [d, k] = mp_are_phases_from_dft(p, deg, thr)
% MP_ARE_PHASES_FROM_DFT Tests if phase shifts are from DFT in N points.
%   d = MP_ARE_PHASES_FROM_DFT(p) tests if phase shifts p are from DFT in N
%   points. The input p is a cell array or matrix having N rows and K
%   columns.
%
%   Phase shifts in p are from DFT in N points is they match to the k-th
%   basis vector.
%
%   d = MP_ARE_PHASES_FROM_DFT(p, deg), where deg is a string, specifies
%   which units the input phases in p use. deg may be 'deg' for degrees,
%   'rad' for radians, and 'grad' for gradians. If omitted or empty then
%   'deg' is assumed.
%
%   d = MP_ARE_PHASES_FROM_DFT(p, deg, thr) specifies a threshold to detect
%   DFT vector. Threshold is applied to the sum of absolute differenced
%   between supplied and DFT basis vectors. If empty or omitted thr is set
%   to 0.001.
%
%   [d, k] = MP_ARE_PHASES_FROM_DFT(p) also returns DFT vector indices
%   stored in k. The output k is a row vector which has as many elements as
%   there are columns in p. The value in k is -1 if particular column from
%   p is not a DFT basis vector and is the index of the DFT basis vector
%   otherwise. Note that k starts from 0 and ends in N-1.
%
%   See also MP_GET_BASEBAND_SHIFTS.

% $Revision: 1.0 $  $Date: 2017/03/14 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 3);

% Set phase conversion factor to get radians.
if (3 > nargin) || isempty(deg); deg = 'deg'; end;
assert( ischar(deg) );

if strcmpi(deg, 'deg')
    p_to_rad = single(pi / 180);
elseif strcmpi(deg, 'rad')
    p_to_rad = single(1.0);
elseif strcmpi(deg, 'grad')
    p_to_rad = single(pi / 200);
else
    error('Unknown phase units!');
end

% Set threshold.
if (4 > nargin) || isempty(thr); thr = 0.001; end;
assert( isnumeric(thr) && (1 == numel(thr)) && (0 < thr) );

% Convert input to radians.
if iscell(p)
    P = single(cell2mat(p)) * p_to_rad;
else
    if isvector(p)
        P = single(p(:)) * p_to_rad;        
    else
        P = single(p) * p_to_rad;
    end
end

P = mp_get_baseband_shifts(P, 'rad');

% Get DFT size and basis angles.
N = size(P,1);
PB = angle(ifft(eye(N)));
PB(PB < 0) = PB(PB < 0) + 2 * pi;

% Get indices.
M = size(P,2);
k = zeros(1, M);
for i = 1 : M
    diffs = bsxfun(@minus, PB, P(:,i));
    dst = sum(abs(diffs));
    [dst_best, idx] = min(dst);
    if dst_best < thr
        k(i) = idx;
    end
end
k = k - 1;

if all( (0 <= k) & (k < N) )
    d = true;
else
    d = false;
end
