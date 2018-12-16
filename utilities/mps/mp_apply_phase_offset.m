function pb = mp_apply_phase_offset(p, offset, deg)
% MP_APPLY_PHASE_OFFSET Applies phase offset and returns baseband shifts.
%   pb = MP_APPLY_PHASE_OFFSET(p, offset) computes p + offset and then
%   rewinds offseted phase shifts into the base interval of [0,360]
%   degrees. Both p and offset must be in degrees.
%
%   pb = MP_APPLY_PHASE_OFFSET(p, offset, deg), where deg is a string,
%   specifies which units the input phases in p use. deg may be 'deg' for
%   degrees, 'rad' for radians, and 'grad' for gradians. If omitted or
%   empty then 'deg' is assumed.
%
%   See also MP_GET_BASEBAND_SHIFTS.

% $Revision: 1.0 $  $Date: 2017/03/14 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 3);

assert( isnumeric(offset) && (1 == numel(offset)) );

% Apply offset.
if iscell(p)
    po = cell(size(p));
    for i = 1 : numel(p)
        assert( isnumeric(p{i}) );
        po{i} = p{i} + offset;
    end    
else
    assert( isnumeric(p) );
    po = p + offset;    
end

% Shift offsets to baseband.
if (3 > nargin)
    pb = mp_get_baseband_shifts(po);
else
    pb = mp_get_baseband_shifts(po, deg);
end