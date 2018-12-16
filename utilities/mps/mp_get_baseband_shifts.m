function pb = mp_get_baseband_shifts(p, deg)
% MP_GET_BASEBAND_SHIFTS Returns baseband phase shifts.
%   pb = MP_GET_BASEBAND_SHIFTS(p) rewinds phase shifts p in degrees into
%   the base interval [0,360] degrees.
%
%   pb = MP_GET_BASEBAND_SHIFTS(p, deg), where deg is a string, specifies
%   which units the input phases in p use. deg may be 'deg' for degrees,
%   'rad' for radians, and 'grad' for gradians. If omitted or empty then
%   'deg' is assumed.
%
%   See also MP_ARE_PHASE_STEPS_UNIQUE.

% $Revision: 1.0 $  $Date: 2017/03/14 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);

% Set phase conversion factor to get radians.
if (2 > nargin) || isempty(deg); deg = 'deg'; end;
assert( ischar(deg) );

if strcmpi(deg, 'deg')
    mod_value = 360;
elseif strcmpi(deg, 'rad')
    mod_value = 2 * pi;
elseif strcmpi(deg, 'grad')
    mod_value = 400;
else
    error('Unknown phase units!');
end

if iscell(p)
    pb = cell(size(p));
    for i = 1 : numel(p)
        assert( isnumeric(p{i}) );
        pb{i} = mod(p{i}, mod_value);
    end
else
    assert( isnumeric(p) );
    pb = mod(p, mod_value);
end

