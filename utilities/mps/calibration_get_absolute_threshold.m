function [abs_thr, scale] = calibration_get_absolute_threshold(L, rel_thr)
% CALIBRATION_GET_ABSOLUTE_THRESHOLD Returns absolute threshold.
%   abs_thr = CALIBRATION_GET_ABSOLUTE_THRESHOLD(L, rel_thr) returns
%   absolute threshold abs_thr that is computed from the relative threshold
%   rel_thr and image type in L. rel_thr must be a number in [0,1]
%   interval. L is a cell array where every element contains either image
%   data or image filename.
%
%   If all images stored in cell array are not of the same type then
%   abs_thr will be set to NaN.
%
%   [abs_thr, scale] = CALIBRATION_GET_ABSOLUTE_THRESHOLD(L, rel_thr) also
%   returns a number scale so abs_thr = scale * rel_thr;
%
%   See also CALIBRATION_SELECT_PIXELS.

% $Revision: 1.0 $  $Date: 2015/11/25 $
% $Author(s): Tomislav Petkovic $

narginchk(2, 2);
nargoutchk(1, 2);

assert( iscell(L) && (0 < numel(L)) );
assert( (1 == numel(rel_thr)) && (0 <= rel_thr) && (rel_thr <= 1) );

% Threshold scale depends on the input image datatype.
if ischar(L{1})
    img = calibration_load_image_as_gray(L{1});
    type_L = class(img);
else
    type_L = class(L{1});
end

all_valid = true;
for i = 2 : numel(L)
    if ischar(L{i})
        img = calibration_load_image_as_gray(L{1});
        type_i = class(img);
    else
        type_i = class(L{i});
    end
    if ~strcmpi(type_L, type_i)
        all_valid = false;
        break;
    end
end

if ~all_valid
    warning('Input images are NOT of the same type!');
    scale = NaN;
else
    scale = 1.0;
    if strcmpi(type_L, 'uint8'); scale = 2^8 - 1; end;
    if strcmpi(type_L, 'uint16'); scale = 2^16 - 1; end;
    if strcmpi(type_L, 'uint32'); scale = 2^32 - 1; end;
    if strcmpi(type_L, 'int8'); scale = 2^7 - 1; end;
    if strcmpi(type_L, 'int16'); scale = 2^16 - 1; end;
    if strcmpi(type_L, 'int32'); scale = 2^31 - 1; end;
end

abs_thr = scale * rel_thr;