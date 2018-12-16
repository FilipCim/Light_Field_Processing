function [D, D_type, Lmax, Lmin] = calibration_dynamic_range_raw(L, L_type, C)
% CALIBRATION_DYNAMIC_RANGE_RAW Estimate dynamic range of image set.
%   [D, D_type] = CALIBRATION_DYNAMIC_RANGE_RAW(L, L_type) returns
%   estimated dynamic range for the input image set. The input L is a cell
%   array where every element of the cell array contains image data or
%   image filename. The input L_type is a cell array of image formats.
%   The output D is an image matrix that contains dynamic range for every
%   pixel.
%
%   [D, D_type, Lmax, Lmin] = CALIBRATION_DYNAMIC_RANGE_RAW(L, L_type) also
%   returns minimal and maximal value in L. Note that D = Lmax - Lmin.
%
%   [D, D_type] = CALIBRATION_DYNAMIC_RANGE_RAW(L, L_type, C) multiplies
%   input images values with coefficients stored in the cell array C.
%   Possible application of premultiplication are to compensate for
%   different exposure times of input images or to compensate for different
%   number of color levels of input images.
%
%   See also CALIBRATION_DYNAMIC_RANGE, CALIBRATION_SELECT_PIXELS,
%   CALIBRATION_COMBINE_DYNAMIC_RANGES.

% $Revision: 1.0 $  $Date: 2017/07/20 $
% $Author(s): Tomislav Petkovic $

narginchk(2, 3);
nargoutchk(2, 4);

assert( iscell(L) && iscell(L_type) && all(size(L) == size(L_type)) );

if ( (3 > nargin) || isempty(C) ); C = cell(size(L)); end;
assert( all(size(L) == size(C)) );

% Get number of images and assert that we have the minimum number of images
% required for dynamic range estimation.
N = numel(L);
assert( 2 <= N );

% Initialize maximum and minimum.
Lmax = single(-Inf);
Lmin = single(Inf);

% Sequentially load all images and adjust maximum and minimum value as
% needed.
N1 = [];
N2 = [];
N3 = [];
D_type = [];
for i = 1 : N
    
    % Load image and convert it to double.
    if ischar(L{i})
        [raw, raw_type] = calibration_load_image_as_gray(L{i});
    else
        raw = L{i};
        raw_type = L_type{i};
    end
    assert( ~isempty(raw) );
    assert( ~isempty(raw_type) && ischar(raw_type) );
    
    % Get image information and verify buffers for consistency. We
    % currently support identical buffers only.
    if isempty(N1); N1 = size(raw, 1); else assert(size(raw, 1) == N1); end;
    if isempty(N2); N2 = size(raw, 2); else assert(size(raw, 2) == N2); end;
    if isempty(N3); N3 = size(raw, 3); else assert(size(raw, 3) == N3); end;
    if isempty(D_type); D_type = raw_type; else assert(strcmpi(D_type, raw_type)); end;
    
    % Premultiply if additional coefficients are available. Note that we
    % multiply coefficients instead of image values to reduce number of
    % multiplication operations.
    raw = single(raw);
    if ~isempty( C{i} )
        raw = raw * C{i};
    end
    
    % Adjust minimum and maximum value.
    Lmax = max(Lmax, raw);
    Lmin = min(Lmin, raw);
end

% Estimate dynamic range.
D = Lmax - Lmin;
if 1 < N3
    D = min(D, [], 3);
end
D_type = ps_raw_format_from_matrix(D, D_type);