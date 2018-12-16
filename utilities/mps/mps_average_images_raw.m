function [T, T_type] = mps_average_images_raw(I, I_type, S, options)
% MPS_AVERAGE_IMAGES_RAW Averages MPS images to create texture.
%   [T, T_type] = MPS_AVERAGE_IMAGES_RAW(I, I_type) averages all images in
%   I to create a texture. The input I may be a cell array of strings where
%   each string is a valid image file name or may be a cell array of images
%   where every cell is a greyscale or RGB image. The input I_type is
%   simply a cell array of strings describing image datatype; if some
%   element of I contains a filename then correspoding I_type is ignored.
%   The output texture T is scaled to the full range of the datatype of the
%   first image in I.
%
%   [T, T_type] = MPS_AVERAGE_IMAGES_RAW(I, I_type, S) uses string S as
%   output data type.
%
%   [T, T_type] = MPS_AVERAGE_IMAGES_RAW(I, I_type, [], 'rescale') always
%   rescales average so it takes full dynamic range of the input datatype.
%
%   See also MPS_AVERAGE_IMAGES, CALIBRATION_LOAD_IMAGE_AS_RAW.

% $Revision: 1.0 $  $Date: 2017/07/20 $
% $Author(s): Tomislav Petkovic $

narginchk(2, 4);
nargoutchk(2, 2);

% Sanity check.
assert( iscell(I) && iscell(I_type) && all(size(I) == size(I_type)) );

% Set default datatype.
if (3 > nargin) || isempty(S); S = ''; end;
assert( ischar(S) );

if (4 > nargin) || isempty(options); options = ''; end;
assert( ischar(options) );

% Fetch images and sum them together.
N1 = [];
N2 = [];
N3 = [];
N4 = numel(I);
acc = single([]);
T_type = [];
for i = 1 : N4
    
    % Fetch image.
    if ischar(I{i})
        [raw, raw_type] = calibration_load_image_as_raw(I{i});
    else
        raw = I{i};
        raw_type = I_type{i};
    end
    assert( ~isempty(raw) );
    assert( ~isempty(raw_type) && ischar(raw_type) );
    
    % Get image information and verify buffers for consistency. We
    % currently support accumulation of identical buffers only.
    if isempty(S); S = class(raw); end;
    if isempty(N1); N1 = size(raw, 1); else assert(size(raw, 1) == N1); end;
    if isempty(N2); N2 = size(raw, 2); else assert(size(raw, 2) == N2); end;
    if isempty(N3); N3 = size(raw, 3); else assert(size(raw, 3) == N3); end;
    if isempty(T_type); T_type = raw_type; else assert(strcmpi(T_type, raw_type)); end;
    
    % Preallocate accumulator.
    if isempty(acc); acc = zeros([N1 N2 N3], 'single'); end;
    
    % Add image to accumulator.
    acc = acc + single(raw);
end

if strcmpi(options, 'rescale')
    
    % Scale amplitudes to [0,1] interval.
    acc_max = max(acc(:));
    acc_min = min(acc(:));
    acc_span = acc_max - acc_min;
    if 0 ~= acc_span
        scl = 1.0 / acc_span;
        acc_scaled = (acc - acc_min) * scl;
    else
        if 0 ~= acc_max
            scl = 1.0 / acc_max;
            acc_scaled = acc * scl;
        else
            acc_scaled = acc - acc_min;
        end
    end
    rescaled = true;
    
else
    
    % Scale by number of images.
    scl = 1.0 / N4;
    acc_scaled = acc * scl;
    rescaled = false;
    
end

% Convert to output type.
if strcmpi(S, 'uint8')
    if rescaled
        T = uint8( (pow2(8) - 1) * acc_scaled );
    else
        T = uint8( acc_scaled );
    end
    T_type = ps_raw_format_from_matrix(T, T_type);
elseif strcmpi(S, 'uint16')
    if rescaled
        T = uint16( (pow2(16) - 1) * acc_scaled );
    else
        T = uint16( acc_scaled );
    end
    T_type = ps_raw_format_from_matrix(T, T_type);
elseif strcmpi(S, 'uint32')
    if rescaled
        T = uint32( (pow2(31) - 1) * acc_scaled );
    else
        T = uint32( acc_scaled );
    end
    T_type = ps_raw_format_from_matrix(T, T_type);
elseif strcmpi(S, 'single')
    T = single(acc_scaled);
    T_type = ps_raw_format_from_matrix(T, T_type);
elseif strcmpi(S, 'double')
    T = double(acc_scaled);
    T_type = ps_raw_format_from_matrix(T, T_type);
else
    warning('Unsupported input data type; texture may be improperly scaled.');
    T = acc_scaled;
    T_type = 'IDT_UNKNOWN';
end
