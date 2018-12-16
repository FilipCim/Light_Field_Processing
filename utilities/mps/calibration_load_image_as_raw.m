function [raw, raw_format] = calibration_load_image_as_raw(filename)
% CALIBRATION_LOAD_IMAGE_AS_RAW Loads image.
%   [raw, raw_format] = CALIBRATION_LOAD_IMAGE_AS_RAW(filename) returns
%   image stored in the file specified. Input image may be in RAW format or
%   in any image format that is recognized by IMREAD. No spatial operations
%   are performed on the input, i.e. if the input RAW format is Bayer then
%   the output is also Bayer.
%
%   See also PS_READ_RAW.

% $Revision: 1.0 $  $Date: 2017/07/20 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 1);
nargoutchk(1, 2);

% Assume loading failed.
raw = [];
raw_format = [];

try
    
    % Get file extension.
    [~, ~, ext] = fileparts(filename);
    
    % Try to load RAW image.
    if strcmpi(ext, '.raw') || strcmpi(ext, '.xml')
        [raw, raw_format, buffer, buffer_format, buffer_stride, width, height] = ps_read_raw(filename);
    else
        raw = [];
        raw_format = [];
    end
    
    % If RAW load failed then fallback to imread.
    if isempty(raw)
        [img, map] = imread(filename);
        if ~isempty(map)
            raw = ind2rgb(img, map);
        else
            raw = img;
        end
        raw_format = ps_raw_format_from_matrix(raw);
        
        if ~isempty(raw) && ( isempty(raw_format) || strcmpi(raw_format, 'IDT_UNKNOWN') )
            error('Unsupported input datatype.');
        end
    elseif isempty(raw_format) || strcmpi(raw_format, 'IDT_UNKNOWN')
        error('Unsupported input datatype.');
    end
    
catch
    
    warning(['Invalid image file: ' filename]);
    
end