function raw_format = ps_raw_format_from_matrix(img, img_format)
% PS_RAW_FORMAT_FROM_MATRIX Returns description string for input matrix.
%   raw_format = PS_RAW_FORMAT_FROM_MATRIX(img) returns string which
%   describes the memory format of the image img.
%
%   raw_format = PS_RAW_FORMAT_FROM_MATRIX(img, img_format) also takes
%   img_format into consideration which is needed to distinguish between
%   one-channel formats such as grayscale and Bayer formats.
%
%   See also PS_READ_RAW, PS_WRITE_RAW, PS_DATA_TYPE_FROM_MATRIX.

% $Revision: 1.0 $  $Date: 2017/07/20 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);

% Set default input.
if (2 > nargin) || isempty(img_format); img_format = ''; end;
assert( ischar(img_format) );

% Get machine endianness.
[~, ~, endian] = computer;

% Infere format depending on the datatype.
raw_format = 'IDT_UNKNOWN';
if ismatrix(img)
    
    if ps_is_bayer(img_format)
        data_type = ps_data_type_from_matrix(img);
        [~, bayer_suffix] = ps_bayer_tile_from_format(img_format);
        if ~isempty(data_type) && ~isempty(bayer_suffix)
            raw_format = cat(2, 'IDT', '_', data_type, '_', bayer_suffix);
            if 'B' == endian
                raw_format = cat(2, raw_format, '_', 'BigEndian');
            end
        end
    elseif isa(img, 'logical')
        raw_format = 'IDT_8U_BINARY';
    elseif isa(img, 'uint8') || isa(img, 'int8')
        data_type = ps_data_type_from_matrix(img);
        if ~isempty(data_type)
            raw_format = cat(2, 'IDT', '_', data_type, '_', 'GRAY');
        end
    else
        data_type = ps_data_type_from_matrix(img);
        if ~isempty(data_type)
            raw_format = cat(2, 'IDT', '_', data_type, '_', 'GRAY');
            if 'B' == endian
                raw_format = cat(2, raw_format, '_', 'BigEndian');
            end
        end
    end
    
elseif (3 == ndims(img)) && (3 == size(img,3))
    
    assert( ~ps_is_bayer(img_format) );
    data_type = ps_data_type_from_matrix(img);
    if ~isempty(data_type)
        raw_format = cat(2, 'IDT', '_', data_type, '_', 'RGB', '_', 'Planar');
        if 'B' == endian
            raw_format = cat(2, raw_format, '_', 'BigEndian');
        end
    end
    
end

if strcmpi(raw_format, 'IDT_UNKNOWN')
    warning('No known raw format matches the input.');
end