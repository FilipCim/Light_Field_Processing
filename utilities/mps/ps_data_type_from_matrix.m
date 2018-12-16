function data_type = ps_data_type_from_matrix(img)
% PS_DATA_TYPE_FROM_MATRIX Returns data type description string.
%   data_type = PS_DATA_TYPE_FROM_MATRIX(img) returns description string
%   which describes the data type of matrix img.
%
%   See also PS_RAW_FORMAT_FROM_MATRIX, PS_DATA_TYPE_FROM_FORMAT.

% $Revision: 1.0 $  $Date: 2017/07/20 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 1);

if isa(img, 'logical')
    data_type = '8U';
elseif isa(img, 'uint8')
    data_type = '8U';
elseif isa(img, 'uint16')
    data_type = '16U';
elseif isa(img, 'uint32')
    data_type = '32U';
elseif isa(img, 'int8')
    data_type = '8S';
elseif isa(img, 'int16')
    data_type = '16S';
elseif isa(img, 'int32')
    data_type = '32S';
elseif isa(img, 'double')
    data_type = '64F';
elseif isa(img, 'single')
    data_type = '32F';
else
    error('Unsupported input datatype.');
end