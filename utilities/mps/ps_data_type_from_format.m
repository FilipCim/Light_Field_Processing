function data_type = ps_data_type_from_format(format)
% PS_DATA_TYPE_FROM_FORMAT Returns data type description string.
%   data_type = PS_DATA_TYPE_FROM_FORMAT(format) returns description string
%   which describes the data type indicated by the string format.
%
%   See also PS_DATA_TYPE_FROM_MATRIX.

% $Revision: 1.0 $  $Date: 2017/09/08 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 1);

assert( ischar(format) );

if have_one_instance_of_substring(format, '64F')
    data_type = '64F';
elseif have_one_instance_of_substring(format, '32F')
    data_type = '32F';
elseif have_one_instance_of_substring(format, '8U')
    data_type = '8U';
elseif have_one_instance_of_substring(format, '16U')
    data_type = '16U';
elseif have_one_instance_of_substring(format, '32U')
    data_type = '32U';
elseif have_one_instance_of_substring(format, '64U')
    data_type = '64U';
elseif have_one_instance_of_substring(format, '8S')
    data_type = '8S';
elseif have_one_instance_of_substring(format, '16S')
    data_type = '16S';
elseif have_one_instance_of_substring(format, '32S')
    data_type = '32S';
elseif have_one_instance_of_substring(format, '64S')
    data_type = '64S';
else
    error('Unsupported input datatype.');
end

end


function TF = have_one_instance_of_substring(str, substr)
% HAVE_ONE_INSTANCE_OF_SUBSTRING Checks if substring appears only one in a string.
%   TF = HAVE_ONE_INSTANCE_OF_SUBSTRING(str, substr) returns true if a
%   substring substr appears exactly once in a string str and false
%   otherwise. Comparison is not case sensitive.

assert( ischar(str) && ischar(substr) );
K = strfind(upper(str), upper(substr));
TF = (1 == numel(K));

end