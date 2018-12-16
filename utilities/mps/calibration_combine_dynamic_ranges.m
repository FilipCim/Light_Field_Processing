function R = calibration_combine_dynamic_ranges(A, B, method)
% CALIBRATION_COMBINE_DYNAMIC_RANGES Combines multiple dynamic ranges maps.
%   R = CALIBRATION_COMBINE_DYNAMIC_RANGES(A, B) returns single dynamic
%   range image R which is a combination of dynamic ranges images A and B.
%   The inputs A and B must have the same size. Then output R is simply
%   min(A,B).
%
%   Each of inputs A or B may also be a cell array of images. Then each
%   pixel of the output R is a minimum of all corresponding input pixels.
%
%   If any of input is an empty array it is ignored in the computation. The
%   output will be an empty array only if all inputs are empty.
%
%   R = CALIBRATION_COMBINE_DYNAMIC_RANGES(A, B, method) also accepts the
%   range selection method. The input method is a string which may be:
%     a) min - R is computed as the minimum of input values,
%     b) max - R is computed as the maximum of input values.
%
%   See also CALIBRATION_DYNAMIC_RANGE.

% $Revision: 1.0 $  $Date: 2017/06/28 $
% $Author(s): Tomislav Petkovic $

narginchk(2, 3);

% Set defaults.
if (3 > nargin) || isempty(method); method = 'min'; end;
assert( ischar(method) );

% Get combination method.
if strcmpi(method, 'min')
    process = @(x,y) min(x,y);
elseif strcmpi(method, 'max')
    process = @(x,y) max(x,y);
else
    error('Unknown combination method!');
end

% Process data.
if isnumeric(A) && isnumeric(B)
    
    if ~isempty(A) && ~isempty(B)
        
        % Sanity check.
        assert( all( size(A) == size(B) ) );
        
        % Convert to double and combine.
        a = double(A(:));
        b = double(B(:));
        R = process(a, b);
        R = reshape(R, size(A));
        
    elseif ~isempty(A) && isempty(B)
        
        % Set range to A as B is empty.
        R = double(A);
        
    elseif isempty(A) && ~isempty(B)
        
        % Set range to B as A is empty.
        R = double(B);
        
    else
        
        % Both inputs are empty matrices.
        R = [];
        
    end
    
elseif isnumeric(A) && iscell(B)
    
    % Set range to A and then sequentially update it using the data in B.
    R = double(A);
    for i = 1 : numel(B)
        R = calibration_combine_dynamic_ranges(R, B{i}, method);
    end
    
elseif iscell(A) && isnumeric(B)
    
    % Set range to B and then sequentially update it using the data in A.
    R = double(B);
    for i = 1 : numel(A)
        R = calibration_combine_dynamic_ranges(R, A{i}, method);
    end
    
elseif iscell(A) && iscell(B)
    
    % Set range to first element of A and then process the rest.
    if ~isempty(A); R = double(A{1}); else R = []; end;
    for i = 2 : numel(A)
        R = calibration_combine_dynamic_ranges(R, A{i}, method);
    end
    R = calibration_combine_dynamic_ranges(R, B, method);
    
end