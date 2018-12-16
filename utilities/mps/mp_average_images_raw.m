function [A, A_type] = mp_average_images_raw(X, X_type, Y, Y_type)
% MP_AVERAGE_IMAGES_RAW Combines multiple raw images into one.
%   [A, A_type] = MP_AVERAGE_IMAGES_RAW(X, X_type, Y, Y_type) returns the
%   average A of X and Y.
%
%   If the inputs X and Y are images of the same size then A is (X+Y)/2.
%
%   If any of the input X or Y is a cell array then that cell array must
%   contain images of the same size. The output A is then the average of
%   all input images in X and Y.
%
%   This function is intended to combine DC components of the MPS into one
%   image which may be used as a texture for reconstruction. DC components
%   may be obtained from MP_WRAPPED_PHASE or MPS_WRAPPED_PHASE functions.
%
%   [A, A_type] = MP_AVERAGE_IMAGES_RAW(X, X_type) averages images in cell
%   array X only.
%
%   See also MP_AVERAGE_IMAGES, MP_WRAPPED_PHASE, MPS_WRAPPED_PHASE.

% $Revision: 1.0 $  $Date: 2017/07/20 $
% $Author(s): Tomislav Petkovic $

narginchk(2, 4);
nargoutchk(2, 2);

% Start with empty output.
A = [];
A_type = [];
N = 0;

% Accumulate inputs.
[A, A_type, N] = accumulate(X, X_type, A, A_type, N);
if 2 < nargin; [A, A_type, N] = accumulate(Y, Y_type, A, A_type, N); end;
A_type = ps_raw_format_from_matrix(A, A_type);

% Scale output.
if 0 < N
    A = (1 / N) * A;
else
    A = [];
end

end


function [A, A_type, N] = accumulate(X, X_type, A, A_type, N)
% ACCUMULATE Accumulates values in X.
%   [A, A_type, N, sz] = ACCUMULATE(X, X_type, A, A_type, N) accumulates
%   values from X into the accumulator A. The counter N is increased for
%   the number of items in X. Pixel data types are checked for consistency.

% Initialize accumulator to zero or get accumulator size.
if isempty(A)
    A = 0;
    A_type = [];
    sz = [];
else
    sz = size(A);
    assert( ischar(A_type) );
end

if iscell(X)
    assert(iscell(X_type));
    
    % Accumulate each item from a cell array.
    for i = 1 : numel(X)
        assert( isnumeric(X{i}) && ischar(X_type{i}) );
        if ~isempty(X{i})
            % Validate image size.
            if ~isempty(sz)
                assert( all(size(X{i}) == sz) );
            else
                sz = size(X{i});
            end
            
            % Validate image type.
            if ~isempty(A_type)
                assert( strcmpi(A_type, X_type{i}) );
            else
                A_type = X_type{i};
            end
            
            % Accumulate.
            A = A + double(X{i});
            N = N + 1;
        end
    end
elseif isnumeric(X)
    assert(ischar(X_type));
    
    % Add X into accumulator.
    if ~isempty(X)
        % Validate image size.
        if ~isempty(sz); assert( all(size(X) == sz) ); end;
        
        % Validate image type.
        if ~isempty(A_type)
            assert( strcmpi(A_type, X_type{i}) );
        else
            A_type = X_type{i};
        end
        
        % Accumulated.
        A = A + double(X);
        N = N + 1;
    end
else
    % Abort execution.
    error('Unsupported input datatype.');
end

end