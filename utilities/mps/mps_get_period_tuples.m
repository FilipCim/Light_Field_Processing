function [k, width, w_min, w_max] = mps_get_period_tuples(lambda, width)
% MPS_GET_PERIOD_TUPLES Returns all valid period tuples.
%   k = MPS_GET_PERIOD_TUPLES(lambda) returns all valid period tuples k for
%   input periods lambda. All elements of lambda must be whole numbers. If
%   there are M elements in lambda then k is N x M matrix whose each row
%   contains one of M-tuples of periods.
%
%   k = MPS_GET_PERIOD_TUPLES(lambda, width) uses width as the maximal
%   allowed unwrapped value and returns only tuples k that correspond to
%   unwrapped values in the [0,width] interval. If width is ommited or
%   empty then maximal allowable width is used. If width is supplied then
%   elements of lambda may also be real numbers.
%
%   [k, width] = MPS_GET_PERIOD_TUPLES(lambda) also returns used width.
%
%   [k, width, w_min, w_max] = MPS_GET_PERIOD_TUPLES(lambda) also returns
%   bounaries for each of generated period tuples.
%
%   See also MPS_GET_LINES, MPS_NUMBER_OF_TUPLES.

% $Revision: 1.2 $  $Date: 2016/11/23 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);
nargoutchk(0, 4);

lambda = reshape(lambda, 1, numel(lambda));

% Set default width.
if all( round(lambda) == lambda )
    
    % Check if input periods (wavelengts) are relatively prime integers and
    % reduce them if they are not.
    gcm_lambda = mps_gcd(lambda);
    if 1 < gcm_lambda
        warning('Periods are not relatively prime integers.');
        lambda = lambda ./ gcm_lambda;
    end    
    
    width_max = mps_lcm(lambda);
    if (2 > nargin) || isempty(width); width = width_max; end;
    assert( isnumeric(width) && (1 == numel(width)) );
end

% Wavelengths should normally be whole numbers so number-theoretical
% generator may be used to generate period-tuples. If the wavelengths are
% not whole numbers then we must use a different generator for period
% tuples. Additionally, the number-theoretical method of generating period
% tuples is quite slow if many wavelengths are used so we use this code
% only for two wavelengths.
if all( round(lambda) == lambda ) && (round(width) == width) && (numel(lambda) < 3)
    
    if width > width_max
        error(['Invalid input! Supplied width must be in [0,' num2str(width_max) '] interval']);
    end
    assert( (0 <= width) && (width <= width_max) );
    
    % Pre-allocate output.
    N = mps_number_of_tuples(lambda, width);
    k = zeros(N, numel(lambda), 'int32');
    
    % Create valid tuples. Note we use step of 1 which is always valid; to
    % speed up the generation of tuples the step of the for loop may be
    % increased from 1 to the greatest common divisor of all lambda.
    j = 2;
    for i = min(lambda) : width - 1
        inc = ( 0 == rem(i, lambda) );
        if any(inc)
            k(j, :) = k(j-1, :) + int32(inc);
            j = j + 1;
        end
    end
    assert( N + 1 == j );
    
else
    
    % Check width.
    if (0 == exist('width','var')); error('Width must be supplied if lambdas are not whole numbers!'); end;
    assert( isnumeric(width) && (1 == numel(width)) && (0 < width) );
    
    % Create valid tuples.
    k = zeros(1, numel(lambda), 'int32');
    j = 2;
    next_boundary = lambda;
    x = min(next_boundary);
    inc = next_boundary <= x;
    while (x < width) && any(inc)
        next_boundary(inc) = next_boundary(inc) + lambda(inc);
        k(j, :) = k(j - 1, :) + int32(inc);
        j = j + 1;
        x = min(next_boundary);
        inc = next_boundary <= x;
    end
    
    % Get number of period tuples.
    N = size(k, 1);
    
end

% Compute interval boundaries.
if (2 < nargout)
    stops = double(k);
    scale = double(lambda);
    for i = 1 : N; stops(i, :) = stops(i, :) .* scale; end;
    w_min = max(stops, [], 2);
    if (3 < nargout)
        w_max = w_min(2:end, :);
        w_max(N, :) = width;
    end
end