function [k, width, idx] = mps_get_wrapped_tuples(lambda, width, method)
% MPS_GET_WRAPPED_TUPLES Returns valid wrapped period tuples.
%   k = MPS_GET_WRAPPED_TUPLES(lambda) returns all valid wrapped period
%   tuples k for input periods lambda. All elements of lambda must be whole
%   numbers. If there are M elements in lambda then k is N x M matrix whose
%   each row contains one of M-tuples of wrapped periods.
%
%   k = MPS_GET_WRAPPED_TUPLES(lambda, width) uses width as the maximal
%   allowed unwrapped value and returns only wrapped tuples k that
%   correspond to unwrapped values in the [0,width] interval. If width is
%   ommited or empty then maximal allowable width is used. If width is
%   supplied then elements of lambda may also be real numbers.
%
%   [k, width] = MPS_GET_WRAPPED_TUPLES(lambda) also returns used width.
%
%   [k, width, idx] = MPS_GET_WRAPPED_TUPLES(lambda) also returns the index
%   idx of the unwrapped period tuple as returned by MPS_GET_PERIOD_TUPLES.
%
%   k = MPS_GET_WRAPPED_TUPLES(lambda, width, method) uses selected method
%   to generate wrapped period order tuples. The input method may be
%   'generative' or 'vectorized'. Generative method uses a loop to generate
%   all wrapped tuples and is verbatim to the procedure BOUNDARYVECTORS
%   from the article T. Petkoviæ, T. Pribaniæ, M. Ðonliæ "Temporal phase
%   unwrapping using orthographic projection" doi:j.optlaseng.2016.09.006.
%   Vectorized method (default) is optimized for Matlab.
%
%   See also MPS_GET_PERIOD_TUPLES.

% $Revision: 1.1 $  $Date: 2016/11/25 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 3);
nargoutchk(0, 3);

% Set default method.
if (3 > nargin) || isempty(method); method = 'vectorized'; end;
assert( ischar(method) );

lambda = reshape(lambda, 1, numel(lambda));

% Set default width.
if all( round(lambda) == lambda )
    width_max = mps_lcm(lambda);
    if (2 > nargin) || isempty(width); width = width_max; end;
    assert( (0 <= width) && (width <= width_max) );
end
assert( isnumeric(width) && (1 == numel(width)) && (0 < width) );


% Create outputs.
k = zeros(0, numel(lambda), 'int32');
j = 1;
idx = zeros(0, 1, 'int32');

% Create first wrapped tuple. The first wrapped tuple covers all corners of
% the hyper-cube except ones whose all coordinates are 0 or 2*pi.
k_tmp = zeros(1, numel(lambda), 'int32');
j_tmp = 1;
I = 1 : numel(lambda);
s = numel(I);
if strcmpi(method, 'generative')
    
    % Generative method using a loop.
    assert(s < 32);
    b = 1;
    while b < 2^s-1
        k(j, :) = k_tmp(1, :);
        idx(j, 1) = 1;
        for i = 1 : s
            if 1 == bitget(uint32(b), i)
                k(j, I(i)) = k(j, I(i)) - 1;
            end
        end
        j = j + 1;
        b = b + 1;
    end
    
    use_vectorized_code = false;
elseif strcmpi(method, 'vectorized') || strcmpi(method, 'vectorised')
    
    % Vectorized method.
    steps = -rem( floor((1:(2^s-2)).' * pow2(0:-1:-(s-1))), 2);
    k_steps = bsxfun(@plus, steps, double(k_tmp(1,:)));
    n_steps = size(steps, 1);
    k((1:n_steps)+j-1, :) = int32(k_steps);
    idx((1:n_steps)+j-1) = 1;
    j = j + n_steps;
    
    use_vectorized_code = true;
    
else
    error('Unknown method to generate wrapped tuples!');
end


% Create remaining wrapped tuples.
j_tmp = j_tmp + 1;
next_boundary = lambda;
x = min(next_boundary);
inc = next_boundary <= x;
while (x < width)
    k_tmp(j_tmp, :) = k_tmp(j_tmp - 1, :) + int32(inc);
    next_boundary(inc) = next_boundary(inc) + lambda(inc);
    if 1 < sum(inc)
        I = find(inc);
        s = numel(I);
        if (false == use_vectorized_code)
            assert(s < 32);
            b = 1;
            while b < 2^s-1
                k(j, :) = k_tmp(j_tmp, :);
                idx(j, 1) = j_tmp;
                for i = 1 : s
                    if 1 == bitget(uint32(b), i)
                        k(j, I(i)) = k(j, I(i)) - 1;
                    end
                end
                j = j + 1;
                b = b + 1;
            end
        else
            steps = -rem( floor((1:(2^s-2)).' * pow2(0:-1:-(s-1))), 2);
            n_steps = size(steps, 1);
            k_steps = repmat(double(k_tmp(j_tmp, :)), n_steps, 1);
            k_steps(:,I) = bsxfun(@plus, steps, k_steps(:,I));            
            k((1:n_steps)+j-1, :) = int32(k_steps);
            idx((1:n_steps)+j-1) = j_tmp;
            j = j + n_steps;
        end
    end
    j_tmp = j_tmp + 1;
    x = min(next_boundary);
    inc = next_boundary <= x;
    
    % Enforce memory limit of 500MB on generated wrapped tuples.
    if 4 * numel(k) > 500 * pow2(20)
        fprintf(2, 'Generated %d wrapped tuples.\nProcessed %.2f %% of the unwrapping interval.\n', size(k,1), x/width);
        error('Hardcoded memory limit of 500MB exceeded. Terminating!');
    end
end