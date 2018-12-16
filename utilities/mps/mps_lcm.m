function lcm = mps_lcm(a, b)
% MPS_LCM Computes the least common multiple of two whole numbers.
%   lcm = MPS_LCM(a, b) returns the least common multiple of two whole
%   numbers a and b.
%
%   lcm = MPS_LCM(x) where x is a vector of whole numbers returns the least
%   common multiple of all numbers in x. Note that all elements in x must
%   be whole numbers.
%
%   See also MPS_GCD, MPS_EXTENDED_GCD.

% $Revision: 1.1 $  $Date: 2016/11/23 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);
nargoutchk(0, 1);

% If there are two inputs then compute the least common multiple directly.
% Otherwise, if there is only one input argument with more than two
% elements call MPS_LCM iteratively as the least common multiple may be
% computed recursively.
if 2 == nargin
    
    assert( isnumeric(a) && (1 == numel(a)) );
    assert( isnumeric(b) && (1 == numel(b)) );
    
    % Test if inputs are whole numbers.
    a_is_whole = (round(a) == a);
    b_is_whole = (round(b) == b);
    if (false == a_is_whole); a = round(a); end;
    if (false == b_is_whole); b = round(b); end;
    if (false == a_is_whole) || (false == b_is_whole); warning('Rounding input(s) to the closest whole number!'); end;    

    % Test for a possibility of integer overflow.
    assert( isa(a, class(b)) );    
    if isa(a, 'double')
        limit = pow2(53); % 53 binary digits
        assert( (a < limit) && (b < limit) && (a * b <= limit) );
    elseif isa(a, 'single')
        limit = pow2(24); % 24 binary digits
        assert( (a < limit) && (b < limit) && (double(a) * double(b) <= limit) );        
    elseif isa(a, 'int32')
        limit = pow2(31) - 1;
        assert( (a < limit) && (b < limit) && (double(a) * double(b) <= limit) );
    elseif isa(a, 'int16')
        limit = pow2(15) - 1;
        assert( (a < limit) && (b < limit) && (double(a) * double(b) <= limit) );        
    elseif isa(a, 'int8')
        limit = pow2(8) - 1;
        assert( (a < limit) && (b < limit) && (double(a) * double(b) <= limit) );
    else
        warning('Unsupported data type! Skipping integer overflow tests.');
    end
    
    gcd = mps_gcd(a, b);
    lcm = abs(a * b) / gcd;
    
else
    
    assert( isnumeric(a) && all(a(:) == round(a(:))) );
    
    if (1 < numel(a))
        lcm = mps_lcm(a(1), a(2));
        for i = 3 : numel(a); lcm = mps_lcm(lcm, a(i)); end;
    elseif 1 == numel(a)
        lcm = a(1);
    else
        lcm = [];
    end
    
end