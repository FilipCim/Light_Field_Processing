function gcd = mps_gcd(a, b)
% MPS_GCD Computes the greatest common divisor of two whole numbers.
%   gcd = MPS_GCD(a, b) returns the greatest common divisor of two whole
%   numbers a and b. Function uses Euclidean algorithm.
%
%   gcd = MPS_GCD(x) where x is a vector of whole numbers returns the
%   greatest common divisor of all numbers in x. All elements in x must be
%   whole numbers.
%
%   See also MPS_LCM, MPS_EXTENDED_GCD.

% $Revision: 1.3 $  $Date: 2017/04/20 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);
nargoutchk(0, 1);

% If there are two inputs then compute the greatest common divisor
% directly. Otherwise, if there is only one input argument with more than
% two elements call MPS_GCD iteratively as the greatest common divisor may
% be computed recursively.
if 2 == nargin
    
    assert( isnumeric(a) && (1 == numel(a)) );
    assert( isnumeric(b) && (1 == numel(b)) );    
    
    % Test if inputs are whole numbers.
    a_is_whole = (round(a) == a);
    b_is_whole = (round(b) == b);
    if (false == a_is_whole); a = round(a); end;
    if (false == b_is_whole); b = round(b); end;
    if (false == a_is_whole) || (false == b_is_whole); warning('Rounding input(s) to the closest whole number (s)!'); end;
    
    % Test if input numbers have all digits.
    assert( isa(a, class(b)) );    
    if isa(a, 'double')
        limit = pow2(53); % 53 binary digits
        assert( (a < limit) && (b < limit) );
    elseif isa(a, 'single')
        limit = pow2(24); % 24 binary digits
        assert( (a < limit) && (b < limit) );
    end    
    
    % The second input must be non-zero to avoid infinite loop. Note that
    % the first input may be zero in which case the while loop terminates
    % on its first iteration.
    assert( 0 ~= b );
    
    % Euclidean algorithm.
    r0 = a;
    r1 = b;
    r = 1;
    
    while 0 ~= r
        r = rem(r0, r1);
        r0 = r1;
        r1 = r;
    end
    
    gcd = r0;
    assert( round(gcd) == gcd );

else
    
    assert( isnumeric(a) && all(a(:) == round(a(:))) );
    assert( all(a(:) ~= 0) );
    
    if 1 < numel(a)
        gcd = mps_gcd(a(1), a(2));
        for i = 3 : numel(a); gcd = mps_gcd(gcd, a(i)); end;
    elseif 1 == numel(a)
        gcd = a(1);
    else
        gcd = [];
    end
    
end