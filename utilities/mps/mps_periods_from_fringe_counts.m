function [lambda, width] = mps_periods_from_fringe_counts(count, width)
% MPS_PERIODS_FROM_FRINGE_COUNTS Returns periods from fringe counts.
%   lambda = MPS_PERIODS_FROM_FRINGE_COUNTS(count) returns an array of
%   whole numbers lambda which holds periods for fringes with fringe counts
%   stored in array count. All numbers stored in count must be whole.
%
%   For MPS methods each fringe pattern has a number lambda which describes
%   in how many intervals one fringe period has to be divided. For MPS
%   structured light strategies sometimes it is simpler to define fringe
%   patterns by the number of fringes that fit into one projector screen,
%   which is fringe count, instead of using lambda, which is a virutal
%   number that defines number of segments in one period. This functions
%   converts the number of fringes into number of segments per period.
%
%   [lambda, width] = MPS_PERIODS_FROM_FRINGE_COUNTS(count) also returns
%   maximal (virtual) width of the projector screen. When computing the
%   projector coordinate unwrapped value may be divided by width and
%   multiplied by true projector width in px to obtain required coordinate
%   in projector coordinate system. Note that maximal allowable width is
%   the least common multiple of all fringe counts.
%
%   [lambda, width] = MPS_PERIODS_FROM_FRINGE_COUNTS(count, width) uses
%   supplied width to compute periods. If width is supplied then both count
%   and lambda may be real numbers.
%
%   See also MPS_FRINGE_COUNTS_FROM_PERIOD, MPS_LCM.

% $Revision: 1.2 $  $Date: 2017/06/19 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);
nargoutchk(0, 2);

% Check inputs.
all_integers = all( round(count) == count );
have_width = (1 < nargin) && ~isempty(width);

if all_integers
    
    % Check if inputs are relatively prime integers and reduce them if they
    % are not.
    gcm_count = mps_gcd(count);
    if 1 < gcm_count
        warning('Counts are not relatively prime integers.');
        count = count ./ gcm_count;
    end
    
    % For integer counts the maximal allowed width is the least common
    % multiple of numbers stored in count.
    lcm_count = mps_lcm(count);
    if ~have_width; width = lcm_count;  end;
    
    assert( isnumeric(width) && (1 == numel(width)) && (0 < width) && (width <= lcm_count) );
    
else
    
    % For real counts the maximal allowed width is not well defined
    % therefore we warn the user that width should be selected. If the
    % width is not supplied then compute one as least common multiple of
    % integer counts only.
    if (1 == nargin) || ~have_width
        warning('Width should be defined for non-integer fringe counts!');
        
        lcm_count = mps_lcm( count(round(count) == count) );
        width = lcm_count;
    end
    
    assert( isnumeric(width) && (1 == numel(width)) && (0 < width) );
    
end

% Compute periods.
lambda = width ./ count;

% Reduce periods.
if all_integers
    
    gcd_lambda = mps_gcd(lambda);
    if 1 < gcd_lambda
        warning('Periods are not relatively prime integers!');
        lambda = lambda ./ gcm_lambda;
    end
    
    lcm_lambda = mps_lcm(lambda);
    if width > lcm_lambda
        warning('Reducing width to LCM of periods!');
        width = lcm_lambda;
    end
        
    assert(lcm_lambda == lcm_count);   
    assert( all( round(lambda) == lambda ) );
end