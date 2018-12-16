function [V, X0, k, width, w_min, w_max] = mps_get_lines(lambda, width)
% MPS_GET_LINES Returns line equations for chosen set of frequencies.
%   [V, X0] = MPS_GET_LINES(lambda) returns set of lines for input periods
%   lambda. All elements of lambda must be whole numbers. Output lines are
%   defined by their direction vector V and points X0. If there are M
%   elements in lambda then both V is a N element vector and X0 is a M x N
%   matrix whose each column contains a M dimensional point, i.e. all
%   points on i-th output line are given by the equation
%      X0(:, i) + w * V,
%   where the real number w is line parameter. Note the line direction
%   vector V is not normalized.
%
%   [V, X0] = MPS_GET_LINES(lambda, width) uses supplied width and returns
%   only lines that correspond to period tuples in the interval
%   [0, width], i.e. only lines whose parameter w is in [0, width]
%   interval for wrapped phases in the hypercube [-pi, pi]^M. If width is
%   ommited or empty then the maximal allowable width is used. When width
%   is supplied then values is lambda may be real numbers.
%
%   [V, X0, k] = MPS_GET_LINES(lambda) also returns period tuples k that
%   correspond to each line.
%
%   [V, X0, k, width, w_min, w_max] = MPS_GET_LINES(lambda) also returns
%   used width and limits w_min and w_max on line parameter w for each of
%   returned lines. Limit values may be used to visualize the paritions
%   centers in M-dimensional space as line segments, i.e. for i-th line
%   compute and visualize the segment X0(:, i) + w * V for parameter w in
%   [w_min(i), w_max(i)] interval.
%
%   See also MPS_GET_PERIOD_TUPLES, MPS_NUMBER_OF_TUPLES.

% $Revision: 1.0 $  $Date: 2016/03/30 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);
nargoutchk(2, 6);

% Get period tuples.
if 1 == nargin; width = []; end;
[k, width, w_min, w_max] = mps_get_period_tuples(lambda, width);
N = size(k, 1);
M = size(k, 2);
assert( numel(lambda) == M );

% Compute line equations.
tau = 2 * pi;
X0 = - tau * double(k.');
V = tau ./ double(lambda(:));