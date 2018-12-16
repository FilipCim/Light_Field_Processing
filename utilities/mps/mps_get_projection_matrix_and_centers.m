function [O, Xk, Xv] = mps_get_projection_matrix_and_centers(lambda, width)
% MPS_GET_PROJECTION_MATRIX_AND_CENTERS Returns ortographic projection and projected tuple center points.
%   [O, Xk, Xv] = MPS_GET_PROJECTION_MATRIX_AND_CENTERS(lambda, width)
%   returns ortographic projection matrix O, tuple center points Xk, and
%   points and wrapped period-order numbers stored in Xv that are required
%   for MPS unwrapping.
%
%   The orthographic projection matrix O is (N - 1) x N matrix, where N is
%   the number of fringe periods in lambda. The matrix O transforms wrapped
%   phases from the space of dimension N to the space of co-dimension 1.
%
%   The output Xk is (N - 1) x K matrix whose each column stores one
%   (N - 1) dimensional point which is a center of a corresponding period
%   tuple as returned by MPS_GET_PERIOD_TUPLES function, i.e. columns of Xk
%   correspond to rows of MPS_GET_PERIOD_TUPLES output.
%
%   The output Xv is a three-column cell array which stores wrapped-around
%   vertices and corresponding wrapped period-order numbers. Wrapped-around
%   vertices are stored in the first column of Xv and are located exactly
%   on the edges of the phase hypercube (except for two vertices having
%   coordinates all 0 and all 2pi, which correspond to the first and the
%   last regular period tuple). Wrapped period-order numbers are stored in
%   the second column of Xv. Third column stores the index of the next
%   true period-order number.
%
%   All elements of lambda must be whole numbers (unless width is
%   specified). The input width may be omitted or empty; if so then default
%   maximal allowed value is used. If width is given the the resulting
%   line/point constellations are ones limited to [0,width] interval.
%
%   Note there are two strategies to unwrap the wrapped phase using MPS:
%   a) using slicing hyperplanes which divide the wrapped phase hypercube
%      into disjoint partitions around central lines, and
%   b) using the property that all central lines are parallel which allows
%      projecting the hypercube onto hyperplane of co-dimension 1 and then
%      unwrapping the phase by using the nearest-neighbour method.
%   The first and more cumbersome approach is implemented in functions
%   MPS_GET_DECODING_TABLES_AND_SLICING_HYPERPLANES and
%   MPS_UNWRAP_PHASE_HYPERPLANES. The second and more elegant approach is
%   implemented in functions MPS_GET_PROJECTION_MATRIX_AND_CENTERS and
%   MPS_UNWRAP_PHASE_NN.
%
%   Note that functions MPS_UNWRAP_PHASE_L1, MPS_GET_SLICING_EQUATIONS, and
%   MPS_GET_CENTRAL_EQUATIONS are specialised implementation of the second
%   method b) for two frequency MPS.
%
%   See also MPS_UNWRAP_PHASE_NN, MPS_GET_LINES, MPS_GET_PERIOD_TUPLES.

% $Revision: 1.1 $  $Date: 2016/06/14 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);
nargoutchk(3, 3);

% Get lines.
if 1 == nargin; width = []; end;
[V, X0, ~, width] = mps_get_lines(lambda, width);

% Get projection matrix.
O = null(V.').';

% Project lines.
Xk = O * X0;

% Generate wrapped-around phase vectors.
tau = 2 * pi;
[k, ~, idx] = mps_get_wrapped_tuples(lambda, width);
X1 = - tau * double(k.');
Xw = O * X1;

M = size(Xw, 2);

% Create last output.
Xv = cell(M, 3);
for i = 1 : M
    Xv{i, 1} = Xw(:, i);
    Xv{i, 2} = k(i, :);
    Xv{i, 3} = idx(i);
end