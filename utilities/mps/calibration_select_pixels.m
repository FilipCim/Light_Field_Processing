function [row, col, lut] = calibration_select_pixels(rng, thr, step_row, step_col, and_mask)
% CALIBRATION_SELECT_PIXELS Selects pixels and extracts their coordinates.
%   [row, col] = CALIBRATION_SELECT_PIXELS(rng, thr) returns two vectors
%   that hold coordinates of pixels in image rng that have amplitude larger
%   than thr.
%
%   Function behaves differently depending on the input rng:
%   1) If rng is a numeric matrix then all matrix elements are compared
%      against the threshold thr. In this case function outputs coordinates
%      of elements which are larger than thr.
%   2) If rng is a logical matrix then thr must be set to empty matrix.
%      In this case function outputs coordinates of elements of rng which
%      are set to true.
%   3) If rng is a cell array then function outputs coordinates of pixels
%      which satisfy the larger-than-threshold condition for all images in
%      the rng cell array.
%
%   [row, col, lut] = CALIBRATION_SELECT_PIXELS(rng, thr) also returns
%   reverse index stored in lut. lut is a matrix of the same size as rng
%   that at pixel (row, col) stores the index where that exact point is
%   located in the output vectors row and col. All points that are not in
%   the output vectors will have -1 as the value of reverse index lut.
%
%   [row, col, lut] = CALIBRATION_SELECT_PIXELS(rng, thr, step_row, step_col)
%   applies spatial decimation on the points in rng so every step_row row
%   and every step_col column are considered. Inputs step_row and step_col
%   must be whole positive numbers or empty matrices. If particular step
%   size is empty then decimation is not performed.
%
%   [row, col, lut] = CALIBRATION_SELECT_PIXELS(rng, thr, step_row, step_col, and_mask)
%   uses and_mask as additional mask for point selection. The input
%   and_mask must be a logical matrix of the same size as rng.
%
%   See also CALIBRATION_GET_ABSOLUTE_THRESHOLD.

% $Revision: 1.2 $  $Date: 2016/11/25 $
% $Author(s): Tomislav Petkovic $

narginchk(2, 5);
nargoutchk(2, 3);

% Define decimation factors.
if (3 > nargin) || isempty(step_row); step_row = 1; end;
assert( (1 == numel(step_row)) && (0 < step_row) && (round(step_row) == step_row) );
step_row = int32(step_row);

if (4 > nargin) || isempty(step_col); step_col = 1; end;
assert( (1 == numel(step_col)) && (0 < step_col) && (round(step_col) == step_col) );
step_col = int32(step_col);

% Threshold the input image and decimate.
if iscell(rng)
    N = numel(rng);
    size_rng = size(rng{1});
    mask = true(size_rng);
    assert( numel(mask) < 2^31 - 1 );
    
    if ~iscell(thr)
        assert( 1 == numel(thr) );
        thr_val = thr;
        thr = cell(size(rng));
        for i = 1 : N; thr{i} = thr_val; end;
    else
        assert( all(size(rng) == size(thr)) );
    end
    
    for i = 1 : N
        assert( ismatrix(rng{i}) && all(size(rng{i}) == size_rng) );
        assert( numel(rng{i}) < 2^31 - 1 );
        
        if islogical(rng{i})
            mask = mask & rng{i};
        else
            assert( isnumeric(rng{i}) && (1 == numel(thr{i})) );
            mask = mask & (rng{i} > thr{i});
        end
    end
    
else
    assert( ismatrix(rng) );
    size_rng = size(rng);
    assert( numel(rng) < 2^31 - 1 );
    
    if islogical(rng)
        assert( isempty(thr) );
        mask = rng;
    else
        assert( isnumeric(rng) && (1 == numel(thr)) );
        mask = (rng > thr);
    end
end

% Apply and mask (if any).
if (5 <= nargin) && ~isempty(and_mask)
    assert( islogical(and_mask) && all(size(and_mask) == size_rng) );
    mask = mask & and_mask;
end

% Decimate final result.
mask_decimated = mask(1:step_row:end, 1:step_col:end);

% Create output coordinates.
[row, col] = find(mask_decimated);

assert( numel(row) < 2^31 - 1 );
row = (int32(row) - 1) * step_row + 1;
col = (int32(col) - 1) * step_col + 1;

% Create reverse index look-up-table.
idx = row + (col - 1) * size(mask, 1);
lut = repmat(int32(-1), size_rng);
lut(idx) = 1 : numel(row);