function I = mp_get_illuminated(A_col, thr_col, A_row, thr_row)
% MP_GET_ILLUMINATED Returns logical array indicating illuminated pixels.
%   I = MP_GET_ILLUMINATED(A_col, thr_col, A_row, thr_row) returns and
%   image I which indicates illuminated pixels.
%
%   Inputs A_col and A_row are amplitudes of MPS signals for colum and row
%   code respectively; each input may be an image or a cell array
%   containing multiple images. Inputs thr_col and thr_row are thresholds
%   to apply the the amplitude images; if A_col and A_row are single images
%   then thr_col and thr_row are simple numbers and if A_col and A_row are
%   cell arrays then thr_col and thr_row must also be cell arrays.
%
%   Pixels in the output mask I are set to true only if all input
%   amplitudes are above given thresholds.
%
%   I = MP_GET_ILLUMINATED(A, thr) does the same for a single input only.
%
%   See also MP_WRAPPED_PHASE, MPS_WRAPPED_PHASE.

% $Revision: 1.0 $  $Date: 2017/06/28 $
% $Author(s): Tomislav Petkovic $

narginchk(2, 4);
nargoutchk(1, 1);

% Start with empty output.
I = [];

% Accumulate inputs.
I = update_mask(A_col, thr_col, I);
if 2 < nargin; I = update_mask(A_row, thr_row, I); end;

end


function I = update_mask(X, thr, I)
% UPDATE_MASK Updates input mask I depending on comparison result.
%   I = UPDATE_MASK(X, thr, I) returns updates mask I which is

% Initialize mask or get its size.
if isempty(I)
    I = true;
    sz = [];
else
    sz = size(I);
end

if iscell(X)
    
    assert( iscell(thr) );
    assert( numel(X) == numel(thr) );
    
    % Sequentially update mask by processing each item from a cell array.
    for i = 1 : numel(X)
        if ~isempty(X{i})
            if ~isempty(sz)
                assert( all(size(X{i}) == sz) );
            else
                sz = size(X{i});
            end
            assert( isnumeric(thr{i}) && (1 == numel(thr{i})) );
            I = I & (X{i} > thr{i});
        else
            assert( isempty(thr{i}) );
        end
    end
    
elseif isnumeric(X)
    
    if ~isempty(X)
        if ~isempty(sz); assert( all(size(X) == sz) ); end;
        assert( isnumeric(thr) && (1 == numel(thr)) );
        I = I & (X > thr);
    else
        assert( isempty(thr) );
    end
    
else
    
    % Abort execution.
    error('Unsupported input datatype.');
    
end

end
