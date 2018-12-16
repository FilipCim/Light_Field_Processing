function [I, I_type] = calibration_load_images_as_raw(L)
% CALIBRATION_LOAD_IMAGES_AS_RAW Loads specified images as raw.
%   [I, I_type] = CALIBRATION_LOAD_IMAGES_AS_RAW(L) returns a cell array I
%   containing raw image data and a cell array I_type which indicates raw
%   image type. The input L is a cell array containing a list of filenames.
%   Input images may be in RAW format or in any image format that is
%   recognized by IMREAD.
%
%   See also PS_READ_RAW.

% $Revision: 1.0 $  $Date: 2017/07/20 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 1);
nargoutchk(2, 2);

assert( iscell(L) && ~isempty(L) );

% Preallocate file list.
I = cell(size(L));
I_type = cell(size(L));

% Cycle through L and load all images.
for i = 1 : numel(L)
    [I{i}, I_type{i}] = calibration_load_image_as_raw(L{i});
end