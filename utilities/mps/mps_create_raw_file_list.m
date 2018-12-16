function [L, n, p, prj, ID] = mps_create_raw_file_list(D, tag)
% MPS_CREATE_RAW_FILE_LIST Creates a list of all raw images in directory.
%   L = MPS_CREATE_RAW_FILE_LIST(D) returns a cell array L of all raw
%   images found in a directory D that have a corresponding XML metadata
%   file. Elements of L are absolute file names (including the path).
%
%   [L, n, p] = MPS_CREATE_RAW_FILE_LIST(D) also parses filenames
%   and tries to extract values of n and p from the filenames, if they are
%   present. All are output as cell arrays of same size as the file list L.
%   An empty matrix indicates that value was not present in the filename.
%
%   The expected format of a filename that encodes n and p values is:
%
%       image_n?_p?.raw
%
%   where ? denotes a decimal number. Note that: a) instead of image any
%   word may be used, b) the order of appearance of _n? and _p? tokens in a
%   filename may change, c) token parsing is case insensitive, and d)
%   value will be parsed only if it appears exactly once. If any of n or p
%   values cannot be parsed from the filename than the corresponding cell
%   array element will be an empty matrix.
%
%   Function recognizes raw image by their extension. A raw image must be
%   accompanied by an xml file of the same name that contains image
%   metadata. Note that this function does not parse the xml file nor
%   checks if the raw file contains any image data.
%
%   [L, n, p, prj, ID] = MPS_CREATE_RAW_FILE_LIST(D) additionally pareses
%   filenames for the values of projector name prj and file index ID.
%
%   For multi-projector setup an image filename encodes additional values:
%
%       prj_?_image_ID?_n?_p?.raw
%
%   where the first ? denotes any combination of uppercase letters and
%   numbers while other ? denote decimal numbers. Token _ID? is parsed in
%   the same way as _n? and _p? tokens whiel token prj_?_ must appear at
%   the beginning of the filename.
%
%   L = MPS_CREATE_RAW_FILE_LIST(D, tag) returns only items whose filenames
%   contains the specified tag, e.g. set tag to 'row' to retrieve only
%   filenames which contain 'row'. If omitted or empty then no filtering is
%   performed.
%
%   See also MPS_CREATE_FILE_LIST, CALIBRATION_LOAD_IMAGES,
%   PS_CREATE_RAW_DATA_MATRIX, PS_CREATE_DATA_MATRIX.

% $Revision: 1.4 $  $Date: 2018/05/21 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);
nargoutchk(1, 5);

assert( (0 < numel(D)) && (true == ischar(D)) );

% Get directory content.
files = dir(D);
if 0 == numel(files)
    error('Directory not found.')
else
    cur_D = pwd;
    cd(D);
    pathstr = pwd;
    cd(cur_D);
end

% Valid image extensions recognized by imread.
extensions_data = {'.raw'};
extensions_metadata = {'.raw'};

% Parse filenames for multiple-phase shift parameters.
parse_parameters = (1 < nargout);

% Filter filenames if input tag is specified.
filter_filenames = (1 < nargin) && ~isempty(tag);
if filter_filenames
    assert( ischar(tag) && (0 < numel(tag)) );
end

% Parse found files and retain images.
L = {};
if parse_parameters
    n = {};
    p = {};
    prj = {};
    ID = {};
end
count = 0;
for i = 1 : numel(files)
    
    % Get file extension.
    [~, name, ext] = fileparts(files(i).name);
    if isempty(ext); continue; end;
    
    % Check file extension againts known image extensions.
    is_image = false;
    has_metadata = false;
    for j = 1 : numel(extensions_data)
        if strcmpi(ext, extensions_data{j})
            is_image = true;
            
            % Check if found image has metadata file.
            for k = 1 : numel(extensions_metadata)
                if exist([pathstr filesep name extensions_metadata{k}], 'file') == 2
                    has_metadata = true;
                    break;
                end
            end
            
            break;
        end
    end
    
    % Check if input filename matches the requested tag.
    tag_matches = true;
    if filter_filenames
        tag_matches = ~isempty( strfind(name, tag) );
    end
    
    % Add file to output list and parse it to find parametes.
    if is_image && has_metadata && tag_matches
        filename = [pathstr filesep files(i).name];
        assert( 2 == exist(filename, 'file') );
        
        count = count + 1;
        L{count} = filename;
        
        if parse_parameters
            [startIndex, endIndex] = regexp(filename, '_[nN][0-9]+(\.[0-9]*)?');
            if (1 == numel(startIndex)) && (1 == numel(endIndex))
                n{count} = str2double(filename(startIndex+2 : endIndex));
            else
                n{count} = [];
            end
            
            [startIndex, endIndex] = regexp(filename, '_[pP][0-9]+(\.[0-9]*)?');
            if (1 == numel(startIndex)) && (1 == numel(endIndex))
                p{count} = str2double(filename(startIndex+2 : endIndex));
            else
                p{count} = [];
            end
            
            [startIndex, endIndex] = regexp(files(i).name, '^prj_[A-Z0-9]_');
            if (1 == numel(startIndex)) && (1 == numel(endIndex))
                prj{count} = files(i).name(startIndex+4 : endIndex-1);
            else
                prj{count} = '';
            end
            
            [startIndex, endIndex] = regexp(filename, '_[iI][dD][0-9]+(\.[0-9]*)?');
            if (1 == numel(startIndex)) && (1 == numel(endIndex))
                ID{count} = str2double(filename(startIndex+3 : endIndex));
            else
                ID{count} = [];
            end
        end
    end
    
end