function [raw, raw_format, buffer, buffer_format, buffer_stride, width, height] = ps_read_raw(filename)
% PS_READ_RAW Reads raw image data.
%   [raw, raw_format, buffer, buffer_format, buffer_stride, width, height] = PS_READ_RAW(filename)
%   reads raw image data and tries to interpret it.
%
%   A raw image is comprised of two files, one with the extension .raw that
%   stores the raw image buffer data and one with the extension .xml that
%   stores the image metadata. Name of any of these two may be given as
%   input.
%
%   Outputs are as follows:
%   raw is an image matrix with the precision not smaller than the input
%   data that contains decoded image data, or is an empty matrix if the raw
%   image data cannot be interpreted.
%   raw_format is a string describing the image data pixel format, see
%   BatchAcquisition.h header for details of all available formats.
%   buffer contains the raw image data (bytes) stored as uint8.
%   buffer_format is a string describing the image data pixel format, see
%   BatchAcquisition.h header for details of all available formats.
%   buffer_stride is length of one image row in bytes for the matrix buffer.
%   width is the image width in pixels.
%   height is the image height in pixels.
%
%   If input raw image data is Bayer then demosaicing is NOT performed. To
%   decode the Bayer pattern use function DEMOSAIC from the Image
%   Processing Toolbox.
%
%   If input raw image daty type contains RGB(A)/BGR(A) values in any order
%   then the output raw is always converted to a planar RGB image.
%
%   Function always corrects input endianness to the machine endianness
%   therefore the output string raw_format will never be a format which
%   ends with _BigEndian. There is no need to correct the data endianness.
%
%   See also PS_CREATE_RAW_FILE_LIST, PS_CREATE_RAW_DATA_MATRIX, DEMOSAIC.

% $Revision: 1.1 $  $Date: 2017/01/10 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 1);
nargoutchk(7, 7);

% Create default outputs.
raw = [];
raw_format = [];
buffer = [];
buffer_format = [];
buffer_stride = [];
width = [];
height = [];

% Decompose input filename.
[pathstr, name, ~] = fileparts(filename);

% We need to read the metadata stored in .xml file and the raw data stored
% in .raw file.
filename_xml = [pathstr filesep name '.xml'];
filename_raw = [pathstr filesep name '.raw'];

% Check if file exist.
xml_exist = (2 == exist(filename_xml, 'file'));
raw_exist = (2 == exist(filename_xml, 'file'));
if ~xml_exist; warning('Cannot find .xml metadata file.'); end;
if ~raw_exist; warning('Cannot find .raw image file.'); end;
if ~xml_exist || ~raw_exist; return; end;

% Read and parse .xml file.
xDoc = xmlread(filename_xml);

% Get image metadata.
buffer_size = xDoc.getDocumentElement.getElementsByTagName('BufferSize').item(0).getTextContent;
buffer_format = xDoc.getDocumentElement.getElementsByTagName('PixelFormat').item(0).getTextContent;
width = xDoc.getDocumentElement.getElementsByTagName('Width').item(0).getTextContent;
height = xDoc.getDocumentElement.getElementsByTagName('Height').item(0).getTextContent;
buffer_stride = xDoc.getDocumentElement.getElementsByTagName('Stride').item(0).getTextContent;

% Convert to Matlab format.
buffer_size = str2double(buffer_size);
buffer_format = char(buffer_format);
width = str2double(width);
height = str2double(height);
buffer_stride = str2double(buffer_stride);

% Load RAW data into uint8 buffer. Note that this buffer is row-major.
FID = fopen(filename_raw);
assert(-1 ~= FID);
if (-1 == FID)
    return;
else
    buffer = fread(FID, Inf, 'uint8=>uint8');
    fclose(FID);
    FID = -1;
end
assert(numel(buffer) == buffer_size);

% NOTES ON DATA CONVERSION
%
% Note that MATLAB uses column-major storage while the buffer is row-major
% so rows and columns must be swapped during or after the copy operation.
% Profiling showed that copy operation into a temporary matrix followed by
% transpose operation is faster than direct row-colum swap in the body of
% the for loop which copies and decodes the data (probably due to better
% cache utilization and cache-efficient implementation of the transpose).
%
% All copy operations are done using uint8 datatype which is then
% converted to the desired datatype by using typecast function.
%
% Regarding machine endianness to avoid code duplication default
% implementations are for little endian machines. For big endian machines
% if possible we change the input buffer format to one required for
% specific data conversion.

% Get machine endianness.
[~, ~, endian] = computer;
if 'B' == endian
    
    % This switch statement contains only conversions which cannot be
    % implemented by little endian conversion code; these mostly include
    % packed formats.
    switch buffer_format
        
        case 'IDT_12U_GRAY_Packed'
            img_8u = zeros(width, 2*height, 'uint8');
            idx_odd = int32(1:2:width);
            idx_even = int32(2:2:width);
            idx_1 = int32(1:3:3*width/2);
            idx_2 = int32(2:3:3*width/2);
            idx_3 = int32(3:3:3*width/2);
            for i = 1 : height
                offset = int32( (i-1) * buffer_stride );
                byte0 = buffer( idx_1 + offset );
                byte1 = buffer( idx_2 + offset );
                byte2 = buffer( idx_3 + offset );
                img_8u(idx_odd,  2*i  ) = bitand(byte1, 240);
                img_8u(idx_odd,  2*i-1) = byte0;
                img_8u(idx_even, 2*i  ) = bitshift(bitand(byte1, 15), 4);
                img_8u(idx_even, 2*i-1) = byte2;
            end
            img_8u = img_8u.';
            raw = typecast(img_8u(:), 'uint16');
            raw = reshape(raw, height, width);
            raw_format = 'IDT_16U_GRAY';
            return;
            
        case 'IDT_16U_GRAY'
            buffer_format = 'IDT_16U_GRAY_BigEndian';
            
        case 'IDT_16U_GRAY_BigEndian'
            buffer_format = 'IDT_16U_GRAY';
            
        case 'IDT_32U_GRAY'
            img_8u = zeros(width, 4*height, 'uint8');
            idx_1 = int32(1:4:4*width);
            idx_2 = int32(2:4:4*width);
            idx_3 = int32(3:4:4*width);
            idx_4 = int32(4:4:4*width);
            for i = 1 : height
                offset = int32( (i-1) * buffer_stride );
                img_8u(:, 4*i-3) = buffer( idx_4 + offset );
                img_8u(:, 4*i-2) = buffer( idx_3 + offset );
                img_8u(:, 4*i-1) = buffer( idx_2 + offset );
                img_8u(:, 4*i  ) = buffer( idx_1 + offset );
            end
            img_8u = img_8u.';
            raw = typecast(img_8u(:), 'uint32');
            raw = reshape(raw, height, width);
            raw_format = 'IDT_32U_GRAY';
            return;
            
        case 'IDT_16S_GRAY'
            buffer_format = 'IDT_16S_GRAY_BigEndian';
            
        case 'IDT_16S_GRAY_BigEndian'
            buffer_format = 'IDT_16S_GRAY';
            
        case 'IDT_32S_GRAY'
            img_8u = zeros(width, 4*height, 'uint8');
            idx_1 = int32(1:4:4*width);
            idx_2 = int32(2:4:4*width);
            idx_3 = int32(3:4:4*width);
            idx_4 = int32(4:4:4*width);
            for i = 1 : height
                offset = int32( (i-1) * buffer_stride );
                img_8u(:, 4*i-3) = buffer( idx_4 + offset );
                img_8u(:, 4*i-2) = buffer( idx_3 + offset );
                img_8u(:, 4*i-1) = buffer( idx_2 + offset );
                img_8u(:, 4*i  ) = buffer( idx_1 + offset );
            end
            img_8u = img_8u.';
            raw = typecast(img_8u(:), 'int32');
            raw = reshape(raw, height, width);
            raw_format = 'IDT_32S_GRAY';
            return;
            
        case 'IDT_12U_BayerGR_Packed'
            img_8u = zeros(width, 2*height, 'uint8');
            idx_odd = int32(1:2:width);
            idx_even = int32(2:2:width);
            idx_1 = int32(1:3:3*width/2);
            idx_2 = int32(2:3:3*width/2);
            idx_3 = int32(3:3:3*width/2);
            for i = 1 : height
                offset = int32( (i-1) * buffer_stride );
                byte0 = buffer( idx_1 + offset );
                byte1 = buffer( idx_2 + offset );
                byte2 = buffer( idx_3 + offset );
                img_8u(idx_odd,  2*i  ) = bitand(byte1, 240);
                img_8u(idx_odd,  2*i-1) = byte0;
                img_8u(idx_even, 2*i  ) = bitshift(bitand(byte1, 15), 4);
                img_8u(idx_even, 2*i-1) = byte2;
            end
            img_8u = img_8u.';
            raw = typecast(img_8u(:), 'uint16');
            raw = reshape(raw, height, width);
            raw_format = 'IDT_16U_BayerGR';
            return;
            
        case 'IDT_12U_BayerRG_Packed'
            img_8u = zeros(width, 2*height, 'uint8');
            idx_odd = int32(1:2:width);
            idx_even = int32(2:2:width);
            idx_1 = int32(1:3:3*width/2);
            idx_2 = int32(2:3:3*width/2);
            idx_3 = int32(3:3:3*width/2);
            for i = 1 : height
                offset = int32( (i-1) * buffer_stride );
                byte0 = buffer( idx_1 + offset );
                byte1 = buffer( idx_2 + offset );
                byte2 = buffer( idx_3 + offset );
                img_8u(idx_odd,  2*i  ) = bitand(byte1, 240);
                img_8u(idx_odd,  2*i-1) = byte0;
                img_8u(idx_even, 2*i  ) = bitshift(bitand(byte1, 15), 4);
                img_8u(idx_even, 2*i-1) = byte2;
            end
            img_8u = img_8u.';
            raw = typecast(img_8u(:), 'uint16');
            raw = reshape(raw, height, width);
            raw_format = 'IDT_16U_BayerRG';
            return;
            
        case 'IDT_12U_BayerGB_Packed'
            img_8u = zeros(width, 2*height, 'uint8');
            idx_odd = int32(1:2:width);
            idx_even = int32(2:2:width);
            idx_1 = int32(1:3:3*width/2);
            idx_2 = int32(2:3:3*width/2);
            idx_3 = int32(3:3:3*width/2);
            for i = 1 : height
                offset = int32( (i-1) * buffer_stride );
                byte0 = buffer( idx_1 + offset );
                byte1 = buffer( idx_2 + offset );
                byte2 = buffer( idx_3 + offset );
                img_8u(idx_odd,  2*i  ) = bitand(byte1, 240);
                img_8u(idx_odd,  2*i-1) = byte0;
                img_8u(idx_even, 2*i  ) = bitshift(bitand(byte1, 15), 4);
                img_8u(idx_even, 2*i-1) = byte2;
            end
            img_8u = img_8u.';
            raw = typecast(img_8u(:), 'uint16');
            raw = reshape(raw, height, width);
            raw_format = 'IDT_16U_BayerGB';
            return;
            
        case 'IDT_12U_BayerBG_Packed'
            img_8u = zeros(width, 2*height, 'uint8');
            idx_odd = int32(1:2:width);
            idx_even = int32(2:2:width);
            idx_1 = int32(1:3:3*width/2);
            idx_2 = int32(2:3:3*width/2);
            idx_3 = int32(3:3:3*width/2);
            for i = 1 : height
                offset = int32( (i-1) * buffer_stride );
                byte0 = buffer( idx_1 + offset );
                byte1 = buffer( idx_2 + offset );
                byte2 = buffer( idx_3 + offset );
                img_8u(idx_odd,  2*i  ) = bitand(byte1, 240);
                img_8u(idx_odd,  2*i-1) = byte0;
                img_8u(idx_even, 2*i  ) = bitshift(bitand(byte1, 15), 4);
                img_8u(idx_even, 2*i-1) = byte2;
            end
            img_8u = img_8u.';
            raw = typecast(img_8u(:), 'uint16');
            raw = reshape(raw, height, width);
            raw_format = 'IDT_16U_BayerBG';
            return;
            
        case 'IDT_16U_BayerGR'
            buffer_format = 'IDT_16U_BayerGR_BigEndian';
            
        case 'IDT_16U_BayerRG'
            buffer_format = 'IDT_16U_BayerRG_BigEndian';
            
        case 'IDT_16U_BayerGB'
            buffer_format = 'IDT_16U_BayerGB_BigEndian';
            
        case 'IDT_16U_BayerBG'
            buffer_format = 'IDT_16U_BayerBG_BigEndian';
            
        case 'IDT_16U_BayerGR_BigEndian'
            buffer_format = 'IDT_16U_BayerGR';
            
        case 'IDT_16U_BayerRG_BigEndian'
            buffer_format = 'IDT_16U_BayerRG';
            
        case 'IDT_16U_BayerGB_BigEndian'
            buffer_format = 'IDT_16U_BayerGB';
            
        case 'IDT_16U_BayerBG_BigEndian'
            buffer_format = 'IDT_16U_BayerBG';
            
    end
    
elseif 'L' == endian
    % Nothing to do.
else
    error('Unknown machine endianness!');
end

% Convert data for little endian machine.
switch buffer_format
    
    case 'IDT_UNKNOWN'
        error('Unknown image format!');
        
    case 'IDT_8U_BINARY'
        raw = zeros(width, height, 'uint8');
        idx = int32(1:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            raw(:, i) = buffer( idx + offset );
        end
        raw = raw.';
        raw_format = 'IDT_8U_BINARY';
        return;
        
    case 'IDT_8U_GRAY'
        raw = zeros(width, height, 'uint8');
        idx = int32(1:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            raw(:, i) = buffer( idx + offset );
        end
        raw = raw.';
        raw_format = 'IDT_8U_GRAY';
        return;
        
    case 'IDT_12U_GRAY_Packed'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        idx_1 = int32(1:3:3*width/2);
        idx_2 = int32(2:3:3*width/2);
        idx_3 = int32(3:3:3*width/2);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            byte0 = buffer( idx_1 + offset );
            byte1 = buffer( idx_2 + offset );
            byte2 = buffer( idx_3 + offset );
            img_8u(idx_odd,  2*i-1) = bitand(byte1, 240);
            img_8u(idx_odd,  2*i  ) = byte0;
            img_8u(idx_even, 2*i-1) = bitshift(bitand(byte1, 15), 4);
            img_8u(idx_even, 2*i  ) = byte2;
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_GRAY';
        return;
        
    case 'IDT_16U_GRAY'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_odd + offset );
            img_8u(:, 2*i  ) = buffer( idx_even + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_GRAY';
        return;
        
    case 'IDT_16U_GRAY_BigEndian'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_even + offset );
            img_8u(:, 2*i  ) = buffer( idx_odd + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_GRAY';
        return;
        
    case 'IDT_32U_GRAY'
        img_8u = zeros(width, 4*height, 'uint8');
        idx_1 = int32(1:4:4*width);
        idx_2 = int32(2:4:4*width);
        idx_3 = int32(3:4:4*width);
        idx_4 = int32(4:4:4*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 4*i-3) = buffer( idx_1 + offset );
            img_8u(:, 4*i-2) = buffer( idx_2 + offset );
            img_8u(:, 4*i-1) = buffer( idx_3 + offset );
            img_8u(:, 4*i  ) = buffer( idx_4 + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint32');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_32U_GRAY';
        return;    
        
    case 'IDT_8S_GRAY'
        img_8u = zeros(width, height, 'uint8');
        idx = int32(1:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, i) = buffer( idx + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'int8');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_8S_GRAY';
        return;
        
    case 'IDT_16S_GRAY'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_odd + offset );
            img_8u(:, 2*i  ) = buffer( idx_even + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'int16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16S_GRAY';
        return;
        
    case 'IDT_16S_GRAY_BigEndian'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_even + offset );
            img_8u(:, 2*i  ) = buffer( idx_odd + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'int16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16S_GRAY';
        return;
        
    case 'IDT_32S_GRAY'
        img_8u = zeros(width, 4*height, 'uint8');
        idx_1 = int32(1:4:4*width);
        idx_2 = int32(2:4:4*width);
        idx_3 = int32(3:4:4*width);
        idx_4 = int32(4:4:4*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 4*i-3) = buffer( idx_1 + offset );
            img_8u(:, 4*i-2) = buffer( idx_2 + offset );
            img_8u(:, 4*i-1) = buffer( idx_3 + offset );
            img_8u(:, 4*i  ) = buffer( idx_4 + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'int32');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_32S_GRAY';
        return;
        
    case 'IDT_32F_GRAY'
        img_8u = zeros(width, 4*height, 'uint8');
        idx_1 = int32(1:4:4*width);
        idx_2 = int32(2:4:4*width);
        idx_3 = int32(3:4:4*width);
        idx_4 = int32(4:4:4*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 4*i-3) = buffer( idx_1 + offset );
            img_8u(:, 4*i-2) = buffer( idx_2 + offset );
            img_8u(:, 4*i-1) = buffer( idx_3 + offset );
            img_8u(:, 4*i  ) = buffer( idx_4 + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'single');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_32F_GRAY';
        return;
        
    case 'IDT_64F_GRAY'
        img_8u = zeros(width, 8*height, 'uint8');
        idx_1 = int32(1:8:8*width);
        idx_2 = int32(2:8:8*width);
        idx_3 = int32(3:8:8*width);
        idx_4 = int32(4:8:8*width);
        idx_5 = int32(5:8:8*width);
        idx_6 = int32(6:8:8*width);
        idx_7 = int32(7:8:8*width);
        idx_8 = int32(8:8:8*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 8*i-7) = buffer( idx_1 + offset );
            img_8u(:, 8*i-6) = buffer( idx_2 + offset );
            img_8u(:, 8*i-5) = buffer( idx_3 + offset );
            img_8u(:, 8*i-4) = buffer( idx_4 + offset );
            img_8u(:, 8*i-3) = buffer( idx_5 + offset );
            img_8u(:, 8*i-2) = buffer( idx_6 + offset );
            img_8u(:, 8*i-1) = buffer( idx_7 + offset );
            img_8u(:, 8*i  ) = buffer( idx_8 + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'double');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_64F_GRAY';
        return;
        
    case 'IDT_8U_BayerGR'
        img_8u = zeros(width, height, 'uint8');
        idx = int32(1:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, i) = buffer( idx + offset );
        end
        raw = img_8u.';
        raw_format = 'IDT_8U_BayerGR';
        return;
        
    case 'IDT_8U_BayerRG'
        img_8u = zeros(width, height, 'uint8');
        idx = int32(1:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, i) = buffer( idx + offset );
        end
        raw = img_8u.';
        raw_format = 'IDT_8U_BayerRG';
        return;
        
    case 'IDT_8U_BayerGB'
        img_8u = zeros(width, height, 'uint8');
        idx = int32(1:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, i) = buffer( idx + offset );
        end
        raw = img_8u.';
        raw_format = 'IDT_8U_BayerGB';
        return;
        
    case 'IDT_8U_BayerBG'
        img_8u = zeros(width, height, 'uint8');
        idx = int32(1:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, i) = buffer( idx + offset );
        end
        raw = img_8u.';
        raw_format = 'IDT_8U_BayerBG';
        return;
        
    case 'IDT_12U_BayerGR_Packed'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        idx_1 = int32(1:3:3*width/2);
        idx_2 = int32(2:3:3*width/2);
        idx_3 = int32(3:3:3*width/2);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            byte0 = buffer( idx_1 + offset );
            byte1 = buffer( idx_2 + offset );
            byte2 = buffer( idx_3 + offset );
            img_8u(idx_odd,  2*i-1) = bitand(byte1, 240);
            img_8u(idx_odd,  2*i  ) = byte0;
            img_8u(idx_even, 2*i-1) = bitshift(bitand(byte1, 15), 4);
            img_8u(idx_even, 2*i  ) = byte2;
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerGR';
        return;
        
    case 'IDT_12U_BayerRG_Packed'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        idx_1 = int32(1:3:3*width/2);
        idx_2 = int32(2:3:3*width/2);
        idx_3 = int32(3:3:3*width/2);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            byte0 = buffer( idx_1 + offset );
            byte1 = buffer( idx_2 + offset );
            byte2 = buffer( idx_3 + offset );
            img_8u(idx_odd,  2*i-1) = bitand(byte1, 240);
            img_8u(idx_odd,  2*i  ) = byte0;
            img_8u(idx_even, 2*i-1) = bitshift(bitand(byte1, 15), 4);
            img_8u(idx_even, 2*i  ) = byte2;
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerRG';
        return;
        
    case 'IDT_12U_BayerGB_Packed'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        idx_1 = int32(1:3:3*width/2);
        idx_2 = int32(2:3:3*width/2);
        idx_3 = int32(3:3:3*width/2);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            byte0 = buffer( idx_1 + offset );
            byte1 = buffer( idx_2 + offset );
            byte2 = buffer( idx_3 + offset );
            img_8u(idx_odd,  2*i-1) = bitand(byte1, 240);
            img_8u(idx_odd,  2*i  ) = byte0;
            img_8u(idx_even, 2*i-1) = bitshift(bitand(byte1, 15), 4);
            img_8u(idx_even, 2*i  ) = byte2;
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerGB';
        return;
        
    case 'IDT_12U_BayerBG_Packed'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        idx_1 = int32(1:3:3*width/2);
        idx_2 = int32(2:3:3*width/2);
        idx_3 = int32(3:3:3*width/2);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            byte0 = buffer( idx_1 + offset );
            byte1 = buffer( idx_2 + offset );
            byte2 = buffer( idx_3 + offset );
            img_8u(idx_odd,  2*i-1) = bitand(byte1, 240);
            img_8u(idx_odd,  2*i  ) = byte0;
            img_8u(idx_even, 2*i-1) = bitshift(bitand(byte1, 15), 4);
            img_8u(idx_even, 2*i  ) = byte2;
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerBG';
        return;
        
    case 'IDT_16U_BayerGR'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:2*width);
        idx_even = int32(2:2:2*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_odd + offset );
            img_8u(:, 2*i  ) = buffer( idx_even + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerGR';
        return;
        
    case 'IDT_16U_BayerRG'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:2*width);
        idx_even = int32(2:2:2*width);
        for i = 1 : height
            offset = int32((i-1) * buffer_stride);
            img_8u(:, 2*i-1) = buffer( idx_odd + offset );
            img_8u(:, 2*i  ) = buffer( idx_even + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerRG';
        return;
        
    case 'IDT_16U_BayerGB'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:2*width);
        idx_even = int32(2:2:2*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_odd + offset );
            img_8u(:, 2*i  ) = buffer( idx_even + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerGB';
        return;
        
    case 'IDT_16U_BayerBG'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:2*width);
        idx_even = int32(2:2:2*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_odd + offset );
            img_8u(:, 2*i  ) = buffer( idx_even + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerBG';
        return;
        
    case 'IDT_16U_BayerGR_BigEndian'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:2*width);
        idx_even = int32(2:2:2*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_even + offset );
            img_8u(:, 2*i  ) = buffer( idx_odd + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerGR';
        return;
        
    case 'IDT_16U_BayerRG_BigEndian'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:2*width);
        idx_even = int32(2:2:2*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_even + offset );
            img_8u(:, 2*i  ) = buffer( idx_odd + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerRG';
        return;
        
    case 'IDT_16U_BayerGB_BigEndian'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:2*width);
        idx_even = int32(2:2:2*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_even + offset );
            img_8u(:, 2*i  ) = buffer( idx_odd + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerGB';
        return;
        
    case 'IDT_16U_BayerBG_BigEndian'
        img_8u = zeros(width, 2*height, 'uint8');
        idx_odd = int32(1:2:2*width);
        idx_even = int32(2:2:2*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 2*i-1) = buffer( idx_even + offset );
            img_8u(:, 2*i  ) = buffer( idx_odd + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'uint16');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_16U_BayerBG';
        return;
        
    case 'IDT_64F_BayerRG'
        img_8u = zeros(width, 8*height, 'uint8');
        idx_1 = int32(1:8:8*width);
        idx_2 = int32(2:8:8*width);
        idx_3 = int32(3:8:8*width);
        idx_4 = int32(4:8:8*width);
        idx_5 = int32(5:8:8*width);
        idx_6 = int32(6:8:8*width);
        idx_7 = int32(7:8:8*width);
        idx_8 = int32(8:8:8*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, 8*i-7) = buffer( idx_1 + offset );
            img_8u(:, 8*i-6) = buffer( idx_2 + offset );
            img_8u(:, 8*i-5) = buffer( idx_3 + offset );
            img_8u(:, 8*i-4) = buffer( idx_4 + offset );
            img_8u(:, 8*i-3) = buffer( idx_5 + offset );
            img_8u(:, 8*i-2) = buffer( idx_6 + offset );
            img_8u(:, 8*i-1) = buffer( idx_7 + offset );
            img_8u(:, 8*i  ) = buffer( idx_8 + offset );
        end
        img_8u = img_8u.';
        raw = typecast(img_8u(:), 'double');
        raw = reshape(raw, height, width);
        raw_format = 'IDT_64F_BayerRG';
        return;
        
    case 'IDT_8U_RGB'
        img_8u = zeros(width, height, 3, 'uint8');
        idx_1 = int32(1:3:3*width);
        idx_2 = int32(2:3:3*width);
        idx_3 = int32(3:3:3*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, i, 1) = buffer( idx_1 + offset );
            img_8u(:, i, 2) = buffer( idx_2 + offset );
            img_8u(:, i, 3) = buffer( idx_3 + offset );
        end
        raw = zeros(height, width, 3, 'uint8');
        raw(:, :, 1) = img_8u(:, :, 1).';
        raw(:, :, 2) = img_8u(:, :, 2).';
        raw(:, :, 3) = img_8u(:, :, 3).';
        raw_format = 'IDT_8U_RGB_Planar';
        return;
        
    case 'IDT_8U_RGB_Planar'
        error('Converter for planar RGB is not implemented!');
        
    case 'IDT_8U_RGBA'
        img_8u = zeros(width, height, 3, 'uint8');
        idx_1 = int32(1:4:4*width);
        idx_2 = int32(2:4:4*width);
        idx_3 = int32(3:4:4*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, i, 1) = buffer( idx_1 + offset );
            img_8u(:, i, 2) = buffer( idx_2 + offset );
            img_8u(:, i, 3) = buffer( idx_3 + offset );
        end
        raw = zeros(height, width, 3, 'uint8');
        raw(:, :, 1) = img_8u(:, :, 1).';
        raw(:, :, 2) = img_8u(:, :, 2).';
        raw(:, :, 3) = img_8u(:, :, 3).';
        raw_format = 'IDT_8U_RGB_Planar';
        return;
        
    case 'IDT_8U_BGR'
        img_8u = zeros(width, height, 3, 'uint8');
        idx_1 = int32(1:3:3*width);
        idx_2 = int32(2:3:3*width);
        idx_3 = int32(3:3:3*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, i, 1) = buffer( idx_3 + offset );
            img_8u(:, i, 2) = buffer( idx_2 + offset );
            img_8u(:, i, 3) = buffer( idx_1 + offset );
        end
        raw = zeros(height, width, 3, 'uint8');
        raw(:, :, 1) = img_8u(:, :, 1).';
        raw(:, :, 2) = img_8u(:, :, 2).';
        raw(:, :, 3) = img_8u(:, :, 3).';
        raw_format = 'IDT_8U_RGB_Planar';
        return;
        
    case 'IDT_8U_BGRA'
        img_8u = zeros(width, height, 3, 'uint8');
        idx_1 = int32(1:4:4*width);
        idx_2 = int32(2:4:4*width);
        idx_3 = int32(3:4:4*width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            img_8u(:, i, 1) = buffer( idx_3 + offset );
            img_8u(:, i, 2) = buffer( idx_2 + offset );
            img_8u(:, i, 3) = buffer( idx_1 + offset );
        end
        raw = zeros(height, width, 3, 'uint8');
        raw(:, :, 1) = img_8u(:, :, 1).';
        raw(:, :, 2) = img_8u(:, :, 2).';
        raw(:, :, 3) = img_8u(:, :, 3).';
        raw_format = 'IDT_8U_RGB_Planar';
        return;
        
    case 'IDT_8U_YUV411'
        error('Converter for YUV411 is not implemented!');
        
    case 'IDT_8U_YUV422'
        img_8u = zeros(width, height, 3, 'uint8');
        idx_1 = int32(1:4:2*width);
        idx_2 = int32(2:4:2*width);
        idx_3 = int32(3:4:2*width);
        idx_4 = int32(4:4:2*width);
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            U  = double( buffer( idx_1 + offset ) ) - 128.0;
            Y1 = double( buffer( idx_2 + offset ) );
            V  = double( buffer( idx_3 + offset ) ) - 128.0;
            Y2 = double( buffer( idx_4 + offset ) );
            
            UVB = 1.77200 * U;
            UVG = 0.34414 * U + 0.71414 * V;
            UVR =               1.40200 * V;
            
            img_8u(idx_odd, i, 1) = uint8( Y1 + UVR );
            img_8u(idx_odd, i, 2) = uint8( Y1 - UVG );
            img_8u(idx_odd, i, 3) = uint8( Y1 + UVB );
            
            img_8u(idx_even, i, 1) = uint8( Y2 + UVR );
            img_8u(idx_even, i, 2) = uint8( Y2 - UVG );
            img_8u(idx_even, i, 3) = uint8( Y2 + UVG );
        end
        raw = zeros(height, width, 3, 'uint8');
        raw(:, :, 1) = img_8u(:, :, 1).';
        raw(:, :, 2) = img_8u(:, :, 2).';
        raw(:, :, 3) = img_8u(:, :, 3).';
        raw_format = 'IDT_8U_RGB_Planar';
        return;
        
    case 'IDT_8U_YUV422_BT601'
        img_8u = zeros(width, height, 3, 'uint8');
        idx_1 = int32(1:4:2*width);
        idx_2 = int32(2:4:2*width);
        idx_3 = int32(3:4:2*width);
        idx_4 = int32(4:4:2*width);
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            U  =           double( buffer( idx_1 + offset ) ) - 128.0;
            Y1 = 1.16438 * double( buffer( idx_2 + offset ) ) -  16.0;
            V  =           double( buffer( idx_3 + offset ) ) - 128.0;
            Y2 = 1.16438 * double( buffer( idx_4 + offset ) ) -  16.0;
            
            UVB = 2.01723 * U;
            UVG = 0.39176 * U + 0.81297 * V;
            UVR =               1.59603 * V;
            
            img_8u(idx_odd, i, 1) = uint8( Y1 + UVR );
            img_8u(idx_odd, i, 2) = uint8( Y1 - UVG );
            img_8u(idx_odd, i, 3) = uint8( Y1 + UVB );
            
            img_8u(idx_even, i, 1) = uint8( Y2 + UVR );
            img_8u(idx_even, i, 2) = uint8( Y2 - UVG );
            img_8u(idx_even, i, 3) = uint8( Y2 + UVG );
        end
        raw = zeros(height, width, 3, 'uint8');
        raw(:, :, 1) = img_8u(:, :, 1).';
        raw(:, :, 2) = img_8u(:, :, 2).';
        raw(:, :, 3) = img_8u(:, :, 3).';
        raw_format = 'IDT_8U_RGB_Planar';
        return;
        
    case 'IDT_8U_YUV422_BT709'
        img_8u = zeros(width, height, 3, 'uint8');
        idx_1 = int32(1:4:2*width);
        idx_2 = int32(2:4:2*width);
        idx_3 = int32(3:4:2*width);
        idx_4 = int32(4:4:2*width);
        idx_odd = int32(1:2:width);
        idx_even = int32(2:2:width);
        for i = 1 : height
            offset = int32( (i-1) * buffer_stride );
            U  =           double( buffer( idx_1 + offset ) ) - 128.0;
            Y1 = 1.16438 * double( buffer( idx_2 + offset ) ) -  16.0;
            V  =           double( buffer( idx_3 + offset ) ) - 128.0;
            Y2 = 1.16438 * double( buffer( idx_4 + offset ) ) -  16.0;
            
            UVB = 2.11240 * U;
            UVG = 0.21325 * U + 0.53291 * V;
            UVR =               1.79274 * V;
            
            img_8u(idx_odd, i, 1) = uint8( Y1 + UVR );
            img_8u(idx_odd, i, 2) = uint8( Y1 - UVG );
            img_8u(idx_odd, i, 3) = uint8( Y1 + UVB );
            
            img_8u(idx_even, i, 1) = uint8( Y2 + UVR );
            img_8u(idx_even, i, 2) = uint8( Y2 - UVG );
            img_8u(idx_even, i, 3) = uint8( Y2 + UVG );
        end
        raw = zeros(height, width, 3, 'uint8');
        raw(:, :, 1) = img_8u(:, :, 1).';
        raw(:, :, 2) = img_8u(:, :, 2).';
        raw(:, :, 3) = img_8u(:, :, 3).';
        raw_format = 'IDT_8U_RGB_Planar';
        return;
        
    case 'IDT_8U_YUV444'
        error('Converter for YUV444 is not implemented!');
        
    case 'IDT_8U_UYV444'
        error('Converter for YUV444 is not implemented!');
        
    otherwise
        error('Unknown image format!');
end