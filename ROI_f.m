function [roi] = ROI_f (img, debug)

% Finding ROI in image - finds a ROI that is calculated form binary picture
% In this instance, white projected picture has been used to estimate the
% surface of the screen - ROI is then the maximum surface area that is
% above the preset threshold

% After finding ROI, it concatenates it into a TIFF image of the ROI
% regions on images

% variable debug for displaying images, 1 displays, 0 doesn't
% =========================================================================
    if debug == 1
        figure()
        imshow(img);
        title('Original Image');
    end
     
    [counts,x] = imhist(img,1024); 
    T = otsuthresh(counts); 
    imgBin = imbinarize(img, 0.75*T);
        
    % find area of all contours and corresponding bboxes; define the ROI as the 
    % bbox of the largest area contour 
    s = regionprops(imgBin, 'area', 'boundingBox'); 
    [~, maxAreaIdx] = max(cell2mat({s.Area}));
    roi = s(maxAreaIdx).BoundingBox; 
    
    if debug == 1
        figure()
        imshow(imgBin);
        hold on;
        rectangle('Position', [roi(1),roi(2),roi(3),roi(4)],'EdgeColor','r','LineWidth',2 )
        hold off;
    end
    end
