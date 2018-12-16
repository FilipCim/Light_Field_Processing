function [synth_mask, synth_mask_number_of_ones, phi]=loadMeasurementMasks(directoryPath, fileExtension, varargin)
%loadMeasurementMasks - function loads synthetic masks from input directory, crops whole mask to
%only one submask and plots the loaded masks;
%it outputs submasks and number of ones in each of them, and also it
%outputs matrix phi

p=inputParser;

p.addRequired('directoryPath');
p.addRequired('fileExtension');
p.addParameter('Plot', 0, @isnumeric);

p.parse(directoryPath, fileExtension, varargin{:});

files = dir(fullfile(directoryPath , strcat('*', fileExtension)));

for fileNo = 1:length(files)
    fileName = files(fileNo).name;

    synth_mask{fileNo} = double(imread(strcat(directoryPath, fileName)))./256;
    synth_mask{fileNo} = double(synth_mask{fileNo} > 0);
    % this is used to check if used masks consist of 50% black and white
    % pixels
    synth_mask_number_of_ones(fileNo)=sum(synth_mask{fileNo}(:)); 
end

for fileNo=1:length(synth_mask)
    % create measurement matrix from different measurement masks
    measurement_vector{fileNo}=synth_mask{fileNo}(:);
    
    if(fileNo==1)
        phi=measurement_vector{1};
    else
        phi=[phi measurement_vector{fileNo}];
    end
    
    if(p.Results.Plot)
        % plot measurement masks
        figure(100), imagesc(reshape(phi(:,fileNo),[8 8])), colormap gray, title('Measurement mask'), drawnow
    end
end

phi=phi';

if(p.Results.Plot)
    % plot measurement matrix phi
    figure, imagesc(phi), colormap gray, title('Measurement matrix - phi')
end
