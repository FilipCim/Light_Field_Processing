
directoryPath = '/home/ivan/Development/Codes/SlideshowCapture/imagesProjector/hadamardMeasurementMaskImages';
hResolution = 1920;
vResolution = 1080;
plotBool = 0;

if(exist(directoryPath,'dir')~=7)
    mkdir(directoryPath)
end

cd(directoryPath)


% subdirectory for mask tiles
mkdir('Mask Tiles')
mkdir('Mask Tiles_8x8')

% cd 'Mask Tiles'


noOfPhases=4;
blockSize=8;
noOfMasks=64;

for maskNo=1:noOfMasks
    
    % generate binary mask with exactly 50% of ones
    %     randomMask{maskNo} = mod( reshape(randperm(mask_size*mask_size), mask_size, mask_size), 2 );
    
    %     randomMask{maskNo} = false(1,64) ;
    %     randomMask{maskNo}(1:32) = true ;
    %     randomMask{maskNo}(maskNo) = true ;
    
    %     randomMask{maskNo} = randomMask{maskNo}(randperm(numel(randomMask{maskNo})));
    %     randomMask{maskNo}=double(reshape(randomMask{maskNo},[8 8]));
    
    phi = walsh(blockSize^2);
    phi = phi(1:noOfMasks,:);
    
    randomMask{maskNo} = phi(maskNo,:);
    randomMask{maskNo}=double(reshape(randomMask{maskNo},[8 8]));
    
    cd 'Mask Tiles_8x8'
    fileNameString2=sprintf('maskTile_8x8_%02d.png',maskNo);
    imwrite(randomMask{maskNo}, fileNameString2);
    
    cd ..
    
    zerosMask=zeros(8,8);
    
    randomMask{maskNo}=imresize(randomMask{maskNo}, 2, 'nearest');
    zerosMask=imresize(zerosMask, 2, 'nearest');
    
    
    
    if(plotBool==1)
        % plot one mask
        figure(101)
        imagesc(randomMask{maskNo});
        drawnow
        colormap gray
        title('Single Mask Preview')
    end
    
    
    % generate mask with frequencyOfOnes*100% of white pixels in a mask
    %     randomMask{maskNo}=binornd(ones(mask_size,mask_size), frequencyOfOnes);
    
    
    % generate all white pixels mask
    %     randomMask{maskNo}=ones(mask_size,mask_size);
    
    
    for phaseNo=1:noOfPhases
        
        switch phaseNo
            case 1
                maskTile{phaseNo}=[randomMask{maskNo} zerosMask
                    zerosMask zerosMask];
            case 2
                maskTile{phaseNo}=[zerosMask randomMask{maskNo}
                    zerosMask zerosMask];
            case 3
                maskTile{phaseNo}=[zerosMask zerosMask
                    randomMask{maskNo} zerosMask];
            case 4
                maskTile{phaseNo}=[zerosMask zerosMask
                    zerosMask randomMask{maskNo}];
        end
        
        
        wholeMask{maskNo}{phaseNo}=repmat(maskTile{phaseNo}, [round(vResolution/(blockSize*2)), round(hResolution/(blockSize*2))]);
        wholeMask{maskNo}{phaseNo}=wholeMask{maskNo}{phaseNo}(1:vResolution,1:hResolution);
        
        wholeMask{maskNo}{phaseNo}(vResolution-7:vResolution,1:hResolution)=0;
        
        %         wholeMask{maskNo}{phaseNo}(1:16,1:hResolution)=0;
        %         wholeMask{maskNo}{phaseNo}(1:vResolution,1:16)=0;
        %         wholeMask{maskNo}{phaseNo}(1:vResolution,hResolution-16:hResolution)=0;
        
        
        % draw all masks
        if(plotBool==1)
            figure(102)
            colormap gray
            imagesc(wholeMask{maskNo}{phaseNo})
            title(['Whole Mask - maskNo: ', num2str(maskNo), '- phaseNo: ', num2str(phaseNo)])
            drawnow
        end
        
        cd(directoryPath)
        
        fileNameString1=sprintf('measurement_mask_%d_%02d.png',phaseNo, maskNo);
        imwrite(wholeMask{maskNo}{phaseNo}, fileNameString1);
        
        
        cd 'Mask Tiles'
        fileNameString2=sprintf('maskTile_%02d.png',maskNo);
        imwrite(randomMask{maskNo}, fileNameString2);
        
        cd ..
    end
    
end
    
