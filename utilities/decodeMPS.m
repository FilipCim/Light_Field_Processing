function [rows, cols, cameraPoints, projectorPoints, lut] = decodeMPS(L_directory, params)
%% SET PARAMETERS
% Set input directories and other script parameters.

% Base name for PLY outputs.
basename = 'reconstruct_mps';

% If defined tag is appended at the end of output .PLY filenames.
basename_tag = '';

% Image data directories. Note that all parameters are shared accross
% all directories indicated here.
data_dirs = { ...
    L_directory, ... % 1
    };

% Camera subdirectories; these must exist in the input directory.
L_camera = { ...
    '.', ... % Camera 0
    };
num_cam = numel(L_camera);

% Unique camera IDs; these must exist in XML calibration file.
% For our scanners top camera always has higher/larger UID.
camera_names = { ...
    '16069296', ... % Camera 0
    };
assert(numel(camera_names) == num_cam);

% Projector name; must exist in XML calibration file.
projector_name = '\\.\DISPLAY2';

% Calibration data.
% calibration_xml_file = 'e:\Output\2017-06-27 Canon LV-WX310ST + Fujinon HF12.5SA-1\calibration_canon_20170627.xml';
%calibration_xml_file = 'e:\Output\2017-06-27 Canon LV-WX310ST + Fujinon HF12.5SA-1\calibration_canon_20170627_zhang.xml';
%calibration_xml_file = 'e:\Output\2017-06-27 Canon LV-WX310ST + Fujinon HF12.5SA-1\calibration_canon_20170627_bunadj8.xml';

% Required minimal dynamic range for input data.
rel_thr = 0.01;

% Phase offset used when the sequence was generated; if any of
% MP_pattern_*.m scripts was used then phase offset is 90 degrees.
phase_offset_deg = 90;

% Input tag for row code.
tag_rows = 'row';

% Input tag for column code.
tag_cols = 'col';

% Bounding sphere for the 3D data; set to empty to skip clipping.
radius = 5000;

% Distance for filtering out bad triangulation points; set to positive
% infinity to skip clipping.
ray_dst_thr = 50;

% Threshold on unwrapped normalized phase deviation; set to positive
% infinity to skip clipping. Note that normalized phase is in [0,1] range.
phs_dev_thr = 0.005;
phs_dev_col = 5; % Width of block in which deviation is computed.
phs_dev_row = 5; % Height of block in which deviation is computed.

% Decimation factors; both factors must be whole numbers. Set to 1 for
% pixel-wise reconstruction; if set to a number k lager than 1 then every
% k-th pixel will be reconstructed.
step_row = 1;
step_col = 1;

% Flag and parameters that control normal estimation. Set the flag to false
% if speed is required as normal estimation is currently somewhat slow.
% Processing window is of MxM size where M must be odd.
estimate_normals = false;
nrm_M = 7; % Must be odd.

% Color channel gains in RGB order.
rgb_gains = 1.2 * [1.39 1.00 1.25];
assert( all(size(rgb_gains) == [1 3]) );


%% PROCESS ALL DATA DIRECTORIES
% Sequentially process all data directories indicated by the user.

for idx_dir = 1 : numel(data_dirs)
    
    % Fetch directory name.
    data_dir = data_dirs{idx_dir};
    
    disp(' ');
    disp(['Processing directory ' data_dir '.']);
    
    
    %% CREATE INPUT AND OUTPUT DIRECTORY NAMES
    % Create lists of input and output directory names and check if they exist.
    % After this operation a cell array L_input should contain all input camera
    % subdirectories and a cell array L_output should contain output
    % subdirectory for each camera. As there may be multiple cameras per one
    % output directory we also set camera tags in L_output_tag cell array so
    % filenames may be different.
    
    disp('Creating lists of input and output subdirectories.');
    
    L_input = cell(num_cam, 1);
    L_output = cell(num_cam, 1);
    L_output_tag = cell(num_cam, 1);
    
    for j = 1 : num_cam
        L_input{j} = cat(2, data_dir, filesep, L_camera{j});
        assert( 7 == exist(L_input{j}, 'dir') );
        
        L_output{j} = data_dir;
        assert( 7 == exist(L_output{j}, 'dir') );
        
        L_output_tag{j} = cat(2, '_CAM', num2str(j-1));
    end
    
    
    %% CREATE INPUT FILENAMES
    % For each input camera subdirectory create a list of filenames for row and
    % for column MPS code. Note that input phases and fringe counts are
    % automatically extracted from filenames. Also note we use specific tags in
    % the filenames to separate row from column code.
    
    disp('Creating lists of input and output files.');
    
    L_row = cell(num_cam, 1);
    Ln_row = cell(num_cam, 1);
    Lp_row = cell(num_cam, 1);
    
    L_col = cell(num_cam, 1);
    Ln_col = cell(num_cam, 1);
    Lp_col = cell(num_cam, 1);
    
    for j = 1 : num_cam
        
        pathstr = L_input{j};
        
        % Create file list for columns. Only files identified by tag_cols
        % are selected.
        [L_tmp, Ln_tmp, Lp_tmp] = mps_create_raw_file_list(pathstr, tag_cols);
        if isempty(L_tmp); [L_tmp, Ln_tmp, Lp_tmp] = mps_create_file_list(pathstr, {}, tag_cols); end;
        
        L_col{j} = L_tmp;
        Ln_col{j} = Ln_tmp;
        Lp_col{j} = Lp_tmp;
        
        % Create file list for rows. Only files identified by tag_rows are
        % selected.
        [L_tmp, Ln_tmp, Lp_tmp] = mps_create_raw_file_list(pathstr, tag_rows);
        if isempty(L_tmp); [L_tmp, Ln_tmp, Lp_tmp] = mps_create_file_list(pathstr, {}, tag_rows); end;
        
        L_row{j} = L_tmp;
        Ln_row{j} = Ln_tmp;
        Lp_row{j} = Lp_tmp;
        
        % We must have either row or column data.
        assert( ~isempty(L_col{j}) || ~isempty(L_row{j}) );
    end
    
    
    %% 3D RECONSTRUCTION
    % Now for each input directory first perform MPS decoding to find row and
    % column projector coordinate for each pixel. The compute the texture and
    % identifiy illuminated pixels. Finally, perform triangulation to get point
    % cloud.
    
    % Sanity check.
    assert( numel(L_row) == num_cam );
    assert( numel(L_col) == num_cam );
    assert( numel(L_input) == num_cam );
    assert( numel(L_output) == num_cam );
    assert( numel(L_output_tag) == num_cam );
    
    % Load calibration data.
    %     [cam, prj] = gc_load_geometry(calibration_xml_file);
    
    % Sequentially process each camera.
    for i = 1 : num_cam
        
        %% PROCESS COLUMNS
        % Process projector columns if they exist.
        
        % Get unique fringe counts.
        [cnt_col, ia] = unique(cell2mat(Ln_col{i}));
        [~, ib] = sort(ia);
        count = cnt_col(ib);
        
        if 1 < numel(count)
            
            disp(['Processing column data for CAM' L_camera{i} '.']);
            
            % Prepare parameters for MPS fringe encoding.
            [lambda, W] = mps_periods_from_fringe_counts(count);
            k = mps_get_period_tuples(lambda, W);
            [O, Xk, Xv] = mps_get_projection_matrix_and_centers(lambda, W);
            wgt = mps_get_weights(lambda);
            
            % Process each wavelength.
            M = numel(L_col{i});
            K = numel(count);
            abs_thr_col = cell(K, 1);
            scl_col = cell(K, 1);
            rng_col = cell(K, 1);
            WP_col = cell(K, 1);
            A_col = cell(K, 1);
            DC_col = cell(K, 1);
            texture_col = cell(K, 1);
            texture_col_type = cell(K, 1);
            range_col = cell(K, 1);
            range_col_type = cell(K, 1);
            for j = 1 : K
                
                % Select images for current wavelength.
                idx = [];
                for m = 1 : M; if (count(j) == Ln_col{i}{m}); idx = [idx m]; end; end;
                
                % Load images into memory.
                [I, I_type] = calibration_load_images_as_raw(L_col{i}(idx));
                
                % Get threshold.
                [abs_thr_col{j}, scl_col{j}] = calibration_get_absolute_threshold(I, rel_thr);
                
                % Compute texture.
                [texture_col{j}, texture_col_type{j}] = mps_average_images_raw(I, I_type, 'single');
                
                % Compute dynamic range.
                [range_col{j}, range_col_type{j}] = calibration_dynamic_range_raw(I, I_type);
                
                % Compute wrapped phase.
                phases_deg = cell2mat( Lp_col{i}(idx) ).';
                [WP_col{j}, A_col{j}, DC_col{j}] = mp_wrapped_phase(I, phases_deg, -phase_offset_deg);
                
            end
            
            % Unwrap phase.
            dst_phase_min_col = mps_get_minimal_distance(Xk, Xv);
            [P_col, idxk_col, dst_phase_col] = mps_unwrap_phase_NN(WP_col, k, O, Xk, Xv, count, wgt);
            
        else
            
            disp(['Missing column data for CAM' L_camera{i} '.']);
            P_col = [];
            
        end
        
        
        %% PROCESS ROWS
        % Process projector rows if they exist.
        
        % Get unique fringe counts.
        [cnt_row, ia] = unique(cell2mat(Ln_row{i}));
        [~, ib] = sort(ia);
        count = cnt_row(ib);
        
        if 1 < numel(count)
            
            disp(['Processing row data for CAM' L_camera{i} '.']);
            
            % Prepare parameters for MPS fringe encoding.
            [lambda, W] = mps_periods_from_fringe_counts(count);
            k = mps_get_period_tuples(lambda, W);
            [O, Xk, Xv] = mps_get_projection_matrix_and_centers(lambda, W);
            wgt = mps_get_weights(lambda);
            
            % Process each wavelength.
            M = numel(L_row{i});
            K = numel(count);
            abs_thr_row = cell(K, 1);
            scl_row = cell(K, 1);
            rng_row = cell(K, 1);
            WP_row = cell(K, 1);
            A_row = cell(K, 1);
            DC_row = cell(K, 1);
            texture_row = cell(K, 1);
            texture_row_type = cell(K, 1);
            range_row = cell(K, 1);
            range_row_type = cell(K, 1);
            for j = 1 : K
                
                % Select images for current wavelength.
                idx = [];
                for m = 1 : M; if (count(j) == Ln_row{i}{m}); idx = [idx m]; end; end;
                
                % Load images into memory.
                [I, I_type] = calibration_load_images_as_raw(L_row{i}(idx));
                
                % Get threshold and scale.
                [abs_thr_row{j}, scl_row{j}] = calibration_get_absolute_threshold(I, rel_thr);
                
                % Compute texture.
                [texture_row{j}, texture_row_type{j}] = mps_average_images_raw(I, I_type, 'single');
                
                % Compute dynamic range.
                [range_row{j}, range_row_type{j}] = calibration_dynamic_range_raw(I, I_type);
                
                % Compute wrapped phase.
                phases_deg = cell2mat( Lp_row{i}(idx) ).';
                [WP_row{j}, A_row{j}, DC_row{j}] = mp_wrapped_phase(I, phases_deg, -phase_offset_deg);
                
            end
            
            % Unwrap phase.
            dst_phase_min_row = mps_get_minimal_distance(Xk, Xv);
            [P_row, idxk_row, dst_phase_row] = mps_unwrap_phase_NN(WP_row, k, O, Xk, Xv, count, wgt);
            
        else
            
            disp(['Missing row data for CAM' L_camera{i} '.']);
            P_row = [];
            
        end
        
        
        %% FIND TEXTURE AND ILLUMINATED AREA
        % Illuminated area may determined by tresholding pixel dynamic ranges
        % or the amplitudes of sinusoid for each pixel. Only pixels which are
        % illuminated in all input images are selected as valid.
        % Texture is computed by averaging all DC components which are results
        % of wrapped phase computations. For this procedure to be effective the
        % number of phase shifts must be large enough so there is no fine
        % ripple visible in the resulting texture image.
        
        have_cols = ~isempty(P_col);
        have_rows = ~isempty(P_row);
        
        if have_rows && have_cols
            illuminated = mp_get_illuminated(range_col, abs_thr_col, range_row, abs_thr_row);
            [texture_raw, texture_type] = mp_average_images_raw(texture_col, texture_col_type, texture_row, texture_row_type);
            dst_phase = calibration_combine_dynamic_ranges(dst_phase_row, dst_phase_col, 'max');
        elseif have_rows && ~have_cols
            illuminated = mp_get_illuminated(range_row, abs_thr_row);
            [texture_raw, texture_type] = mp_average_images(texture_row, texture_row_type);
            dst_phase = dst_phase_row;
        elseif ~have_rows && have_cols
            illuminated = mp_get_illuminated(range_col, abs_thr_col);
            [texture_raw, texture_type] = mp_average_images(texture_col, texture_col_type);
            dst_phase = dst_phase_col;
        else
            illuminated = [];
            texture_raw = [];
            texture_type = '';
            dst_phase = [];
        end
        
           
        
        %% SELECT PIXELS FOR TRIANGULATION
        % Create coordinate arrays required for the 3D reconstruction. Note we
        % only select points that are illuminated.
        
        % Select camera pixels.
        [camera_row, camera_col, lut] = calibration_select_pixels(illuminated, [], step_row, step_col);
        
        % Select projector pixels.
        if ~isempty(P_col)
            projector_col = P_col(-1 ~= lut) * params.screenWidth;
        else
            projector_col = [];
        end
        
        if ~isempty(P_row)
            projector_row = P_row(-1 ~= lut) * params.screenHeight;
        else
            projector_row = [];
        end
        %%
        
        cameraPoints = double([camera_col, camera_row]);
        projectorPoints = double([projector_col, projector_row]);
        
        rows = P_row*params.screenHeight;
        cols = P_col*params.screenWidth;
        
%         if(params.cropFlag)
%             rows = imcrop(rows, round([params.roi(1), 0, params.roi(3), size(rows, 1)]));
%             cols = imcrop(cols, round([params.roi(1), 0, params.roi(3), size(cols, 1)]));
%         end
        
        figure
        subplot(121), imagesc(rows)
        subplot(122), imagesc(cols)
        
        
        
    end
    
end