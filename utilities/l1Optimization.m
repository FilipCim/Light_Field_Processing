function [ image_est ] = l1Optimization( y, psi_inv, theta, reconstructionOptions )
%l1Optimization Summary of this function goes here
%   Detailed explanation goes here

[M, N] = size(theta);
y = reshape(y, M, 1);

y_m = y;
y_mk = y_m;


if(strcmp(lower(reconstructionOptions.Package), 'sedumi'))
    %% SeDuMi SOLVER
    
    % standard dual form: data conditioning for minimum L1 - SeDuMi algorithm
    % input data preparation for the problem
    
    
    
    b = [ spalloc(N,1,0); -sparse(ones(N,1)) ];
    
    At = [ -sparse(theta)   ,     spalloc(M,N,0)    ;...
        sparse(theta)   ,     spalloc(M,N,0)    ;...
        speye(N)        ,    -speye(N)          ;...
        -speye(N)        ,    -speye(N)          ;...
        spalloc(N,N,0)  ,    -speye(N)          ];
    
    % SEDUMI OPTIMIZATION
    
    
    % Standard dual form: data conditioning for minimum L1
    
    c = [ -sparse(y); sparse(y); spalloc(3*N,1,0) ];
    
    % Optimization
    pars.fid=0; % suppress output
    K.l = max(size(At));
    [~,s_est]=sedumi(At,b,c,K,pars); % SeDuMi
    
    % Output data processing
    s_est=s_est(:);
    s_est=s_est(1:N);
    
    image_est = (psi_inv * s_est).';
    
elseif(strcmp(lower(reconstructionOptions.Package), 'sedumi2'))
    
    tic
    
    [s_est, status] = cs_sr06(y, theta);
    
    status.numerr;
    % Output data processing
    s_est=s_est(:);
    s_est=s_est(1:N);
    
    
    image_est = (psi_inv * s_est).';
    
    toc
elseif(strcmp(lower(reconstructionOptions.Package), 'sedumi2_noise'))
    tic
    [s_est, status] = cs_sr07(y, theta, 0.4, 0);
    
    status.numerr;
    % Output data processing
    s_est=s_est(:);
    s_est=s_est(1:N);
    
    image_est = (psi_inv * s_est).';
    toc
    
elseif(strcmp(lower(reconstructionOptions.Package), 'mex'))
%                 theta=theta./repmat(sqrt(sum(theta.^2)),[size(theta,1) 1]);
                %                 param.pos=true;
%                 param.lambda2=0.01
%                 param.ols=true;
                param.mode=2;
%                 param.pos=1;
                param.cholesky=true;
                param.lambda=reconstructionOptions.lambda; % not more than 20 non-zeros coefficients
%                 param.L=32; % not more than 10 non-zeros coefficients
                param.eps=1e-5; % squared norm of the residual should be less than 0.1
                param.numThreads=-1; % number of processors/cores to use; the default choice is -1
                % and uses all the cores of the machine
%                 W=rand(size(theta,2),1);
%                 W=1:size(theta,2);
                W=size(theta,2):-1:1;
%                 W=ones(size(theta,2),1)';
                W=W';
                
%                 s_est = mexLassoWeighted(y_m, full(theta), W, param);
                s_est = mexLasso(y_m, full(theta), param);

                %                 s_est = wmpalg('omp', y ,psi)
                
                image_est = (psi_inv * s_est);    
    
    
    
    
elseif(strcmp(lower(reconstructionOptions.Package), 'cvx'))
    
    %% CVX SOLVER - alternative to SeDuMi
    
    image_est = [];
    epsilon = max(max(std(y)));
    
    
    % Reconstructing initial signal
    %     cvx_solver(reconstructionOptions.Algorithm)
    
    cvx_begin quiet
    cvx_precision high
    
    variable s_est(64, 1);
    minimize(norm(s_est, 1));
    subject to
    norm(theta*s_est-y,1)<=epsilon
    %     theta * s_est == y;
    cvx_end
    
    image_est = (psi_inv * s_est).';
    
elseif(strcmp(lower(reconstructionOptions.Package), 'sparsa'))
    tic
    
    for i = 1:1
        
        tau = 0.001;
        
        [s_est, ~, objective] = SpaRSA(y_mk, theta, tau,...
            'Debias',0,...
            'Initialization',0,...
            'StopCriterion',1,...
            'ToleranceA',1e-10,...
            'Verbose', 0,...
            'Continuation', 0);
        
        image_est = (psi_inv * s_est).';
        
        %         image_est(k:k+block_size-1, l:l+block_size-1)= reshape(signal_est,[block_size block_size]);
        
        % --- bk update ---
        b_Axk = y_m - theta*s_est;
        
        test1(:,i)=b_Axk;
        
        if(norm(b_Axk) < 1e-10)
            disp('kraj');
            break;
        end
        
        y_mk = y_mk + b_Axk;
        
    end
    
    toc
    
    %     figure(100)
    %     imshow(image_est, 'InitialMagnification', 'fit'), title('Image Reconstruction'), colormap gray, axis image
    %     drawnow
    
elseif(strcmp(lower(reconstructionOptions.Package), 'sparsa2'))
    tau = 0.0000001;
    noOfBregIter = 5;
    
    A = @(x) theta*x;
    AT = @(x) theta'*x;
    
    % denoising function
    tv_iters = 10;
    
    %     Psi = @(x,th) tvdenoise(x, 2/th ,tv_iters);
    Psi = @(x, th) soft(x, th);
    %     Psi = @(x,th) hard(x, th);
    
    %     Psi = @(x,th) tvdenoise_sitcm(x, 2/th, tv_iters);
    
    % regularizer
    %     Phi = @(x) TVnorm_sitcm(x);
    %     Phi = @(x) l0norm(x);
    Phi = @(x) norm(x, 1);
    %       Phi = @(x)
    
    s_est = bregman_sparsa(y, theta, tau, noOfBregIter, A, AT, Psi, Phi);
    
    image_est = (psi_inv * double(s_est)).';
    
end

end


