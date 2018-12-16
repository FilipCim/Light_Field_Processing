function [WP, a, DC] = mp_wrapped_phase(L, p, offset, deg)
% MP_WRAPPED_PHASE Computes wrapped phase from images.
%   WP = MP_WRAPPED_PHASE(L, p) computes wrapped phase WP from images
%   specified by the cell array L. p is a cell array having the same first
%   dimension as L. For more details see function MPS_WRAPPED_PHASE.
%
%   This function differs from MPS_WRAPPED_PHASE in the following way:
%   The function first checks if all input phases in p are a match to
%   the DFT basis. If phases do not match then the call is redirected to
%   MPS_WRAPPED_PHASE. If phases do match than a more efficient phase
%   computation method which utilizes FFT is used.
%   
%   WP = MP_WRAPPED_PHASE(L, p, offset) shifts phases for the specified
%   offest. The offset is one used during the sequence generation, i.e. if
%   the fringe pattern was generated using the phase offset of 90 degrees
%   then an offset of -90 degrees will be supplied to this function as
%   offset.
%
%   WP = MP_WRAPPED_PHASE(L, p, offset, deg), where deg is a string,
%   specifies which units the input phases in p use. deg may be 'deg' for
%   degrees, 'rad' for radians, and 'grad' for gradians. If omitted or
%   empty then 'deg' is assumed. Output wrapped phase is always in radians.
%
%   [WP, a, DC] = MPS_WRAPPED_PHASE(L, p) also returns amplitudes a and the
%   direct component DC.
%
%   See also MPS_WRAPPED_PHASE, MP_APPLY_PHASE_OFFSET,
%   MP_ARE_PHASES_FROM_DFT.

narginchk(2, 4);
nargoutchk(1, 3);

assert( iscell(L) );

% Set default parameters.
if (3 > nargin) || isempty(offset); offset = 0; end;
assert( isnumeric(offset) && (1 == numel(offset)) );

if (4 > nargin) || isempty(deg); deg = ''; end;
assert( ischar(deg) );

% Get phase offset in radians.
if strcmpi(deg, 'deg') || isempty(deg)
    offset_rad = single(offset * pi / 180);
elseif strcmpi(deg, 'rad')
    offset_rad = single(offset);
elseif strcmpi(deg, 'grad')
    offset_rad = single(offset * pi / 200);
else
    error('Unknown phase units!');
end

% Apply forward phase offset. The phase offset is required as shifted
% phases cannot be matched to DFT basis.
po = mp_apply_phase_offset(p, offset, deg);

% Check if input phases are from DFT.
[use_fft, k] = mp_are_phases_from_dft(po, deg);

% Redirect call to MPS_WRAPPED_PHASE if inputs are not from DFT basis.
% Note that phase offset is not required for MPS_WRAPPED_PHASE as it is
% more general.
if false == use_fft
    warning('Supplied phases do not match a DFT basis.');    
    [WP, a, DC] = mps_wrapped_phase(L, p, deg);    
    return;
end

% Load all images and construct data matrix.
N1 = [];
N2 = [];
N3 = numel(L);
D_gray = single([]);
for i = 1 : N3
    
    % Load image and convert it to grayscale single.
    if ischar(L{i})
        img = calibration_load_image_as_gray(L{i});
    else
        img = L{i};
    end
    img = single(img);
    if 3 == size(img, 3); img = ps_rgb2gray(img); end;
    assert( ~isempty(img) && (1 == size(img,3)) );
    
    if isempty(N1); N1 = size(img, 1); else assert(size(img, 1) == N1); end;
    if isempty(N2); N2 = size(img, 2); else assert(size(img, 2) == N2); end;
    if isempty(D_gray); D_gray = zeros([N1*N2 N3], 'single'); end;
        
    D_gray(:, i) = reshape(img, [N1*N2 1] );
end
D_gray = shiftdim(D_gray, 1);
D_gray = single(D_gray);

% Compute DFT in N points.
D_spectrum = fft(D_gray);

% Get DC component.
scl = single(1 / N3);
if 2 < nargout
    DC = scl * real(D_spectrum(1,:));
    DC = reshape(DC, [N1 N2]);
end

% Get other components.
K = mp_get_baseband_steps(N3, k);
Np = numel(K);
WP = cell(Np, 1);
if 1 < nargout; a = cell(Np, 1); end;
for i = 1 : Np
    k1 = K(i) + 1;
    k2 = N3 - K(i) + 1;
    WP{i} = angle( D_spectrum(k2,:) );
    WP{i} = reshape(WP{i}, [N1 N2]);
    WP{i} = mod(WP{i} + (pi - offset_rad), 2*pi);
    
    if 1 < nargout
        a{i} = scl * 2.0 * abs(D_spectrum(k1,:));
        a{i} = reshape(a{i}, [N1 N2]);
    end
end

if (1 == numel(WP)); WP = WP{1}; end;
if (1 == numel(a)); a = a{1}; end;