%% DBFH_all_segments_simultaneously_noisy.m
% Four captures for double-blind holography with tiles 1 & 2:
% 1) only odd tiles,  2) only even tiles
% 3) All tiles, even overlapping odd uing NOMINAL tilts only,
% 4) All tiles, even overlapping odd uing NOMINAL tilts only with pi/2 phase shift on even,

clear all; clc;close all;
winsomnia(true);
parula_with_nan_white = [1 1 1; parula(256)];
font_size = 16;

%% Parameters
img_res = 2*2^10;
seg_px  = 140;     % flat-to-flat pixels per hex TODO: why only 140 works?
f0_m    = 120;     % focal length [m]
fft_res = img_res;    % PSF FFT size
seg_flat_diam_m = 6; % mirror size [m]
lambda = 550e-9;
k = 2*pi/lambda;
Rm = seg_flat_diam_m/sqrt(3);

%% noise specs - https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-detector-overview/nircam-detector-performance#gsc.tab=0
num_frames = 100;
E_peak = 75e3; % JWST well capacity ~ 100K e-
t_exp   = 1;        % seconds
nGroups = 50;       % up-the-ramp groups
dark_rate_e = 0.001; % NIRCam-like dark current
read_noise_e = 13;  % CDS read noise
bg_rate_e = 0; % optional background

%% Build diffraction-limited commands (uses your make_segments)
segments = make_segments(img_res, seg_px, f0_m);

[~, ~, U37] = render_selected_tiles(segments, [1:37]);
[~, I37]    = pupil_fft2(U37, fft_res);

% fig_unstacked = figure;imagesc((I37));  axis image  ij ; colorbar; colormap gray;
% title('scrambled and unstacked','Color','w');

%% add nontrivial phases to the segments
% pre-phasing tolerances:
sigma_opd_nm = 200; % nm
sigma_tilt_nrad = 50; % nrad
sigma_defocus_um = 5*50; % um

[piston_std, tilt_std, defocus_std] = stds_for_scramble(lambda, Rm, f0_m, ...
    sigma_opd_nm, sigma_tilt_nrad, sigma_defocus_um);

% scramble_segments expects tolerances in rads
segments = scramble_segments(segments, ...
    'piston_std', piston_std, 'tilt_std_x', tilt_std, 'tilt_std_y', tilt_std,'defocus_std', defocus_std,'seed', 4);

[phi37_scrambled, mask37_scrambled, ~] = render_selected_tiles(segments, [1:37]);
mask37_scrambled_nans = double(mask37_scrambled);
mask37_scrambled_nans(mask37_scrambled_nans==0) = nan;

figure; imagesc(phi37_scrambled.*mask37_scrambled_nans);  axis image ij off; colorbar; colormap(parula_with_nan_white);
set(gca,'FontSize',font_size);
exportgraphics(gcf,'figures\scrambled_phases.png');


%% --- only odd tiles ---
[phi_odd, M_odd, U_odd] = render_selected_tiles(segments, [1:2:37]);
[~, I_odd]    = pupil_fft2(U_odd, fft_res);


% --- only even tiles ---
[phi_even, M_even, U_even] = render_selected_tiles(segments, [2:2:37]);
[~, I_even]    = pupil_fft2(U_even, fft_res);


% --- odd+even ---
[phi37, ~, U37] = render_selected_tiles(segments, [1:37]);
[~, I37]    = pupil_fft2(U37, fft_res);


% --- odd+even(+pi/2)---
segments_pi         = segments;
segments_pi.pistons([2:2:37]) = segments_pi.pistons([2:2:37]) + pi/2;
[phi37p, ~, U37p] = render_selected_tiles(segments_pi, [1:37]);
[~, I37p]    = pupil_fft2(U37p, fft_res);


%% add poisson shot-noise
sim_to_e = E_peak / max(I37(:)); % scaling factor to move from simulation counts to realistic electron count-per-pixel

I_odd_e = sim_to_e*I_odd;
I_even_e = sim_to_e*I_even;
I37_e = sim_to_e*I37;
I37p_e = sim_to_e*I37p;

I_odd_e_noisy = zeros(size(I_odd_e));
I_even_e_noisy = zeros(size(I_even_e));
I37_e_noisy = zeros(size(I37_e));
I37p_e_noisy = zeros(size(I37p_e));

for ind_frame = 1:num_frames % add noisy frames iteratively
    I_odd_e_noisy = I_odd_e_noisy + add_jwst_noise(I_odd_e, t_exp,...
        'dark_rate_e',dark_rate_e, ...     % NIRCam-like dark current
        'read_noise_e', read_noise_e,  ...     % CDS read noise
        'n_groups', nGroups, ...
        'bg_rate_e', bg_rate_e);            % optional background

    I_even_e_noisy = I_even_e_noisy + add_jwst_noise(I_even_e, t_exp,...
        'dark_rate_e',dark_rate_e, ...     % NIRCam-like dark current
        'read_noise_e', read_noise_e,  ...     % CDS read noise
        'n_groups', nGroups, ...
        'bg_rate_e', bg_rate_e);            % optional background

    I37_e_noisy = I37_e_noisy + add_jwst_noise(I37_e, t_exp,...
        'dark_rate_e',dark_rate_e, ...     % NIRCam-like dark current
        'read_noise_e', read_noise_e,  ...     % CDS read noise
        'n_groups', nGroups, ...
        'bg_rate_e', bg_rate_e);            % optional background

    I37p_e_noisy = I37p_e_noisy + add_jwst_noise(I37p_e, t_exp,...
        'dark_rate_e',dark_rate_e, ...     % NIRCam-like dark current
        'read_noise_e', read_noise_e,  ...     % CDS read noise
        'n_groups', nGroups, ...
        'bg_rate_e', bg_rate_e);            % optional background
end


%% figure captures

figure;tiledlayout(1,4, "TileSpacing", "compact", "Padding", "tight"); % Adjust TileSpacing and Padding
ax1 = nexttile; 
imagesc(log(I_odd_e_noisy));  axis image  ij ; colorbar; colormap gray;
xlim(img_res/12*[-1,1]+img_res/2);
ylim(img_res/12*[-1,1]+img_res/2);
title('only odd');

ax2 = nexttile; imagesc(log(I_even_e_noisy));  axis image  ij ; colorbar; colormap gray;
xlim(img_res/12*[-1,1]+img_res/2);
ylim(img_res/12*[-1,1]+img_res/2);
title('only even');

ax3 = nexttile; imagesc(log(I37_e_noisy));  axis image  ij ; colorbar; colormap gray;
xlim(img_res/12*[-1,1]+img_res/2);
ylim(img_res/12*[-1,1]+img_res/2);
title('odd+even');

ax4 = nexttile; imagesc(log(I37p_e_noisy));  axis image  ij ; colorbar; colormap gray;
xlim(img_res/12*[-1,1]+img_res/2);
ylim(img_res/12*[-1,1]+img_res/2);
title('odd+even(+pi/2)');

% insets
scale = 3;
inset_width = ax1.Position(3)/scale;
inset_height = ax1.Position(4)/scale;
inset_offset_x = ax1.Position(3)*(scale-1)/(scale);
inset_offset_y = ax1.Position(4)*(scale-1)/(scale);

ax1_inset = axes('Position',[ax1.Position(1:2),0,0]+[inset_offset_x,inset_offset_y,inset_width,inset_height]);
imagesc(log(I_odd_e_noisy));  axis image  ij ;
ax1_inset.XColor = 'w'; 
ax1_inset.YColor = 'w'; 
ax1_inset.FontSize = 4;
ax1_inset.XTick = 0:500:img_res;
ax1_inset.XTickLabelRotation = 90;
ax1_inset.YTick = 0:500:img_res;


ax2_inset = axes('Position',[ax2.Position(1:2),0,0]+[inset_offset_x,inset_offset_y,inset_width,inset_height]);
imagesc(log(I_even_e_noisy));  axis image  ij ;
ax2_inset.XColor = 'w'; 
ax2_inset.YColor = 'w'; 
ax2_inset.FontSize = 4;
ax2_inset.XTick = 0:500:img_res;
ax2_inset.XTickLabelRotation = 90;
ax2_inset.YTick = 0:500:img_res;

ax3_inset = axes('Position',[ax3.Position(1:2),0,0]+[inset_offset_x,inset_offset_y,inset_width,inset_height]);
imagesc(log(I37_e_noisy));  axis image  ij ;
ax3_inset.XColor = 'w'; 
ax3_inset.YColor = 'w'; 
ax3_inset.FontSize = 4;
ax3_inset.XTick = 0:500:img_res;
ax3_inset.XTickLabelRotation = 90;
ax3_inset.YTick = 0:500:img_res;

ax4_inset = axes('Position',[ax4.Position(1:2),0,0]+[inset_offset_x,inset_offset_y,inset_width,inset_height]);
imagesc(log(I37p_e_noisy));  axis image  ij ;
ax4_inset.XColor = 'w'; 
ax4_inset.YColor = 'w'; 
ax4_inset.FontSize = 4;
ax4_inset.XTick = 0:500:img_res;
ax4_inset.XTickLabelRotation = 90;
ax4_inset.YTick = 0:500:img_res;

exportgraphics(gcf,'figures\inputs_intensities_to_DBFH_all_segments.png');


%% prepare DBFH data

magA_c = sqrt(I_odd_e_noisy);
magB_c = sqrt(I_even_e_noisy);
C0   = I37_e_noisy - I_odd_e_noisy - I_even_e_noisy;
C90  = I37p_e_noisy - I_odd_e_noisy - I_even_e_noisy;
S_c = 0.5*(C0+1i*C90);

magA = ifftshift(magA_c);
magB = ifftshift(magB_c);
S = ifftshift(S_c);

fprintf('rel. RMS of |S|-|A||B|: %.3g\n',rms(abs(S(:))-magA(:).*magB(:)) / mean(magA(:).*magB(:)));

figure;
subplot(1,3,1);imagesc(log(magA_c));title('log(magA)');axis image ij off;
subplot(1,3,2);imagesc(log(magB_c));title('log(magB)');axis image ij off;
subplot(1,3,3);imagesc(log(abs(S_c)));title('log(abs((S)))');axis image ij off;

suppA = M_odd;
suppB = M_even;


%% --- solve ---

% --- Build the linear system as operators (no huge explicit matrix) ---
sys = dbh_prepare_system(magA, magB, S, suppA, suppB, 'lambda', 1e-3);

% Real-stacked wrapper for lsqr (adapts complex A/A' to lsqr's API)
Afun = @(x,tf) Afun_lsqr(x, tf, sys.Amul, sys.ATmul, sys.n, sys.m);

% RHS and initial guess (real-stacked)
bR = [real(sys.b); imag(sys.b)];
x0 = zeros(2*sys.n,1);

assert(numel(bR) == 2*sys.m, 'bR must be length 2*m');
assert(numel(x0) == 2*sys.n, 'x0 must be length 2*n');

% --- Solve ---
% lsqr
tic;
[xR,flag,relres,iter] = lsqr(Afun, bR, 1e-10, 1e5, [], [], x0);
solve_time_lsqr = toc;

% Back to complex unknown on overlap K (phasor of B in frequency domain)
zB = xR(1:sys.n) + 1i*xR(sys.n+1:end);
zB = zB ./ max(abs(zB), 1e-12);   % (optional) unit-modulus projection

% Recover zA on the overlap via S = A * conj(B): zA = (S / (|A||B|)) zB
ZB = zeros(size(S));  ZB(sys.K) = zB;
ZA = zeros(size(S));
den = max(magA.*magB, 1e-15);
ZA(sys.K) = (S(sys.K) ./ den(sys.K)) .* zB;

%% --- Reconstruct pupil fields for the two tiles (DFT layout -> iFFT) ---
Ahat_rec = magA .* ZA;      % Â_rec
Bhat_rec = magB .* ZB;      % B̂_rec
a_rec    = ifft2(Ahat_rec); % pupil field of tile 1 (complex)
b_rec    = ifft2(Bhat_rec); % pupil field of tile 2 (complex)

phi1_rec = angle(a_rec);    % recovered pupil phase (tile 1)
phi2_rec = angle(b_rec);    % recovered pupil phase (tile 2)

% originals (your known pupil phases for tiles 1,2)
phi1_0  = rm_piston_phi(phi_odd,  suppA);
phi2_0  = rm_piston_phi(phi_even,  suppB);

% recovered (from DBH)
phi1r_0 = rm_piston_phi(phi1_rec, suppA);
phi2r_0 = rm_piston_phi(phi2_rec, suppB);

% differences & RMS (on their own masks)
d1   = angle(exp(1i*(phi1r_0 - phi1_0))) .* suppA;
d2   = angle(exp(1i*(phi2r_0 - phi2_0))) .* suppB;
rms1 = sqrt(mean(d1(logical(suppA)).^2,'omitnan'));
rms2 = sqrt(mean(d2(logical(suppB)).^2,'omitnan'));
fprintf('RMS phase error (tile 1): %.3g rad\n', rms1);
fprintf('RMS phase error (tile 2): %.3g rad\n', rms2);

if 1
    %%
    [y, x] = find(suppA);
    xc_tile1 = mean(x);
    yc_tile1 = mean(y);
    [y, x] = find(suppB);
    xc_tile2 = mean(x);
    yc_tile2 = mean(y);

    nan_mask_tile1 = double(suppA);
    nan_mask_tile1(nan_mask_tile1==0) = nan;

    nan_mask_tile2 = double(suppB);
    nan_mask_tile2(nan_mask_tile2==0) = nan;

    cl1 = [-pi pi];
    cl2 = 2e-1*[-pi pi];

    font_size = 16;

    figure;
    tiledlayout(2,3, "TileSpacing", "compact", "Padding", "tight"); % Adjust TileSpacing and Padding
    ax1 = nexttile; imagesc(phi1_0.*nan_mask_tile1);  axis image ij off; clim(cl1); colorbar; colormap(parula_with_nan_white);
    title(['odd tiles - original']);set(gca,'FontSize',font_size);
    ax2 = nexttile; imagesc(phi1r_0.*nan_mask_tile1); axis image ij off; clim(cl1); colorbar;
    title(['odd tiles – recovered']);set(gca,'FontSize',font_size);
    ax3 = nexttile; imagesc(d1.*nan_mask_tile1);      axis image ij off;clim(cl2);colorbar;
    title(sprintf('odd tiles – diff (RMS=%.3g)',rms1));set(gca,'FontSize',font_size);
    ax4 = nexttile; imagesc(phi2_0.*nan_mask_tile2);  axis image ij off; clim(cl1); colorbar;
    title(['even tiles – original']);set(gca,'FontSize',font_size);
    ax5 = nexttile; imagesc(phi2r_0.*nan_mask_tile2); axis image ij off; clim(cl1); colorbar;
    title(['even tiles – recovered']);set(gca,'FontSize',font_size);
    ax6 = nexttile; imagesc(d2.*nan_mask_tile2);      axis image ij off; clim(cl2);colorbar;
    title(sprintf('even tiles – diff (RMS=%.3g)',rms2));set(gca,'FontSize',font_size);

    linkaxes([ax1,ax2,ax3]);
    linkaxes([ax4,ax5,ax6]);

    xlim(ax1,img_res/2+4*seg_px*[-1,1]);
    ylim(ax1,img_res/2+4*seg_px*[-1,1]);
    xlim(ax6,img_res/2+4*seg_px*[-1,1]);
    ylim(ax6,img_res/2+4*seg_px*[-1,1]);
    exportgraphics(gcf,'figures\DBFH_odd_even.png');

end


%% --- end scripts ---
winsomnia(false);

%% ----------------- helpers (local functions) -----------------
function [segments_overlap,dtx,dty] = overlap_by_nominal_tilt_unstacked(unstacked_segments, ind1, ind2,unstack_ratio, fftN, signNom)
%OVERLAP_BY_NOMINAL_TILT_UNSTACKED  Use ONLY nominal tilts to overlap PSFs in UNSTACKED space.
%   segments_overlap = overlap_by_nominal_tilt_unstacked(unstacked_segments, ind1, ind2, fftN, signNom)
%
% Inputs
%   unstacked_segments : struct with residual (UNSTACKED) tilts
%   ind1, ind2         : segment indices (PSF of ind2 moves to ind1)
%   fftN               : optional FFT size for debug print (default = meta.img_res)
%   signNom            : optional sign for applying nominal delta; default = -1
%                        (i.e., tx2 <- tx2 - (tx_nom1 - tx_nom2))
%
% Output
%   segments_overlap   : same struct (still UNSTACKED); only ind2 tilts adjusted.

if nargin < 5 || isempty(fftN)
    fftN = unstacked_segments.meta.img_res;
end
if nargin < 6 || isempty(signNom)
    signNom = -1;   % your observed convention needs a minus
end

segments_overlap = unstacked_segments;
N = numel(unstacked_segments.tilt_x);
assert(all([ind1>=1 , ind1<=N , ind2>=1 , ind2<=N]), 'ind1/ind2 out of range');

% --- geometry / nominal model ---
Rpx  = unstacked_segments.meta.seg_flat_diam_px / sqrt(3);
aX   = 1.5*Rpx;  aY = sqrt(3)*Rpx;
axial = generate_axial_37();
q = axial(:,1); r = axial(:,2);
Xc = aX*q; Yc = aY*(r + q/2);
u0 = Xc/Rpx; v0 = Yc/Rpx;

f0   = unstacked_segments.meta.focal_length_m;
lam  = unstacked_segments.meta.lambda_m;
pixm = unstacked_segments.meta.pixel_pitch_m;
if isinf(f0) || f0==0
    c0 = 0;
else
    Rm = Rpx * pixm;
    c0 = -(pi/lam) * (Rm^2) * (1/f0);   % phi = c0 (u^2+v^2)
end

tx_nom = unstack_ratio*2*c0*u0;   ty_nom = unstack_ratio*2*c0*v0;

% nominal delta: (1) - (2)
dtx = tx_nom(ind1) - tx_nom(ind2);
dty = ty_nom(ind1) - ty_nom(ind2);

% apply with chosen sign (default: minus)
tx = segments_overlap.tilt_x(:);
ty = segments_overlap.tilt_y(:);

% Predicted PSF shift (debug) before
kx_b = (tx/(2*pi*Rpx))*fftN;  ky_b = (ty/(2*pi*Rpx))*fftN;

tx(ind2) = tx(ind2) + signNom * dtx;
ty(ind2) = ty(ind2) + signNom * dty;

% Predicted PSF shift (debug) after
kx_a = (tx/(2*pi*Rpx))*fftN;  ky_a = (ty/(2*pi*Rpx))*fftN;

if numel(tx)==1
    fprintf(['[overlap_nominal_unstacked] sign=%+d, Δt_nom=(%.4g,%.4g) rad  |  ' ...
        'PSF px BEFORE i1=(%.3f,%.3f) i2=(%.3f,%.3f)  ->  AFTER i1=(%.3f,%.3f) i2=(%.3f,%.3f)\n'], ...
        signNom, dtx, dty, ...
        kx_b(ind1), ky_b(ind1), kx_b(ind2), ky_b(ind2), ...
        kx_a(ind1),  ky_a(ind1),  kx_a(ind2),  ky_a(ind2));
end

segments_overlap.tilt_x = reshape(tx, size(unstacked_segments.tilt_x));
segments_overlap.tilt_y = reshape(ty, size(unstacked_segments.tilt_y));
end


function Y = center_padcrop(X, outSz, padval)
if nargin < 3, padval = 0; end
[h,w] = size(X);  H = outSz(1);  W = outSz(2);
Y = cast(padval*ones(H,W), 'like', X);
r0 = floor((H - h)/2);  c0 = floor((W - w)/2);
rs = max(1,1-r0) : min(h, H-r0);
cs = max(1,1-c0) : min(w, W-c0);
rd = (rs + r0);  cd = (cs + c0);
Y(rd, cd) = X(rs, cs);
end

function y = Afun_lsqr(x, transpFlag, Amul, ATmul, n, m)
% Real-stacked adapter for LSQR.
% - If transpFlag = 'notransp' (or 0): y = [Re;Im]( A * (Re+ i Im) )
% - If transpFlag = 'transp'   (or 1): y = [Re;Im]( A' * (Re+ i Im) )

if (ischar(transpFlag) && strcmpi(transpFlag,'transp')) || isequal(transpFlag,1)
    xr = x(1:m);  xi = x(m+1:end);
    yc = ATmul(xr + 1i*xi);       % complex m→n
    y  = [real(yc); imag(yc)];
else
    xr = x(1:n);  xi = x(n+1:end);
    yc = Amul(xr + 1i*xi);        % complex n→m
    y  = [real(yc); imag(yc)];
end
end

function phi0 = rm_piston_phi(phi, mask)
% Remove piston using circular mean over masked pixels.
% mask can be logical or numeric (0/1); outputs zero outside mask.
m = logical(mask);
if ~any(m,'all')
    phi0 = zeros(size(phi), 'like', phi);
    return
end
% circular mean phase on the mask
phibar = angle(mean(exp(1i*phi(m)), 'omitnan'));
% subtract and wrap to [-pi,pi), then zero outside mask
phi0 = angle(exp(1i*(phi - phibar))) .* m;
end

function dtilt2 = nominal_pair_dtilt(unstacked_segments, ind1, ind2)
% Return the NOMINAL tilt to move tile ind2's *unstacked* PSF onto tile ind1's,
% using only the nominal focusing model (no measured/residual terms).

% geometry
Rpx  = unstacked_segments.meta.seg_flat_diam_px / sqrt(3);
aX   = 1.5 * Rpx;
aY   = sqrt(3) * Rpx;
axial = generate_axial_37();
q = axial(:,1); r = axial(:,2);

Xc = aX * q;                   % pixels (mosaic centered)
Yc = aY * (r + q/2);
u  = Xc / Rpx;
v  = Yc / Rpx;

% nominal quadratic coefficient
f0   = unstacked_segments.meta.focal_length_m;
lam  = unstacked_segments.meta.lambda_m;
pixm = unstacked_segments.meta.pixel_pitch_m;
Rm   = Rpx * pixm;             % meters
if isinf(f0) || f0 == 0
    c0 = 0;
else
    c0 = -(pi/lam) * (Rm^2) * (1/f0);   % phi = c0*(u^2+v^2)
end

% delta tilt to move #2 onto #1 in the UNSTACKED plane (nominal only)
du = u(ind2) - u(ind1);
dv = v(ind2) - v(ind1);
dtilt2 = [ 2*c0*du,  2*c0*dv ];         % radians in the u,v basis
end

% % -test-
% function dtilt_nom = nominal_pair_dtilt(unstacked_segments, i1, i2)
%     Rpx  = unstacked_segments.meta.seg_flat_diam_px / sqrt(3);
%     aX   = 1.5*Rpx;  aY = sqrt(3)*Rpx;
%     axial = generate_axial_37();
%     q=axial(:,1); r=axial(:,2);
%     u = (aX*q)/Rpx;  v = (aY*(r+q/2))/Rpx;
%
%     f0   = unstacked_segments.meta.focal_length_m;
%     lam  = unstacked_segments.meta.lambda_m;
%     pixm = unstacked_segments.meta.pixel_pitch_m;
%     Rm   = Rpx * pixm;
%     c0   = (isinf(f0) || f0==0) ? 0 : -(pi/lam) * (Rm^2) * (1/f0);
%     % delta tilt of tile 2 to land on tile 1 in UNSTACKED space (nominal only)
%     dtilt_nom = [ 2*c0*(u(i1)-u(i2)), 2*c0*(v(i1)-v(i2)) ];
% end
