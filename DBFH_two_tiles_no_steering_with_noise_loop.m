%% DBFH_two_tiles_loop_noise_no_steering.m
% Four captures for double-blind holography with tiles 1 & 2:
% 1) only tile 1,  2) only tile 2,
% 3) tiles 1+2 overlapped using NOMINAL tilts only,
% 4) tiles 1+2 overlapped + add piston pi/2 to tile 2.

clear all; clc;close all;
winsomnia(true);
parula_with_nan_white = [1 1 1; parula(256)];
font_size = 16;

%% Parameters
img_res = 2^10; % 2^11
seg_px  = 140/2^2;     % flat-to-flat pixels per hex
f0_m    = 120;     % focal length [m]
fft_res = img_res;    % PSF FFT size
seg_flat_diam_m = 1; % mirror size [m]
lambda = 550e-9;
k = 2*pi/lambda;
Rm = seg_flat_diam_m/sqrt(3);

%% Build diffraction-limited commands
segments = make_segments(img_res, seg_px, f0_m);

%% add nontrivial phases to the segments
% pre-phasing tolerances:
sigma_opd_nm = 200; % nm
sigma_tilt_nrad = 6*50; % nrad
sigma_defocus_um = 36*5*50; % um

[piston_std, tilt_std, defocus_std] = stds_for_scramble(lambda, Rm, f0_m, ...
    sigma_opd_nm, sigma_tilt_nrad, sigma_defocus_um);

% scramble_segments expects tolerances in rads
segments = scramble_segments(segments,'piston_std', piston_std, 'tilt_std_x', tilt_std, 'tilt_std_y', tilt_std,'defocus_std', defocus_std,'seed', 4);

[phi37_scrambled, mask37_scrambled, ~] = render_selected_tiles(segments, [1:37]);
mask37_scrambled_nans = double(mask37_scrambled);
mask37_scrambled_nans(mask37_scrambled_nans==0) = nan;

figure; imagesc(phi37_scrambled.*mask37_scrambled_nans);  axis image ij off; colorbar; colormap(parula_with_nan_white);
set(gca,'FontSize',font_size);

%% ---- Capture 1: only tile 1 ----
tile_1_ind = [3];
[phi1, M1, U1] = render_selected_tiles(segments, tile_1_ind);
[~, I1]    = pupil_fft2(U1, fft_res);

%% ---- Capture 2: only tile 2 ----
tile_2_ind = [8];
[phi2, M2, U2] = render_selected_tiles(segments, tile_2_ind);
[~, I2]    = pupil_fft2(U2, fft_res);

%% ---- Capture 3: tiles 1 & 2 overlapped using NOMINAL tilt only ----
% Align tile 2 spot to tile 1 spot in UNSTACKED tilt space, then re-stack.
[phi12, M12, U12] = render_selected_tiles(segments, [tile_1_ind tile_2_ind]);
[~, I12]     = pupil_fft2(U12, fft_res);

%% ---- Capture 4: tiles 1 & 2 overlapped + piston pi/2 on tile 2 ----
segments_pi2         = segments;
segments_pi2.pistons(tile_2_ind) = segments_pi2.pistons(tile_2_ind) + pi/2;
[phi12p, M12p, U12p] = render_selected_tiles(segments_pi2, [tile_1_ind tile_2_ind]);
[~, I12p]      = pupil_fft2(U12p, fft_res);

%% Display (optional)
hex_grid = draw_hex_grid(segments);

figure;
tiledlayout(3,2,'Padding','tight','TileSpacing','compact');
nexttile; imagesc(log(I1));  axis image ij off; colormap parula; colorbar; title(['Tile ',num2str(tile_1_ind),' only']);set(gca,'FontSize',font_size);
nexttile; imagesc(log(I2));  axis image ij off; colormap parula; colorbar; title(['Tile ',num2str(tile_2_ind),' only']);set(gca,'FontSize',font_size);
nexttile; imagesc(log(I12)); axis image ij off; colormap parula; colorbar; title(['Tiles ',num2str(tile_1_ind),'+',num2str(tile_2_ind),' overlapped (nominal)']);set(gca,'FontSize',font_size);
nexttile; imagesc(log(I12p));axis image ij off; colormap parula; colorbar; title(['Tiles ',num2str(tile_1_ind),'+',num2str(tile_2_ind),' overlapped, piston \pi/2 on latter']);set(gca,'FontSize',font_size);
nexttile; imagesc(M1);hold on;imagesc(hex_grid,'AlphaData',0.3);axis image ij off; colormap parula; colorbar; title(['CS tile ',num2str(tile_1_ind)]);set(gca,'FontSize',font_size);
nexttile; imagesc(M2);hold on;imagesc(hex_grid,'AlphaData',0.3);axis image ij off; colormap parula; colorbar; title(['CS tile ',num2str(tile_2_ind)]);set(gca,'FontSize',font_size);

%% solve
Rpx = segments.meta.seg_flat_diam_px / sqrt(3);

g = gpuDevice;
I1g = gpuArray(I1); I2g = gpuArray(I2); I12g = gpuArray(I12); I12pg = gpuArray(I12p);

peak_intensity_vec = 10.^[2:9];
sim_to_e = peak_intensity_vec ./ max(I12(:)); % scaling factor to move from simulation counts to realistic electron count-per-pixel

for ind_peak_instenisty = 1:numel(peak_intensity_vec)

    this_sim_to_e = sim_to_e(ind_peak_instenisty);

    % apply adequate noise condition
    I1_noisy = add_jwst_noise(this_sim_to_e*I1g, 1,'dark_rate_e', 0 ,'read_noise_e', 0 ,'n_groups', 1,'bg_rate_e', 0,'seed',2);
    I2_noisy = add_jwst_noise(this_sim_to_e*I2g, 1,'dark_rate_e', 0 ,'read_noise_e', 0 ,'n_groups', 1,'bg_rate_e', 0,'seed',2);
    I12_noisy = add_jwst_noise(this_sim_to_e*I12g, 1,'dark_rate_e', 0 ,'read_noise_e', 0 ,'n_groups', 1,'bg_rate_e', 0,'seed',2);
    I12p_noisy = add_jwst_noise(this_sim_to_e*I12pg, 1,'dark_rate_e', 0 ,'read_noise_e', 0 ,'n_groups', 1,'bg_rate_e', 0,'seed',2);

    prep = dbh_prepare_from_four_centered(I1_noisy, I2_noisy, I12_noisy, I12p_noisy, [0 0], [0 0], fft_res, Rpx);

    fprintf('max imag(I2_al) / max real(I2_al) = %.3g\n', ...
        max(abs(imag(prep.I2_al(:)))) / max(abs(real(prep.I2_al(:)))));
    noise_score = rms(abs(prep.S(:))-prep.magA(:).*prep.magB(:)) / mean(prep.magA(:).*prep.magB(:));
    fprintf('rel. RMS of |S|-|A||B|: %.3g\n',noise_score);
    [~,mean_noise_ind] = min(abs(noise_score - mean(noise_score)));

    % --- DBH inputs in DFT layout (unshifted) ---
    magA = ifftshift(prep.magA);           % |Â|
    magB = ifftshift(prep.magB);           % |B̂|
    S    = ifftshift(prep.S);              % Â·B̂*  (complex)

    figure;
    subplot(1,3,1);imagesc(log(prep.magA));title('log(magA)');axis image ij off;
    subplot(1,3,2);imagesc(log(prep.magB));title('log(magB)');axis image ij off;
    subplot(1,3,3);imagesc(log(abs((prep.S))));title('log(abs((S)))');axis image ij off;

    % supports must be same size as magA/magB/S; center pad if needed
    suppA = gpuArray(center_padcrop(M1, size(magA), false));
    suppB = gpuArray(center_padcrop(M2, size(magA), false));

    % --- Build the linear system as operators (no huge explicit matrix) ---
    sys = dbh_prepare_system(magA, magB, S, suppA, suppB, 'lambda',0*1e-3);

    % Real-stacked wrapper for lsqr (adapts complex A/A' to lsqr's API)
    Afun = @(x,tf) Afun_lsqr(x, tf, sys.Amul, sys.ATmul, sys.n, sys.m);

    % RHS and initial guess (real-stacked)
    bR = [real(sys.b); imag(sys.b)];
    x0 = zeros(2*sys.n,1);

    assert(numel(bR) == 2*sys.m, 'bR must be length 2*m');
    assert(numel(x0) == 2*sys.n, 'x0 must be length 2*n');
    % x0 = xR;
    % --- Solve ---
    tic;
    [xR,flag,relres,iter,resvec] = lsqr(Afun, bR, 1e-12, 3e4, [], [], x0);
    solve_time = toc
    %
    % Back to complex unknown on overlap K (phasor of B in frequency domain)
    zB = xR(1:sys.n) + 1i*xR(sys.n+1:end);
    zB = zB ./ max(abs(zB), 1e-12);   % (optional) unit-modulus projection

    % Recover zA on the overlap via S = A * conj(B): zA = (S / (|A||B|)) zB
    ZB = zeros(size(S));  ZB(sys.K) = zB;
    ZA = zeros(size(S));
    den = max(magA.*magB, 1e-15);
    ZA(sys.K) = (S(sys.K) ./ den(sys.K)) .* zB;

    % --- Reconstruct pupil fields for the two tiles (DFT layout -> iFFT) ---
    Ahat_rec = magA .* ZA;      % Â_rec
    Bhat_rec = magB .* ZB;      % B̂_rec
    a_rec    = ifft2(Ahat_rec); % pupil field of tile 1 (complex)
    b_rec    = ifft2(Bhat_rec); % pupil field of tile 2 (complex)

    phi1_rec = angle(a_rec);    % recovered pupil phase (tile 1)
    phi2_rec = angle(b_rec);    % recovered pupil phase (tile 2)

    % display results (original vs. recovered, piston removed)
    M1 = logical(M1);
    M2 = logical(M2);

    % originals (your known pupil phases for tiles 1,2)
    phi1_0  = rm_piston_phi(phi1,  M1);
    % phi2_0  = rm_piston_phi(phi2,  M2);
    phi2_0  = rm_piston_phi(phi12.*M2,  M2);

    % recovered (from DBH)
    phi1r_0 = rm_piston_phi(phi1_rec, M1);
    phi2r_0 = rm_piston_phi(phi2_rec, M2);

    % differences & RMS (on their own masks)
    d1   = angle(exp(1i*(phi1r_0 - phi1_0))) .* M1;
    d2   = angle(exp(1i*(phi2r_0 - phi2_0))) .* M2;
    rms1(ind_peak_instenisty) = sqrt(mean(d1(M1).^2,'omitnan'));
    rms2(ind_peak_instenisty) = sqrt(mean(d2(M2).^2,'omitnan'));
    fprintf('RMS phase error (tile 1): %.3g rad\n', rms1(ind_peak_instenisty));
    fprintf('RMS phase error (tile 2): %.3g rad\n', rms2(ind_peak_instenisty));

    % --- errors ---
    max_rms(ind_peak_instenisty) = max(rms1(ind_peak_instenisty),rms2(ind_peak_instenisty))
    piston_error_rms_m(ind_peak_instenisty) =  max_rms(ind_peak_instenisty)/sqrt(sum(M1(:))) / (2*pi) * (lambda);
    tilt_error_rms_rads(ind_peak_instenisty) =  2*max_rms(ind_peak_instenisty)*(lambda)/(2*pi)/sqrt(sum(M1(:)))/(seg_flat_diam_m/2);
    defocus_error_rms_m(ind_peak_instenisty) = 4*sqrt(12)*(f0_m/(seg_flat_diam_m/2))^2*(max_rms(ind_peak_instenisty)*(lambda)/(2*pi))/sqrt(sum(M1(:)));

end

save('data\noisy_DBFH_workspace_4.mat',"-v7.3")


%% --- noise figure

gem = orderedcolors("gem");
% max_rms = max([rms1;rms2],[],1);
n_plot = 7;

figure;
ax1 = subplot(4,1,1);loglog(peak_intensity_vec(1:n_plot),max_rms(1:n_plot),'-o','LineWidth',2,'Color',gem(1,:));grid on;
xlabel('peak intensity [e^{-}]');ylabel('\phi_{px} RMSE [rad]');ax1.FontSize = 12;ax1.LineWidth = 1.5;
ax2 = subplot(4,1,2);loglog(peak_intensity_vec(1:n_plot),piston_error_rms_m(1:n_plot),'-o','LineWidth',2,'Color',gem(2,:));grid on;
xlabel('peak intensity [e^{-}]');ylabel('piston RMSE [m]');ax2.FontSize = 12;ax2.LineWidth = 1.5;
ax3 = subplot(4,1,3);loglog(peak_intensity_vec(1:n_plot),tilt_error_rms_rads(1:n_plot),'-o','LineWidth',2,'Color',gem(3,:));grid on;
xlabel('peak intensity [e^{-}]');ylabel('tilt x/y RMSE [rad]');ax3.FontSize = 12;ax3.LineWidth = 1.5;
ax4 = subplot(4,1,4);loglog(peak_intensity_vec(1:n_plot),defocus_error_rms_m(1:n_plot),'-o','LineWidth',2,'Color',gem(4,:));grid on;
xlabel('peak intensity [e^{-}]');ylabel('defocus RMSE [m]');ax4.FontSize = 12;ax4.LineWidth = 1.5;

% exportgraphics(gcf,'figures\errors_vs_peak_intensity.png');


%% --- Figure ---

[y, x] = find(M1);
xc_tile1 = mean(x);
yc_tile1 = mean(y);
[y, x] = find(M2);
xc_tile2 = mean(x);
yc_tile2 = mean(y);

nan_mask_tile1 = double(M1);
nan_mask_tile1(nan_mask_tile1==0) = nan;

nan_mask_tile2 = double(M2);
nan_mask_tile2(nan_mask_tile2==0) = nan;

cl1 = [-pi pi];
cl2 = 5e-2*[-pi pi];

font_size = 16;

global_phase = mean((phi1r_0-phi1_0).*nan_mask_tile1,'all',"omitnan");
global_phase = mean((phi2r_0-phi2_0).*nan_mask_tile2,'all',"omitnan");


figure;
tiledlayout(2,3, "TileSpacing", "compact", "Padding", "tight"); % Adjust TileSpacing and Padding

ax1 = nexttile; imagesc(phi1_0.*nan_mask_tile1);  axis image ij off; clim(cl1); colorbar; colormap(parula_with_nan_white);
title(['Tile  ',num2str(tile_1_ind),'  – original']);set(gca,'FontSize',font_size);


ax2 = nexttile; imagesc(phi1r_0.*nan_mask_tile1); axis image ij off; clim(cl1); colorbar;
title(['Tile ',num2str(tile_1_ind),' – recovered']);set(gca,'FontSize',font_size);

ax3 = nexttile; imagesc(d1.*nan_mask_tile1);      axis image ij off;clim(cl2);colorbar;
title(sprintf('Tile %d – diff (RMS=%.3g)',tile_1_ind,rms1));set(gca,'FontSize',font_size);

ax4 = nexttile; imagesc(phi2_0.*nan_mask_tile2);  axis image ij off; clim(cl1); colorbar;
title(['Tile  ',num2str(tile_2_ind),'  – original']);set(gca,'FontSize',font_size);

ax5 = nexttile; imagesc(phi2r_0.*nan_mask_tile2); axis image ij off; clim(cl1); colorbar;
title(['Tile  ',num2str(tile_2_ind),'  – recovered']);set(gca,'FontSize',font_size);

ax6 = nexttile; imagesc(d2.*nan_mask_tile2);      axis image ij off; clim(cl2);colorbar;
title(sprintf('Tile %d – diff (RMS=%.3g)',tile_2_ind,rms2));set(gca,'FontSize',font_size);


linkaxes([ax1,ax2,ax3]);
linkaxes([ax4,ax5,ax6]);

xlim(ax1,xc_tile1+0.8*seg_px*[-1,1]);
ylim(ax1,yc_tile1+0.8*seg_px*[-1,1]);

xlim(ax6,xc_tile2+0.8*seg_px*[-1,1]);
ylim(ax6,yc_tile2+0.8*seg_px*[-1,1]);

% exportgraphics(gcf,'figures\DBFH_tiles_3_8.png');

% save('data\noisy_DBFH_workspace.mat',"-v7.3")


%% --- end scripts ---
winsomnia(false);


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
