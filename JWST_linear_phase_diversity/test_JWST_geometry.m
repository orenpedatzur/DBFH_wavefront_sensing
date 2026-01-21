close all;clear all;clc;
winsomnia(true);
rng(2);

pupil_amp = fitsread('jwst_pupil_amp.fits');        % double, 0/1
seg_index = fitsread('jwst_pupil_segments.fits');   % int16, 0..18

% zero pad
pad_n = 256;
pupil_amp = padarray(pupil_amp,[pad_n pad_n]);
seg_index = padarray(seg_index,[pad_n pad_n]);

N = size(pupil_amp, 1);
D = 6.603464;                      % JWST pupil diameter (m)
[x,y] = meshgrid(linspace(-D/2,D/2,N));

lambda = 2e-6;
k = 2*pi/lambda;
W = zeros(N);

% JWST-like 1-sigma values (wavefront, meters)
sigma_piston = 30e-9;              % 20 nm rms piston per segment
sigma_tip    = 10*10e-9 / (D/2);      % tip/tilt give ~10 nm across half aperture
sigma_tilt   = sigma_tip;
sigma_ast    = 50*10e-9 / (D/2)^2;    % astig terms ~10 nm over full segment

for s = 1:18
    mask = (seg_index == s);
    if ~any(mask,'all'), continue; end

    % segment centroid
    xs = mean(x(mask));
    ys = mean(y(mask));
    dx = x - xs;
    dy = y - ys;

    % draw random coefficients (Gaussian, zero mean, given stds)
    piston = sigma_piston * randn;
    tip    = sigma_tip    * randn;   % x-slope
    tilt   = sigma_tilt   * randn;   % y-slope
    astx   = sigma_ast    * randn;   % (x^2 - y^2)
    asty   = sigma_ast    * randn;   % (2xy)

    W_s = piston ...
        + tip  .* dx + tilt .* dy ...
        + astx .* (dx.^2 - dy.^2) ...
        + asty .* (2*dx.*dy);

    W(mask) = W(mask) + W_s(mask);
end

P = pupil_amp .* exp(1i * k * W);   % complex pupil with randomized segment WFE

[I_foc, I_div, U_foc, U_div, meta] = jwst_psf_with_diversity(P);

tic;
[phi_est, Xf_phase, Xd_phase,info] = solve_linear_diversity_full_fp(I_foc, I_div, pupil_amp, meta);
toc

% tic;
% [phi_est, Xf_phase, info] = solve_linear_diversity_lsqr(I_foc, I_div, pupil_amp, meta);
% toc

figure;
subplot(1,2,1);imagesc(log(I_foc));
subplot(1,2,2);imagesc(log(I_div));

figure;
subplot(1,2,1);imagesc(angle(P));colorbar;
subplot(1,2,2);imagesc(phi_est);colorbar;   


winsomnia(false);


