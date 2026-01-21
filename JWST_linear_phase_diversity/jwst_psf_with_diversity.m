function [I_focus, I_div, U_focus, U_div, meta] = ...
    jwst_psf_with_diversity(P, lambda, D, f, diversity_waves)

% jwst_psf_with_diversity
%   Compute JWST-like focused and defocused (phase-diversity) PSFs
%   from a complex pupil field P.
%
% Inputs:
%   P               - complex pupil field on an N x N grid
%   lambda          - wavelength [m]       (default: 2e-6)
%   D               - primary diameter [m] (default: 6.5)
%   f               - effective focal length [m] (default: 131.4)
%   diversity_waves - defocus PV at pupil edge, in waves (default: 8)
%
% Outputs:
%   I_focus, I_div  - intensity PSFs (normalized)
%   U_focus, U_div  - complex fields in image plane
%   meta            - struct with sampling info, etc.

    if nargin < 2 || isempty(lambda),        lambda = 2e-6;   end
    if nargin < 3 || isempty(D),             D = 6.5;         end
    if nargin < 4 || isempty(f),             f = 131.4;       end
    if nargin < 5 || isempty(diversity_waves), diversity_waves = 8; end

    [N, M] = size(P);
    if N ~= M
        error('P must be square.');
    end

    % Pupil coordinates in meters (assume P spans the primary diameter D)
    [x, y] = meshgrid(linspace(-D/2, D/2, N));
    r = sqrt(x.^2 + y.^2);
    rho = r / (D/2);   % normalized radius, rho = 1 at outer edge

    k = 2*pi / lambda;

    % -------- Helpers --------
    fft2c  = @(u) fftshift(fft2(ifftshift(u)));
    ifft2c = @(u) fftshift(ifft2(ifftshift(u)));

    % -------- Focused PSF --------
    U_focus = fft2c(P);
    % normalize so total energy = 1
    % U_focus = U_focus / sqrt(sum(abs(U_focus(:)).^2));
    I_focus = abs(U_focus).^2;

    % -------- Diversity defocus (Zernike-like) --------
    % Defocus OPD with PV = diversity_waves * lambda across pupil
    % W_def(ρ) = (diversity_waves * λ / 2) * (2ρ^2 - 1)
    W_def = (diversity_waves * lambda / 2) * (2*rho.^2 - 1);

    % Only apply defocus where P has aperture (to avoid NaNs outside pupil)
    P_div = P .* exp(1i * k * W_def);

    % Diversity PSF
    U_div = fft2c(P_div);
    % U_div = U_div / sqrt(sum(abs(U_div(:)).^2));
    I_div = abs(U_div).^2;

    % -------- Meta info --------
    % Angular/pixel sampling: Δθ ≈ λ / (N * Δx), Δx = D / N
    delta_x = D / N;
    delta_theta = lambda / (N * delta_x);      % radians per pixel
    plate_scale = f * delta_theta;            % meters per pixel at focal plane

    meta = struct();
    meta.lambda       = lambda;
    meta.D            = D;
    meta.f            = f;
    meta.div_waves    = diversity_waves;
    meta.delta_theta  = delta_theta;
    meta.plate_scale  = plate_scale;  % focal plane meters per pixel
end
