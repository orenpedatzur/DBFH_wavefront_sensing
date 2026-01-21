function [phi_pupil, Xf_phase, Xd_phase, info] = ...
    solve_linear_diversity_full_fp(I_focus, I_div, pupil_amp, meta, opts)
% solve_linear_diversity_full_fp
% Linear phase-diversity solver with full [X_f; X_d] unknown vector and
% defocus relation enforced in the focal plane.
%
% Inputs:
%   I_focus, I_div : N x N intensity images
%   pupil_amp      : N x N pupil support (0/1)
%   meta           : struct with fields: .lambda, .D, .div_waves
%   opts           : optional struct: .tol, .max_eigs_iter, .verbose
%
% Outputs:
%   phi_pupil      : estimated pupil phase [rad], N x N
%   Xf_phase       : N x N focal-plane phasor for focus
%   Xd_phase       : N x N focal-plane phasor for diversity
%   info           : diagnostics

    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'tol'),           opts.tol = 1e-10; end
    if ~isfield(opts,'max_eigs_iter'), opts.max_eigs_iter = 1e3; end
    if ~isfield(opts,'verbose'),       opts.verbose = true; end

    [N,M] = size(I_focus);
    if N ~= M, error('I_focus must be square'); end
    if any(size(I_div) ~= [N N]), error('I_div must match I_focus'); end

    lambda = meta.lambda;
    D      = meta.D;
    waves  = meta.div_waves;
    k      = 2*pi / lambda;

    % centered FFT helpers
    fft2c  = @(u) fftshift(fft2(ifftshift(u)));
    ifft2c = @(u) fftshift(ifft2(ifftshift(u)));

    % image-plane amplitudes
    A_f = sqrt(max(I_focus,0));
    A_d = sqrt(max(I_div,0));

    % defocus phase in pupil plane
    [xp,yp] = meshgrid(linspace(-D/2, D/2, N));
    r   = sqrt(xp.^2 + yp.^2);
    rho = r / (D/2);
    W_def    = (waves*lambda/2) * (2*rho.^2 - 1);
    phaseDef = exp(1i * k * W_def);

    % compact support in pupil
    CS_mask    = (pupil_amp ~= 0);
    outside_CS = ~CS_mask;
    C  = nnz(outside_CS);
    S  = nnz(CS_mask);

    nX        = N^2;
    nUnknowns = 2 * nX;
    nEqs      = 2*C + nX;   % Uf(out), Ud(out), focal-plane defocus (all pixels)

    overdet  = nEqs >  nUnknowns;
    underdet = nEqs <  nUnknowns;

    if opts.verbose
        fprintf('solve_linear_diversity_full_fp:\n');
        fprintf('  unknowns   = %d\n', nUnknowns);
        fprintf('  equations  = %d\n', nEqs);
        if overdet
            fprintf('  system is OVERdetermined.\n');
        elseif underdet
            fprintf('  system is UNDERdetermined (consider zero-padding / extra diversities).\n');
        else
            fprintf('  system is square (nEqs == nUnknowns).\n');
        end
    end

    % ===== Linear operator L : [X_f; X_d] -> stacked constraints =====
    function y = L_forward(z)
        Xf_vec = z(1:nX);
        Xd_vec = z(nX+1:end);
        Xf = reshape(Xf_vec, N, N);
        Xd = reshape(Xd_vec, N, N);

        % pupil fields
        Uf = ifft2c(A_f .* Xf);
        Ud = ifft2c(A_d .* Xd);

        % 1) CS focus: Uf(outside_CS) = 0
        c1 = Uf(outside_CS);

        % 2) CS diversity: Ud(outside_CS) = 0
        c2 = Ud(outside_CS);

        % 3) focal-plane defocus relation:
        %    A_d X_d - FFT{ phaseDef .* Uf } = 0   at all pixels
        rhs = fft2c(phaseDef .* Uf);
        c3  = A_d .* Xd - rhs;        % N x N
        c3  = c3(:);

        y = [c1; c2; c3];
    end

    % ===== Adjoint L^H =====
    function z = L_adj(y)
        idx1 = C;
        idx2 = 2*C;
        c1 = y(1:idx1);
        c2 = y(idx1+1:idx2);
        c3 = y(idx2+1:end);

        % reshape pieces
        c3_mat = reshape(c3, N, N);

        % gradients in pupil plane
        Uf_grad = zeros(N);
        Ud_grad = zeros(N);

        % from c1 = Uf(outside)
        Uf_grad(outside_CS) = Uf_grad(outside_CS) + c1;

        % from c2 = Ud(outside)
        Ud_grad(outside_CS) = Ud_grad(outside_CS) + c2;

        % from c3 = A_d X_d - FFT{ phaseDef .* Uf }
        % w.r.t X_d: grad_Xd += conj(A_d) .* c3_mat
        Zd_grad = conj(A_d) .* c3_mat;

        % w.r.t Uf: term is - FFT{ phaseDef .* Uf }
        % adjoint: grad_Uf += - conj(phaseDef) .* IFFT(c3_mat)
        Uf_grad = Uf_grad - conj(phaseDef) .* ifft2c(c3_mat);

        % back to image plane from Uf_grad, Ud_grad
        Zf_grad = fft2c(Uf_grad);
        Xf_grad = conj(A_f) .* Zf_grad;

        % Ud_grad -> Xd extra contribution (CS only)
        Zd_extra = fft2c(Ud_grad);
        Xd_extra = conj(A_d) .* Zd_extra;

        Xd_grad = Zd_grad + Xd_extra;

        z = [Xf_grad(:); Xd_grad(:)];
    end

    % ===== Operator for A^H A =====
    function y = ATA(z)
        y = L_adj(L_forward(z));
    end

    Afun = @(z) ATA(z);

    % eigs options
    eigs_opts.isreal = false;
    eigs_opts.issym  = true;
    eigs_opts.tol    = opts.tol;
    eigs_opts.maxit  = opts.max_eigs_iter;

    [v, dmin,exit_flag] = eigs(Afun, nUnknowns, 1, 'sm', eigs_opts);

    if opts.verbose
        fprintf('  smallest eigenvalue of A^H A ≈ %g\n', real(dmin));
    end

    Xf_est = reshape(v(1:nX),      N, N);
    Xd_est = reshape(v(nX+1:end),  N, N);

    Xf_phase = exp(1i * angle(Xf_est));
    Xd_phase = exp(1i * angle(Xd_est));

    Uf_est = ifft2c(A_f .* Xf_phase);
    Uf_est(~CS_mask) = 0;
    phi_pupil = angle(Uf_est) .* CS_mask;

    info = struct();
    info.nUnknowns      = nUnknowns;
    info.nEqs           = nEqs;
    info.overdetermined = overdet;
    info.underdetermined = underdet;
    info.lambda_min     = dmin;
    info.exit_flag     = exit_flag;

end
