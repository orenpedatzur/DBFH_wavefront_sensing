function [phi_pupil, Xf_phase, info] = solve_linear_diversity( ...
            I_focus, I_div, pupil_amp, meta, opts)
% solve_linear_diversity
% Solve the linear phase-diversity system from eqs (12–14) in the memo.
% Unknown: focal-plane phasor X_f (N x N).
%
% Inputs:
%   I_focus, I_div : N x N intensity images
%   pupil_amp      : N x N pupil support (0/1, defines CS in pupil)
%   meta           : struct with at least:
%                    .lambda, .D, .div_waves
%   opts           : optional struct
%                    .max_eigs_iter, .tol, .verbose
%
% Outputs:
%   phi_pupil      : estimated pupil phase (rad), N x N
%   Xf_phase       : estimated focal-plane phasor, N x N (|X|=1)
%   info           : struct with equation counts, eigenvalue, etc.

    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'max_eigs_iter'), opts.max_eigs_iter = 200; end
    if ~isfield(opts,'tol'),           opts.tol           = 1e-7; end
    if ~isfield(opts,'verbose'),       opts.verbose       = true; end

    [N, M] = size(I_focus);
    if N ~= M
        error('I_focus must be square.');
    end
    if any(size(I_div) ~= [N N])
        error('I_div must match I_focus size.');
    end

    lambda = meta.lambda;
    D      = meta.D;
    waves  = meta.div_waves;
    k      = 2*pi / lambda;

    % --- Amplitudes in image plane ---
    A_f = sqrt(max(I_focus, 0));
    % A_d not needed explicitly after eliminating X_d

    % --- Defocus kernel in focal plane: K = F{exp(i k W_def)} ---
    [xp, yp] = meshgrid(linspace(-D/2, D/2, N));
    r   = sqrt(xp.^2 + yp.^2);
    rho = r / (D/2);
    W_def = (waves * lambda / 2) * (2*rho.^2 - 1);
    phase_def = exp(1i * k * W_def);
    K_def = fftshift(fft2(fftshift(phase_def)));
    % K_def = K_def / sqrt(mean(abs(K_def(:)).^2));   % unit RMS
    K_def = K_def / max(abs(K_def(:)));            % max |K| = 1

    % --- Compact support masks in pupil plane ---
    CS_mask    = (pupil_amp ~= 0);
    outside_CS = ~CS_mask;
    C = nnz(outside_CS);
    n_unknowns = N^2;
    n_eqs      = 2 * C;

    if opts.verbose
        fprintf('solve_linear_diversity: unknowns = %d, eqs = %d\n', ...
                n_unknowns, n_eqs);
        if n_eqs > n_unknowns
            fprintf('  System is overdetermined (good).\n');
        else
            fprintf('  System is UNDERdetermined (consider zero-padding).\n');
        end
    end

    % --- Helper: forward operator L(X_f) ---
    function y = L_forward(xvec)
        X = reshape(xvec, N, N);

        % Focus constraint: Uf = F^-1[A_f X]
        Zf = A_f .* X;
        Uf = fftshift(ifft2(ifftshift(Zf)));
        c1 = Uf(outside_CS);

        % Diversity constraint: Ud = F^-1[(A_f X) * K_def]
        Y  = conv2(Zf, K_def, 'same');
        Ud = fftshift(ifft2(ifftshift(Y)));
        c2 = Ud(outside_CS);

        y = [c1; c2];
    end

    % --- Helper: adjoint operator L^H ---
    function xvec = L_adj(zvec)
        len1 = C;
        z1 = zvec(1:len1);
        z2 = zvec(len1+1:end);

        % Back to full N x N arrays, zeros inside CS
        Uf_back = zeros(N);
        Uf_back(outside_CS) = z1;

        Ud_back = zeros(N);
        Ud_back(outside_CS) = z2;

        % Adjoint of Uf = ifft2(Zf) is Zf_grad = fft2(Uf_back)
        Zf_grad1 =  fftshift(fft2( ifftshift(Uf_back)));
        grad1 = conj(A_f) .* Zf_grad1;

        % For diversity: Ud = ifft2(conv2(Zf, K_def, 'same'))
        % Adjoint: Zf_grad2 = conv2(fft2(Ud_back), rot90(conj(K_def),2), 'same')
        T = fftshift(fft2(ifftshift(Ud_back)));
        Yadj = conv2(T, rot90(conj(K_def),2), 'same');
        grad2 = conj(A_f) .* Yadj;

        X_grad = grad1 + grad2;
        xvec = X_grad(:);
    end

    % --- Operator for A^H A ---
    function y = ATA(xvec)
        y = L_adj(L_forward(xvec));
    end

    % --- Find smallest singular vector of L via eigs on A^H A ---
    Afun = @(x) ATA(x);
    eigs_opts.isreal = false;
    eigs_opts.issym  = true;
    eigs_opts.maxit  = opts.max_eigs_iter;
    eigs_opts.tol    = opts.tol;

    [v, dmin,flag] = eigs(Afun, n_unknowns, 1, 'sm', eigs_opts);

    if opts.verbose
        fprintf('  Smallest eigenvalue of A^H A ~ %g\n', real(dmin));
    end

    Xf_est = reshape(v, N, N);

    % Project to unit modulus to get phasor
    Xf_phase = exp(1i * angle(Xf_est));

    % Pupil field & phase
    Uf_est = fftshift(ifft2(ifftshift(A_f .* Xf_phase)));
    % Uf_est(~CS_mask) = 0;
    phi_pupil = angle(Uf_est) .* CS_mask;

    % --- info struct ---
    info = struct();
    info.n_unknowns = n_unknowns;
    info.n_eqs      = n_eqs;
    info.overdetermined = (n_eqs > n_unknowns);
    info.lambda_min = dmin;
    info.convergence_flag = flag;
end
