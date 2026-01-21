function [phi_pupil, Xf_phase, info] = solve_linear_diversity_lsqr( ...
                I_focus, I_div, pupil_amp, meta, opts)
% solve_linear_diversity_lsqr
% LSQR-based solver for the linear phase-diversity system in X_f.
%
% Inputs
%   I_focus, I_div : N x N intensity images
%   pupil_amp      : N x N pupil support (0/1)
%   meta           : struct with .lambda, .D, .div_waves
%   opts           : optional: .tol, .maxit, .verbose
%
% Outputs
%   phi_pupil      : estimated pupil phase [rad]
%   Xf_phase       : N x N focal-plane phasor (|X|=1)
%   info           : struct with LSQR diagnostics

if nargin < 5, opts = struct(); end
if ~isfield(opts,'tol'),    opts.tol    = 1e-6;  end
if ~isfield(opts,'maxit'),  opts.maxit  = 200;   end
if ~isfield(opts,'verbose'),opts.verbose = true; end

[N,M] = size(I_focus);
if N ~= M, error('I_focus must be square'); end
if any(size(I_div) ~= [N N]), error('I_div size mismatch'); end

lambda = meta.lambda;
D      = meta.D;
waves  = meta.div_waves;
k      = 2*pi/lambda;

% -------- amplitudes in image plane --------
A_f = sqrt(max(I_focus,0));

% -------- defocus kernel K_def in image plane --------
[xp,yp] = meshgrid(linspace(-D/2, D/2, N));
r   = sqrt(xp.^2 + yp.^2);
rho = r / (D/2);

W_def    = (waves*lambda/2) * (2*rho.^2 - 1);
phaseDef = exp(1i*k*W_def) .* pupil_amp;      % only inside pupil
K_def    = fftshift(fft2(ifftshift(phaseDef)));

CS_mask    = (pupil_amp ~= 0);
outside_CS = ~CS_mask;
C          = nnz(outside_CS);

nUnknowns  = N^2;
nEqs       = 2*C;

if opts.verbose
    fprintf('LSQR solver: unknowns = %d, eqs = %d\n', nUnknowns, nEqs);
end

% Flat index for gauge pixel (fix global phase) - take central pixel
j0 = sub2ind([N N], ceil(N/2), ceil(N/2));
all_idx   = (1:nUnknowns).';
var_idx   = all_idx(all_idx ~= j0);   % unknowns except gauge
nVar      = numel(var_idx);

% ---------- Linear operator A*x = constraints ----------
% L_forward: full operator mapping X_f -> stacked constraints
    function y = L_forward(xvec)
        X = reshape(xvec, N, N);

        % Focus CS: Uf = F^-1[A_f X]
        Zf = A_f .* X;
        Uf = fftshift(ifft2(ifftshift(Zf)));
        c1 = Uf(outside_CS);

        % Diversity CS: Ud = F^-1[(A_f X) * K_def]
        Y  = conv2(Zf, K_def, 'same');
        Ud = fftshift(ifft2(ifftshift(Y)));
        c2 = Ud(outside_CS);

        y = [c1; c2];
    end

% Column a0 = A * e_j0  (effect of the gauge pixel alone)
e_j0 = zeros(nUnknowns,1); e_j0(j0) = 1;
a0   = L_forward(e_j0);
b    = -a0;  % move gauge column to RHS

% Afun represents A2 (all columns except j0) for LSQR
    function y = Afun(z,transp_flag)
        if strcmp(transp_flag,'notransp')
            % z is length nVar; build full x with gauge=0
            xfull        = zeros(nUnknowns,1);
            xfull(var_idx) = z;
            y = L_forward(xfull);         % = A2*z
        else
            % adjoint: A2^H * r
            r = z;   % 2C x 1
            % Back-propagate through adjoint of L_forward
            len1 = C;
            r1 = r(1:len1);
            r2 = r(len1+1:end);

            Uf_back = zeros(N); Uf_back(outside_CS) = r1;
            Ud_back = zeros(N); Ud_back(outside_CS) = r2;

            Zf_grad1 = fftshift(fft2(ifftshift(Uf_back)));
            grad1    = conj(A_f) .* Zf_grad1;

            T        = fftshift(fft2(ifftshift(Ud_back)));
            Yadj     = conv2(T, rot90(conj(K_def),2), 'same');
            grad2    = conj(A_f) .* Yadj;

            X_grad = grad1 + grad2;          % full N x N
            gfull  = X_grad(:);              % length nUnknowns
            y      = gfull(var_idx);         % drop gauge component
        end
    end

% ---------- LSQR solve A2*z = b ----------
[ z,flag,relres,iter ] = lsqr(@Afun, b, opts.tol, opts.maxit);

if opts.verbose
    fprintf('LSQR flag=%d, relres=%.3e, iter=%d\n', flag, relres, iter);
end

% reconstruct full X_f vector with gauge pixel fixed to 1
x_full          = zeros(nUnknowns,1);
x_full(j0)      = 1;
x_full(var_idx) = z;
Xf_est          = reshape(x_full, N, N);

% impose unit-modulus phasor
Xf_phase = exp(1i*angle(Xf_est));

% pupil field and phase
Uf_est = fftshift(ifft2(ifftshift(A_f .* Xf_phase)));
Uf_est(~CS_mask) = 0;
phi_pupil = angle(Uf_est) .* CS_mask;

info = struct();
info.flag   = flag;
info.relres = relres;
info.iter   = iter;
info.nEqs   = nEqs;
info.nUnknowns = nUnknowns;
end
