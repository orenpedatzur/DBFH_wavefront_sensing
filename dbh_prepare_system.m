function sys = dbh_prepare_system(magA, magB, S, suppA, suppB, opts)
%DBH_PREPARE_SYSTEM  Build the linear DBFH system A x = b as operators.
%   sys = dbh_prepare_system(magA,magB,S,suppA,suppB,opts)
%
% Unknown is the complex phasor zB(k) = exp(i*phiB(k)) on the frequency
% overlap K = (|A|>0 & |B|>0). Equations:
%   ifft2( (S./|B|) .* ZB ) = 0  outside suppA
%   ifft2(  |B|      .* ZB ) = 0  outside suppB
% plus one anchor row that fixes global phase of ZB at one overlap pixel.
%
% Inputs
%   magA, magB : NxN double, nonnegative (|Â|, |B̂|)
%   S          : NxN complex, cross term (Â · B̂*)
%   suppA,B    : NxN logical, pupil supports (true = inside)
%   opts (all optional):
%       .lambda    (default 0)   Tikhonov weight (adds λ·I rows)
%       .anchor    (default [])  linear index in K to anchor; auto if empty
%       .anchorVal (default 1+0i) complex value at the anchor
%
% Output: struct sys with
%   sys.Amul(z)    -> y = A*z        (complex operator)
%   sys.ATmul(y)   -> x = A'*y       (complex adjoint)
%   sys.b          -> RHS (complex)
%   sys.m, sys.n   -> sizes (#rows, #cols)
%   sys.K          -> logical overlap mask (unknown support)
%   sys.idxA_out   -> indices of pixels outside suppA (pupil plane)
%   sys.idxB_out   -> indices of pixels outside suppB
%
% Helper to use with real-stacked LSQR is shown after the function.

    arguments
        magA double
        magB double
        S    double
        suppA logical
        suppB logical
        opts.lambda (1,1) double {mustBeNonnegative} = 0
        opts.anchor = []
        opts.anchorVal (1,1) double = 1+0i
    end

    % ---- sanity / shared sizes ----
    [H,W] = size(magA);
    assert(isequal(size(magA),size(magB),size(S),size(suppA),size(suppB)), ...
        'All inputs must be HxW and of equal size.');

    % Overlap (where unknown lives)
    epsM = 1e-15;
    K = (magA > epsM) & (magB > epsM) & isfinite(S);
    % central_fraction = 100;
    % thrA = max(magA(round(H/2-H/central_fraction):round(H/2+H/central_fraction),round(W/2-W/central_fraction):round(W/2+W/central_fraction)),[],'all');
    % thrB = max(magB(round(H/2-H/central_fraction):round(H/2+H/central_fraction),round(W/2-W/central_fraction):round(W/2+W/central_fraction)),[],'all');
    % K = (magA > thrA) & (magB > thrB) & isfinite(S);

    nK = nnz(K);
    if nK==0, error('No frequency overlap between |A| and |B|.'); end

    % Amplitudes for the two constraint blocks
    W_A = S ./ max(magB, epsM);   % zB -> spectrum used for A-support constraint
    W_B = magB;                   % zB -> spectrum used for B-support constraint

    % Rows correspond to pixels OUTSIDE supports (in pupil plane)
    idxA_out = find(~suppA);
    idxB_out = find(~suppB);
    nA = numel(idxA_out);  nB = numel(idxB_out);

    % Anchor (fix global phase of ZB at a single overlap pixel)
    if isempty(opts.anchor)
        [rK,cK] = find(K);
        [~,mid] = min( (rK - median(rK)).^2 + (cK - median(cK)).^2 );
        anchor_idx = sub2ind([H,W], rK(mid), cK(mid));
    else
        anchor_idx = opts.anchor;
        assert(K(anchor_idx), 'opts.anchor must lie in the overlap mask K.');
    end
    anchor_mask = false(H,W); anchor_mask(anchor_idx) = true;

    % Sizes of A
    m = nA + nB + 1 + (opts.lambda>0)*nK;   % rows
    n = nK;                                 % cols (unknowns)

    % Embedding/extraction helpers
    toFull = @(zvec) embed_on_mask(zvec, K, H, W);   % n->H×W (zeros elsewhere)
    toVec  = @(X) X(K);                               % H×W -> n (on K)

    scale = H*W;   % adjoint scaling: (ifft2)^H = (1/N) fft2

    % ----- define A and A' as function handles (complex operators) -----
    Amul = @(zvec) local_A_times(zvec, H,W, idxA_out,idxB_out, ...
                                 W_A,W_B, anchor_mask, opts.anchorVal, ...
                                 opts.lambda, toFull);

    ATmul = @(yvec) local_AH_times(yvec, H,W, idxA_out,idxB_out, ...
                                   W_A,W_B, anchor_mask, opts.lambda, ...
                                   toVec, scale);

    % RHS b: zeros except the anchor row
    b = zeros(m,1);  b(nA+nB+1) = opts.anchorVal;

    % Package
    sys = struct('Amul',Amul, 'ATmul',ATmul, ...
                 'b',b, 'm',m, 'n',n, ...
                 'K',K, 'idxA_out',idxA_out, 'idxB_out',idxB_out);
end

% ====== internal helpers ======

function y = local_A_times(zvec, H,W, idxA_out,idxB_out, W_A,W_B, anchor_mask, anchor_val, lambda, toFull)
    ZB = toFull(zvec);                         % H×W, zeros outside K
    a = ifft2(W_A .* ZB);  yA = a(idxA_out);   % A-support violation
    b = ifft2(W_B .* ZB);  yB = b(idxB_out);   % B-support violation
    yAnchor = sum(ZB(anchor_mask),'all');      % picks that pixel
    y = [yA; yB; yAnchor];
    if lambda > 0, y = [y; lambda*zvec]; end   % Tikhonov rows
end

function x = local_AH_times(yvec, H,W, idxA_out,idxB_out, W_A,W_B, anchor_mask, lambda, toVec, scale)
    nA = numel(idxA_out); nB = numel(idxB_out);
    if lambda > 0
        yA = yvec(1:nA); yB = yvec(nA+1:nA+nB);
        yAnchor = yvec(nA+nB+1); yLam = yvec(nA+nB+2:end);
    else
        yA = yvec(1:nA); yB = yvec(nA+1:nA+nB);
        yAnchor = yvec(nA+nB+1);
    end

    Aout = zeros(H,W); Aout(idxA_out) = yA;
    Bout = zeros(H,W); Bout(idxB_out) = yB;

    % adjoint of a = ifft2(Q .* ZB)  is  (1/N) Q^* .* fft2(a_out)
    g1 = conj(W_A) .* (fft2(Aout)/scale);
    g2 = conj(W_B) .* (fft2(Bout)/scale);
    g3 = zeros(H,W); g3(anchor_mask) = yAnchor;

    G = g1 + g2 + g3;          % H×W
    x = toVec(G);              % n×1 (on K)
    if lambda > 0, x = x + lambda*yLam; end
end

function X = embed_on_mask(x, K, H, W)
    X = zeros(H,W); X(K) = x;
end
