function Fvec = dftmtx2d(sz_x)
%refer to sec 2.6 of Linear Algebra, Gilbert Strang

% image vector needs to be rehsaped like so:
% x_vec = reshape(x,[],1);
% Apply fft2 like so:
% f_x_vec = Fvec * x_vec; %input,output are vectors
% convert like so:
% f_x = reshape(f_x_vec,size(x,1),size(x,2));

m = sz_x(1);
n = sz_x(2);

Fvec = zeros(m*n,m*n);
for i=1:m*n
    inpvec = zeros(m*n,1);
    inpvec(i) = 1;
    inpvec = reshape(inpvec,m,n);
    %finding out one col of the transformation. since i/p and o/p basis are
    %canonical, nothing further is required
    opvec = fft2(inpvec); 
    opvec = reshape(opvec,[],1);
    Fvec(:,i) = opvec;
end