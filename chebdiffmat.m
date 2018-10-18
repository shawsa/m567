function [x, DM] = chebdiffmat(N,M,a,b)
%CHEBDIFFMAT  Chebyshev pseudospectral differentiation matrices over [a,b]
%  The function [x, DM] =  chebdiffmat(N,M,a,b) computes the differentiation
%  matrices D1, D2, ..., DM on Chebyshev nodes spaced over the interval
%  [a,b].
%
%  Input:
%  N:        Size of differentiation matrix.
%  M:        Number of derivatives required (integer).
%  [a,b]:    Interval for the approximation.
%  Note:     0 < M <= N-1.
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
%
%  The code implements two strategies for enhanced accuracy suggested by
%          W.Don and S. Solomonoff in SIAM J. Sci. Comp. Vol. 6, pp.
%          1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric identities
%  to avoid the computation of differences x(k)-x(j) and (b) the use of the
%  "flipping trick" which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
%
%  Note added May 2003:  It may, in fact, be slightly better not to
%  implement the strategies (a) and (b).   Please consult the following
%  paper for details:
%          R. Baltensperger and M.R. Trummer. Spectral Differencing with a
%          Twist, by , to appear in SIAM J. Sci. Comp.
%
%  Original author: J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified
%  by JACW, May 2003.
%
%  Modified by G.B. Wright 2018 to scale to any interval, and a few minor
%  improvements for later versions of MATLAB.

% Generate the DMs over [-1,1].
[x,DM] = chebdm(N,M);

% Flip nodes so the are sorted from -1 to 1.
x = flipud(x);
% Change interval to [a,b].
x = 0.5*((x+1)*b-(x-1)*a);
temp = 1;
% Adjust the DMs by the appropriate scale factors.
for j=1:M
    temp = temp*(2/(b-a));
    DM(:,:,j) = rot90(DM(:,:,j),2)*temp;
end

end

function [x,DM] = chebdm(N,M)
L = logical(speye(N));                % Logical identity matrix.
ov = ones(N,1);                       % Vector of ones.

n1 = floor(N/2); n2  = ceil(N/2);     % Indices used for flipping trick.

k = (0:N-1).';                        % Compute theta vector.
th = k*pi/(N-1);

x = sin(pi*(N-1:-2:1-N)'/(2*(N-1)));  % Compute Chebyshev points.

T = repmat(th/2,1,N);
DX = 2*sin(T'+T).*sin(T'-T);               % Trigonometric identity.
DX = [DX(1:n1,:); -rot90(DX(1:n2,:),2)];   % Flipping trick.
DX(L) = ov;                                % Put 1's on the main diagonal of DX.

C = toeplitz((-1).^k);                     % C is the matrix with
C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2;      % entries c(k)/c(j)
C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;

Z = 1./DX;                           % Z contains entries 1/(x(k)-x(j))
Z(L) = zeros(N,1);                   % with zeros on the diagonal.

D = eye(N);                          % D contains diff. matrices.
DM = zeros(N,N,M);
for ell = 1:M
    D = ell*Z.*(C.*(diag(D)*ov.') - D);      % Off-diagonals
    D(L) = -D*ov;                            % Correct main diagonal of D
    DM(:,:,ell) = D;                         % Store current D in DM
end
end
