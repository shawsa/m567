N = 10;
M = 2;

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