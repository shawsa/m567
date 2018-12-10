% Computes the stream function for the velocity field (u,v) by solving the 
% Poisson equation
% laplacian(psi) = -curl((u,v)) = u_y - v_x 
% with homogeneous Dirichlet boundary conditions.
%
% Values of the stream function live at the cell orders of the MAC grid.
function stream = streamfunction(u,v,grid)

%   IDEAS FOR IMPROVEMENT
%   Speed:
%   a. Use a fast sine transform solver to solve the Poisson equation.
%   b. Use multigrid

% Vorticity lives at the corners of the MAC grid.
vort = vorticity(u,v,grid);

% Use homogeneous boundary conditions.

mx = size(vort,2)-2;  % Number of interior points
my = size(vort,1)-2;

% Construct the Laplacian operator
D2x = toeplitz(spdiags([-2 1],[0 1],1,mx))/grid.dx^2;
D2y = toeplitz(spdiags([-2 1],[0 1],1,my))/grid.dy^2;

% Laplacian
L = kron(D2x,speye(my)) + kron(speye(mx),D2y);

stream = zeros(size(vort));
stream(2:end-1,2:end-1) = reshape(L\reshape(vort(2:end-1,2:end-1),[],1),mx,my);

end
