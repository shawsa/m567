% Sets up the Marker-and-Cell (MAC) grid for the Navier-Stokes equation
% solver with Dirichlet bounddary conditions.  Grid is assumed to be
% square with length (L) and m interior cell centers.
% Arrangement for variables on the MAC grid is
% 
%    -----v----- -----v----- 
%   |           |           |    u = horizontal velocity component
%   |           |           |        (vertical cell edges)
%   u     p     u     p     u    
%   |           |           |    v = vertical velocity component
%   |           |           |        (horizontal cell edges)
%    -----v----- -----v-----     
%   |           |           |    p = pressure (cell centers)
%   |           |           |
%   u     p     u     p     u 
%   |           |           |
%   |           |           |
%    -----v----- -----v----- 
%
% The grid for each variable is assumed to contain ghost cells where that
% variable does not have a given boundary value.  This gives the following
% dimensions for the grid variables (y,x)
% u : (m+2)-by-(m+1)  --> u(1,:) and u(m+2,:) are ghost vals
% v : (m+1)-by-(m+2)  --> v(:,1) and v(:,m+2) are ghost vals
% p : (m+2)-by-(m+2)  --> p(:,1), p(:,m+2), p(:,1), p(:,m+2) are ghost vals
%
function grid = setupGrid(L,m,cfl,Re)

% Total number of the variables (including ghost cells) in each direction
grid.mux = m+1;  % Total number of u velocity values x-direction
grid.muy = m+2;  % Total number of u velocity values y-direction
grid.mvx = m+2;  % Total number of v velocity values x-direction
grid.mvy = m+1;  % Total number of v velocity values y-direction
grid.mpx = m+2;  % Total number of p values x-direction
grid.mpy = m+2;  % Total number of p values y-direction

% Setup the interior indicies for the Dirichlet problem
grid.jux = 2:grid.mux-1; % u velocity unknowns x-direction
grid.juy = 2:grid.muy-1; % u velocity unknowns y-direction

grid.jvx = 2:grid.mvx-1; % v velocity unknowns x-direction
grid.jvy = 2:grid.mvy-1; % v velocity unknowns y-direction

grid.jpx = 2:grid.mpx-1; % p pressure unknowns x-direction
grid.jpy = 2:grid.mpy-1; % p pressure unknowns y-direction

grid.L = L;
grid.dx = L/m;
grid.dy = L/m;
grid.dt = cfl*min(grid.dx,grid.dy);

% Laplacian for the pressure step (interior nodes).
mx = grid.mpx-2;
my = grid.mpy-2;
D2x = toeplitz(spdiags([-2 1],[0 1],1,mx))/grid.dx^2;
D2y = toeplitz(spdiags([-2 1],[0 1],1,my))/grid.dy^2;
% Modify for Neumann boundary conditions
D2x(1,1) = -D2x(1,2); D2x(mx,mx) = -D2x(mx,mx-1);
D2y(1,1) = -D2y(1,2); D2y(my,my) = -D2y(my,my-1);
grid.Lp = kron(D2x,speye(my)) + kron(speye(mx),D2y);
% Add Lagrange multiplier to Laplacian to make the system invertible. This
% is necessary because the pressure only unique up to an additive constant.
grid.Lp = [[grid.Lp ones(mx*my,1)];[ones(1,mx*my) 0]];

% Laplacian for the diffusion step: u variable

% Number of interior points in y and x directions
mx = grid.mux-2;
my = grid.muy-2;
D2x = (1/grid.dx^2)*toeplitz(spdiags([-2 1],[0 1],1,mx));
D2y = (1/grid.dy^2)*toeplitz(spdiags([-2 1],[0 1],1,my));
D2y(1,1) = -3*(1/grid.dy^2);  % Account for boundary extrapolation.
D2y(my,my) = -3*(1/grid.dy^2);
grid.Lu = kron(D2x,speye(my)) + kron(speye(mx),D2y);

% Laplacian for the diffusion step: v variable

% Number of interior points in y and x directions
mx = grid.mvx-2;
my = grid.mvy-2;
D2x = (1/grid.dx^2)*toeplitz(spdiags([-2 1],[0 1],1,mx));
D2y = (1/grid.dy^2)*toeplitz(spdiags([-2 1],[0 1],1,my));
D2x(1,1) = -3*(1/grid.dx^2);  % Account for boundary extrapolation.
D2x(mx,mx) = -3*(1/grid.dx^2);
grid.Lv = kron(D2x,speye(my)) + kron(speye(mx),D2y);

end
