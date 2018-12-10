%
% Solves the advection-diffusion step using upwinding for the advection
% term and Crank-Nicolson for the diffusion term:
%
% (I - alpha*laplacian)ustar = (I + alpha*laplacian)*u - dt*((u,v).grad(u))
% (I - alpha*laplacian)vstar = (I + alpha*laplacian)*v - dt*((u,v).grad(v))
%
% where alpha = 0.5*dt/Re.  Laplacian(u) is computed using second-order
% centered finite-differences, while (u,v).grad(u) and (u,v).grad(v) are
% computed using upwinding.
function [ustar,vstar] = advectDiffuse(u,v,p,dt,grid,Re)

%
% Compute advective derivatives (u,v).grad(u) and (u,v).grad(v)
%

% This uses first-order upwinding. You will get a big improvement by
% changing this to something that is more accurate.
%   IDEAS FOR IMPROVEMENT
%   a. Use Lax-Wendroff to compute uadv,vadv at the t+dt/2
%   b. Use 2nd-order centered-differences
%   c. Use 4th-order compact differencing to compute grad(u) and grad(v)
%   d. Use either (b) and (c) and an extrapolation routine to approximate 
%      the advective derivatives at t + dt.  This will require storing the 
%      values of advDerivU and advDerivV from the previous time-step.

%[advDerivU,advDerivV] = advectUW(u,v,dt,grid);
[advDerivU,advDerivV] = advect_O2_centered(u,v,grid);


%
% Do the diffusion part implicitly with Crank-Nicolson.
%

%   The code below uses a direct method for solving the system.  You will
%   get an improvement in speed by switching to something better.  
%   IDEAS FOR IMPROVEMENT
%   a. Use ADI
%   b. Use a fast direct solver based on the fast sine transform.
%   c. Use multigrid
%   d. Use conjugate gradient method (look at the pcg function in MATLAB).
%   e. Pre-compute Cholesky factorization of the matrices and use them at
%      each time-step to solve system (not the best improvement).
%
%   You can also get better accuracy results by switching to better methods.  
%   Some ideas would be
%   a. Use the TR-BDF2 method to advance the system in time.
%   b. Include the gradient of the pressure in the time update
%   c. Use a 4th-order compact scheme for the Laplacian.

% Update the boundary conditions and ghost cell values for ustar.
[ustar,vstar] = updateGhostCellsDiffuse(u,v,p,dt,grid,0);

% Make code slightly easier to read by defining local variables
dx = grid.dx;
dy = grid.dy;

alpha = 0.5*dt/Re;

%
% Solve u equation first
% 

% Number of interior points in y and x directions
mx = grid.mux-2;
my = grid.muy-2;

% Identity matrix in the Crank-Nicolson method
I = speye(mx*my);

% Compute right hand-side (I+L)*u on the grid.
f = u(grid.juy,grid.jux);
f = f + alpha*laplacianu(u,grid) - dt*advDerivU;

% Incorporate boundary conditions from ustar
f(1,:) = f(1,:) + (alpha/dy^2)*ustar(1,grid.jux);
f(my,:) = f(my,:) + (alpha/dy^2)*ustar(grid.muy,grid.jux);
f(:,1) = f(:,1) + (alpha/dx^2)*ustar(grid.juy,1);
f(:,mx) = f(:,mx) + (alpha/dx^2)*ustar(grid.juy,grid.mux);

% Solve using direct solver
ustar(grid.juy,grid.jux) = reshape((I - alpha*grid.Lu)\f(:),my,mx);

%
% Solve v equation second
% 

% Number of interior points in y and x directions
mx = grid.mvx-2;
my = grid.mvy-2;

% Identity matrix in the Crank-Nicolson method
I = speye(mx*my);

% Compute right hand-side (I+L)*v on the grid.
f = v(grid.jvy,grid.jvx);
f = f + alpha*laplacianv(v,grid) - dt*advDerivV;

% Incorporate boundary conditions from vstar
f(1,:) = f(1,:) + (alpha/dy^2)*vstar(1,grid.jvx);
f(my,:) = f(my,:) + (alpha/dy^2)*vstar(grid.mvy,grid.jvx);
f(:,1) = f(:,1) + (alpha/dx^2)*vstar(grid.jvy,1);
f(:,mx) = f(:,mx) + (alpha/dx^2)*vstar(grid.jvy,grid.mvx);

% Solve using direct solver
vstar(grid.jvy,grid.jvx) = reshape((I - alpha*grid.Lv)\f(:),my,mx);

[ustar,vstar] = updateGhostCellsDiffuse(ustar,vstar,p,dt,grid,1);

end

function lap = laplacianu(u,grid)

lap = (-2*u(grid.juy,grid.jux) + u(grid.juy+1,grid.jux) + u(grid.juy-1,grid.jux))/grid.dy^2 + ...
      (-2*u(grid.juy,grid.jux) + u(grid.juy,grid.jux+1) + u(grid.juy,grid.jux-1))/grid.dx^2;
end

function lap = laplacianv(v,grid)

lap = (-2*v(grid.jvy,grid.jvx) + v(grid.jvy+1,grid.jvx) + v(grid.jvy-1,grid.jvx))/grid.dy^2 + ...
      (-2*v(grid.jvy,grid.jvx) + v(grid.jvy,grid.jvx+1) + v(grid.jvy,grid.jvx-1))/grid.dx^2;
end


