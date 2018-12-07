% Computes the pressure by solving the pressure Poisson equation
% Laplacian(p) = div(u,v)/dt
% with Neumman boundary conditions.
%
% Current implementation uses a direct solver.

 function [p,lagmult] = pressure(ustar,vstar,dt,grid)

div = divergence(ustar,vstar,grid)/dt;

%   IDEAS FOR IMPROVEMENT
%   Speed:
%   a. Use a fast cosine transform solver
%   b. Use an iterative method like multigrid
%   c. Pre-compute LU factorization of laplacian matrix and use it
%      each time-step to solve system (not the best improvement).
%
%   Better method:
%   a.  Use 4th order compact scheme

[mx,my] = size(div);

% Solve for the pressure (in column re-order form)
pcr = grid.Lp\[div(:);0];

% The last entry of p is the value of the Lagrange multiplier and gives an
% idea of how far off the divergence is from satisfying the compatibility
% condition that its integral be zero.
lagmult = pcr(end);

% Make pcr into a matrix
p = zeros(grid.mpy,grid.mpx);
p(grid.jpy,grid.jpx) = reshape(pcr(1:end-1),my,mx);

% Fill in the ghost values (using Neumann BC)
p(1,:) = p(2,:); p(:,1) = p(:,2);
p(grid.mpy,:) = p(grid.mpy-1,:); p(:,grid.mpx) = p(:,grid.mpx-1);

end
