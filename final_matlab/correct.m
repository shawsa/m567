% Correct the velocity field (ustar,vstar) accoriding the Helmholtz
% decomposition:
% u = ustar - dt*dp/dx
% v = vstar - dt*dp/dy
function [u,v] = correct(u,v,p,dt,grid)

% IDEAS FOR IMPROVEMENT:
% Better method:
% a. Use compact scheme for computing the gradient

% Compute x-component of the gradient at vertical cell-edges
px = (p(:,3:grid.mpx-1)-p(:,2:grid.mpx-2))/grid.dx;
% Correct u
u(:,grid.jux) = u(:,grid.jux) - dt*px;

% Compute y-component of the gradient at horizontal cell-edges
py = (p(3:grid.mpy-1,:)-p(2:grid.mpx-2,:))/grid.dy;
% Correct v
v(grid.jvy,:) = v(grid.jvy,:) - dt*py;

end