% Updates the ghost cells and boundary values for the diffusion step. These
% values are based on the lid driven cavity problem with unit speed on the
% top wall.  The flag extrap is used for controlling wether the interior
% points should be extrapolated, with 1 meaning extrap, 0 meaning no
% extrap.  This is helpful for implementing BC in direct implicit solvers
% where one may not want to extrapolate the interior points.
function [u,v] = updateGhostCellsDiffuse(u,v,p,dt,grid,extrap)

if nargin == 5
    extrap = 1;
end

%
% Do u-component first
%

% Compute the dp/dx along the bottom and top walls at each cell edge
pxb = diff(p(2,:))/grid.dx;
pxt = diff(p(grid.mpy-1,:))/grid.dx;

% u = 0 on vertical boundaries
u(:,1) = 0;
u(:,grid.mux) = 0;

% Linearly extrapolate to get the top and bottom ghost cells using the fact
% that we know what u is these boundaries.
u(1,:) = -extrap*u(2,:) + 2*dt*pxb;
u(grid.muy,:) = -extrap*u(grid.muy-1,:) + 2*(1 + dt*pxt);

% Average to the corners
u(1,1) = 0.5*u(1,1);
u(1,grid.mux) = 0.5*u(1,grid.mux);
u(grid.muy,1) = 0.5*u(grid.muy,1);
u(grid.muy,grid.mux) = 0.5*u(grid.muy,grid.mux);

%
% Do v-component second
%

% Compute the dp/dy along the vertical walls at each cell edge
pyw = diff(p(:,2))/grid.dy;
pye = diff(p(:,grid.mpx-1))/grid.dy;

% v = 0 on horizontal boundaries, so this gives that the average of ghost
% value and edge value are zero on the boundary.
v(1,:) = 0;
v(grid.mvy,:) = 0;

% Linearly extrapolate to get the left and right ghost cells using the fact
% that we know what v is on these boundaries.
v(:,1) = -extrap*v(:,2) + 2*dt*pyw;
v(:,grid.mvx) = -extrap*v(:,grid.mvx-1) + 2*dt*pye;

% Average to the corners
v(1,1) = 0.5*v(1,1);
v(1,grid.mvx) = 0.5*v(1,grid.mvx);
v(grid.mvy,1) = 0.5*v(grid.mvy,1);
v(grid.mvy,grid.mvx) = 0.5*v(grid.mvy,grid.mvx);


end



