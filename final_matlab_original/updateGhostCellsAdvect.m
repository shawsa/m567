% Updates the ghost cells and boundary values for the advective step. These
% values are based on the lid driven cavity problem with unit speed on the
% top wall.  Specifically, the boundary conditions are
% u = 1 on top wall, u = 0 on all other walls,
% v = 0 on all walls
function [u,v] = updateGhostCellsAdvect(u,v,grid)

% u on vertical walls (no extrapolation necessary)
u(:,1) = 0; 
u(:,grid.mux) = 0;

% u on bottom wall
u(1,:) = -u(2,:);

% u on top wall
u(grid.muy,:) = 2 - u(grid.muy-1,:);

% v on horizontal walls (no extrapolation necessary)
v(1,:) = 0; 
v(grid.mvy,:) = 0;

% v on vertical walls
v(:,1) = -v(:,2);
v(:,grid.mvx) = -v(:,grid.mvx-1);

end



