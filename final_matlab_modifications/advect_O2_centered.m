function [uadv,vadv] = advect_O2_centered(u,v,grid)
%ADVECT_O2_CENTERED Summary of this function goes here
%   Detailed explanation goes here

%update ghost nodes
[u,v] = updateGhostCellsAdvect(u,v,grid);

h = grid.dx;
% v at us
V = avgYEdgToXEdg(v(:,grid.jvx));
% du/dx at us
F = (u(grid.juy,grid.jux+1) - u(grid.juy,grid.jux-1))/(2*h);
% du/dy at us
G = (u(grid.juy+1,grid.jux) - u(grid.juy-1,grid.jux))/(2*h);
uadv = u(grid.juy,grid.jux).*F + V.*G;

% u at vs
U = avgXEdgToYEdg(u(grid.juy,:));
% dv/dx at us
F = (v(grid.jvy,grid.jvx+1) - v(grid.jvy,grid.jvx-1))/(2*h);
% dv/dy at us
G = (v(grid.jvy+1,grid.jvx) - v(grid.jvy-1,grid.jvx))/(2*h);
vadv = U.*F + v(grid.jvy,grid.jvx).*G;

end

