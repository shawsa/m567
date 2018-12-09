% Computes the advective derivatives (u,v).grad(u) and (u,v).grad(v) in the
% Navier-Stokes equations.
%
% Scheme is based on first-order upwinding (donor-cell).  Values are
% assumed to reside on a MAC (marker-and-cell) or staggered grid. Advective
% derivatives are computed at the interior cell edges.
%
% Input:
%    u - vertical edge value to be advected
%    v - horizontal edge value to be advected
%    dt    - time-step
%    grid  - structure for the MAC grid
%
% It is assumed that all variables contain ghost cells.  The values at
% these cells will be set in this function.
function [uadv,vadv] = advectUW(u,v,dt,grid)

%
% Compute (u,v).grad(u) component first
%

% Average u to the cell centers and v to the corners
U = avgXEdgToCntr(u(grid.juy,:));
V = avgYEdgToCrnr(v(:,grid.jvx));

F = zeros(size(U));   % flux at vertical edges
G = zeros(size(V));   % flux at horizontal edges

[Ny,Nx] = size(F);

% First the vertical edges
id = (U > 0);
F = id.*u(2:Ny+1,1:Nx) + ~id.*u(2:Ny+1,2:Nx+1);

[Ny,Nx] = size(G);

% Second the horizontal edges
id = V > 0;
G = id.*u(1:Ny,2:Nx+1) + ~id.*u(2:Ny+1,2:Nx+1);

% Need v on the vertical edge to compute (u,v) dot grad (u)
V = avgYEdgToXEdg(v(:,grid.jvx));

uadv = (1/grid.dx)*u(grid.juy,grid.jux).*diff(F,1,2) + ...
       (1/grid.dy)*V.*diff(G,1,1);

%
% Compute (u,v).grad(v) component second
%

% Average v to the cell centers and u to the corners
V = avgYEdgToCntr(v(:,grid.jvx));
U = avgXEdgToCrnr(u(grid.juy,:));

F = zeros(size(U));   % flux at vertical edges
G = zeros(size(V));   % flux at horizontal edges

[Ny,Nx] = size(G);

% First the horizontal edges
id = (V > 0);
G = id.*v(1:Ny,2:Nx+1) + ~id.*v(2:Ny+1,2:Nx+1);

[Ny,Nx] = size(F);

% Second the vertical edges
id = (U > 0);
F = id.*v(2:Ny+1,1:Nx) + ~id.*v(2:Ny+1,2:Nx+1);

% Need u on the horizontal edge to compute (u,v) dot grad (v)
U = avgXEdgToYEdg(u(grid.juy,:));

vadv = (1/grid.dy)*v(grid.jvy,grid.jvx).*diff(G,1,1) + ...
       (1/grid.dx)*U.*diff(F,1,2);

end
