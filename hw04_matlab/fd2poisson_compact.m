function [u,x,y] = fd2poisson_compact(ffun,gfun,a,b,m)
h = (b-a)/(m+1);
[x,y] = meshgrid(a:h:b);   % Uniform mesh, including boundary points.
idx = 2:m+1;
idy = 2:m+1;
% Compute boundary terms, south, north, east, west
ubs = feval(gfun,x(1,1:m+2),y(1,1:m+2));     % Include corners
ubn = feval(gfun,x(m+2,1:m+2),y(m+2,1:m+2)); % Include corners
ube = feval(gfun,x(idy,m+2),y(idy,m+2));     % No corners
ubw = feval(gfun,x(idy,1),y(idy,1));         % No corners 

% Evaluate the RHS of Poisson's equation at the interior points.
f = feval(ffun,x(idy,idx),y(idy,idx));
f = reshape(f,m*m,1);

% generate left hand side matrix
main_diag = ones(m,1)*[4 -20 4];
sub_diag = ones(m,1)*[1 4 1];
sup_eye = spdiags(ones(m), 1, m, m);

lhs = kron(speye(m), spdiags(main_diag, [-1 0 1], m,m));
lhs = lhs + kron(sup_eye, spdiags(sub_diag, [-1 0 1], m,m));
lhs = lhs + kron(sup_eye', spdiags(sub_diag, [-1 0 1], m,m));
lhs = 1/(6*h^2) * lhs;

% Generate right hand side matrix
z = ones(m,1)*[1 4 1];
F2 = spdiags(z, [-1 0 1], m, m);
F2x = kron(F2,speye(m));
F2y = kron(speye(m),F2);
rhs = F2x+F2y;

rhs = 1/12 * rhs*f;

%constructing this will simplify the boundary conditions
u = [ubs;[ubw,zeros(m,m),ube];ubn];
% Adjust rhs for boundary terms
boundary_adjust = zeros(m,m);
boundary_adjust(1,:) = boundary_adjust(1,:) + 4*u(1,2:end-1);
boundary_adjust(end,:) = boundary_adjust(end,:) + 4*u(end,2:end-1);
boundary_adjust(:,1) = boundary_adjust(:,1) + 4*u(2:end-1,1);
boundary_adjust(:,end) = boundary_adjust(:,end) + 4*u(2:end-1,end);
% boundary_adjust(:,1) = - ubw;             % West
% boundary_adjust(:,m) = - ube;             % East
% boundary_adjust(1,1:m) = - ubs(idx);    % South
% boundary_adjust(m,1:m) = - ubn(idx);    % North
boundary_adjust = -1/(6*h^2) * boundary_adjust;

% include boundary for f
boundary_adjust(1,:) = boundary_adjust(1,:) + 1/12*feval(ffun,x(1, idx),y(1,idy));
boundary_adjust(end,:) = boundary_adjust(end,:) + 1/12*feval(ffun,x(end, idx),y(end,idy));
boundary_adjust(:,1) = boundary_adjust(:,1) + 1/12*feval(ffun,x(idx, 1),y(idx,1));
boundary_adjust(:,end) = boundary_adjust(:,end) + 1/12*feval(ffun,x(idx, end),y(idx,end));

boundary_adjust = reshape(boundary_adjust, m^2, 1);

rhs = rhs + boundary_adjust;

% Solve the system
u = lhs\rhs;

% Convert u from a column vector to a matrix to make it easier to work with
% for plotting.
u = reshape(u,m,m);

% Append on to u the boundary values from the Dirichlet condition.
u = [ubs;[ubw,u,ube];ubn];
 
end



