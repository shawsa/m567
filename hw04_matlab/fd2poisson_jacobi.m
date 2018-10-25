% Numerical approximation to Poisson's equation over the square [a,b]x[a,b] with
% Dirichlet boundary conditions.  Uses a uniform mesh with (n+2)x(n+2) total
% points (i.e, n interior grid points).
% Input:
%     ffun : the RHS of poisson equation (i.e. the Laplacian of u).
%     gfun : the boundary function representing the Dirichlet B.C.
%      a,b : the interval defining the square
%        m : m+2 is the number of points in either direction of the mesh.
% Ouput:
%        u : the numerical solution of Poisson equation at the mesh points.
%      x,y : the uniform mesh.
%
function [u,x,y] = fd2poissonsor(ffun,gfun,a,b,m)

tol = 10^-8;
max_iter = 10^4;

h = (b-a)/(m+1);   % Mesh spacing

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

r_tol = tol * norm(f);

rhs = h^2 * f;

% Adjust f for boundary terms
% f(:,1) = f(:,1) - ubw/h^2;             % West
% f(:,m) = f(:,m) - ube/h^2;             % East
% f(1,1:m) = f(1,1:m) - ubs(idx)/h^2;    % South
% f(m,1:m) = f(m,1:m) - ubn(idx)/h^2;    % North

u = ones(m,m);
u = [ubs;[ubw,u,ube];ubn];

%omega = 2/(1+sin(pi*h));
omega = 1;

for i=1:max_iter
    % Calculate Residual
    res = rhs + 4.* u(2:end-1, 2:end-1);
    res = res - u(1:end-2, 2:end-1);
    res = res - u(3:end, 2:end-1);
    res = res - u(2:end-1, 1:end-2);
    res = res - u(2:end-1, 3:end);
    %res = res - h^2 * f;
    if norm(res) <= r_tol
       break; 
    end
    u(2:end-1, 2:end-1) = u(2:end-1, 2:end-1) + omega*(-1/4)*res;
end
if i==max_iter
   warning('Max iterations reached in SOR.'); 
end
 
end



