%% Set up the Lid-Driven Cavity problem.

% Driver for the Lid-Driven Cavity problem.  Comparisons of the results
% will be made with those reported in 
% Ghia, Ghia, & Shen. High-Re Solutions for Incomprossible Flow Using the
% Navier-Stokes Equations and a Multigrid Method.  Journal of Computational
% Physics, 48, 387-411 (1982).
% Which we abbreviate as GGS82

% Parameters (Reynolds number). Comparisons should be made for Re=100 and
% Re=1000.
Re = 100; 

% Contour levels for plotting the stream function to compare GGS82.
cntrs = [-1e-10 -1e-7 -1e-5 -1e-4 -1e-2 -3e-2 -5e-2 -7e-2 -9e-2 -0.1 -0.11 -0.115 -0.1175];
cntrs = [cntrs 1e-8 1e-7 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];

% Max CFL number to run at
cfl = 0.5; 

% Grid - only square grids are handled
L = 1;    % Unit box
m = 2^5;  % Number of cell-centers on the grid.
grid = setupGrid(L,m,cfl,Re);
dt    = grid.dt;

% Time parameters
time  = 0; 
tfinal = 40;
tstep = tfinal/dt;  % number of time-steps 

% Create the grid of values for plotting the stream function
x = linspace(0,L,grid.mpx-1); y = linspace(0,L,grid.mpy-1);
[xx,yy] = meshgrid(x,y);

% Initialize velocity and pressure, include ghost cells.
u = zeros(grid.muy,grid.mux); 
v = zeros(grid.mvy,grid.mvx); 
p = zeros(grid.mpy,grid.mpx);

% Initialize the stream function 
sf = zeros(grid.mpy-1,grid.mpx-1);
% Variables for measuring steady-state
sfold = sf;
uold = u;
vold = v;

%% Solve the Navier-Stokes equations
for k=1:tstep    
    % Iterate to get the ustar correct.  3 times is typically sufficient.
    for j = 1:3
        % Advection-diffusion step
        [ustar,vstar] = advectDiffuse(u,v,p,dt,grid,Re);
        % Pressure step
        [p,lgmult] = pressure(ustar,vstar,dt,grid);
    end
    
    % Correct the velocity field
    [u,v] = correct(ustar,vstar,p,dt,grid);
    
    % Update the time
    time = k*dt;

    % Print and plot diagnositics every 100 time-steps
    if (mod(k,100)==0)
        fprintf('time %f, step %i out of %i\n',time,k,tstep)

        figure(1);        
        sf = streamfunction(u,v,grid);
        contour(xx,yy,sf,cntrs);
        axis([0 1 0 1]), axis image; drawnow
        pause(0.1);
        
        % Measure how close we are to steady state by comparing the
        % relative two-norm of the difference between the new stream
        % function and the previous stream function and the previous
        % velocities
        convSf = norm(sf(:)-sfold(:),2)/norm(sf(:),2);
        convU = max(norm(u(:)-uold(:),inf),norm(u(:)-uold(:),inf));
        sfold = sf;
        uold = u;
        vold = v;
        
        % Compute the max-velocity
        ucntr = avgXEdgToCntr(u(2:end-1,:)); 
        vcntr = avgYEdgToCntr(v(:,2:end-1));
        maxvel = max(max(sqrt(ucntr(2:end-1,2:end-1).^2 + vcntr(2:end-1,2:end-1).^2)));
        
        % Compute how far off the corrected velocity is from satisfying the
        % correct boundary condition for the flow.
        errUBndry = max(abs(0.5*(u(end,2:end-1)+u(end-1,2:end-1))-1));
        
        % Print out some diagnostics
        fprintf('Max velocity=%1.4e, Error at boundary=%1.4e, Convergence=%1.4e, %1.4e\n',maxvel,errUBndry,convSf,convU);
    end
end

%%  Compare to the GGS82 data

% Generate the grid of values (with ghost points) for u-velocity
y = [-grid.dy/2 linspace(grid.dy,grid.L-grid.dy/2,grid.muy-2) 1+grid.dy/2];
x = linspace(0,grid.L,grid.mux);
[xx_u,yy_u] = meshgrid(x,y);

% Generate the grid of values (with ghost points) for v-velocity
x = [-grid.dx/2 linspace(grid.dx/2,grid.L-grid.dx/2,grid.mvx-2) 1+grid.dx/2];
y = linspace(0,grid.L,grid.mvy);
[xx_v,yy_v] = meshgrid(x,y);

% Locations to evaluate the u component of the velocity from GGS82
yiu = [1 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547 0];
xiu = 0.5 + 0*yiu;

% Locations to evaluate the v component of the velocity from GGS82
xiv = [1 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0];
yiv = 0.5 + 0*xiv;

% Determine the values for comparison to GS82
if Re == 100
    % Values of u at (xiu,yiu)
    uGGS82 = [1 0.84123 0.78871 0.73722 0.68717 0.23151 0.00332 -0.13641 -0.20581 -0.21090 -0.15662 -0.10150 -0.06434 -0.04775 -0.04192 -0.03717 0];
    % Values of v at (xiv,yiv)
    vGGS82 = [0 -0.05906 -0.07391 -0.08864 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.10890 0.10091 0.09233 0];
elseif Re == 1000
    % Values of u at (xiu,yiu)
    uGGS82 = [1 0.65928 0.57492 0.51117 0.46604 0.33304 0.18719 0.05702 -0.06080 -0.10648 -0.27805 -0.38289 -0.29370 -0.22220 -0.20196 -0.18109 0];
    % Values of v at (xiv,yiv)
    vGGS82 = [0 -0.21388 -0.27669 -0.33714 -0.39188 -0.51550 -0.42665 -0.31966 0.02526 0.32235 0.33075 0.37095 0.32627 0.30353 0.29012 0.27485 0];
else
    error('No GGS82 Data for this Re.  Data only avaialble for Re=100 and Re=1000');
end

% Interpolate the computed velocities to the positions the GGS82 data is
% given.
ui = interp2(xx_u,yy_u,u,xiu,yiu,'linear');
vi = interp2(xx_v,yy_v,v,xiv,yiv,'linear');

figure(2)
subplot(1,2,1), plot(yiu,uGGS82,'o-',yiu,ui,'x-')
legend('GGS82','Current'), xlabel('y'), ylabel('u'), title(sprintf('Re=%d',Re))
subplot(1,2,2), plot(xiv,vGGS82,'o-',xiv,vi,'x-'), 
legend('GGS82','Current'), xlabel('x'), ylabel('v'), title(sprintf('Re=%d',Re))
