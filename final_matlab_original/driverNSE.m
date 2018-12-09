%% Set up the Lid-Driven Cavity problem.

% Driver for the Lid-Driven Cavity problem.  Comparisons of the results
% will be made with those reported in 
% Ghia, Ghia, & Shen. High-Re Solutions for Incomprossible Flow Using the
% Navier-Stokes Equations and a Multigrid Method.  Journal of Computational
% Physics, 48, 387-411 (1982).
% Which we abbreviate as GGS82

% Reynolds number.  Re=100, 400, 1000 are the options for comparisons.
Re = 100; 

% Contour levels for plotting the stream function to compare GGS82.
cntrs = [-1e-10 -1e-7 -1e-5 -1e-4 -1e-2 -3e-2 -5e-2 -7e-2 -9e-2 -0.1 -0.11 -0.115 -0.1175];
cntrs = [cntrs 1e-8 1e-7 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];

% Max CFL number to run at
cfl = 0.5; 

% Grid - only square grids are handled
L = 1;    % Unit box
m = 2^6;  % Number of cell-centers on the grid.
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

    % Advection-diffusion step
    [ustar,vstar] = advectDiffuse(u,v,p,dt,grid,Re);
    % Pressure step
    [p,lgmult] = pressure(ustar,vstar,dt,grid);
    
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

%%  Compare to the GGS82 and MSA09 results

% Generate the grid of values (with ghost points) for u-velocity
y = [-grid.dy/2 linspace(grid.dy/2,grid.L-grid.dy/2,grid.muy-2) 1+grid.dy/2];
x = linspace(0,grid.L,grid.mux);
[xx_u,yy_u] = meshgrid(x,y);

% Generate the grid of values (with ghost points) for v-velocity
x = [-grid.dx/2 linspace(grid.dx/2,grid.L-grid.dx/2,grid.mvx-2) 1+grid.dx/2];
y = linspace(0,grid.L,grid.mvy);
[xx_v,yy_v] = meshgrid(x,y);

% Locations to evaluate the u component of the velocity from GGS82
yiu = linspace(0,1,129);
idy = [129 126 125 124 123 110 95 80 65 59 37 23 14 10 9 8 1];
yiu = yiu(idy);
xiu = 0.5 + 0*yiu;

% Locations to evaluate the u component of the velocity from MSA09
yiu2 = linspace(0,1,17); yiu2 = yiu2(2:end-1);
xiu2 = 0.5 + 0*yiu2;

% Locations to evaluate the v component of the velocity from GGS82
xiv = linspace(0,1,129);
idx = [129 125 124 123 122 117 111 104 65 31 30 21 13 11 10 9 1];
xiv = xiv(idx);
yiv = 0.5 + 0*xiv;

% Locations to evaluate the v component of the velocity from MSA09
xiv2 = linspace(0,1,17); xiv2 = xiv2(2:end-1);
yiv2 = 0.5 + 0*xiv2;

% Determine the values for comparison to GS82
if Re == 100
    % Values of u at (xiu,yiu)
    uGGS82 = [1 0.84123 0.78871 0.73722 0.68717 0.23151 0.00332 -0.13641 -0.20581 -0.21090 -0.15662 -0.10150 -0.06434 -0.04775 -0.04192 -0.03717 0];
    uMSA09 = [-0.041974991 -0.077125399 -0.109816214 -0.141930064 -0.172712391 -0.198470859 -0.212962392 -0.209149142 -0.182080595 -0.131256301 -0.060245594 0.027874448 0.140425325 0.310557090 0.597466694];
    % Values of v at (xiv,yiv)
    vGGS82 = [0 -0.05906 -0.07391 -0.08864 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.10890 0.10091 0.09233 0];
    vMSA09 = [0.094807616 0.149243000 0.174342933 0.179243328 0.169132064 0.145730201 0.108775865 0.057536559 -0.007748504 -0.084066715 -0.163010143 -0.227827313 -0.253768577 -0.218690812 -0.123318170];
elseif Re == 400
    % Values of u at (xiu,yiu)
    uGGS82 = [1 0.75837 0.68439 0.61756 0.55892 0.29093 0.16256 0.02135 -0.11477 -0.17119 -0.32726 -0.24299 -0.14612 -0.10338 -0.09266 -0.08186 0];
    uMSA09 = [-0.092599260 -0.178748051 -0.263917200 -0.321229080 -0.320251090 -0.266306350 -0.190730560 -0.115053628 -0.042568947 0.030243020 0.105456010 0.181306850 0.252203840 0.316829690 0.469580199];
    % Values of v at (xiv,yiv)
    vGGS82 = [0 -0.12146 -0.15663 -0.19254 -0.22847 -0.23827 -0.44993 -0.38598 0.05186 0.30174 0.30203 0.28124 0.22965 0.20920 0.19713 0.18360 0];
    vMSA09 = [0.185132290 0.262251260 0.297479230 0.300960030 0.268310960 0.206571390 0.130571694 0.052058082 -0.024714514 -0.100884164 -0.182109238 -0.280990219 -0.400042350 -0.449011850 -0.270354943];
elseif Re == 1000
    % Values of u at (xiu,yiu)
    uGGS82 = [1 0.65928 0.57492 0.51117 0.46604 0.33304 0.18719 0.05702 -0.06080 -0.10648 -0.27805 -0.38289 -0.29370 -0.22220 -0.20196 -0.18109 0];
    uMSA09 = [-0.202330048 -0.347845100 -0.384409400 -0.318946100 -0.245693700 -0.183732100 -0.123410460 -0.062056130 0.000561800 0.065248742 0.133572570 0.207914610 0.288442400 0.362545400 0.422932100];
    % Values of v at (xiv,yiv)
    vGGS82 = [0 -0.21388 -0.27669 -0.33714 -0.39188 -0.51550 -0.42665 -0.31966 0.02526 0.32235 0.33075 0.37095 0.32627 0.30353 0.29012 0.27485 0];
    vMSA09 = [0.28070570 0.36504180 0.36785270 0.30710428 0.23126839 0.16056422 0.09296931 0.02579946 -0.04184068 -0.11079830 -0.18167970 -0.25338150 -0.33156670 -0.46777560 -0.45615254];
else
    error('No GGS82 Data for this Re.  Data only avaialble for Re=100, 400, and 1000');
end

% Interpolate the computed velocities to the positions the GGS82 data is
% given.
ui = interp2(xx_u,yy_u,u,xiu,yiu,'linear');
vi = interp2(xx_v,yy_v,v,xiv,yiv,'linear');
ui2 = interp2(xx_u,yy_u,u,xiu2,yiu2,'linear');
vi2 = interp2(xx_v,yy_v,v,xiv2,yiv2,'linear');

% Plot the reference solution with the computed solution for the u,v
% components.
figure(2)
[tempy,id] = sort([yiu yiu2]);
tempu = [ui ui2]; tempu = tempu(id);
[tempx,id] = sort([xiv xiv2]);
tempv = [vi vi2]; tempv = tempv(id);
subplot(1,2,1), plot(yiu,uGGS82,'ko',yiu2,uMSA09,'b^',yy_u(:,m/2+1),u(:,m/2+1),'r-')
legend('GGS82','MSA09','Current'), xlabel('y'), ylabel('u'), title(sprintf('Re=%d',Re))
subplot(1,2,2), plot(xiv,vGGS82,'ko',xiv2,vMSA09,'b^',xx_v(m/2+1,:),v(m/2+1,:),'r-')
legend('GGS82','MSA09','Current'), xlabel('x'), ylabel('v'), title(sprintf('Re=%d',Re))

% Plot the error in the computed solution as measured against the MSA09
% results, which are more accurate that GGS82.
figure(3)
subplot(1,2,1), plot(yiu2,uMSA09-ui2,'x-'), hold on
subplot(1,2,2), plot(xiv2,vMSA09-vi2,'x-'), hold on
