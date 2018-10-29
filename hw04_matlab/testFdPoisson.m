% Script for testing fd2poisson over the square [a,b]x[a,b]
a = 0; 
b = 1;
k = 6;
m = 2^k-1;  % Number of interior grid points in one direction
h = (b-a)/(m+1);

% % homework problem
% % Laplacian(u) = f
% f = @(x,y) 10*pi^2*(1+cos(4*pi*(x+2*y))-2*sin(2*pi*(x+2*y))).*exp(sin(2*pi*(x+2*y)));  
% % u = g on Boundary
% g = @(x,y) exp(sin(2*pi*(x+2*y)));       

% worksheet problem
% Laplacian(u) = f
f = @(x,y) -5*pi^2*sin(pi*x).*cos(2*pi*y);
% u = g on Boundary
g = @(x,y) sin(pi*x).*cos(2*pi*y);

% Exact solution is g.
uexact = @(x,y) g(x,y);                     

% % Compute and time the solution
% tic
% [u,x,y] = fd2poisson(f,g,a,b,m);
% gedirect = toc;
% fprintf('GE takes \t\t%d s\n', gedirect);
% 
% % Compute and time the solution
% tic
% [u,x,y] = fd2poissonsp(f,g,a,b,m);
% gesparse = toc;
% fprintf('Sparse GE takes \t%d s\n', gesparse);


% Compute and time the solution
tic
[u,x,y] = fd2poissonsor(f,g,a,b,m);
sor = toc;
fprintf('SOR takes \t\t%d s\n', sor);

% % Compute and time the solution
% tic
% [u,x,y] = fd2poissondst(f,g,a,b,m);
% dstTime = toc;
% fprintf('DST takes \t\t%d s\n', dstTime);
% 
% % Compute and time the solution
% tic
% [u,x,y] = fd2poissonmg(f,g,a,b,m);
% mgTime = toc;
% fprintf('Multigrid takes \t%d s\n', mgTime);



fprintf('\n');
%% Test compact
a = 0; 
b = 1;
k = 5;
m = 2^k-1;  % Number of interior grid points in one direction
h = (b-a)/(m+1);

% Laplacian(u) = f
f = @(x,y) 10*pi^2*(1+cos(4*pi*(x+2*y))-2*sin(2*pi*(x+2*y))).*exp(sin(2*pi*(x+2*y)));  
% u = g on Boundary
g = @(x,y) exp(sin(2*pi*(x+2*y)));            

% Exact solution is g.
uexact = @(x,y) g(x,y);                     

% % Compute and time the solution
[u,x,y] = fd2poisson_compact(f,g,a,b,m);

%% Plot solution
%my_surf = surf(x,y,u)
%set(my_surf, 'linestyle','none');
figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
surf(x,y,u,'linestyle','none'), xlabel('x'), ylabel('y'), zlabel('u(x,y)'),
title(strcat('Numerical Solution to Poisson Equation, h=',num2str(h)));

% Plot error
figure, set(gcf,'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]), 
surf(x,y,u-uexact(x,y), 'linestyle','none'),xlabel('x'),ylabel('y'), zlabel('Error'), 
title(strcat('Error, h=',num2str(h)));
%light
%shading interp

%% Timings

runs = 10;
kmin=3;
kmax=10;
ge_cutoff = 6;

% parameters
a = 0; 
b = 1;

% Laplacian(u) = f
f = @(x,y) 10*pi^2*(1+cos(4*pi*(x+2*y))-2*sin(2*pi*(x+2*y))).*exp(sin(2*pi*(x+2*y)));  
% u = g on Boundary
g = @(x,y) exp(sin(2*pi*(x+2*y)));            
% Exact solution is g.
uexact = @(x,y) g(x,y);                     


direct_times = zeros(ge_cutoff-kmin+1,runs);
sor_times = zeros(kmax-kmin+1,runs);
sparse_times = zeros(kmax-kmin+1,runs);
dst_times = zeros(kmax-kmin+1,runs);
mg_times = zeros(kmax-kmin+1,runs);

for run=(1:runs)
    for k=(kmin:kmax)
        fprintf('run %d, k=%d started...', run, k);
        m = 2^k-1;
        h = (b-a)/(m+1);

        if k <= ge_cutoff
            tic
            [u,x,y] = fd2poisson(f,g,a,b,m);
            direct_times(k-kmin+1,run) = toc;
        end

        tic
        [u,x,y] = fd2poissonsor(f,g,a,b,m);
        sor_times(k-kmin+1,run) = toc;

        tic
        [u,x,y] = fd2poissonsp(f,g,a,b,m);
        sparse_times(k-kmin+1,run) = toc;



        % Compute and time the solution
        tic
        [u,x,y] = fd2poissondst(f,g,a,b,m);
        dst_times(k-kmin+1,run) = toc;
        % Compute and time the solution
        tic
        [u,x,y] = fd2poissonmg(f,g,a,b,m);
        mg_times(k-kmin+1,run) = toc;
        
        fprintf('complete\n');
    end
end
%%

direct_times_ave = mean(direct_times,2);
sor_times_ave = mean(sor_times,2);
sparse_times_ave = mean(sparse_times,2);
dst_times_ave = mean(dst_times,2);
mg_times_ave = mean(mg_times,2);

%cubic regression on direct method
p = polyfit((2.^(kmin:ge_cutoff)'-1).^2, direct_times_ave, 3);
ge_times_extrapolated = polyval(p, (2.^(ge_cutoff+1:kmax)'-1).^2);
direct_times_ave_all = [direct_times_ave; ge_times_extrapolated]

table = [(kmin:kmax)' (2.^(kmin:kmax)'-1).^2 direct_times_ave_all sor_times_ave sparse_times_ave dst_times_ave mg_times_ave];
save('timings.mat', 'table');


