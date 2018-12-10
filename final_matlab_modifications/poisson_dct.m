function [pcr] = poisson_dct(div, grid)
%POISSON_DCT Summary of this function goes here
%   Detailed explanation goes here
m = grid.m;
h = grid.dx;
p_hat = idct2(dct2(div,2),1);
cos_vec = cos( (0:m-1)  * pi/m);
denom = 2* (cos_vec' + cos_vec - 2);
denom(1,1) = 1;
u_hat = h^2 * p_hat./denom;
u_hat(1,1) = 0;
pcr = dct2(idct2(u_hat,2),1);
end

