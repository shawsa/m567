% Computes the averages of values at horizontal edges of a MAC grid (i.e. the 
% edges where the y-direction velocity is defined) to the cell corners.
% avg = 0  corresponds to harmonic averaging.
% avg = 1 corresponds to arithmetic averaging.
function crnravg = avgYEdgToCrnr(edg,avg)

Nx = size(edg,2);

if nargin <= 3
    avg = 0;
end

if avg == 1    % harmonic averaging
   crnravg = 1./edg(:,1:Nx-1) + 1./edg(:,2:Nx);
   crnravg = 2./crnravg;
else           % arithmetic averaging
   crnravg = (edg(:,1:Nx-1) + edg(:,2:Nx))/2;
end

