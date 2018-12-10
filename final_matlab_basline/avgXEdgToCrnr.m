% Computes the averages of values at vertical edges of a MAC grid (i.e. the 
% edges where the x-direction velocity is defined) to the cell corners.
% avg = 0  corresponds to harmonic averaging.
% avg = 1 corresponds to arithmetic averaging.
function crnravg = avgXEdgToCrnr(edg,avg)

Ny = size(edg,1);

if nargin <= 3
    avg = 0;
end

if avg == 1    % harmonic averaging
   crnravg = 1./edg(1:Ny-1,:) + 1./edg(2:Ny,:);
   crnravg = 2./crnravg;
else           % arithmetic averaging
   crnravg = (edg(1:Ny-1,:) + edg(2:Ny,:))/2;
end

