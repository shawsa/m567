% Computes the averages of values at vertical edges of a MAC grid (i.e. the 
% edges where the x-direction velocity is defined) to the cell centers.
% avg = 0  corresponds to harmonic averaging.
% avg = 1 corresponds to arithmetic averaging.
function cntravg = avgXEdgToCntr(edg,avg)

Nx = size(edg,2);

if nargin == 1
    avg = 0;
end

if avg == 1    % harmonic averaging
   cntravg = 1./edg(:,1:Nx-1) + 1./edg(:,2:Nx);
   cntravg = 2./cntravg;
else           % arithmetic averaging
   cntravg = (edg(:,1:Nx-1) + edg(:,2:Nx))/2;
end

end
