% Computes the averages of values at vertical edges of a MAC grid (i.e. the 
% edges where the x-direction velocity is defined) to the horizontal edges
% (i.e. where the y-direction velocity is defined).
% avg = 0  corresponds to harmonic averaging.
% avg = 1 corresponds to arithmetic averaging.
function cntravg = avgXEdgToYEdg(edg,avg)

[Ny,Nx] = size(edg);

if nargin == 1
    avg = 0;
end

if avg == 1    % harmonic averaging
   cntravg = 1./edg(1:Ny-1,:) + 1./edg(2:Ny,:);
   cntravg = 4./(cntravg(:,1:Nx-1) + cntravg(:,2:Nx));
else           % arithmetic averaging
   cntravg = (edg(1:Ny-1,:) + edg(2:Ny,:));
   cntravg = (cntravg(:,1:Nx-1) + cntravg(:,2:Nx))/4;
end

end