% Computes the averages of values at horizontal edges of a MAC grid (i.e. the 
% edges where the y-direction velocity is defined) to the vertical edges
% (i.e. where the x-direction velocity is defined).
% avg = 0  corresponds to harmonic averaging.
% avg = 1 corresponds to arithmetic averaging.
function cntravg = avgYEdgToXEdg(edg,avg)

[Ny,Nx] = size(edg);

if nargin == 1
    avg = 0;
end

% Use four-point averages
if avg == 1    % harmonic averaging
   cntravg = 1./edg(:,1:Nx-1) + 1./edg(:,2:Nx);
   cntravg = 4./(cntravg(1:Ny-1,:) + cntravg(2:Ny,:));
else           % arithmetic averaging
   cntravg = (edg(:,1:Nx-1) + edg(:,2:Nx));
   cntravg = (cntravg(1:Ny-1,:) + cntravg(2:Ny,:))/4;
end

end