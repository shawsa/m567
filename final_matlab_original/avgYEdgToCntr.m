% Computes the averages of values at horizontal edges of a MAC grid (i.e. the 
% edges where the y-direction velocity is defined) to the cell centers.
% avg = 0  corresponds to harmonic averaging.
% avg = 1 corresponds to arithmetic averaging.
function cntravg = avgYEdgToCntr(edg,avg)

Ny = size(edg,1);

if nargin == 1
    avg = 0;
end

% Use two-point averages

if avg == 1    % harmonic averaging
   cntravg = 1./edg(1:Ny-1,:) + 1./edg(2:Ny,:);
   cntravg = 2./cntravg;
else           % arithmetic averaging
   cntravg = (edg(1:Ny-1,:) + edg(2:Ny,:))/2;
end

end
