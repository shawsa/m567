% Computes the 2D vorticity of the vector field (u,v).  The value of the
% vorticity is computed at the cell-corners of the MAC grid.
function vort = vorticity(u,v,grid)

% Compute u_y - v_x at the interior points. The vorticity will live at the
% cell-corners.
vort = (u(2:grid.muy,:)-u(1:grid.muy-1,:))/grid.dy - ...
      (v(:,2:grid.mvx)-v(:,1:grid.mvx-1))/grid.dx;
  
end