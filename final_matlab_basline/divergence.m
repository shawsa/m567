% Computes the divergence of the vector field (u,v).  The value of the
% divergence is computed at the cell-centers of the MAC grid.
function div = divergence(u,v,grid)

% Compute u_x + v_y at the interior points. The divergence will live at the
% cell-centers.
div = (u(grid.juy,2:grid.mux)-u(grid.juy,1:grid.mux-1))/grid.dx + ...
      (v(2:grid.mvy,grid.jvx)-v(1:grid.mvy-1,grid.jvx))/grid.dy;
  
end