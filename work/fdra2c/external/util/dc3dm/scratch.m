  function [x y] = mTensorMesh (rs)
  % Get x and y so that every element of the mesh is represented.
    c = mData(rs);
    x = linspace(c.xlim(1) + 0.5*c.dx, c.xlim(2) - 0.5*c.dx,...
                 round(diff(c.xlim)/c.dx));
    y = linspace(c.ylim(1) + 0.5*c.dy, c.ylim(2) - 0.5*c.dy,...
                 round(diff(c.ylim)/c.dy));
  end
  
