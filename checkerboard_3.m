function C = checkerboard_3(m, n)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

  % generate the parity map
  px = mod(1: m, 2);
  p = mod(1 : n, 2);
  % pass the xor operator, a column and a row vector
  % containing the parity data
  C = bsxfun(@xor, px', p); 
  
  figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
  % Get rid of tool bar and pulldown menus that are along top of figure.
  set(gcf, 'Toolbar', 'none', 'Menu', 'none'); imshow(C);
end
