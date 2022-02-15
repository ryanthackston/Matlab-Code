function [C, K] = checkerb(m, n)
%UNTITLED5 Summary of this function goes here
% m = 19; n = 21;
%   m = 150;
%   n = 250;

  % generate the parity map
  px = mod(1: m, 2);
  p = mod(1 : n, 2);
  % pass the xor operator, a column and a row vector
  % containing the parity data
  C = bsxfun(@xor, px', p); 
  
  figure(1); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
  % Get rid of tool bar and pulldown menus that are along top of figure.
  set(gcf, 'Toolbar', 'none', 'Menu', 'none'); 
  imshow(C, 'Border', 'tight');
  
  K = flip(C);
  
end
