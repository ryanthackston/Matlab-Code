function [ a ] = checkerboard(m,n)
a = zeros(n);
for i = 1:m
    for j = 1:n
        if (i == j)
            a (i, j) = 0;
        elseif (mod(j, 2) == 0) && (mod(i,2) == 0)
             a(i,j) = 0;
        elseif (mod(j, 2) == 0) || (mod(i,2) == 0)
            a(i,j) = 1;
        end
    end
end

  figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
  % Get rid of tool bar and pulldown menus that are along top of figure.
  set(gcf, 'Toolbar', 'none', 'Menu', 'none'); imshow(a);

end

