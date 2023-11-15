function dX = FHN(t,X)
dX = X; x = X(1,:); y = X(2,:);

dX(1,:) = y;
dX(2,:) = - x + 4*y - y.^3;
end