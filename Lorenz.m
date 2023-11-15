function dX = Lorenz(t,X)
dX = X; 
%p = [1, 31, 25];
% x = X(1,:); y = X(2,:); z = X(3,:);
p = [8/3, 28, 10];
% dX(1,:) = y.*z - p(1).*x;
% dX(2,:) = -x.*z + p(2).*z - y;
% dX(3,:) = -p(3).*z + p(3).*y;

s = 10;
r = 28;
b = 8/3;

x = X(1,:); y = X(2,:); z = X(3,:);

dX(1,:) = -s*x + s*y;
dX(2,:) = -z.*x + r*x - y;
dX(3,:) = x.*y - b*z;
end