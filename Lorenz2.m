function dX = Lorenz2(t,X)
dX = X; 
s = 10;
r = 28;
b = 8/3;

%Lorenz original
% x = X(1,:); y = X(2,:); z = X(3,:);
% 
% dX(1,:) = -s*x + s*y;
% dX(2,:) = -z.*x + r*x - y;
% dX(3,:) = x.*y - b*z;


%Lorenz transformed
x = X(1,:); u = X(2,:); v = X(3,:);
y = u/s + x;
znum = (- v/s - y - u + r*x);
z = znum./x;
dz = (x.*y - b*z);
dy = (r*x - znum - y);

dX(1,:) = u;
dX(2,:) = v;
dX(3,:) = s*(-v + r*u - dz.*x - z.*u - dy);

end