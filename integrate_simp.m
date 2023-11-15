function [y, ind] = integrate_simp(datarow,y0,h)
%[y, ind] = integrate_simp(datarow,y0,h)
% y is integral
% ind is a vector of new indices
    x = datarow;
    N = length(x);
    ind = 1:2:N;
    y = zeros(1,length(ind));
    y(1) = y0;
    h3 = h/3;
    for i = 3:2:N %Adams-Moulton 3
        y(i) = y(i-2) + h3*(x(i) + 4*x(i-1) + x(i-2)); %simps
    end
    y = y(ind);
end