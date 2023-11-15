function [y, ind] = integrate_bool(datarow,y0,h)
%[y, ind] = integrate_simp(datarow,y0,h)
% y is integral
% ind is a vector of new indices
    x = datarow;
    N = length(x);
    ind = 1:4:N;
    y = zeros(1,length(ind));
    y(1) = y0;
    h45 = 2*h/45;
    for i = 5:4:N %bool 5
        y(i) = y(i-4) + h45*(7*x(i) + 32*x(i-1) + 12*x(i-2) + 32*x(i-3) + 7*x(i - 4));
    end
    y = y(ind);
end