function y = integrate_Euler(datarow,y0,h)
    x = datarow;
    N = length(x);
    y = zeros(1,N);
    y(1) = y0;
    for i = 2:N
        y(i) = y(i-1) + h*x(i-1);
    end
end