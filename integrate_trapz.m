function y = integrate_trapz(datarow,y0,htime)
    N = length(datarow);
    if length(htime) == 1
        h = htime;
        htime = (0:(N - 1))*h;
    end
    x = datarow;
    N = length(x);
    y = zeros(1,N);
    y(1) = y0;
    for i = 2:N
        h = htime(i) - htime(i-1);
        y(i) = y(i-1) + 0.5*h*(x(i-1) + x(i));
    end
end