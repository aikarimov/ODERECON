function y = diff2ord2(datarow)
    x = datarow;
    [N, M] = size(x);
    y = zeros(N,M);
    for i = 1:1
        y(i,:) = 2*x(i,:) - 5*x(i+1,:) + 4*x(i+2,:) - 1*x(i+3,:); %2 ord
    end
    for i = 2:N - 1 %2
        y(i,:) = x(i + 1,:) - 2*x(i,:) + x(i - 1,:);  
    end
    for i = N:N
        y(i,:) = 2*x(i,:) - 5*x(i-1,:) + 4*x(i-2,:) - 1*x(i-3,:); %2 ord
    end
end