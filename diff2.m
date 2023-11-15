function y = diff2(datarow)
    x = datarow;
    [N, M] = size(x);
    y = zeros(N,M);
    y(1,:) = -3/2*x(1,:) + 2*x(2,:) - 1/2*x(3,:);
    for i = 2:N - 1
        y(i,:) = (x(i+1,:) - x(i-1,:))/2;
    end
    y(N,:) = 3/2*x(N,:) - 2*x(N-1,:) + 1/2*x(N-2,:);
end