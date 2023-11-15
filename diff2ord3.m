function y = diff2ord3(datarow)
    x = datarow;
    p2 = 1/2;
    [N, M] = size(x);
    y = zeros(N,M);
    for i = 1:3 %2 ord
        y(i,:) = -5/2*x(i,:) + 9*x(i+1,:) - 12*x(i+2,:) + 7*x(i+3,:) - 3/2*x(i+4,:); 
    end
    for i = 4:N - 3 %4 ord
        %-1/12	16/12	 -30/12	16/12	-1/12
        y(i,:) = p2*(-x(i - 2,:) + 2*x(i - 1,:) - 2*x(i + 1,:) + x(i + 2,:));  
    end
    for i = N - 2:N %4 ord
        y(i,:) = 5/2*x(i,:) - 9*x(i-1,:) + 12*x(i-2,:) - 7*x(i-3,:) + 3/2*x(i-4,:); 
    end
end