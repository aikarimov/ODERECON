function y = diff4ord2(datarow)
    x = datarow;
    p12 = 1/12;
    [N, M] = size(x);
    y = zeros(N,M);
    for i = 1:2
        y(i,:) = 15/4*x(i,:) - 77/6*x(i+1,:) + 107/6*x(i+2,:) - 13*x(i+3,:) + 61/12*x(i+4,:) - 5/6*x(i+5,:); %4 ord
    end
    for i = 3:N - 2 %4 ord
        %-1/12	16/12	 -30/12	16/12	-1/12
        y(i,:) = p12*(-x(i + 2,:) + 16*x(i + 1,:) - 30*x(i,:) + 16*x(i - 1,:) - x(i - 2,:));  
    end
    for i = N-1:N
        y(i,:) = 15/4*x(i,:) - 77/6*x(i-1,:) + 107/6*x(i-2,:) - 13*x(i-3,:) + 61/12*x(i-4,:) - 5/6*x(i-5,:); %4 ord
    end
end