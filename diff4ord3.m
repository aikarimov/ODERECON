function y = diff4ord3(datarow)
    x = datarow;
    p8 = 1/8;
    [N, M] = size(x);
    y = zeros(N,M);
    for i = 1:3 %2 ord
        y(i,:) = -49/8*x(i,:) + 29*x(i+1,:) - 461/8*x(i+2,:) + 62*x(i+3,:) - 307/8*x(i+4,:) + 13 *x(i+5,:) - 15/8*x(i+6,:); 
    end
    for i = 4:N - 3 %4 ord
        %-1/12	16/12	 -30/12	16/12	-1/12
        y(i,:) = -p8*(x(i + 3,:) - 8*x(i + 2,:) + 13*x(i + 1,:) -13*x(i - 1,:) + 8*x(i - 2,:) - x(i - 3,:));  
    end
    for i = N - 2:N %4 ord
        y(i,:) = 49/8*x(i,:) - 29*x(i-1,:) + 461/8*x(i-2,:) - 62*x(i-3,:) + 307/8*x(i-4,:) - 13 *x(i-5,:) + 15/8*x(i-6,:); 
    end
end