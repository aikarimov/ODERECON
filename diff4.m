function y = diff4(varargin)
% DIFF4 find numerical derivative with 4 order accuracy
%   Y = DIFF4(X) differentiates the X with stepsize 1, X is N x DIM array, N >= 5 is the number of observations, DIM is deminsion
%
%   Y = DIFF4(X, T) differentiates X in times defined in T, T is 1 x N array  

x = varargin{1,1};
[N, M] = size(x);
y = zeros(N,M);

if nargin == 2
    t = varargin{1,2};
    for i = 1:N
         %generate time stepping coefficients
         B = zeros(5,1);
         
         switch i
             case {1,2} %right differences
                 irow = 1:5;
                 h = t(2) - t(1);
             case {N - 1, N} %left differences
                 irow = N - 4:N;
                 h = t(N) - t(N-1);
             otherwise %center differences
                 irow = i - 2:i + 2;
                 h = t(i+1) - t(i);
         end
         n = (t(irow) - t(i))/h;
         
         Mtr = rot90(vander(n)); %Vandermonde matrix
         B(2) = 1;
         A = Mtr\B;
         y(i,:) = A'*x(irow,:)/h; %calculate new value of y
    end
else 
    p12 = 1/12;
    
    y(1,:) = -25/12*x(1,:) + 4*x(2,:) - 3*x(3,:) + 4/3*x(4,:) - 1/4*x(5,:);
    y(2,:) = -1/4*x(1,:) - 5/6*x(2,:) + 3/2*x(3,:) - 1/2*x(4,:) + 1/12*x(5,:);

    for i = 3:N - 2
        %-1/12	2/3	 0	-2/3	1/12
        y(i,:) = p12*(-x(i + 2,:) + 8*x(i + 1,:) - 8*x(i - 1,:) + x(i - 2,:));  
    end

    y(N-1,:) = 1/4*x(N,:) + 5/6*x(N-1,:) - 3/2*x(N-2,:) + 1/2*x(N-3,:) - 1/12*x(N-4,:);
    y(N,:)   = 25/12*x(N,:) - 4*x(N-1,:) + 3*x(N-2,:) - 4/3*x(N-3,:) + 1/4*x(N-4,:);
end
   
end