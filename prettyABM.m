function prettyABM(H,T)
%PRETTYABM prints results of ABM - LSM reconstruction in a simple
%print-like form
%   PRETTYABM(H,T)
%   H is a cell array of size 1 x M, where M is dimension, containing lists
%   of coefficients by terms
%   T is a cell array of size 1 x M, where M is dimension, containing basis
%   terms obtained by delMinorTerms
[~, M] = size(T);
for i = 1:M %loop by number of functions
    str = ['f' , num2str(i), ' = ']; %string for entries
    h = H{1,i};
    t = T{1,i};
    [N, ~] = size(h); %number of terms
    [~, Mt] = size(t);
    for j = 1:N %loop by number of terms
        flag1 = 0; %flag, if 1 then do not draw * before monomial
        if j > 1 %formatting of the line
            if h(j) > 0
                if ~isequal(num2str(h(j)),'1')
                    str = [str, ' + ', num2str(h(j))];
                else
                    str = [str, ' +'];
                    flag1 = 1;
                end
            else
                if ~isequal(num2str(h(j)),'-1')
                    str = [str, ' - ', num2str(-h(j))];
                else
                    str = [str, ' -'];
                    flag1 = 1;
                end
            end
        else
            if isequal(num2str(h(j)),'-1')
                str = [str, '-'];
                flag1 = 1;
            else
                if ~isequal(num2str(h(j)),'1')
                    str = [str, num2str(h(j))];
                else
                    flag1 = 1;
                end
            end
        end
        for k = 1:Mt %in each term, display entries
            if t(j,k) ~= 0
                if flag1
                    str = [str, ' x', num2str(k)]; %x1, x2, x3 ...
                    flag1 = 0;
                else
                    str = [str, '*x', num2str(k)]; %x1, x2, x3 ...
                end
                if t(j,k) ~= 1
                    str = [str, '^', num2str(t(j,k))];% x1^2, x2^3
                end
            end
        end
    end
    disp(str);
end
end