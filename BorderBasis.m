function [L, l]= BorderBasis(O, d, sigma)
%  [L, l] = BORDERBASIS(O, sigma, d) returns L - all border elements of order ideal
%  O of order d, where O is of size [K x M], M is dimension, K is number of elements in O,
%  l is number of elements in L
%  sigma - l x M order of monomials in a form (e.g. for 2-dimensional data):
%   1    x    y    x^2   xy  y^2  x^2y  xy^2  x^2y^2
%  [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1;  1 2;  2 2] 

[K, M] = size(O);
l = 0; %number of elements in a border basis
L = [];
order = sum(O, 2);
for i = 1:K %for every element in O
    if( order(i) == d - 1) %if its order is the highest for border
        for j = 1:M %try to append to each element by dimension
            s = zeros(1,M); s(j) = 1;
            cand = O(i,:) + s; %candidate monomial, check if it is not contained in L or O
            if ~(israwcontained(cand, L) && israwcontained(cand, L))
                L = [cand; L];
                l = l + 1;
            end
        end
    end
end

%order L w.r.t sigma
ind = zeros(1,l);
for i = 1:l
    [~, ind(i)] = israwcontained(L(i,:), sigma); %find positions of L(i,:) in sigma
end
[~,I] = sort(ind); %sort positions
L = L(I,:); %re-order L

end