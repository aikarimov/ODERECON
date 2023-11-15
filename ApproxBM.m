function [G, O] = ApproxBM(X, eps, sigma)
%APPROXBM finds an order ideal O and a border basis G for the input data X. 
%The algorithm ensures that any linear combination of monomials in O does
%not vanish on X
%  [G O] = APPROXBM(X, V, eps, sigma) evaluates approximate border basis G and
%  order ideal O for a data point set X and reference values V.
%  O is of size K x M where M is dimension, K is number of elements in an
%  order ideal
%  G is cell array of size 1 x Q, where Q is number of elements in a border
%  basis
%  X is of size N x M, where N is number of data points
%  eps - error tolerance for linear independency
%  sigma - monomial ordering: l x M order of monomials in a form (e.g. for 2-dimensional data):
%   1    x    y    x^2   xy  y^2  x^2y  xy^2  x^2y^2
%  [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1;  1 2;  2 2] 

[N, M] = size(X);

d = 1;
O = zeros(1,M);
G = cell(0);

M = ones(N,1);

%initialize L - set of monomials of order 1
L = [];
ords = sum(sigma,2); %orders of monomials
maxord = max(ords);
lenS = length(ords);

l = 0; k = 1;
for i = 1:lenS %collect all linear terms
    %while ords(i) == 1 % collect in L all linear terms
    if ords(i) == 1 % collect in L all linear terms
        L = [L; sigma(i,:)];
        l = l + 1;
        if lenS < i %exit to ensure that ords(i) exists
            break;
        end
    end
end

gctr = 0;

while (~isempty(L)) && (d <= maxord) % limit: max degree of elements in basis is maxord
    for i = 1:l
        t = zeros(l,1); 
        t(i) = 1; %extract i-th monomial
        A = [EvalPoly(t,X,L), M];
        B = transpose(A)*A;
        [Vlam,Dlam] = eig(B);
        [lmin,I] = min(abs(diag(Dlam)));
        if sqrt(lmin) <= eps
            s = transpose(Vlam(:,I(1)));
            gctr = gctr + 1;
            G{1, gctr} = s;
        else
            O = [L(i,:); O];
            k = k + 1;
            M = A;
        end
    end
    d = d + 1;
    [L, l] = BorderBasis(O, d, sigma);
end

%order O w.r.t. sigma
ind = zeros(1,k);
for i = 1:k
    [~, ind(i)] = israwcontained(O(i,:), sigma); %find positions of L(i,:) in sigma
end
[~,I] = sort(ind); %sort positions
O = O(I,:); %re-order O

Ored = []; %reduced O with excluded elements absent in sigma
for i = 1:k
    if(ind(I(i)) ~= 0)
        Ored = [Ored; O(i,:)];
    end
end
O = Ored;
end
