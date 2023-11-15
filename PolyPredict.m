function P = PolyPredict(X, H, Tau)
%  P = POLYPREDICT(X, H, Tau) evaluates polynomial in datapoints X of size N x M, M is
%  dimension, N is number of data points made up of monomials
%  p - values of polynomes in X, N x Q
%  H{1,j} - coefficients by terms t, where t[i] is a monomial of the
%  corresponding order, H{1,j} size is L x Q
%  Tau{1,j} - L x M ordered monomials (e.g., w.r.t. degree-lexicographic order) like this, for 2-dimensional data:
%   1    x    y    x^2   xy  y^2  x^2y  xy^2  x^2y^2
%  [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1;  1 2;  2 2] 

%get dimension
[~,R] = size(H);

[N,~] = size(X);

P = zeros(N,R);
for j = 1:R
    P(:,j) = EvalPoly(H{1,j},X,Tau{1,j});
end

end