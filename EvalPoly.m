function p = EvalPoly(h, X, T)
%  p = EVALPOLY(h, X, T) evaluates polynomial in datapoints X of size N x M, M is
%  dimension, N is number of data points made up of monomials
%  p - values of polynomes in X, N x Q
%  h - coefficients by terms t, where t[i] is a monomial of the
%  corresponding order, h size is L x Q
%  T - L x M ordered monomials (e.g., w.r.t. degree-lexicographic order) like this, for 2-dimensional data:
%   1    x    y    x^2   xy  y^2  x^2y  xy^2  x^2y^2
%  [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1;  1 2;  2 2] 

%get dimension
[N, M] = size(X);
[L, Q] = size(h);

p = zeros(N,Q); %evaluated polynomial

for i = 1:L
   p = p + h(i,:).*prod(X.^repmat(T(i,:),N,1),2);
end

end