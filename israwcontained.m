function [b, I] = israwcontained(V, W)
% ISRAWCONTAINED determines whether a row vector V is contained in W and
% returns its position
%[b, I] = ISRAWCONTAINED(V, W) returns logical b, true if raw V is contained in a matrix W,
%otherwise false, and I - the position of V in W (or 0 if V is not contained in W)
%   V is a raw vector of size 1 x M
%   W is a matrix of size N x M

b = 0;
I = 0;
[N, M] = size(W);

for k = 1:N
    if ( isequal(V , W(k,:))) % if V is contained in W
        b = 1;
        I = k;
    end
end

b = logical(b);