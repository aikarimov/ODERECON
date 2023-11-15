function [h , T] = delMinorTerms(varargin)
%DELMINORTERMS constructs and interpolating polynomial on the order ideal O
% [h, T] = DELMINORTERMS(X,V,O,eta)
%
% [h, T] = DELMINORTERMS(X,V,O,eta, h0)
%
% [h, T] = DELMINORTERMS(X,V,O,eta, h0, deleteminor)
% 
% [h, T] = DELMINORTERMS(X,V,O,eta, h0, deleteminor,alpha)
% 
% [h, T] = DELMINORTERMS(X,V,O,eta, h0, deleteminor,alpha, irls)
% 
% returns h - set of coefficients by ideal
%monomials in O, T is the reduced ideal
% X is N x M matrix, where N is the number of data points, M is input dimension
% V is N x 1 matrix
% O is a sorted order ideal of size L x M, e.g. in 2-dimesional case it may be
%   1    x    y    x^2   xy  y^2  x^2y  xy^2  x^2y^2
%  [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1;  1 2;  2 2] 
% eta is a tolerance
% h0 is an initial value of h
% set deleteminor to 0 if you want to preserve the same size of T as O, 1
% otherwise (default)
% alpha is regularization term for L1-regularization
% irls is a flag whether we use IRLS (irls = 1) or OLS (irls = 0)

alp = 0.0; %regularization term (for L1 regularization)
irls = 0;

deleteminor = 1;

if nargin >= 4
    X = varargin{1,1}; [N, ~] = size(X);
    V = varargin{1,2};
    O = varargin{1,3}; [L, M] = size(O);
    eta = varargin{1,4};
end

if nargin >= 5
    h = varargin{1,5};
end

if nargin >= 8
    irls = varargin{1,8};
end

if nargin == 4
    E = EvalPoly(eye(L),X,O);
    if irls
        h = IRLS(E,V,X,O,eta,alp);
    else
        h = (E'*E)\(E'*V);
    end
end

if nargin == 6
    deleteminor = varargin{1,6};
end

if nargin == 7
    alp = varargin{1,7};
end




T = O;
htmp = h;
Ttmp = T;

zeroingT = ones(L,1);
reindex = 1:L;
L0 = L;

while 1/N*norm(V - EvalPoly(htmp,X,Ttmp)) <= eta && L > 1
    h = htmp;
    T = Ttmp;
    
    %find a minimal norm of the monomial in interpolation
    minval = inf;
    mink = 1;
    for k = 1:L
        t = zeros(L,1);
        t(k) = h(k); %extract k-th monomial
        tmp = norm(EvalPoly(t,X,T));
        if tmp < minval
            minval = tmp;
            mink = k;
        end
    end
    
    %exclude mink-th element from T
    if L > 1
        zeroingT(reindex(mink)) = 0;
        
        if mink > 1 && mink < L
            Ttmp = T([1:mink - 1,mink + 1:end],:);
            reindex = reindex([1:mink - 1,mink + 1:end]);
        end
        if mink == L
            Ttmp = T(1:mink - 1,:);
            reindex = reindex(1:mink - 1);
        end
        if mink == 1
            Ttmp = T(mink + 1:end,:);
            reindex = reindex(mink + 1:end);
        end
        L = L - 1;
       
    end
    E = EvalPoly(eye(L),X,Ttmp);
    
    if irls
        htmp = IRLS(E,V,X,Ttmp,eta,alp);
    else
        %htmp = (E'*E + eta)\(E'*V);
        [Q,R] = qr(E);
        Q1 = Q(:,1:L);
        R1 = R(1:L,1:L);
        htmp = R1\(Q1'*V);


        if alp ~= 0
            opts = optimoptions('fminunc','Display','none');
            htmp = fminunc(@(h)objective1(h,V, X,Ttmp,alp),htmp,opts);
        end
    end
end

if 1/N*norm(V - EvalPoly(htmp,X,Ttmp)) <= eta && L == 1
    h = htmp;
    T = Ttmp;
end

if ~deleteminor
    htmp = zeros(L0,1);
    ctr = 1;
    for i = 1:L0
        if( zeroingT(i) ~= 0) && (ctr <= length(h))
            htmp(i) = h(ctr);
            ctr = ctr + 1;
        end
    end
    h = htmp;
    T = O;
end
end

function fun = objective1(h,V,X,Ttmp,alpha)
%objective for L1-regularized LSM
    fun = (norm(V - EvalPoly(h,X,Ttmp)))^2 + alpha*sum(abs(h));
end