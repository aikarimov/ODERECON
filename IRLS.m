function h = IRLS(E,V,X,T,tol,alpha)
%IRLS iteratively reweighted least squares
% h = IRLS(E,V,X,T,tol, alpha)
% where E are monomial values, V are regression values, X are observations, T are monomials, tol is tolerance, alpha
% is the L1-regularization peramenter
% h is the vector of coefficients by monomials

[N,M] = size(E);
kmax = 100; %maximal number of iterations
wgs = ones(N,1); %weights
del = 1e-4; %minimal inverse weight
I = ones(N,1);
flag = 1; %flag for running
ctr = 1; %iteration count
h0 = zeros(M,1); %allocate memory for h
while flag && ctr < kmax
    if ctr == 1 || alpha == 0
        h = (E'*(wgs.*E))\((wgs.*E)'*V);
    else
        %optimization
        wcol = transpose(max(transpose(wgs)));
        if alpha ~= 0
            opts = optimoptions('fminunc','Display','none');
            h = fminunc(@(hx)(sum(wcol.*(V - EvalPoly(hx,X,T)).^2) + alpha*sum(abs(hx))),h,opts);
        else
            W = diag(wcol);
            h =(E'*W*E)\(E'*W*V);
        end
    end
    %update weights
    wgs = I./(max(del,abs(V - EvalPoly(h,X,T))));
    %stopping criterion:
    if norm(h - h0) < tol
        flag = false;
    end
    h0 = h;
    ctr = ctr + 1;
end
end




