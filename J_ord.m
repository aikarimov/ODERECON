function J = J_ord(hi, Y1, tau, deltaY, ord)
%partial derivatives with order in a point Y1
%ord = 2
%ord = 4

[~,M] = size(Y1);

dY = zeros(1,M);
dY(:,1) = deltaY;
J = zeros(M);

if (ord == 2)
    for i = 1:M
        fp = EvalPoly(hi,Y1 + dY,tau);
        fm = EvalPoly(hi,Y1 - dY,tau);
        
        J(:,i) = 0.5*transpose(fp - fm)/deltaY; %ord 2
        dY = circshift(dY,1);
    end
end

if (ord == 4)
    for i = 1:M
        fp2 = EvalPoly(hi,Y1 + 2*dY,tau);
        fp = EvalPoly(hi,Y1 + dY,tau);
        fm = EvalPoly(hi,Y1 - dY,tau);
        fm2 = EvalPoly(hi,Y1 - 2*dY,tau);
        %1/12	2/3	 0	2/3	1/12
        J(:,i) = transpose(-1/12*fp2 + 2/3*fp - 2/3*fm + 1/12*fm2)/deltaY; %ord 4
        dY = circshift(dY,1);
    end
end

end