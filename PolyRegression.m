function [H, Tau, err] = PolyRegression(varargin)
%POLYREGRESSION reconstructs an M-dimensional system using Kera-Hasegawa process
%get uniformly distributed points from data
%   Syntax:
%
%   [H, Tau, err] = POLYREGRESSION(Y, W, dmax)
%
%   [H, Tau, err] = POLYREGRESSION(Y, W, O)
%
%   [H, Tau, err] = POLYREGRESSION(Y, W, dmin, dmax)
%
%   [H, Tau, err] = POLYREGRESSION(Y, W, dmin, dmax, eta)
%
%   [H, Tau, err] = POLYREGRESSION(Y, W, dmin, dmax, eta, eps)
%
%   [H, Tau, err] = POLYREGRESSION(Y, W, dmin, dmax, eta, eps, deleteminor)
%
%   [H, Tau, err] = POLYREGRESSION(Y, W, dmin, dmax, eta, eps, deleteminor, useirls)
%
%   Input:
%
%   Y is the N x M vector set of system independent variables
%   W is the N x M vector set of system dependent variables (for ODE, derivatives)
%   O is L x M matrix if order ideal: set of monomials
%   dmin - is min order of monomials in the system, dmin <= 0
%   dmax - max order of monomials in the system
%   eta - tolerance for delMinorTerms algorithm
%   eps - tolerance for ApproxBM algorithm
%   deleteminor - should we delete minor terms and reduce order
%   useirls - should we use IRLS or ordirary LSM
%
%   Output:
%
%   H - cell array of 1 x M dimension, each cell contains column vector of coefficients by monomials;
%   Tau - cell array of 1 x M dimension, each cell contains T x M matrix
%   where each row corresponts to the monomial, T is th number of monomials
%   err - error in polynomial estimation
deleteminor = 1;
useirls = 0;

eps = 1e-7;
eta = 1e-5;

O = [];

if nargin == 3
    dmin = 0;
    Y = varargin{1,1};
    W = varargin{1,2};
    if numel(varargin{1,3}) == 1
        dmax = varargin{1,3};
    else
        O = varargin{1,3};
    end
end

if nargin == 4
    Y = varargin{1,1};
    W = varargin{1,2};
    dmin = varargin{1,3};
    dmax = varargin{1,4};
end

if nargin == 5
    Y = varargin{1,1};
    W = varargin{1,2};
    dmin = varargin{1,3};
    dmax = varargin{1,4};
    eta = varargin{1,5};
end

if nargin >= 6
    Y = varargin{1,1};
    W = varargin{1,2};
    dmin = varargin{1,3};
    dmax = varargin{1,4};
    eta = varargin{1,5};
    eps = varargin{1,6};
end

if nargin >= 7
    deleteminor = varargin{1,7};
end

if nargin >= 8
    useirls = varargin{1,8};
end

rng default;
[N,M] = size(Y);

%reconstruct order ideal

if numel(O) == 0
    sigma = deglexord(dmin,dmax,M);
    [~, O] = ApproxBM(Y, eps, sigma); %use approximate Buchberger-Moller algorithm
end

%Use LSM for fitting the equations with the proper coefficients

H = cell(1,M);
Tau = cell(1,M);

err = 0; fsum = 0;
Vs = W;
%reconstruct each equation
for i = 1:M
    [L, ~] = size(O);
    cscl = max(abs(Y));
    E = EvalPoly(eye(L),Y./cscl,O);   

    V = W(:,i);
    cscl2 = EvalPoly(eye(L),cscl,O); 
    if ~useirls
        %hi = (E'*E)\(E'*V)./cscl2'; %initial approximation
        [Q,R] = qr(E);
        Q1 = Q(:,1:L);
        R1 = R(1:L,1:L);
        hi = R1\(Q1'*V)./cscl2';
    else
        hi = IRLS(E,V,Y,O,eta,0);
        hi = hi./cscl2';
    end

    tau = O;
    if deleteminor
        [hi,tau] = delMinorTerms(Y,W(:,i),O,eta,hi,1,0,useirls); %get equation and basis
    end
    Vs(:,i) = EvalPoly(hi,Y,tau);
    H{1,i} = hi;
    Tau{1,i} = tau;
end

for i = 1:N
    err = err + norm(Vs(i,:) - W(i,:)); %by point Y[i], estimate (dY_e - dY_actual)
    fsum = fsum + norm(Vs(i,:));
end
err = err/fsum;
end